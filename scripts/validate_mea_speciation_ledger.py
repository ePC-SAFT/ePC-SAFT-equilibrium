from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LEDGER = REPO_ROOT / "data" / "reference" / "mea_speciation_input_ledger.json"
ALLOWED_OBSERVATION_ROLES = {
    "aggregate",
    "balance-inferred",
    "calibration",
    "contextual",
    "direct",
    "excluded",
    "held-out_validation",
}


class LedgerError(ValueError):
    pass


def require(record: dict[str, Any], field: str, path: str) -> Any:
    if field not in record or record[field] in (None, "", []):
        raise LedgerError(f"missing {path}.{field}")
    return record[field]


def exact_rank(matrix: list[list[int | float]]) -> int:
    rows = [[Fraction(value) for value in row] for row in matrix]
    if not rows:
        return 0
    rank = 0
    column_count = len(rows[0])
    for column in range(column_count):
        pivot = next((row for row in range(rank, len(rows)) if rows[row][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][column]
        rows[rank] = [value / pivot_value for value in rows[rank]]
        for row in range(len(rows)):
            if row == rank or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                value - factor * pivot_entry
                for value, pivot_entry in zip(rows[row], rows[rank], strict=True)
            ]
        rank += 1
        if rank == len(rows):
            break
    return rank


def validate_records(document: dict[str, Any]) -> None:
    for field in ("identity", "status", "scope", "authority"):
        require(document, field, "ledger")

    sections = {
        "species": ("identity", "name", "formula", "charge", "unit", "basis", "source"),
        "reactions": ("identity", "equation", "stoichiometry", "unit", "basis", "source"),
        "equilibrium_constant_candidates": (
            "identity",
            "reaction",
            "coefficients",
            "temperature_range_K",
            "unit",
            "basis",
            "standard_state",
            "source",
            "role",
            "status",
        ),
        "observation_sets": ("identity", "temperature_K", "unit", "basis", "source", "role"),
        "historical_model_evidence": ("identity", "source", "role", "classification"),
        "sources": (
            "identity",
            "doi",
            "zotero_parent_key",
            "zotero_attachment_key",
            "file_sha256",
            "role",
        ),
        "blockers": ("identity", "scope", "resolution"),
    }
    for section, fields in sections.items():
        records = require(document, section, "ledger")
        identities: set[str] = set()
        for index, record in enumerate(records):
            path = f"{section}[{index}]"
            for field in fields:
                require(record, field, path)
            identity = record["identity"]
            if identity in identities:
                raise LedgerError(f"duplicate {path}.identity: {identity}")
            identities.add(identity)

    for field in ("identity", "unit", "basis", "source"):
        require(document["balances"], field, "balances")
        require(document["feed"], field, "feed")
        require(document["provider"], field, "provider")
    for field in (
        "identity",
        "temperature_unit",
        "temperature_basis",
        "pressure_unit",
        "pressure_basis",
        "source",
    ):
        require(document["state_grid"], field, "state_grid")
    for field in (
        "identity",
        "expression",
        "temperature_unit",
        "result_unit",
        "basis",
        "standard_state",
        "source",
    ):
        require(document["equilibrium_constant_formula"], field, "equilibrium_constant_formula")


def validate_chemistry(document: dict[str, Any]) -> tuple[int, int, int]:
    species = document["species"]
    names = [record["name"] for record in species]
    if names != document["species_order"]:
        raise LedgerError("species order does not match the ordered species records")
    species_count = len(species)
    if species_count != 9:
        raise LedgerError("the MEA ledger must contain exactly nine true species")

    balances = document["balances"]["matrix"]
    reactions = [record["stoichiometry"] for record in document["reactions"]]
    if any(len(row) != species_count for row in balances + reactions):
        raise LedgerError("balance and reaction rows must match species_count")

    balance_rank = exact_rank(balances)
    reaction_rank = exact_rank(reactions)
    if balance_rank != document["balances"]["declared_rank"]:
        raise LedgerError("declared balance rank is incorrect")
    if reaction_rank != document["declared_reaction_rank"]:
        raise LedgerError("declared reaction rank is incorrect")
    if balance_rank + reaction_rank != species_count:
        raise LedgerError("balance and reaction ranks do not span the species space")

    elements = document["balances"]["row_order"]
    formula_matrix = [
        [record["formula"][element] for record in species]
        for element in elements
    ]
    if balances != formula_matrix:
        raise LedgerError("balance matrix does not match the species formula records")

    for reaction_index, reaction in enumerate(reactions):
        for balance_index, balance in enumerate(balances):
            if sum(a * b for a, b in zip(balance, reaction, strict=True)) != 0:
                raise LedgerError(
                    f"reactions[{reaction_index}] violates balance row {balance_index}"
                )
        charge_change = sum(
            record["charge"] * coefficient
            for record, coefficient in zip(species, reaction, strict=True)
        )
        if charge_change != 0:
            raise LedgerError(f"reactions[{reaction_index}] violates charge conservation")

    if document["feed"]["initial_charge_equivalents"] != 0:
        raise LedgerError("feed is not electroneutral")
    return species_count, balance_rank, reaction_rank


def validate_roles(document: dict[str, Any]) -> None:
    for index, observation in enumerate(document["observation_sets"]):
        roles = {observation["role"], *observation.get("observable_roles", {}).values()}
        unknown = roles - ALLOWED_OBSERVATION_ROLES
        if unknown:
            raise LedgerError(f"observation_sets[{index}] has unsupported roles: {sorted(unknown)}")
    for index, evidence in enumerate(document["historical_model_evidence"]):
        if evidence["role"] != "model_evidence":
            raise LedgerError(f"historical_model_evidence[{index}] is not labeled model evidence")


def validate_blockers(document: dict[str, Any]) -> None:
    blockers = {record["identity"] for record in document["blockers"]}
    referenced: set[str] = set()
    for species in document["species"]:
        if blocker := species.get("blocker"):
            referenced.add(blocker)
    for record in (
        document["feed"],
        document["state_grid"],
        document["equilibrium_constant_formula"],
        document["provider"],
    ):
        if blocker := record.get("blocker"):
            referenced.add(blocker)
    for observation in document["observation_sets"]:
        referenced.update(observation.get("blockers", []))
    unknown = referenced - blockers
    if unknown:
        raise LedgerError(f"records refer to undeclared blockers: {sorted(unknown)}")
    if document["executable"] or document["status"] != "blocked":
        raise LedgerError("this evidence ledger must remain blocked and non-executable")


def validate(document: dict[str, Any]) -> dict[str, int | bool | str]:
    validate_records(document)
    species_count, balance_rank, reaction_rank = validate_chemistry(document)
    validate_roles(document)
    validate_blockers(document)
    return {
        "balance_rank": balance_rank,
        "blockers": len(document["blockers"]),
        "executable": document["executable"],
        "reaction_rank": reaction_rank,
        "species_count": species_count,
        "status": document["status"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the private MEA input evidence ledger")
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--require-executable", action="store_true")
    args = parser.parse_args()
    try:
        document = json.loads(args.ledger.read_text(encoding="utf-8"))
        report = validate(document)
    except (OSError, json.JSONDecodeError, LedgerError) as error:
        print(str(error), file=sys.stderr)
        return 1
    if args.require_executable and not report["executable"]:
        print("ledger is not executable: unresolved scientific blockers remain", file=sys.stderr)
        return 2
    print(json.dumps(report, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
