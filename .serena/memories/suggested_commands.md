# Suggested commands

- Configure/build existing native tree: `cmake --build build -j2` (after its environment has a non-editable Provider wheel).
- Compact tests: `pytest -q`.
- Focused live HELD2 trace: `pytest -s tests/test_perdomo_held2_trace.py::<test-name> --held2-live`.
- Pure-saturation validation anchor: `python scripts/validate_saturation.py`.
- Ruff: `ruff check .`; strict typing: `mypy src/epcsaft_equilibrium`.
- Final cleanup audit after mutations/processes: `bash "$HOME/.codex/hooks/codex-cleanup.sh" --repo-root .`.
- Prefer bounded focused tests before the full suite; pin BLAS/OpenMP threads for subprocess-heavy numerical runs when needed.