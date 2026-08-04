from __future__ import annotations

METHANE_ETHANE_PARAMETERS = {
    "schema": "epcsaft.parameters",
    "schema_version": 1,
    "components": ("methane", "ethane"),
    "parameters": {
        "mw": (0.016043, 0.03007),
        "m": (1.0, 1.6069),
        "s": (3.7039, 3.5206),
        "e": (150.03, 191.42),
        "k_ij": ((0.0, 0.0), (0.0, 0.0)),
    },
    "options": {"permittivity_model": "none"},
    "validity": {"kind": "unknown"},
}


METHANE_PARAMETERS = {
    "schema": "epcsaft.parameters",
    "schema_version": 1,
    "components": ("methane",),
    "parameters": {
        "mw": (0.016043,),
        "m": (1.0,),
        "s": (3.7039,),
        "e": (150.03,),
        "k_ij": ((0.0,),),
    },
    "options": {"permittivity_model": "none"},
    "validity": {"kind": "unknown"},
}

ETHANE_PARAMETERS = {
    "schema": "epcsaft.parameters",
    "schema_version": 1,
    "components": ("ethane",),
    "parameters": {
        "mw": (0.03007,),
        "m": (1.6069,),
        "s": (3.5206,),
        "e": (191.42,),
        "k_ij": ((0.0,),),
    },
    "options": {"permittivity_model": "none"},
    "validity": {"kind": "unknown"},
}

PROPANE_PARAMETERS = {
    "schema": "epcsaft.parameters",
    "schema_version": 1,
    "components": ("propane",),
    "parameters": {
        "mw": (0.044096,),
        "m": (2.0020,),
        "s": (3.6184,),
        "e": (208.11,),
        "k_ij": ((0.0,),),
    },
    "options": {"permittivity_model": "none"},
    "validity": {"kind": "unknown"},
}

PURE_SATURATION_PARAMETERS = {
    "methane": METHANE_PARAMETERS,
    "ethane": ETHANE_PARAMETERS,
    "propane": PROPANE_PARAMETERS,
}


FIGIEL_REFERENCE_ELECTROLYTE_PARAMETERS = {
    "schema": "epcsaft.parameters",
    "schema_version": 1,
    "components": ("water", "sodium-cation", "chloride-anion"),
    "parameters": {
        "mw": (0.0180153, None, None),
        "m": (1.2047, 1.0, 1.0),
        "s": (None, 2.8232, 2.7560),
        "e": (353.95, 230.0, 170.0),
        "z": (0, 1, -1),
        "d_born": (None, 3.445, 4.100),
        "f_solv": (1.5, 1.0, 1.0),
        "k_ij": (
            (0.0, -0.3, -0.3),
            (-0.3, 0.0, 0.8),
            (-0.3, 0.8, 0.0),
        ),
        "sites": (
            {
                "component_id": "water",
                "site_id": "a",
                "site_class": "donor",
                "multiplicity": 1,
            },
            {
                "component_id": "water",
                "site_id": "b",
                "site_class": "acceptor",
                "multiplicity": 1,
            },
        ),
        "association": (
            {
                "component_id_a": "water",
                "site_id_a": "a",
                "component_id_b": "water",
                "site_id_b": "b",
                "association_energy_over_k": 2425.7,
                "association_volume": 0.04509,
            },
        ),
        "correlations": (
            {
                "component_id": "water",
                "family": "segment_diameter",
                "form": "constant-plus-sum-of-exponentials",
                "constant": 2.7927,
                "terms": (
                    {"amplitude": 10.11, "exponent_coefficient_per_k": -0.01775},
                    {"amplitude": -1.417, "exponent_coefficient_per_k": -0.01146},
                ),
            },
            {
                "component_id": "water",
                "family": "relative_permittivity",
                "form": "constant",
                "constant": 78.09,
            },
        ),
    },
    "options": {
        "epsilon_r_ion": 8.0,
        "a_dh": 7.01,
        "permittivity_model": "ion-fraction-suppression",
    },
    "validity": {
        "kind": "reported-conditions",
        "temperature_min_k": 298.15,
        "temperature_max_k": 298.15,
        "pressure_min_pa": 100000.0,
        "pressure_max_pa": 100000.0,
        "ion_mole_fraction_max": 0.38,
    },
}


KHUDAIDA_HELD2_PARAMETERS = {
    "schema": "epcsaft.parameters",
    "schema_version": 1,
    "components": (
        "water",
        "ethanol",
        "isobutanol",
        "sodium-cation",
        "chloride-anion",
    ),
    "parameters": {
        "molar_mass": (0.01801528, 0.046068, 0.0741216, 0.02298, 0.03545),
        "segment_count": (1.2047, 2.3827, 2.844, 1.0, 1.0),
        "segment_diameter": (None, None, 3.5561, 2.8232, 2.756),
        "dispersion_energy_over_k": (353.95, 198.24, 257.13, 230.0, 170.0),
        "charge_number": (0, 0, 0, 1, -1),
        "born_diameter": (None, None, None, 3.445, 4.1),
        "solvation_factor": (1.5, 1.6, 1.5, 1.0, 1.0),
        "dispersion_interaction": (
            (0.0, -0.06167, -0.0071, -0.3, -0.3),
            (-0.06167, 0.0, 0.00154, 0.05, 0.8),
            (-0.0071, 0.00154, 0.0, 0.0, 0.0),
            (-0.3, 0.05, 0.0, 0.0, 0.317),
            (-0.3, 0.8, 0.0, 0.317, 0.0),
        ),
        "sites": (
            {"component_id": "ethanol", "site_id": "a", "site_class": "a", "multiplicity": 1},
            {"component_id": "ethanol", "site_id": "b", "site_class": "b", "multiplicity": 1},
            {"component_id": "isobutanol", "site_id": "a", "site_class": "a", "multiplicity": 1},
            {"component_id": "isobutanol", "site_id": "b", "site_class": "b", "multiplicity": 1},
            {"component_id": "water", "site_id": "a", "site_class": "a", "multiplicity": 1},
            {"component_id": "water", "site_id": "b", "site_class": "b", "multiplicity": 1},
        ),
        "association": (
            {
                "component_id_a": "ethanol",
                "site_id_a": "a",
                "component_id_b": "isobutanol",
                "site_id_b": "b",
                "association_energy_over_k": 2583.85,
                "association_volume": 0.01119849133244554,
            },
            {
                "component_id_a": "ethanol",
                "site_id_a": "a",
                "component_id_b": "ethanol",
                "site_id_b": "b",
                "association_energy_over_k": 2653.4,
                "association_volume": 0.03238,
            },
            {
                "component_id_a": "ethanol",
                "site_id_a": "b",
                "component_id_b": "isobutanol",
                "site_id_b": "a",
                "association_energy_over_k": 2583.85,
                "association_volume": 0.01119849133244554,
            },
            {
                "component_id_a": "isobutanol",
                "site_id_a": "a",
                "component_id_b": "isobutanol",
                "site_id_b": "b",
                "association_energy_over_k": 2514.3,
                "association_volume": 0.00391,
            },
            {
                "component_id_a": "ethanol",
                "site_id_a": "b",
                "component_id_b": "water",
                "site_id_b": "a",
                "association_energy_over_k": 2539.55,
                "association_volume": 0.03798098222932118,
            },
            {
                "component_id_a": "isobutanol",
                "site_id_a": "b",
                "component_id_b": "water",
                "site_id_b": "a",
                "association_energy_over_k": 2470.0,
                "association_volume": 0.01299623825352431,
            },
            {
                "component_id_a": "water",
                "site_id_a": "a",
                "component_id_b": "water",
                "site_id_b": "b",
                "association_energy_over_k": 2425.7,
                "association_volume": 0.04509,
            },
            {
                "component_id_a": "ethanol",
                "site_id_a": "a",
                "component_id_b": "water",
                "site_id_b": "b",
                "association_energy_over_k": 2539.55,
                "association_volume": 0.03798098222932118,
            },
            {
                "component_id_a": "isobutanol",
                "site_id_a": "a",
                "component_id_b": "water",
                "site_id_b": "b",
                "association_energy_over_k": 2470.0,
                "association_volume": 0.01299623825352431,
            },
        ),
        "correlations": (
            {
                "component_id": "ethanol",
                "family": "relative_permittivity",
                "form": "constant",
                "constant": 24.88,
            },
            {
                "component_id": "ethanol",
                "family": "segment_diameter",
                "form": "constant",
                "constant": 3.1771,
            },
            {
                "component_id": "isobutanol",
                "family": "relative_permittivity",
                "form": "constant",
                "constant": 17.2,
            },
            {
                "component_id": "water",
                "family": "relative_permittivity",
                "form": "constant",
                "constant": 78.09,
            },
            {
                "component_id": "water",
                "family": "segment_diameter",
                "form": "constant-plus-sum-of-exponentials",
                "constant": 2.7927,
                "terms": (
                    {"amplitude": 10.11, "exponent_coefficient_per_k": -0.01775},
                    {"amplitude": -1.417, "exponent_coefficient_per_k": -0.01146},
                ),
            },
        ),
    },
    "options": {
        "ion_fraction_suppression_coefficient": 7.01,
        "ionic_region_relative_permittivity": 8.0,
        "relative_permittivity_formulation": "ion-fraction-suppression",
    },
    "validity": {
        "kind": "reported-conditions",
        "temperature_min_k": 293.15,
        "temperature_max_k": 313.15,
        "pressure_min_pa": 100000.0,
        "pressure_max_pa": 100000.0,
    },
}
