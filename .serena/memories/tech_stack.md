# Tech stack

- Python >=3.13 package built with scikit-build-core and pybind11.
- C++17 native extension `_equilibrium`; CMake >=3.20.
- Installed Provider dependency: `epcsaft==0.2.0.dev0`, non-editable; CMake validates its native SDK header.
- Optimizers: Ipopt via pkg-config; pinned HiGHS 1.15.1 and NLopt 2.11.0 fetched by CMake.
- Python runtime dependency: Pint >=0.25.
- Test/lint/type tooling: pytest, Ruff (100 columns; E/F/I/UP/B/RUF), strict mypy for `epcsaft_equilibrium`.
- Scientific-computing repository profile; exact derivatives and independent physical certification are first-class contracts.
