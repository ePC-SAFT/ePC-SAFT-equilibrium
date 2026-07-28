# Task completion

1. Run the narrowest affected C++/Python scientific tests, then the compact `pytest -q` suite when risk warrants.
2. Rebuild the extension with `cmake --build build -j2` for native changes.
3. Run `ruff check .` and strict mypy for affected Python/package surfaces when changed.
4. For HELD2 changes, retain manufactured root/basin/KKT/Step-6/Stage-III evidence and a truthful opt-in live trace; optimizer success alone is never acceptance.
5. Review `git diff`, preserve unrelated changes, and verify no task-owned scratch artifacts or processes remain.
6. Run `bash "$HOME/.codex/hooks/codex-cleanup.sh" --repo-root .` after changes or launched processes.
7. Governed work must satisfy the active GitHub issue/Project Truss receipt and verification contract before commit/handoff.