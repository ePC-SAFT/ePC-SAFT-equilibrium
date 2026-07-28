# Conventions

- Keep Provider ownership boundaries strict: Equilibrium consumes values, exact derivatives, pressure/domain/packing data through the Provider SDK only.
- Treat optimizer termination, numerical convergence, pressure closure, physical KKT/certification, finite-search completion, and globality as separate evidence.
- Finite HELD/HELD2 searches always report `globality_certificate=not_guaranteed`.
- Fail closed on unresolved Provider/domain/root/certification failures; never replace invalid evaluations with finite penalty objectives.
- Use reduced/scaled thermodynamic residuals and exact coordinate-chain Hessians; do not loosen physical gates to accept a case.
- Tests protect scientific claims and source algorithm steps, not implementation churn. Keep publication-scale/live runs opt-in.
- Preserve unrelated dirty-worktree changes. Avoid legacy compatibility layers and duplicate solver/EOS paths unless explicitly authoritative.