# Results

| Directory | Contents |
|---|---|
| `mouse_pbpk/` | Current figures and CSV summaries regenerated from `model/ode1.0.py` |
| `spatial_pk_1d/` | One-dimensional flow/trapping demonstration |
| `legacy_mouse_pbpk/` | Figures from an earlier mouse-model revision, retained for comparison |

Generated frontend trajectories live in
`sineup-delivery-atlas/public/data/` because the React application consumes
them directly. Current results should never be replaced by files from the
legacy directory without re-running the model and its mass-balance checks.
