# difference-in-differences: a methods explainer

simulation code for the DiD explainer at [pesie.github.io](https://pesie.github.io/projects/did-explainer.html).

## files
- `simulate_data.R` — generates the synthetic panel dataset
- `did_explainer.R` — estimates DiD models and produces figures

## to use
Run `did_explainer.R`. It sources `simulate_data.R` automatically.

## dependencies
`dplyr`, `tidyr`, `ggplot2`, `fixest`