# Difference-in-Differences: A Methods Explainer

Simulation code for the DiD explainer at [pesie.github.io](https://pesie.github.io/projects/did-explainer.html).

## Files
- `simulate_data.R` — generates the synthetic panel dataset
- `did_explainer.R` — estimates DiD models and produces figures

## Usage
Run `did_explainer.R`. It sources `simulate_data.R` automatically.

## Dependencies
`dplyr`, `tidyr`, `ggplot2`, `fixest`