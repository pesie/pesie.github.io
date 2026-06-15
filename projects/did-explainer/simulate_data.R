# simulate_data.R --------------------------------------------------------
# synthetic data for the DiD methods explainer
# sourced by did_explainer.R 
# pesie.github.io/projects/did-explainer
# author: precious esie

library(dplyr)
library(tidyr)

set.seed(54321)

# 1. parameters ----------------------------------------------------------

n_states <- 50
n_treated <- 25
n_time <- 10
treat_at <- 6        # policy takes effect at period 6
true_att <- 8        # true treatment effect
baseline_treat <- 22 # average starting outcome, treated states
baseline_control <- 35 # average starting outcome, control states
common_trend <- 1.5  # outcome change per period, both groups
noise_sd <- 4        # within-state year-to-year variation
state_sd <- 5        # between-state variation in starting levels

# true_att / noise_sd = 2: effect is detectable but not unrealistically clean

# 2. core dataset --------------------------------------------------------
# clean design: parallel trends hold by construction
# both groups share the same time trend

states <- tibble(
  state_id = 1:n_states,
  treated = state_id <= n_treated,
  state_fe = if_else(
    treated,
    rnorm(n_states, mean = baseline_treat,   sd = state_sd),
    rnorm(n_states, mean = baseline_control, sd = state_sd)
  )
)

core <- states %>%
  expand_grid(time = 1:n_time) %>%
  mutate(
    post = time >= treat_at,
    time_trend = common_trend * time,
    effect = if_else(treated & post, true_att, 0),
    outcome = state_fe + time_trend + effect +
      rnorm(n(), mean = 0, sd = noise_sd),
    group = if_else(treated, "Treatment", "Control")
  )
