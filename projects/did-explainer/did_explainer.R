# did_explainer.R --------------------------------------------------------
# simulates data, estimates DiD models, and produces figures
# for the DiD methods explainer
# pesie.github.io/projects/did-explainer
# author: precious esie

library(dplyr)
library(tidyr)
library(ggplot2)
library(fixest)

source("projects/did-explainer/simulate_data.R")

# 1. theme and palette ---------------------------------------------------

col_treat <- "dodgerblue"
col_control <- "gray40"

theme_site <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "#EEEEEE"),
      axis.title = element_text(size = 11),
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "#555555"),
      plot.caption = element_text(size = 9, color = "#AAAAAA"),
      legend.position = "bottom",
      legend.title = element_blank()
    )
}

# 2. group means ---------------------------------------------------------
# summarise to group-time level for plotting
# individual state variation is averaged out — shows the group trend clearly

core_means <- core %>%
  group_by(group, time) %>%
  summarise(mean_outcome = mean(outcome), .groups = "drop")

# 3. figure 1: group means over time ------------------------------------
# the key visual for DiD intuition:
# groups start at different levels (level difference is fine — DiD accounts
# for it) but trend in parallel before the intervention
# the post-intervention divergence is the treatment effect

fig1 <- ggplot(core_means,
               aes(x = time, y = mean_outcome,
                   color = group, group = group)) +
  geom_vline(xintercept = treat_at - 0.5,
             linetype = "dotted", color = "gray20", linewidth = 0.4) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  annotate("text", x = treat_at - 0.4,
           y = 25,
           label = "Intervention",
           size = 3, color = "gray20", hjust = 0) +
  scale_color_manual(values = c("Treatment" = col_treat,
                                "Control" = col_control)) +
  scale_x_continuous(breaks = 1:n_time) +
  scale_y_continuous(limits = c(0, 50)) +
  labs(
    title = "Outcome trends before and after the intervention",
    subtitle = "Groups start at different levels but trend in parallel before the intervention",
    x = "Time", y = "Outcome",
    caption = "Simulated data"
  ) +
  theme_site()

# 4. event study ---------------------------------------------------------
# two-way fixed effects event study using fixest::feols()
# unit fixed effects (state_id) absorb stable between-state differences
# time fixed effects (time) absorb common shocks affecting all states
# i() creates interactions between time_to_treat and treated
# ref = -1 sets the last pre-intervention period as the reference (estimate = 0)
# standard errors clustered at the state level

core_es <- core %>%
  mutate(
    # time relative to intervention: negative = pre, positive = post
    # untreated states set to 0 (absorbed by fixed effects)
    time_to_treat = if_else(treated, time - treat_at, 0),
    state_id = as.integer(state_id)
  )

es_model <- feols(
  outcome ~ i(time_to_treat, treated, ref = -1) | state_id + time,
  data = core_es,
  cluster = ~state_id
)

# extract estimates and confidence intervals
# term names from fixest follow the pattern "time_to_treat::N:treatedTRUE"
# gsub extracts the period number from the term name
es_df <- as.data.frame(coeftable(es_model)) %>%
  tibble::rownames_to_column("term") %>%
  filter(grepl("time_to_treat", term)) %>%
  mutate(
    time = as.integer(gsub("time_to_treat::(.+):treated", "\\1", term)),
    conf.low = Estimate - 1.96 * `Std. Error`,
    conf.high = Estimate + 1.96 * `Std. Error`,
    period = factor(
      if_else(time < 0, "Pre", "Post"),
      levels = c("Pre", "Post")
    )
  ) %>%
  rename(estimate = Estimate)

# add reference period manually,  estimate is zero by construction
ref_row <- data.frame(
  term = "reference",
  estimate = 0,
  conf.low = NA,
  conf.high = NA,
  time = -1,
  period = "Pre"
)

es_df <- bind_rows(es_df, ref_row) %>%
  arrange(time)

# 5. figure 2: event study plot ------------------------------------------
# pre-period estimates should scatter around zero if parallel trends hold
# a visible pre-trend would suggest the groups were already diverging
# before the intervention — a violation of the parallel trends assumption
# post-period estimates show when and how the effect emerged

fig2 <- ggplot(es_df, aes(x = time, y = estimate, color = period)) +
  geom_hline(yintercept = 0, color = "#AAAAAA", linewidth = 0.5) +
  geom_vline(xintercept = -0.5, linetype = "dotted", color = "#AAAAAA") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                width = 0.1, linewidth = 0.5) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  annotate("point", x = -1, y = 0, size = 2, color = "gray60") +
  scale_color_manual(
    breaks = c("Pre", "Post"),
    values = c("Pre" = "gray60", "Post" = col_treat)
  ) +
  scale_x_continuous(breaks = min(es_df$time):max(es_df$time)) +
  labs(
    title = "Event study plot",
    subtitle = "Estimated effect at each period relative to intervention, with 95% confidence intervals",
    x = "Time relative to intervention",
    y = "Estimated effect",
    caption = "Simulated data. Two-way fixed effects estimator."
  ) +
  theme_site()

# 6. DiD summary estimate ------------------------------------------------
# two-way fixed effects regression
# treated:post interaction is the ATT
# treatedTRUE and postTRUE are dropped due to collinearity with fixed effects
# — this is expected behavior, not a misspecification

did_model <- feols(
  outcome ~ treated * post | state_id + time,
  data = core
)

att <- coef(did_model)["treatedTRUE:postTRUE"]
