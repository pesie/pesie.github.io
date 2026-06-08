# HINTS 7 Latent Class Analysis
# Health Information Navigation, Trust and Quality of Care
#
# Research question: Do distinct profiles of health information experience
# predict trust in health institutions and perceived quality of care,
# independent of demographics?
#
# Data: Health Information National Trends Survey (HINTS 7), 2024
# NCI nationally representative survey of U.S. adults
# Survey design: jackknife replicate weights (PERSON_FINWT0-50)

# Set up ----------------------------------------------------------------------
library(tidyverse)
library(poLCA)   # Latent class analysis
library(survey)  # Complex survey design and weighted estimation
library(psych)   # Exploratory factor analysis (tetrachoric correlations)

# Load data
load("hintsLCA/hints7_public.rda")

# User settings ---------------------------------------------------------------

best_k      <- 4       # Selected number of classes (chosen after model comparison)
class_range <- 2:6     # Range of class solutions to evaluate
set.seed(1234)

# Four segmentation variables capturing two dimensions of health information experience:
#   - Perceived inconsistency: do recommendations contradict or change?
#   - Perceived overload: are there too many recommendations?
seg_vars <- c(
  "HealthRecsConflict",
  "HealthRecsChange",
  "TooManyRecommendations",
  "EverythingCauseCancer"
)

# Full variable list to retain from the raw HINTS file
analysis_vars <- c(
  seg_vars,
  "Age",
  "PoliticalViewpoint",
  "Education",
  "IncomeRanges",
  "RaceEthn",
  "GeneralHealth",
  "SocMed_Visited",
  "PERSON_FINWT0",
  "TrustHCSystem",
  "CancerTrustDoctor",
  "CancerTrustFamily",
  "CancerTrustGov",
  "CancerTrustScientists",
  "QualityCare",
  "SeekCancerInfo",
  "MisleadingHealthInfo",
  "ConfidentMedForms",
  "VAR_STRATUM",
  "VAR_CLUSTER",
  "HHID"
)

# HINTS uses text codes for missing/unusable responses; convert all to NA
missing_codes <- c(
  "Missing data (Not Ascertained)",
  "Missing data (Web partial - Question Never Seen)",
  "Unreadable or Non-conforming numeric response",
  "Multiple responses selected in error"
)

# Demographic binary indicators used in the profile table
# Named vector: display label → variable name in dat_out
binary_var_map <- c(
  "Age 18-34"                  = "age_18_34",
  "Age 35-49"                  = "age_35_49",
  "Age 50-64"                  = "age_50_64",
  "Age 65+"                    = "age_65plus",
  "Low income (under $35,000)" = "low_income",
  "College educated"           = "college_edu",
  "Non-Hispanic White"         = "white",
  "Non-Hispanic Black or African American" = "black",
  "Hispanic"                   = "hispanic",
  "Fair/Poor health"           = "fair_poor_health",
  "Frequent social media use"  = "frequent_socmed",
  "Liberal"                    = "liberal",
  "Moderate"                   = "moderate",
  "Conservative"               = "conservative"
)

# Outcome measures: display label → logical expression evaluated within svyby
outcome_map <- c(
  "Trust HC system (A lot)"      = "TrustHCSystem == 'A lot'",
  "Trust doctors (A lot)"        = "CancerTrustDoctor == 'A lot'",
  "Trust scientists (A lot)"     = "CancerTrustScientists == 'A lot'",
  "Trust government (A lot)"     = "CancerTrustGov == 'A lot'",
  "Healthcare quality excellent" = "QualityCare == 'Excellent'"
)

# Segmentation variables for the class profile table
# Uses explicit 0/1 _hi indicators (defined in dat_out) to avoid TRUE/FALSE
# column name ambiguity in svyby output
seg_var_map <- c(
  "Experts give conflicting advice" = "conflict_hi",
  "Recommendations keep changing"   = "change_hi",
  "Too many recommendations"        = "toomany_hi",
  "Everything causes cancer"        = "everything_hi"
)


# Helper functions  -----------------------------------------------------------

# Recode ordinal HINTS responses to binary integers for poLCA
# poLCA requires: 1 = category A, 2 = category B (no 0s)
recode_binary <- function(x, positive_values, negative_values) {
  case_when(
    x %in% positive_values ~ 1L,
    x %in% negative_values ~ 2L,
    TRUE                   ~ NA_integer_
  )
}

# Fit LCA models across a range of class solutions and return fit statistics
fit_lca_models <- function(formula, data, class_range, nrep = 10, maxiter = 3000) {
  models <- map(
    class_range,
    ~ poLCA(
      formula,
      data     = data,
      nclass   = .x,
      nrep     = nrep,       # Multiple random starts to avoid local maxima
      maxiter  = maxiter,
      verbose  = FALSE
    )
  )
  
  names(models) <- as.character(class_range)
  
  fit_stats <- tibble(
    classes        = class_range,
    log_likelihood = map_dbl(models, "llik"),
    AIC            = map_dbl(models, "aic"),
    BIC            = map_dbl(models, "bic")
  ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2)))
  
  list(models = models, fit_stats = fit_stats)
}

# Append class assignments and posterior probabilities to the analytic dataset
add_lca_classes <- function(data, model) {
  posterior <- as_tibble(model$posterior)
  names(posterior) <- paste0("class_prob_", seq_len(ncol(posterior)))
  bind_cols(data %>% mutate(class = model$predclass), posterior)
}

# Compute survey-weighted column percentages for binary indicators by group
# Uses pre-coded 0/1 variables; avoids logical expression parsing issues
weighted_indicator_by_group <- function(design, indicators, group_var, label_name) {
  map_dfr(names(indicators), function(label) {
    result <- svyby(
      as.formula(paste0("~", indicators[[label]])),
      as.formula(paste0("~", group_var)),
      design,
      svymean,
      na.rm    = TRUE,
      vartype  = NULL
    )
    
    estimate_col <- setdiff(names(result), group_var)[1]
    
    result %>%
      as_tibble() %>%
      transmute(
        !!label_name := label,
        !!group_var  := .data[[group_var]],
        pct = round(.data[[estimate_col]] * 100, 1)
      )
  })
}

# Compute survey-weighted percentages for outcomes defined as logical expressions
# Extracts the TRUE column by pattern match rather than constructed name
weighted_expression_by_group <- function(design, expressions, group_var, label_name) {
  map_dfr(names(expressions), function(label) {
    result <- svyby(
      as.formula(paste0("~I(", expressions[[label]], ")")),
      as.formula(paste0("~", group_var)),
      design,
      svymean,
      na.rm    = TRUE,
      vartype  = NULL
    )
    
    true_col <- grep("TRUE$", names(result), value = TRUE)[1]
    
    result %>%
      as_tibble() %>%
      transmute(
        !!label_name := label,
        !!group_var  := .data[[group_var]],
        pct = round(.data[[true_col]] * 100, 1)
      )
  })
}


# Load and prepare data -------------------------------------------------------

# Identify replicate weight columns (PERSON_FINWT1 through PERSON_FINWT50)
rep_weight_vars_public <- names(public) %>%
  str_subset("^PERSON_FINWT[1-9][0-9]*$")

keep_vars <- union(analysis_vars, rep_weight_vars_public)

dat_lca <- public %>%
  dplyr::select(all_of(keep_vars)) %>%
  mutate(
    # Standardize to character before recoding to avoid factor level issues
    across(where(is.factor), as.character),
    across(where(is.character), ~ na_if(.x, "")),
    # Convert HINTS missing codes to NA across all character variables
    across(where(is.character), ~ if_else(.x %in% missing_codes, NA_character_, .x)),
    # Recode segmentation variables to binary integers (1 = high, 2 = low)
    # poLCA does not accept 0/1 or logical inputs
    HealthRecsConflict = recode_binary(
      HealthRecsConflict,
      positive_values = c("Often", "Very Often"),
      negative_values = c("Never", "Rarely")
    ),
    HealthRecsChange = recode_binary(
      HealthRecsChange,
      positive_values = c("Often", "Very Often"),
      negative_values = c("Never", "Rarely")
    ),
    TooManyRecommendations = recode_binary(
      TooManyRecommendations,
      positive_values = c("Strongly agree", "Somewhat agree"),
      negative_values = c("Somewhat disagree", "Strongly disagree")
    ),
    EverythingCauseCancer = recode_binary(
      EverythingCauseCancer,
      positive_values = c("Strongly agree", "Somewhat agree"),
      negative_values = c("Somewhat disagree", "Strongly disagree")
    )
  ) %>%
  # Complete cases on segmentation variables only; outcomes allowed to be missing
  drop_na(all_of(seg_vars))


# Exploratory factor analysis -------------------------------------------------

# Test whether inconsistency and overload items reflect separate dimensions
# before imposing that structure via LCA
# Tetrachoric correlations are appropriate for binary items
efa_items <- dat_lca %>%
  dplyr::select(all_of(seg_vars)) %>%
  mutate(across(everything(), ~ as.integer(.x == 1)))

efa_correlation <- psych::tetrachoric(efa_items)

efa_solution <- psych::fa(
  efa_correlation$rho,
  nfactors = 2,
  rotate   = "oblimin",   # Oblique rotation; factors are likely correlated
  fm       = "minres"
)
# Result: two-factor solution; inconsistency and overload items load separately


# Fit and select latent class models ------------------------------------------

lca_formula <- as.formula(
  paste0("cbind(", paste(seg_vars, collapse = ", "), ") ~ 1")
)

lca_results <- fit_lca_models(
  formula     = lca_formula,
  data        = dat_lca,
  class_range = class_range
)

fit_stats  <- lca_results$fit_stats
lca_models <- lca_results$models
best_model <- lca_models[[as.character(best_k)]]

# BIC plot for model selection appendix
bic_plot <- ggplot(fit_stats, aes(x = classes, y = BIC)) +
  geom_line() +
  geom_point(size = 3) +
  scale_x_continuous(breaks = class_range) +
  labs(title = "Latent class model selection", x = "Number of classes", y = "BIC") +
  theme_minimal()


# Add class assignments and derived variables ---------------------------------

dat_out <- dat_lca %>%
  add_lca_classes(best_model) %>%
  mutate(
    age_numeric = suppressWarnings(as.numeric(Age)),
    
    # Age group indicators
    age_18_34  = as.integer(age_numeric >= 18 & age_numeric <= 34),
    age_35_49  = as.integer(age_numeric >= 35 & age_numeric <= 49),
    age_50_64  = as.integer(age_numeric >= 50 & age_numeric <= 64),
    age_65plus = as.integer(age_numeric >= 65),
    
    # Income: low = under $35,000
    low_income = as.integer(IncomeRanges %in% c(
      "$0 to $9,999",
      "$10,000 to $14,999",
      "$15,000 to $19,999",
      "$20,000 to $34,999"
    )),
    low_income = if_else(is.na(IncomeRanges), NA_integer_, low_income),
    
    college_edu = as.integer(Education %in% c("College graduate", "Postgraduate")),
    
    # Race/ethnicity indicators
    white    = as.integer(RaceEthn == "Non-Hispanic White"),
    black    = as.integer(RaceEthn == "Non-Hispanic Black or African American"),
    hispanic = as.integer(RaceEthn == "Hispanic"),
    
    fair_poor_health = as.integer(GeneralHealth %in% c("Fair", "Poor")),
    frequent_socmed  = as.integer(SocMed_Visited == "Almost everyday"),
    
    # Political viewpoint collapsed to three categories
    liberal      = as.integer(PoliticalViewpoint %in% c("Very Liberal", "Liberal", "Somewhat Liberal")),
    moderate     = as.integer(PoliticalViewpoint == "Moderate"),
    conservative = as.integer(PoliticalViewpoint %in% c("Somewhat Conservative", "Conservative", "Very Conservative")),
    
    # Logical flags used in opening statistics
    recs_change_often   = HealthRecsChange == 1,
    too_many_recs_agree = TooManyRecommendations == 1,
    
    # Explicit 0/1 segmentation indicators for weighted profile tables
    # Avoids TRUE/FALSE column name ambiguity in svyby output
    conflict_hi   = as.integer(HealthRecsConflict == 1),
    change_hi     = as.integer(HealthRecsChange == 1),
    toomany_hi    = as.integer(TooManyRecommendations == 1),
    everything_hi = as.integer(EverythingCauseCancer == 1),
    
    # Propagate NA from source variables to derived indicators
    across(
      c(age_18_34, age_35_49, age_50_64, age_65plus),
      ~ if_else(is.na(age_numeric), NA_integer_, .x)
    ),
    across(
      c(liberal, moderate, conservative),
      ~ if_else(is.na(PoliticalViewpoint), NA_integer_, .x)
    ),
    across(
      c(white, black, hispanic),
      ~ if_else(is.na(RaceEthn), NA_integer_, .x)
    ),
    college_edu      = if_else(is.na(Education),     NA_integer_, college_edu),
    fair_poor_health = if_else(is.na(GeneralHealth),  NA_integer_, fair_poor_health),
    frequent_socmed  = if_else(is.na(SocMed_Visited), NA_integer_, frequent_socmed)
  )


# Survey design ---------------------------------------------------------------

# HINTS uses jackknife replicate weights with 50 replicates
# combined.weights = TRUE because PERSON_FINWT0 is already incorporated into
# the replicate weights (do not divide; see HINTS methodology report)
rep_weight_vars <- names(dat_out) %>%
  str_subset("^PERSON_FINWT[1-9][0-9]*$")

hints_design <- svrepdesign(
  weights          = ~PERSON_FINWT0,
  repweights       = dat_out[, rep_weight_vars],
  type             = "JK1",
  combined.weights = TRUE,
  data             = dat_out
)


# Weighted outputs ------------------------------------------------------------

# Class size distribution (weighted)
weighted_class <- svytable(~class, hints_design) %>%
  prop.table() %>%
  as_tibble(name = "class") %>%
  transmute(
    class = as.integer(class),
    pct   = round(n * 100, 1)
  )

# Segmentation profile: proportion endorsing each item within each class
segmentation_profile <- weighted_indicator_by_group(
  design     = hints_design,
  indicators = seg_var_map,
  group_var  = "class",
  label_name = "variable"
)

segmentation_profile_wide <- segmentation_profile %>%
  pivot_wider(
    names_from   = class,
    values_from  = pct,
    names_prefix = "Class_"
  ) %>%
  arrange(variable)

# Demographic profile: describes classes without defining them
# Demographics are excluded from the LCA to avoid conflating who people are
# with how they experience health information
demographic_profile <- weighted_indicator_by_group(
  design     = hints_design,
  indicators = binary_var_map,
  group_var  = "class",
  label_name = "characteristic"
)

demographic_profile_wide <- demographic_profile %>%
  pivot_wider(
    names_from   = class,
    values_from  = pct,
    names_prefix = "Class_"
  ) %>%
  arrange(characteristic)

# Outcome profile: trust and perceived care quality by class
outcome_profile <- weighted_expression_by_group(
  design      = hints_design,
  expressions = outcome_map,
  group_var   = "class",
  label_name  = "outcome"
)

outcome_profile_wide <- outcome_profile %>%
  pivot_wider(
    names_from   = class,
    values_from  = pct,
    names_prefix = "Class_"
  ) %>%
  arrange(outcome)

# Opening statistics cited in the post introduction
opening_stats <- tibble(
  statistic = c(
    "Health recommendations change often",
    "Too many health recommendations"
  ),
  pct = c(
    coef(svymean(~recs_change_often,   hints_design, na.rm = TRUE))["recs_change_oftenTRUE"],
    coef(svymean(~too_many_recs_agree, hints_design, na.rm = TRUE))["too_many_recs_agreeTRUE"]
  ) * 100
) %>%
  mutate(pct = round(pct, 1))

# Demographics by individual segmentation variable (for supplemental exploration)
# Supplemental: demographics broken out by individual segmentation variable
# rather than by class. Answers a different question than demographic_profile:
# not "who is in each class" but "who tends to endorse each item independently."
# Useful for checking whether item-level variation tracks class-level patterns.

demographics_by_segmentation <- map_dfr(seg_vars, function(seg_var) {
  weighted_indicator_by_group(
    design     = hints_design,
    indicators = binary_var_map,
    group_var  = seg_var,
    label_name = "demographic"
  ) %>%
    dplyr::filter(.data[[seg_var]] == 1) %>%
    mutate(segmentation_variable = seg_var, .before = 1) %>%
    dplyr::select(segmentation_variable, demographic, pct)
})

demographics_by_segmentation_wide <- demographics_by_segmentation %>%
  pivot_wider(
    names_from  = segmentation_variable,
    values_from = pct
  )

# note: no meaningful variation by demographics
