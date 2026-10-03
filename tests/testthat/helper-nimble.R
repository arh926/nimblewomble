# Attach nimble at the top of every test session. nimble's model builder resolves
# internal helpers (e.g., calc_dmnormAltParams) via the search path; without an explicit
# library(nimble) attach, tests that call nimbleModel via nimblewomble get 'R function
# ... does not exist' errors. See test-gp_fit-covariates.R for a fuller explanation of
# the symptom.
suppressPackageStartupMessages(library(nimble))
