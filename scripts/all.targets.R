library(penaltyLearning)
library(furrr)

plan(multiprocess, workers=10)

path_to_input <- "data/"
path_to_output <- "data/"

f <- function(input, output) {
  load(paste0(path_to_input, "/", input,".RData"))
  names(all.modelSelection)[names(all.modelSelection) == "total.errors"] <- "errors"
  all.targets <- targetIntervals(all.modelSelection, c("chunk.name", "sample.id","phi"))
  save(all.targets, file = paste0(path_to_output,"/",output, ".RData"))
}

params <- tibble::tribble(
  ~ input,                                                  ~ output,
  "all_modelSelection_PDPA_negbin_2",                       "all_targets_PDPA_negbin_2",
  "all_modelSelection_PDPA_poisson",                        "all_targets_PDPA_poisson",
  "all_modelSelection_PDPA",                                "all_targets_PDPA",
  "all_modelSelection_PDPA_negbin_2_largest_peak_rule",     "all_targets_PDPA_negbin_2_largest_peak_rule",
  "all_modelSelection_PDPA_negbin_2_smallest_peak_rule",    "all_targets_PDPA_negbin_2_smallest_peak_rule",
  "all_modelSelection_PDPA_poisson_largest_peak_rule",      "all_targets_PDPA_poisson_largest_peak_rule",
  "all_modelSelection_PDPA_poisson_smallest_peak_rule",     "all_targets_PDPA_poisson_smallest_peak_rule",
  "all_modelSelection_PDPA_largest_peak_rule",              "all_targets_PDPA_largest_peak_rule",
  "all_modelSelection_PDPA_smallest_peak_rule",             "all_targets_PDPA_smallest_peak_rule",
  "all_modelSelection_updown_negbin_2",                     "all_targets_updown_negbin_2",
  "all_modelSelection_updown_poisson",                      "all_targets_updown_poisson",
  "all_modelSelection_updown",                              "all_targets_updown"
)

future_pwalk(params, f)