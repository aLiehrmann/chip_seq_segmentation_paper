library(data.table)
library(namedCapture)
library(furrr)
library(purrr)
library(penaltyLearning)

plan(multiprocess, workers=9)
path_to_datasets <- "data/chunks"
path_to_peaks_errors <- "data"
path_to_output <- "data"
pattern <- paste0(path_to_datasets,"/","(?<set_name>.+?)", "/", "(?<chunk_id>[0-9]+)")


f <- function(algo_, peaks_errors, output) {
  models_path <- Sys.glob(file.path(path_to_datasets, "H*", "*",paste0(algo_,".RData")))
  matched <- str_match_named(models_path, pattern)
  params <- data.frame(matched, models_path = models_path)
  models <- do.call(rbind, pmap(params, function(models_path, set_name, chunk_id) {
    load(models_path) 
    dt.models[,chunk.name := paste0(set_name, "/", chunk_id)]
    dt.models
  }))
  models[, algo:=algo_]
  names(models)[names(models) == "sample"] <- "sample.id"
  if(!any(names(models) == "phi")) {
    models[,phi := NA]
  }
  models <- models[, list(chunk.name, sample.id, peaks, changepoints, loss, phi)]
  temp_ls <- ls()
  load(paste0(path_to_peaks_errors, '/', peaks_errors, ".RData"))
  peaks_errors_ <- setdiff(ls(),c('temp_ls',temp_ls))
  peaks_errors_ <- eval(parse(text=peaks_errors_))
  if(!any(names(peaks_errors_) == "phi")) {
    peaks_errors_[,phi := NA]
  }
  all.totals <- peaks_errors_[, list(
    total.fp=sum(fp),
    total.fn=sum(fn),
    total.errors=sum(fp+fn),
    possible.fp=sum(possible.fp),
    possible.fn=sum(possible.tp),
    labels=.N
  ), by=list(chunk.name, sample.id, changepoints, phi)]
  models[,k:=paste0(chunk.name, sample.id, changepoints, phi)]
  all.totals[, k:=paste0(chunk.name, sample.id, changepoints, phi)]
  setkey(models, k)
  setkey(all.totals, k)
  all.loss <- models[all.totals,]
  all.modelSelection <- all.loss[, {
    modelSelection(.SD,
    complexity="changepoints")
  }, by=list(chunk.name, sample.id, phi)]
  all.modelSelection <- all.modelSelection[, list(chunk.name, sample.id, phi, min.log.lambda, max.log.lambda, peaks, changepoints, loss,total.errors, labels)]
  save(all.modelSelection, file = paste0(path_to_output, "/", output, ".RData"))
}


params <- tibble::tribble(
  ~ algo_,            ~ peaks_errors,                                   ~ output,
  "PDPA_negbin_2",   "PDPA_negbin_2_peaks_errors",                      "all_modelSelection_PDPA_negbin_2",
  "PDPA_poisson",    "PDPA_poisson_peaks_errors",                       "all_modelSelection_PDPA_poisson",
  "PDPA",            "PDPA_peaks_errors",                               "all_modelSelection_PDPA",
  "PDPA_negbin_2",   "PDPA_negbin_peaks_errors_2_largest_peak_rule",    "all_modelSelection_PDPA_negbin_2_largest_peak_rule",
  "PDPA_negbin_2",   "PDPA_negbin_peaks_errors_2_smallest_peak_rule",   "all_modelSelection_PDPA_negbin_2_smallest_peak_rule",
  "PDPA_poisson",    "PDPA_poisson_peaks_errors_largest_peak_rule",     "all_modelSelection_PDPA_poisson_largest_peak_rule",
  "PDPA_poisson",    "PDPA_poisson_peaks_errors_smallest_peak_rule",    "all_modelSelection_PDPA_poisson_smallest_peak_rule",
  "PDPA",            "PDPA_peaks_errors_largest_peak_rule",             "all_modelSelection_PDPA_largest_peak_rule",
  "PDPA",            "PDPA_peaks_errors_smallest_peak_rule",            "all_modelSelection_PDPA_smallest_peak_rule",
  "updown_negbin_2", "updown_negbin_2_peaks_errors",                    "all_modelSelection_updown_negbin_2",
  "updown_poisson",  "updown_poisson_peaks_errors",                     "all_modelSelection_updown_poisson",
  "updown",          "updown_peaks_errors",                             "all_modelSelection_updown"
)

future_pwalk(params, f)