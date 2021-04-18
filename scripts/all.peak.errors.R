library(purrr)
library(furrr)
library(data.table)
library(PeakError)
library(stringr)

post_processing_rules <- list(
  maxjump = function(cp, theta, genomic_position) {
    trend <- sign(diff(theta))
    cp_tmp <- cp[trend!=0]
    theta_tmp <- theta[c(trend!=0,TRUE)]
    trend <- trend[trend!=0]
    if (length(trend)>1){
      trend_rle <- rle(trend)
      targets <- which(trend_rle$values[1:(length(trend_rle$values)-1)] == 1 &
      trend_rle$values[2:length(trend_rle$values)] == -1)
      if (length(targets)){
        offset_start <- map_dbl(targets, function(x){ifelse(x-1, sum(trend_rle$length[1:(x-1)]), 0)})
        offset_end <- map_dbl(targets+1, function(x){ifelse(x-1, sum(trend_rle$length[1:(x-1)]), 0)})
        data.table(
          chromStart =  genomic_position[cp_tmp[map2_dbl(offset_start, targets, function(x,y){
            successive_up_changes <- theta_tmp[(x+1):(x+1+trend_rle$length[y])]
            x+which.max(diff(successive_up_changes))
          })]],
          chromEnd =  genomic_position[cp_tmp[map2_dbl(offset_end, targets+1, function(x,y){
            successive_down_changes <- theta_tmp[(x+1):(x+1+trend_rle$length[y])]
            x+which.min(diff(successive_down_changes))
          })]]
        ) 
      } else {
        data.table(chromStart = 0, chromEnd = 0.001)
      }
    } else {
      data.table(chromStart = 0, chromEnd = 0.001)
    }
  }, 
  thinnest_peak = function(cp, theta, genomic_position) {
    trend <- sign(diff(theta))
    cp_tmp <- cp[trend!=0]
    trend <- trend[trend!=0]
    if (length(trend)>1){
      trend_rle <- rle(trend)
      targets <- which(trend_rle$values[1:(length(trend_rle$values)-1)] == 1 &
      trend_rle$values[2:length(trend_rle$values)] == -1)
      if (length(targets)){
        data.table(
          chromStart =  genomic_position[cp_tmp[map_dbl(targets, ~sum(trend_rle$length[1:.x]))]],
          chromEnd =  genomic_position[cp_tmp[map_dbl(targets, ~sum(trend_rle$length[1:(.x+1)])-trend_rle$length[.x+1]+1)]]
        ) 
      } else {
        data.table(chromStart = 0, chromEnd = 0.001)
      }
    } else {
      data.table(chromStart = 0, chromEnd = 0.001)
    }
  },
  largest_peak = function(cp, theta, genomic_position) {
    trend <- sign(diff(theta))
    cp_tmp <- cp[trend!=0]
    trend <- trend[trend!=0]
    if (length(trend)>1){
      trend_rle <- rle(trend)
      targets <- which(trend_rle$values[1:(length(trend_rle$values)-1)] == 1 &
      trend_rle$values[2:length(trend_rle$values)] == -1)
      if (length(targets)){
        data.table(
          chromStart = genomic_position[cp_tmp[map_dbl(targets, ~sum(trend_rle$length[1:.x])-trend_rle$length[.x]+1)]],
          chromEnd =  genomic_position[cp_tmp[map_dbl(targets, ~sum(trend_rle$length[1:(.x+1)]))]]
        ) 
      } else {
        data.table(chromStart = 0, chromEnd = 0.001)
      }
    } else {
      data.table(chromStart = 0, chromEnd = 0.001)
    }
  }
)

peakCalling <- function(post_processing_rule, input_directories, input_model, input_region,
  input_count, output_directory, output) {
  models.peaks.errors <- rbindlist(future_map(input_directories, function(directory_p) {
    load(file.path(directory_p, input_region))
    setDT(regions)
    chunk <- str_match(pattern = "H.*",directory_p)[1]
    load(file.path(directory_p, input_model))
    load(file.path(directory_p, input_count))
    setDT(counts)
    if (is.null(dt.models$phi)){
      dt.models$phi <- 0
    }
    rbindlist(map(unique(dt.models$phi), function(phi_p) {
      rbindlist(map(unique(dt.models$sample), function(sample_p) {
        dt_regions_c <- regions[sample.id==sample_p,]
        dt_counts_c <- counts[sample.id==sample_p,]
        dt_models_c <- dt.models[sample==sample_p & phi==phi_p,]
        rbindlist(map(dt_models_c$changepoints, function(changepoints_p){
          dt_models_cbis <- dt_models_c[changepoints==changepoints_p]
          dt_peaks <- with(dt_models_cbis$fit[[1]],
            post_processing_rules[[post_processing_rule]](
              cp = changepoints[-length(changepoints)],
              theta = parameters,
              genomic_position = dt_counts_c$chromEnd
            )
          )
          dt_peak_errors <- PeakErrorChrom(dt_peaks, dt_regions_c)
          setDT(dt_peak_errors)
          dt_peak_errors[,c("chunk.name","sample.id","rule","peaks","changepoints","phi") := list(chunk, sample_p, post_processing_rule, dt_models_cbis$peaks, changepoints_p, phi_p)]
        }))
      }))
    }))
  }))
  save(models.peaks.errors, file = file.path(output_directory, output))
}


d <- Sys.glob("data/chunks/H*/*")

m <- rbindlist(list(
    list(
    post_processing_rule = "maxjump", 
    input_directories = list(d), 
    input_model = "PDPA.RData",
    input_region = "regions.RData",
    input_count = "counts.RData",
    output_directory = "data", 
    output = "PDPA_peaks_errors.RData"
  ),
  list(
    post_processing_rule = "thinnest_peak", 
    input_directories = list(d), 
    input_model = "PDPA.RData",
    input_region = "regions.RData",
    input_count = "counts.RData",
    output_directory = "data", 
    output = "PDPA_peaks_errors_thinnest_peak_rule.RData"
  ),
  list(
    post_processing_rule = "largest_peak", 
    input_directories = list(d), 
    input_model = "PDPA.RData",
    input_region = "regions.RData",
    input_count = "counts.RData",
    output_directory = "data", 
    output = "PDPA_peaks_errors_largest_peak_rule.RData"
  ),
  list(
    post_processing_rule = "maxjump", 
    input_directories = list(d), 
    input_model = "PDPA_poisson.RData",
    input_region = "regions.RData",
    input_count = "counts.RData",
    output_directory = "data", 
    output = "PDPA_poisson_peaks_errors.RData"
  ),
  list(
    post_processing_rule = "thinnest_peak", 
    input_directories = list(d), 
    input_model = "PDPA_poisson.RData",
    input_region = "regions.RData",
    input_count = "counts.RData",
    output_directory = "data", 
    output = "PDPA_poisson_peaks_errors_thinnest_peak_rule.RData"
  ),
  list(
    post_processing_rule = "largest_peak", 
    input_directories = list(d), 
    input_model = "PDPA_poisson.RData",
    input_region = "regions.RData",
    input_count = "counts.RData",
    output_directory = "data", 
    output = "PDPA_poisson_peaks_errors_largest_peak_rule.RData"
  ),
  list(
    post_processing_rule = "maxjump", 
    input_directories = list(d), 
    input_model = "PDPA_negbin.RData",
    input_region = "regions.RData",
    input_count = "counts.RData",
    output_directory = "data", 
    output = "PDPA_negbin_peaks_errors.RData"
  ),
  list(
    post_processing_rule = "thinnest_peak", 
    input_directories = list(d), 
    input_model = "PDPA_negbin.RData",
    input_region = "regions.RData",
    input_count = "counts.RData",
    output_directory = "data", 
    output = "PDPA_negbin_peaks_errors_thinnest_peak_rule.RData"
  ),
  list(
    post_processing_rule = "largest_peak", 
    input_directories = list(d), 
    input_model = "PDPA_negbin.RData",
    input_region = "regions.RData",
    input_count = "counts.RData",
    output_directory = "data", 
    output = "PDPA_negbin_peaks_errors_largest_peak_rule.RData"
  ),
  list(
    post_processing_rule = "maxjump", 
    input_directories = list(d), 
    input_model = "updown.RData",
    input_region = "regions.RData",
    input_count = "counts.RData",
    output_directory = "data", 
    output = "updown_peaks_errors.RData"
  ),
  list(
    post_processing_rule = "maxjump", 
    input_directories = list(d), 
    input_model = "updown_poisson.RData",
    input_region = "regions.RData",
    input_count = "counts.RData",
    output_directory = "data", 
    output = "updown_poisson_peaks_errors.RData"
  ),
  list(
    post_processing_rule = "maxjump", 
    input_directories = list(d), 
    input_model = "updown_negbin.RData",
    input_region = "regions.RData",
    input_count = "counts.RData",
    output_directory = "data", 
    output = "updown_negbin_peaks_errors.RData"
  )
))

plan(multiprocess, workers = 20)
pwalk(m, peakCalling)
plan(sequential)

# HEURISTICS

files <- Sys.glob("data/chunks/*/*/error/hmcan.broad.trained.RData")
str_extract(files, "H3.*/\\d+")
HMCAN_peaks_errors <- rbindlist(map(files, function(file_p){
  load(file_p)
  setDT(error)
  error[, chunk.name := str_extract(file_p, "H3.*/\\d+")]
  error
}))
save(HMCAN_peaks_errors, file="data/HMCAN_peaks_errors.RData")

files <- Sys.glob("data/chunks/*/*/error/macs.trained.RData")
str_extract(files, "H3.*/\\d+")
MACS_peaks_errors <- rbindlist(map(files, function(file_p){
  load(file_p)
  setDT(error)
  error[, chunk.name := str_extract(file_p, "H3.*/\\d+")]
  error
}))
save(MACS_peaks_errors, file="data/MACS_peaks_errors.RData")