library(data.table)
library(purrr)
library(stringr)
library(ggplot2)
library(furrr)
plan(multisession, workers=15)
options(future.globals.maxSize= 891289600)

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
}


# peak errors 

load("data/HMCAN_peaks_errors.RData")
HMCAN_peaks_errors$sample.id <- levels(HMCAN_peaks_errors$sample.id)[
  HMCAN_peaks_errors$sample.id
]
HMCAN_peaks_errors$param.name <- levels(HMCAN_peaks_errors$param.name)[
  HMCAN_peaks_errors$param.name
]
load("data/MACS_peaks_errors.RData")
MACS_peaks_errors$sample.id <- levels(MACS_peaks_errors$sample.id)[
  MACS_peaks_errors$sample.id
]
MACS_peaks_errors$param.name <- levels(MACS_peaks_errors$param.name)[
  MACS_peaks_errors$param.name
]
load("data/PDPA_peaks_errors.RData")
load("data/PDPA_negbin_2_peaks_errors.RData")
load("data/updown_poisson_peaks_errors.RData")

# predicted models 

load("data_07_04_2021/predicted_models_HMCAN.RData")
HMCAN_predicted_models <- predicted_models
load("data_07_04_2021/predicted_models_MACS.RData")
MACS_predicted_models <- predicted_models
load("data_07_04_2021/predicted_models_PDPA_gaussian_constant_lambda_10.RData")
PDPA_gaussian_constant_lambda_predicted_models <- predicted_models
load("data_07_04_2021/predicted_models_PDPA_gaussian_linear_lambda_10.RData")
PDPA_gaussian_linear_lambda_predicted_models <- predicted_models
load("data_07_04_2021/predicted_models_PDPA_negbin_constant_lambda_constant_phi_10.RData")
PDPA_negbin_constant_lambda_constant_phi_predicted_models <- predicted_models
load("data_07_04_2021/predicted_models_PDPA_negbin_linear_lambda_constant_phi_10.RData")
PDPA_negbin_linear_lambda_constant_phi_predicted_models <- predicted_models
load("data_07_04_2021/predicted_models_updown_poisson_constant_lambda_10.RData")
updown_poisson_constant_lambda_predicted_models <- predicted_models
load("data_07_04_2021/predicted_models_updown_poisson_linear_lambda_10.RData")
updown_poisson_linear_lambda_predicted_models <- predicted_models


files <- Sys.glob("data/chunks/*/*")
future_walk(files, function(files_p){
  load(file.path(files_p, "counts.RData"))
  setDT(counts)
  counts$sample.id <- levels(counts$sample.id)[counts$sample.id]
  chunk_name <- str_extract(files_p, "(?<=data/chunks/).*") 
  is_H3K4me3 <- str_detect(chunk_name, "H3K4me3")
  load(file.path(files_p, "PDPA.RData"))
  PDPA_models <- dt.models
  load(file.path(files_p, "PDPA_negbin_2.RData"))
  PDPA_negbin_models <- dt.models
  load(file.path(files_p, "updown_poisson.RData"))
  updown_poisson_models <- dt.models
  load(file.path(files_p, "hmcan.broad.trained.RData"))
  HMCAN_all_peaks <- peaks
  load(file.path(files_p, "macs.trained.RData"))
  MACS_all_peaks <- peaks
  walk(unique(counts$sample.id), function(sample_p){
    tryCatch({
    # // counts
    counts_c <- counts[sample.id==sample_p,]
    counts_c <- data.table(
      chromStart = counts_c$chromStart,
      chromEnd= counts_c$chromEnd,
      coverage = c(
        rep(counts_c$coverage,3),
        sqrt(counts_c$coverage+3/8)
      ),
      label = factor(
        x=rep(
          x=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          ), 
          each=nrow(counts_c)
        ), 
        levels=c(
          ifelse(is_H3K4me3,"MACS","HMCan"),
          "Up-Down\npoisson",
          "unconstrained\nnegative binomial",
          "unconstrained\ngaussian"
        )
      )
    )
    
    if(is_H3K4me3){
      #// MACS
      MACS_predicted_models_c <- MACS_predicted_models[
        sample.id==sample_p &
        chunk.name==chunk_name
      ]
      MACS_peaks_errors_c <- MACS_peaks_errors[
        sample.id==sample_p &
        chunk.name==chunk_name &  
        param.name == MACS_predicted_models_c$param.name
      ]
      heuristics_peaks_errors_c <- data.table(
        chromStart = MACS_peaks_errors_c$chromStart, 
        chromEnd = MACS_peaks_errors_c$chromEnd, 
        annotation = MACS_peaks_errors_c$annotation, 
        status = MACS_peaks_errors_c$status,
        label = factor(
        x= "MACS", 
        levels=c(
          ifelse(is_H3K4me3,"MACS","HMCan"),
          "Up-Down\npoisson",
          "unconstrained\nnegative binomial",
          "unconstrained\ngaussian"
          )
        )
      )
      tmp_peaks <- MACS_all_peaks[[MACS_predicted_models_c$param.name]]
      tmp_peaks[[1]] <- levels(tmp_peaks[[1]])[tmp_peaks[[1]]]
      heuristics_peaks <- tmp_peaks[
          tmp_peaks[[1]]==sample_p,
      ]
      setDT(heuristics_peaks)
      heuristics_peaks[,sample.id:=NULL]
      heuristics_peaks[,label := factor(
        x= "MACS", 
        levels=c(
          ifelse(is_H3K4me3,"MACS","HMCan"),
          "Up-Down\npoisson",
          "unconstrained\nnegative binomial",
          "unconstrained\ngaussian"
          )
      )]
      # // PDPA gaussian
      PDPA_predicted_models_c <- PDPA_gaussian_linear_lambda_predicted_models[
        sample.id==sample_p &
        chunk.name==chunk_name
      ]
      PDPA_peaks_errors_c <- PDPA_peaks_errors[
        sample.id==sample_p &
        chunk.name==chunk_name &  
        changepoints == PDPA_predicted_models_c$changepoints
      ]
      PDPA_peaks_errors_c <- data.table(
        chromStart = PDPA_peaks_errors_c$chromStart, 
        chromEnd = PDPA_peaks_errors_c$chromEnd, 
        annotation = PDPA_peaks_errors_c$annotation, 
        status = PDPA_peaks_errors_c$status,
        label = factor(
        x= "unconstrained\ngaussian", 
        levels=c(
          ifelse(is_H3K4me3,"MACS","HMCan"),
          "Up-Down\npoisson",
          "unconstrained\nnegative binomial",
          "unconstrained\ngaussian"
          )
        )
      )
      PDPA_models_c <- PDPA_models[
        sample==sample_p & 
        changepoints == PDPA_predicted_models_c$changepoints
      ]$fit[[1]]
      PDPA_changepoints <- data.table(
        pos=counts_c$chromEnd[
          PDPA_models_c$changepoints[-length(PDPA_models_c$changepoints)]
        ],
        label = factor(
          x= "unconstrained\ngaussian", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
        )
      )
      PDPA_segments <- data.table(
        chromStart=counts_c$chromEnd[c(1,
          PDPA_models_c$changepoints[1:(length(PDPA_models_c$changepoints)-1)]
        )],
        chromEnd=counts_c$chromEnd[
          PDPA_models_c$changepoints[1:length(PDPA_models_c$changepoints)]
        ],
        theta = PDPA_models_c$parameters,
        label = factor(
          x= "unconstrained\ngaussian", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
        )
      )
      PDPA_peaks <- maxjump(
        cp = PDPA_models_c$changepoints[
          -length(PDPA_models_c$changepoints)
        ],
        theta = PDPA_models_c$parameters,
        genomic_position=counts_c$chromEnd
      )
      PDPA_peaks[, label := factor(
          x= "unconstrained\ngaussian", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
      )]
      # // PDPA negbin
      PDPA_negbin_predicted_models_c <- PDPA_negbin_linear_lambda_constant_phi_predicted_models[
        sample.id==sample_p &
        chunk.name==chunk_name
      ]
      PDPA_negbin_peaks_errors_c <- PDPA_negbin_2_peaks_errors[
        sample.id==sample_p &
        chunk.name==chunk_name &  
        changepoints == PDPA_negbin_predicted_models_c$changepoints &
        phi == PDPA_negbin_predicted_models_c$phi
      ]
      PDPA_negbin_peaks_errors_c <- data.table(
        chromStart = PDPA_negbin_peaks_errors_c$chromStart, 
        chromEnd = PDPA_negbin_peaks_errors_c$chromEnd, 
        annotation = PDPA_negbin_peaks_errors_c$annotation, 
        status = PDPA_negbin_peaks_errors_c$status,
        label = factor(
        x= "unconstrained\nnegative binomial", 
        levels=c(
          ifelse(is_H3K4me3,"MACS","HMCan"),
          "Up-Down\npoisson",
          "unconstrained\nnegative binomial",
          "unconstrained\ngaussian"
          )
        )
      )
      PDPA_negbin_models_c <- PDPA_negbin_models[
        sample==sample_p & 
        changepoints == PDPA_negbin_predicted_models_c$changepoints
      ]$fit[[1]]
      PDPA_negbin_changepoints <- data.table(
        pos=counts_c$chromEnd[
          PDPA_negbin_models_c$changepoints[-length(PDPA_negbin_models_c$changepoints)]
        ],
        label = factor(
          x= "unconstrained\nnegative binomial", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
        )
      )
      s <- c(1, PDPA_negbin_models_c$changepoints[
        1:(length(PDPA_negbin_models_c$changepoints)-1)
      ])
      e <- PDPA_negbin_models_c$changepoints[
        1:length(PDPA_negbin_models_c$changepoints)
      ]
      PDPA_negbin_segments <- data.table(
        chromStart=counts_c$chromEnd[s],
        chromEnd=counts_c$chromEnd[e],
        theta = pmap_dbl(list(s=s,e=e), function(s,e){
          mean(rep(
            counts_c$coverage[s:e], 
            counts_c$chromEnd[s:e] - counts_c$chromStart[s:e]
          ))}
        ),
        label = factor(
          x= "unconstrained\nnegative binomial", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
        )
      )
      PDPA_negbin_peaks <- maxjump(
        cp = PDPA_negbin_models_c$changepoints[
          -length(PDPA_negbin_models_c$changepoints)
        ],
        theta = PDPA_negbin_segments$theta,
        genomic_position=counts_c$chromEnd
      )
      PDPA_negbin_peaks[, label := factor(
          x= "unconstrained\nnegative binomial", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
      )]  
      # // updoown poisson
      updown_poisson_predicted_models_c <- updown_poisson_linear_lambda_predicted_models[
        sample.id==sample_p &
        chunk.name==chunk_name
      ]
      updown_poisson_peaks_errors_c <- updown_poisson_peaks_errors[
        sample.id==sample_p &
        chunk.name==chunk_name &  
        changepoints == updown_poisson_predicted_models_c$changepoints
      ]
      updown_poisson_peaks_errors_c <- data.table(
        chromStart = updown_poisson_peaks_errors_c$chromStart, 
        chromEnd = updown_poisson_peaks_errors_c$chromEnd, 
        annotation = updown_poisson_peaks_errors_c$annotation, 
        status = updown_poisson_peaks_errors_c$status,
        label = factor(
        x= "Up-Down\npoisson", 
        levels=c(
          ifelse(is_H3K4me3,"MACS","HMCan"),
          "Up-Down\npoisson",
          "unconstrained\nnegative binomial",
          "unconstrained\ngaussian"
          )
        )
      )
      updown_poisson_models_c <- updown_poisson_models[
        sample==sample_p & 
        changepoints == updown_poisson_predicted_models_c$changepoints
      ]$fit[[1]]
      updown_poisson_changepoints <- data.table(
        pos=counts_c$chromEnd[
          updown_poisson_models_c$changepoints[-length(updown_poisson_models_c$changepoints)]
        ],
        label = factor(
          x= "Up-Down\npoisson", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
        )
      )
      updown_poisson_segments <- data.table(
        chromStart=counts_c$chromEnd[c(1,
          updown_poisson_models_c$changepoints[1:(length(updown_poisson_models_c$changepoints)-1)]
        )],
        chromEnd=counts_c$chromEnd[
          updown_poisson_models_c$changepoints[1:length(updown_poisson_models_c$changepoints)]
        ],
        theta = updown_poisson_models_c$parameters,
        label = factor(
          x= "Up-Down\npoisson", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
        )
      ) 
      updown_poisson_peaks <- maxjump(
        cp = updown_poisson_models_c$changepoints[
          -length(updown_poisson_models_c$changepoints)
        ],
        theta = updown_poisson_models_c$parameters,
        genomic_position=counts_c$chromEnd
      )
      updown_poisson_peaks[, label := factor(
          x= "Up-Down\npoisson", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
      )]
    } else {
      HMCAN_predicted_models_c <- HMCAN_predicted_models[
        sample.id==sample_p &
        chunk.name==chunk_name
      ]
      HMCAN_peaks_errors_c <- HMCAN_peaks_errors[
        sample.id==sample_p &
        chunk.name==chunk_name &  
        param.name == HMCAN_predicted_models_c$param.name
      ]
      heuristics_peaks_errors_c <- data.table(
        chromStart = HMCAN_peaks_errors_c$chromStart, 
        chromEnd = HMCAN_peaks_errors_c$chromEnd, 
        annotation = HMCAN_peaks_errors_c$annotation, 
        status = HMCAN_peaks_errors_c$status,
        label = factor(
        x= "HMCan", 
        levels=c(
          ifelse(is_H3K4me3,"MACS","HMCan"),
          "Up-Down\npoisson",
          "unconstrained\nnegative binomial",
          "unconstrained\ngaussian"
          )
        )
      )
      tmp_peaks <- HMCAN_all_peaks[[HMCAN_predicted_models_c$param.name]]
      tmp_peaks[[1]] <- levels(tmp_peaks[[1]])[tmp_peaks[[1]]]
      heuristics_peaks <- tmp_peaks[
          tmp_peaks[[1]]==sample_p,
      ]
      setDT(heuristics_peaks)
      heuristics_peaks[,sample.id:=NULL]
      heuristics_peaks[,label := factor(
        x= "HMCan", 
        levels=c(
          ifelse(is_H3K4me3,"MACS","HMCan"),
          "Up-Down\npoisson",
          "unconstrained\nnegative binomial",
          "unconstrained\ngaussian"
          )
      )]
            # // PDPA gaussian
      PDPA_predicted_models_c <- PDPA_gaussian_constant_lambda_predicted_models[
        sample.id==sample_p &
        chunk.name==chunk_name
      ]
      PDPA_peaks_errors_c <- PDPA_peaks_errors[
        sample.id==sample_p &
        chunk.name==chunk_name &  
        changepoints == PDPA_predicted_models_c$changepoints
      ]
      PDPA_peaks_errors_c <- data.table(
        chromStart = PDPA_peaks_errors_c$chromStart, 
        chromEnd = PDPA_peaks_errors_c$chromEnd, 
        annotation = PDPA_peaks_errors_c$annotation, 
        status = PDPA_peaks_errors_c$status,
        label = factor(
        x= "unconstrained\ngaussian", 
        levels=c(
          ifelse(is_H3K4me3,"MACS","HMCan"),
          "Up-Down\npoisson",
          "unconstrained\nnegative binomial",
          "unconstrained\ngaussian"
          )
        )
      )
      PDPA_models_c <- PDPA_models[
        sample==sample_p & 
        changepoints == PDPA_predicted_models_c$changepoints
      ]$fit[[1]]
      PDPA_changepoints <- data.table(
        pos=counts_c$chromEnd[
          PDPA_models_c$changepoints[-length(PDPA_models_c$changepoints)]
        ],
        label = factor(
          x= "unconstrained\ngaussian", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
        )
      )
      PDPA_segments <- data.table(
        chromStart=counts_c$chromEnd[c(1,
          PDPA_models_c$changepoints[1:(length(PDPA_models_c$changepoints)-1)]
        )],
        chromEnd=counts_c$chromEnd[
          PDPA_models_c$changepoints[1:length(PDPA_models_c$changepoints)]
        ],
        theta = PDPA_models_c$parameters,
        label = factor(
          x= "unconstrained\ngaussian", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
        )
      )
      PDPA_peaks <- maxjump(
        cp = PDPA_models_c$changepoints[
          -length(PDPA_models_c$changepoints)
        ],
        theta = PDPA_models_c$parameters,
        genomic_position=counts_c$chromEnd
      )
      PDPA_peaks[, label := factor(
          x= "unconstrained\ngaussian", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
      )]
      # // PDPA negbin
      PDPA_negbin_predicted_models_c <- PDPA_negbin_constant_lambda_constant_phi_predicted_models[
        sample.id==sample_p &
        chunk.name==chunk_name
      ]
      PDPA_negbin_peaks_errors_c <- PDPA_negbin_2_peaks_errors[
        sample.id==sample_p &
        chunk.name==chunk_name &  
        changepoints == PDPA_negbin_predicted_models_c$changepoints &
        phi == PDPA_negbin_predicted_models_c$phi
      ]
      PDPA_negbin_peaks_errors_c <- data.table(
        chromStart = PDPA_negbin_peaks_errors_c$chromStart, 
        chromEnd = PDPA_negbin_peaks_errors_c$chromEnd, 
        annotation = PDPA_negbin_peaks_errors_c$annotation, 
        status = PDPA_negbin_peaks_errors_c$status,
        label = factor(
        x= "unconstrained\nnegative binomial", 
        levels=c(
          ifelse(is_H3K4me3,"MACS","HMCan"),
          "Up-Down\npoisson",
          "unconstrained\nnegative binomial",
          "unconstrained\ngaussian"
          )
        )
      )
      PDPA_negbin_models_c <- PDPA_negbin_models[
        sample==sample_p & 
        phi == PDPA_negbin_predicted_models_c$phi & 
        changepoints == PDPA_negbin_predicted_models_c$changepoints
      ]$fit[[1]]
      PDPA_negbin_changepoints <- data.table(
        pos=counts_c$chromEnd[
          PDPA_negbin_models_c$changepoints[-length(PDPA_negbin_models_c$changepoints)]
        ],
        label = factor(
          x= "unconstrained\nnegative binomial", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
        )
      )
      s <- c(1, PDPA_negbin_models_c$changepoints[
        1:(length(PDPA_negbin_models_c$changepoints)-1)
      ])
      e <- PDPA_negbin_models_c$changepoints[
        1:length(PDPA_negbin_models_c$changepoints)
      ]
      PDPA_negbin_segments <- data.table(
        chromStart=counts_c$chromEnd[s],
        chromEnd=counts_c$chromEnd[e],
        theta = pmap_dbl(list(s=s,e=e), function(s,e){
          mean(rep(
            counts_c$coverage[s:e], 
            counts_c$chromEnd[s:e] - counts_c$chromStart[s:e]
          ))}
        ),
        label = factor(
          x= "unconstrained\nnegative binomial", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
        )
      )
      PDPA_negbin_peaks <- maxjump(
        cp = PDPA_negbin_models_c$changepoints[
          -length(PDPA_negbin_models_c$changepoints)
        ],
        theta = PDPA_negbin_segments$theta,
        genomic_position=counts_c$chromEnd
      )
      PDPA_negbin_peaks[, label := factor(
          x= "unconstrained\nnegative binomial", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
      )]
      # // updoown poisson
      updown_poisson_predicted_models_c <- updown_poisson_constant_lambda_predicted_models[
        sample.id==sample_p &
        chunk.name==chunk_name
      ]
      updown_poisson_peaks_errors_c <- updown_poisson_peaks_errors[
        sample.id==sample_p &
        chunk.name==chunk_name &  
        changepoints == updown_poisson_predicted_models_c$changepoints
      ]
      updown_poisson_peaks_errors_c <- data.table(
        chromStart = updown_poisson_peaks_errors_c$chromStart, 
        chromEnd = updown_poisson_peaks_errors_c$chromEnd, 
        annotation = updown_poisson_peaks_errors_c$annotation, 
        status = updown_poisson_peaks_errors_c$status,
        label = factor(
        x= "Up-Down\npoisson", 
        levels=c(
          ifelse(is_H3K4me3,"MACS","HMCan"),
          "Up-Down\npoisson",
          "unconstrained\nnegative binomial",
          "unconstrained\ngaussian"
          )
        )
      )
      updown_poisson_models_c <- updown_poisson_models[
        sample==sample_p & 
        changepoints == updown_poisson_predicted_models_c$changepoints
      ]$fit[[1]]
      updown_poisson_changepoints <- data.table(
        pos=counts_c$chromEnd[
          updown_poisson_models_c$changepoints[-length(updown_poisson_models_c$changepoints)]
        ],
        label = factor(
          x= "Up-Down\npoisson", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
        )
      )
      updown_poisson_segments <- data.table(
        chromStart=counts_c$chromEnd[c(1,
          updown_poisson_models_c$changepoints[1:(length(updown_poisson_models_c$changepoints)-1)]
        )],
        chromEnd=counts_c$chromEnd[
          updown_poisson_models_c$changepoints[1:length(updown_poisson_models_c$changepoints)]
        ],
        theta = updown_poisson_models_c$parameters,
        label = factor(
          x= "Up-Down\npoisson", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
        )
      )
      updown_poisson_peaks <- maxjump(
        cp = updown_poisson_models_c$changepoints[
          -length(updown_poisson_models_c$changepoints)
        ],
        theta = updown_poisson_models_c$parameters,
        genomic_position=counts_c$chromEnd
      )
      updown_poisson_peaks[, label := factor(
          x= "Up-Down\npoisson", 
          levels=c(
            ifelse(is_H3K4me3,"MACS","HMCan"),
            "Up-Down\npoisson",
            "unconstrained\nnegative binomial",
            "unconstrained\ngaussian"
          )
      )]
    }

    peaks_errors <- rbind(
      heuristics_peaks_errors_c,
      updown_poisson_peaks_errors_c,
      PDPA_negbin_peaks_errors_c,
      PDPA_peaks_errors_c 
    )
    changepoints <- rbind(
      PDPA_changepoints,
      PDPA_negbin_changepoints,
      updown_poisson_changepoints
    )
    changepoints[,items:=factor("changepoint", levels=c("segment mean", "changepoint", "peak"))]

    segments <- rbind(
      PDPA_segments,
      PDPA_negbin_segments,
      updown_poisson_segments
    )
    segments[,items:=factor("segment mean", levels=c("segment mean", "changepoint", "peak"))]

    peaks <- rbind(
      PDPA_peaks,
      PDPA_negbin_peaks,
      updown_poisson_peaks,
      heuristics_peaks
    )
    peaks[,item:=factor("peak", levels=c("segment mean", "changepoint", "peak"))]
    max_cov_tr <- max(counts_c[label=="unconstrained\ngaussian"]$coverage)
    min_cov_tr <- min(counts_c[label=="unconstrained\ngaussian"]$coverage)
    max_cov <- max(counts_c[label=="Up-Down\npoisson"]$coverage)
    min_cov <- min(counts_c[label=="Up-Down\npoisson"]$coverage)
    peaks[,ymin:=ifelse(
      label == "unconstrained\ngaussian", 
      min_cov_tr - (max_cov_tr - min_cov_tr)*0.1,
      min_cov - (max_cov - min_cov)*0.1
    )]
    peaks[,ymax:=ifelse(
      label == "unconstrained\ngaussian", 
      min_cov_tr ,
      min_cov
    )]

    if(!any(peaks$chromStart<5)){
    g <- ggplot()+
    facet_grid(label~., scales="free")+
    geom_line(
      data = counts_c, 
      aes(
        y=coverage, 
        x=chromStart
      ),
      color="grey60"
    )+
    geom_rect(
      data = peaks_errors,
      aes(
        xmin=chromStart,
        xmax=chromEnd,
        ymin=-Inf,
        ymax=+Inf,
        fill=annotation,
        linetype=status
      ),
      size=0.6,
      alpha=0.3,
      color="black"
    )+
    scale_linetype_manual("status:", values=c(1:length(unique(peaks_errors$status))))+
    guides(linetype = guide_legend(override.aes = list(fill = "white")))+
    scale_fill_manual("labels:",values=c("grey85","Royalblue", "pink", "#a445ee"))+
    geom_vline(
      data = changepoints, 
      aes(
        xintercept=pos,
        color=items
      ),
      size=0.6,
      linetype=2
    )+
    geom_segment(
      data = segments, 
      aes(
        x=chromStart,
        xend=chromEnd, 
        y=theta,
        yend=theta,
        color=items
      )
    )+
    scale_color_manual("segmentation model:", values=c("red","Royalblue"))+
    ggnewscale::new_scale_fill()+
    geom_rect(
      data = peaks,
      aes(
        xmin=chromStart,
        xmax=chromEnd,
        ymin=ymin,
        ymax=ymax,
        fill= item
      ),
      size=0.6,
      color="red"
    )+
    scale_fill_manual("results:",values=c("gold"))+
    xlab("position on chromosome")+
    ylab("normalized coverge")+
    theme_bw()+
    theme(
      strip.background = element_rect(fill="grey95"), 
      text = element_text(size=15),
      legend.position = "bottom",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black")
    )
    ifelse(is_H3K4me3,
    g <- g+ggtitle(paste0(
      chunk_name,"/",sample_p,"\n",
      "MACS, parameter=",MACS_predicted_models_c$param.name,"\n",
      "Up-Down poisson, lambda=",round(exp((updown_poisson_predicted_models_c$max.log.lambda + updown_poisson_predicted_models_c$min.log.lambda)/2),3),"\n",
      "unconstrained negative binomial, lambda=",round(exp((PDPA_negbin_predicted_models_c$max.log.lambda + PDPA_negbin_predicted_models_c$min.log.lambda)/2),3), ", phi=", round(PDPA_negbin_predicted_models_c$phi,3),"\n",
      "unconstrained gaussian, lambda=",round(exp((PDPA_predicted_models_c$max.log.lambda + PDPA_predicted_models_c$min.log.lambda)/2),3),"\n"
    )),
    g<-g+ggtitle(paste0(
      chunk_name,"/",sample_p,"\n",
      "HMCan, parameter=",HMCAN_predicted_models_c$param.name,"\n",
      "Up-Down poisson, lambda=",round(exp((updown_poisson_predicted_models_c$max.log.lambda + updown_poisson_predicted_models_c$min.log.lambda)/2),3),"\n",
      "unconstrained negative binomial, lambda=",round(exp((PDPA_negbin_predicted_models_c$max.log.lambda + PDPA_negbin_predicted_models_c$min.log.lambda)/2),3), ", phi=", round(PDPA_negbin_predicted_models_c$phi,3),"\n",
      "unconstrained gaussian, lambda=",round(exp((PDPA_predicted_models_c$max.log.lambda + PDPA_predicted_models_c$min.log.lambda)/2),3),"\n"
    )))
    ggsave(plot=g, paste0("figures/models/",str_replace(chunk_name,"/","_"),"_",sample_p,".pdf"), width=17, height=9.5)}},
    error = function(e){print(paste0(chunk_name,"/",sample_p))})
  })
})








