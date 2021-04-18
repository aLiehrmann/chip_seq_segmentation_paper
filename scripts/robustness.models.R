library(data.table)
library(purrr)
library(furrr)
library(ggplot2)
library(stringr)
library(aricode)


#HEURISTICS
para <- tibble::tribble(
  ~ model_f, ~ peaks_errors_f, ~ predicted_models_f, ~ models_f,
  "MACS", "MACS_peaks_errors","predicted_models_MACS", "macs.trained",
  "HMCan", "HMCAN_peaks_errors", "predicted_models_HMCAN", "hmcan.broad.trained"
)

signi <- function(p.values) {
  unlist(lapply(p.values, function(p.val){
    if (p.val > 0.05) {
      "ns"
    } else if (p.val > 0.01) {
      "*"
    } else if (p.val > 0.001) {
      "**"
    } else {
      "***"
    }
  }))
}

plan(multisession, workers=20)

classifications_heuristics <- rbindlist(pmap(para, function(model_f, peaks_errors_f, predicted_models_f, models_f){
  load(paste0("data_07_04_2021/", predicted_models_f, ".RData"))
  setkey(predicted_models, chunk.name)
  classification <- rbindlist(future_map(split(predicted_models, by=c("chunk.name"), drop=TRUE), 
    function(predicted_models_p){
      print(predicted_models_p$chunk.name[[1]])
      load(paste0(
        "data/chunks/",
        predicted_models_p$chunk.name[[1]],
        "/",
        models_f,
        ".RData"))
      load(paste0(
        "data/chunks/",
        predicted_models_p$chunk.name[[1]],
        "/",
        "counts.RData"))
      load(paste0(
        "data/chunks/",
        predicted_models_p$chunk.name[[1]],
        "/",
        "regions.RData"))
      setDT(regions)
      lim <- rbindlist(map(split(unique(regions[,list(sample.id, cell.type)]), by=c("cell.type"), drop=TRUE), function(x){
        lim <- rbindlist(pmap(x, function(sample.id, cell.type){
          counts_c <- counts[counts$sample.id == sample.id,]
          data.table(
            chromStart = counts_c$chromStart[[1]],
            chromEnd = counts_c$chromEnd[[length(counts_c$chromEnd)]]
          )
        }))
        data.table(
          chromStart=min(lim$chromStart), 
          chromEnd=max(lim$chromEnd),
          cell.type = x$cell.type[[1]]
        )
      }))
      regions <- data.frame(regions)
      classifications <- rbindlist(pmap(predicted_models_p[,list(sample.id, param.name)], 
        function(sample.id, param.name){
          tryCatch({
          peaks_c <- peaks[[param.name]][
            peaks[[param.name]]$sample.id==sample.id,
          ]
          counts_c <- counts[counts$sample.id == sample.id,]
          regions_c <- regions[regions$sample.id == sample.id,]
          counts_c[1,]$chromStart <- lim[
            cell.type == counts_c$cell.type[[1]]
          ]$chromStart
          counts_c[nrow(counts_c),]$chromEnd <- lim[
            cell.type == counts_c$cell.type[[1]]
          ]$chromEnd
          classification <- rep(0, counts_c[nrow(counts_c),]$chromEnd-counts_c[1,]$chromStart+1)
          start <- counts_c$chromStart[[1]]
          pwalk(peaks_c, function(sample.id, chromStart, chromEnd){
            if ((chromEnd-start)<length(classification) & (chromStart-start)>0) {
            classification[(chromStart-start):(chromEnd-start)] <<-1}
          })
          data.table(
            chunk.name = predicted_models_p$chunk.name[[1]],
            cell.type = regions_c$cell.type[[1]], 
            sample.id = sample.id,
            classification = list(rle(classification))
          )}, error = function(error){data.table()})
        }
      ))
      classifications
    }
  ))
  classification[,model:=model_f]
}))

plan(sequential)  

save(classifications_heuristics, file="data/classifications_heuristics.RData")
load(file="data/classifications_heuristics.RData")

classifications_heuristics <- classifications_heuristics[
  (str_detect(chunk.name, "H3K4me3") & model=="MACS") |
  (str_detect(chunk.name, "H3K36me3") & model=="HMCan") 
]

plan(multisession, workers=25)
library(gtools)
res_heuristics <- rbindlist(future_map(split(classifications_heuristics, by=c("chunk.name", "cell.type", "model"),drop=TRUE), function(x){
  if (nrow(x)>5){
  m <- combinations(nrow(x),2)[sample(nrow(combinations(nrow(x),2)),min(50, nrow(combinations(nrow(x),2)))),]
  print(m)
  data.table(ari = mean(unlist(map(1:nrow(m), function(y){
    NID(#ARI(
      inverse.rle(x[m[y,1],]$classification[[1]]),
      inverse.rle(x[m[y,2],]$classification[[1]])
    )
  }))), chunk.name=x$chunk.name[[1]], cell.type=x$cell.type[[1]], model=x$model[[1]])
  } else (data.table())
}))
plan(sequential)  


res_heuristics[, type_peak:=ifelse(str_detect(chunk.name, "H3K4me3"),"H3K4me3 (sharp peaks)", "H3K36me3 (broad peaks)")]

save(res_heuristics, file="data/res_nid_heuristics.RData")
save(res_heuristics, file="data/res_ari_heuristics.RData")
# load(file="data/res_ari_heuristics.RData")
# load(file="data/res_nid_heuristics.RData")
# res_heuristics[,mean(ari, na.rm=TRUE), by=list(model,type_peak)]
# res_heuristics[,.N, by=list(model,type_peak)]
# ggplot(
#   data=res_heuristics, aes(x=ari,y=factor(model))
# )+
# facet_grid(type_peak~.)+
# geom_boxplot()


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

# updown_poisson_peaks <- maxjump(
#         cp = updown_poisson_models_c$changepoints[
#           -length(updown_poisson_models_c$changepoints)
#         ],
#         theta = updown_poisson_models_c$parameters,
#         genomic_position=counts_c$chromEnd
#       )

para <- tibble::tribble(
  ~ model_f, ~ peaks_errors_f, ~ predicted_models_f, ~ lambda_f, ~ models_f,
  "unconstrained gaussian", "PDPA_peaks_errors","predicted_models_PDPA_gaussian_constant_lambda_10", "constant", "PDPA", 
  "unconstrained negative binomial", "PDPA_negbin_2_peaks_errors", "predicted_models_PDPA_negbin_constant_lambda_constant_phi_10", "constant", "PDPA_negbin_2",
  "unconstrained poisson", "PDPA_poisson_peaks_errors", "predicted_models_PDPA_poisson_constant_lambda_10", "constant", "PDPA_poisson",
  "updown gaussian", "updown_peaks_errors", "predicted_models_updown_gaussian_constant_lambda_10", "constant", "updown",
  "updown poisson", "updown_poisson_peaks_errors", "predicted_models_updown_poisson_constant_lambda_10", "constant", "updown_poisson",
  "updown negative binomial", "updown_negbin_2_peaks_errors", "predicted_models_updown_negbin_constant_lambda_constant_phi_10", "constant", "updown_negbin_2",
  "unconstrained gaussian", "PDPA_peaks_errors","predicted_models_PDPA_gaussian_linear_lambda_10", "linear", "PDPA",
  "unconstrained negative binomial", "PDPA_negbin_2_peaks_errors", "predicted_models_PDPA_negbin_linear_lambda_constant_phi_10", "linear", "PDPA_negbin_2",
  "unconstrained poisson", "PDPA_poisson_peaks_errors", "predicted_models_PDPA_poisson_linear_lambda_10", "linear", "PDPA_poisson",
  "updown gaussian", "updown_peaks_errors", "predicted_models_updown_gaussian_linear_lambda_10", "linear", "updown",
  "updown poisson", "updown_poisson_peaks_errors", "predicted_models_updown_poisson_linear_lambda_10", "linear", "updown_poisson",
  "updown negative binomial", "updown_negbin_2_peaks_errors", "predicted_models_updown_negbin_linear_lambda_constant_phi_10", "linear", "updown_negbin_2"
)


plan(multisession, workers=20)


classifications <- rbindlist(pmap(para, function(model_f, peaks_errors_f, predicted_models_f, lambda_f, models_f){
  load(paste0("data_07_04_2021/", predicted_models_f, ".RData"))
  setkey(predicted_models, chunk.name)

  load(paste0("data_07_04_2021/", predicted_models_f, ".RData"))
  setkey(predicted_models, chunk.name)
  load(paste0("data/",peaks_errors_f,".RData"))
  peaks_errors <- eval(parse(text=peaks_errors_f))
  if (!"phi"%in%names(peaks_errors)) {
    peaks_errors[,phi:=NA]
  }
  classification <- rbindlist(future_map(split(predicted_models, by=c("chunk.name"), drop=TRUE), 
    function(predicted_models_p){
      print(predicted_models_p$chunk.name[[1]])
      load(paste0(
        "data/chunks/",
        predicted_models_p$chunk.name[[1]],
        "/",
        models_f,
        ".RData"))
      load(paste0(
        "data/chunks/",
        predicted_models_p$chunk.name[[1]],
        "/",
        "counts.RData"))
      load(paste0(
        "data/chunks/",
        predicted_models_p$chunk.name[[1]],
        "/",
        "regions.RData"))
      setDT(regions)
      # peut être pas important par type cellulaire, ~ même limites
      lim <- rbindlist(map(split(unique(regions[,list(sample.id, cell.type)]), by=c("cell.type"), drop=TRUE), function(x){
        lim <- rbindlist(pmap(x, function(sample.id, cell.type){
          counts_c <- counts[counts$sample.id == sample.id,]
          data.table(
            chromStart = counts_c$chromStart[[1]],
            chromEnd = counts_c$chromEnd[[length(counts_c$chromEnd)]]
          )
        }))
        data.table(
          chromStart=min(lim$chromStart), 
          chromEnd=max(lim$chromEnd),
          cell.type = x$cell.type[[1]]
        )
      }))
      regions <- data.frame(regions)
      classifications <- rbindlist(pmap(predicted_models_p[,list(sample.id, changepoints, phi)], 
        function(sample.id, changepoints, phi){
          tryCatch({
          changepoints_p <- changepoints
          phi_p <- phi
          if (!"phi"%in%names(dt.models)) {
            fit <- dt.models[sample == sample.id & changepoints==changepoints_p, ]$fit[[1]]
          } else {
            fit <- dt.models[sample == sample.id & phi==phi_p & changepoints==changepoints_p, ]$fit[[1]]
          }
          counts_c <- counts[counts$sample.id == sample.id,]
          regions_c <- regions[regions$sample.id == sample.id,]
          counts_c[1,]$chromStart <- lim[
            cell.type == counts_c$cell.type[[1]]
          ]$chromStart
          counts_c[nrow(counts_c),]$chromEnd <- lim[
            cell.type == counts_c$cell.type[[1]]
          ]$chromEnd
          classification <- rep(0, counts_c[nrow(counts_c),]$chromEnd-counts_c[1,]$chromStart+1)
          start <- counts_c$chromStart[[1]]
          
          peaks_c <- maxjump(
            cp = fit$changepoints[
            -length(fit$changepoints)
            ],
            theta = fit$parameters,
            genomic_position=counts_c$chromEnd
          )

          pwalk(peaks_c, function(chromStart, chromEnd){
            if ((chromEnd-start)<length(classification) & (chromStart-start)>0) {
            classification[(chromStart-start):(chromEnd-start)] <<-1}
          })
          data.table(
            chunk.name = predicted_models_p$chunk.name[[1]],
            cell.type = regions_c$cell.type[[1]], 
            sample.id = sample.id,
            classification = list(rle(classification))
          )}, error = function(error){data.table()})
        }
      ))
      classifications
    }
  ))
  classification[,lambda := lambda_f]
  classification[,model:=model_f]
}))

plan(sequential) 


save(classifications, file = "data/classifications.RData")
load(file = "data/classifications.RData")

classifications <- classifications[
  (str_detect(chunk.name, "H3K4me3") & lambda=="linear") |
  (str_detect(chunk.name, "H3K36me3") & lambda=="constant") 
]
classifications[,lambda:=NULL]


plan(multisession, workers=25)
library(gtools)
res <- rbindlist(future_map(split(classifications, by=c("chunk.name", "cell.type", "model"),drop=TRUE), function(x){
  if (nrow(x)>5){
  m <- combinations(nrow(x),2)[sample(nrow(combinations(nrow(x),2)),min(50, nrow(combinations(nrow(x),2)))),]
  print(m)
  data.table(ari = mean(unlist(map(1:nrow(m), function(y){
    NID(#ARI(
      inverse.rle(x[m[y,1],]$classification[[1]]),
      inverse.rle(x[m[y,2],]$classification[[1]])
    )
  }))), chunk.name=x$chunk.name[[1]], cell.type=x$cell.type[[1]], model=x$model[[1]])
  } else (data.table())
}))
plan(sequential)  


res[, type_peak:=ifelse(str_detect(chunk.name, "H3K4me3"),"H3K4me3 (sharp peaks)", "H3K36me3 (broad peaks)")]

save(res, file = "data/res_nid.RData")
save(res, file = "data/res_ari.RData")
# load(file = "data/res_ari.RData")
# load(file = "data/res_nid.RData")

# res[,mean(ari,na.rm=TRUE), by=list(model,type_peak)]
# res_heuristics[,mean(ari,na.rm=TRUE), by=list(model,type_peak)]
# res[,.N, by=list(model,type_peak)]
# ggplot(
#   data=res, aes(x=ari,y=factor(model))
# )+
# facet_grid(type_peak~.)+
# geom_boxplot()