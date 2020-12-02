library(penaltyLearning)
library(data.table)
library(furrr)
library(purrr)

plan(multiprocess,workers=20)

counts.RData.vec <- Sys.glob(file.path("data/chunks", "H*", "*", "counts.RData"))
chunks <- sub(".*chunks/", "", dirname(counts.RData.vec))

all.features <- do.call(rbind,future_map(seq_along(counts.RData.vec), function(i){
  load(counts.RData.vec[[i]])
  features_dt <- as.data.table(featureMatrix(
    counts, 
    "sample.id", 
    "coverage"
  ), keep.rownames = T)
  names(features_dt)[names(features_dt)=="rn"] <- "sample.id"
  features_dt[,chunk.name:=chunks[[i]]]
}))

save(all.features, file="data/all.features.RData")