library(data.table)
source("scripts/pick.best.index.R")
load("data/dp.peaks.sets.4folds.RData")
load("data/dp.peaks.matrices.RData")

default.params <- c(
  macs.trained="1.30103",
  hmcan.broad.trained="2.30258509299405")

elist <- list()
roc.list <- list()
for(set.name in names(dp.peaks.sets)){
  train.sets <- dp.peaks.sets[[set.name]]
  chunk.list <- dp.peaks.matrices[[set.name]]
  tp.list <- dp.peaks.matrices.tp[[set.name]]
  fp.list <- dp.peaks.matrices.fp[[set.name]]
  for(set.i in seq_along(train.sets)){
    testSet <- paste(set.name, "split", set.i)
    test.chunks <- train.sets[[set.i]]
    train.chunks <- names(chunk.list)[! names(chunk.list) %in% test.chunks]
    baseline.list <- list()
    for(algorithm in c("hmcan.broad.trained", "macs.trained")){
      train.mat.list <- list()
      for(train.chunk in train.chunks){
        train.mat.list[[train.chunk]] <-
          chunk.list[[train.chunk]][[algorithm]]
      }
      n_col <- median(unlist(lapply(train.mat.list, ncol)))
      train.mat <- do.call(rbind, lapply(train.mat.list, function(x){if(ncol(x)==n_col) x}))
      error.curve <- colSums(train.mat)
      error.sorted <- error.curve[order(as.numeric(names(error.curve)))]
      picked <- pick.best.index(error.sorted)
      baseline.list[[algorithm]] <- names(error.sorted)[[picked]]
    }
    cat(sprintf("%d / %d %s\n", set.i, length(train.sets), set.name))
    for(test.chunk in test.chunks){
      test.info <- chunk.list[[test.chunk]]
      test.tp <- tp.list[[test.chunk]]
      test.fp <- fp.list[[test.chunk]]
      for(algorithm in names(default.params)){
        param.list <- list(
          supervised=baseline.list[[algorithm]])
        err.mat <- test.info[[algorithm]]
        tp.mat <- test.tp[[algorithm]]
        fp.mat <- test.fp[[algorithm]]
        for(train.type in names(param.list)){
          param.name <- param.list[[train.type]]
          errors <- err.mat[, param.name]
          sample.id <- names(errors)
          elist[[paste(set.name, set.i, test.chunk, algorithm, train.type)]] <- 
            data.table(set.name, set.i, testSet, test.chunk,
                       algorithm=ifelse(
                         algorithm=="macs.trained", "MACS", "HMCanBroad"),
                       train.type,
                       param.name, sample.id, errors,
                       tp=tp.mat[, param.name],
                       possible.tp=test.tp$possible.tp,
                       fp=fp.mat[, param.name],
                       possible.fp=test.fp$possible.fp,
                       regions=test.info$regions)
        }
        ord <- order(as.numeric(colnames(err.mat)))
        roc.dt <- data.table(
          set.name, set.i, test.chunk, 
          algorithm=ifelse(
            algorithm=="macs.trained", "MACS", "HMCanBroad"),
          param.i=NA_integer_,
          tp=colSums(tp.mat),
          possible.tp=sum(test.tp$possible.tp),
          fp=colSums(fp.mat),
          possible.fp=sum(test.fp$possible.fp))[ord,]
        roc.dt[, param.i := 1:.N]
        roc.list[[paste(set.name, set.i, test.chunk, algorithm)]] <-
          roc.dt
      }#algorithm
    }#test.chunk
  }#set.i
}#set

test.error <- do.call(rbind, elist)
roc <- do.call(rbind, roc.list)

levs <- c(
  MACS="MACS(baseline)",
  HMCanBroad="HMCanBroad(baseline)")
test.error[, algorithm := levs[algorithm] ]
roc[, algorithm := levs[algorithm] ]
test.counts <- test.error[algorithm %in% levs, list(
  errors=sum(errors),
  labels=sum(regions),
  tp=sum(tp),
  possible.tp=sum(possible.tp),
  fp=sum(fp),
  possible.fp=sum(possible.fp)
), by=.(set.name, set.i, algorithm, train.type)]
test.counts[, TPR := tp/possible.tp]
test.counts[, FPR := fp/possible.fp]
possible.counts <- test.counts[algorithm==algorithm[1] & train.type=="supervised", {
  list(
    possible.fn=sum(possible.tp),
    possible.fp=sum(possible.fp)
  )
}, by=list(set.name)]
test.ranges <- test.counts[, list(
  min.labels=min(labels),
  max.labels=max(labels)
), by=.(set.name, set.i)]
test.ranges[, stopifnot(min.labels==max.labels)]
test.counts[, percent.accuracy := (1-errors/labels)*100]
test.mean <- test.counts[, list(
  mean.percent=mean(percent.accuracy)
), by=.(set.name, algorithm, train.type)]
roc.total <- roc[algorithm %in% levs, list(
  tp=sum(tp),
  possible.tp=sum(possible.tp),
  fp=sum(fp),
  possible.fp=sum(possible.fp)
), by=.(set.name, set.i, algorithm, param.i)]
roc.total[, TPR := tp/possible.tp]
roc.total[, FPR := fp/possible.fp]
roc.ord <- roc.total[order(param.i),]
roc.not.cvx <- roc.total[, {
  tpr <- TPR[which.max(FPR)]
  tpr <- 1
  list(
    FPR=c(1, FPR, 0, 1, 1),
    TPR=c(tpr, TPR, 0, 0, tpr)
  )
}, by=.(set.name, set.i, algorithm)]
auc <- roc.not.cvx[, list(
  auc=geometry::polyarea(FPR, TPR)
), by=.(set.name, set.i, algorithm)]

test.counts 
res_cv<- test.counts[, list(fold.i = set.i, accuracy = percent.accuracy, set.name, model=algorithm)]
m <- ifelse(res_cv$model == "MACS(baseline)", "MACS", "HMCan")
res_cv[,model := m]
save(res_cv, file="data/res_cv_heuristics.RData")