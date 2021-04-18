library(data.table)
library(penaltyLearning)
library(furrr)
library(purrr)

plan(multiprocess, workers=20)
path_to_input_data <- "data"
load(paste0(path_to_input_data, "/dp.peaks.sets.RData"))
load(paste0(path_to_input_data, "/all.features.RData"))


constant_lambda <- function(validation_chunks, training_models, i){ #training_chunks, i){
  setkey(all.modelSelection,chunk.name)
  validation_models <- all.modelSelection[validation_chunks]
  training_models <- training_models[!is.na(min.log.lambda),]
  validation_models <- validation_models[!is.na(min.log.lambda),]
  sorted.lambda <- sort(c(training_models$min.log.lambda, training_models$max.log.lambda))
  grid <- rle(sorted.lambda)$values
  errors <- vector(mode="integer", length= length(grid)-1)
  walk(1:nrow(training_models), function(i){
    a <- which(grid>=training_models$min.log.lambda[[i]] & grid<training_models$max.log.lambda[[i]])
    errors[a] <<- errors[a] + training_models$total.errors[[i]]
  })
  i.min <- which.min(errors)
  predicted.lambda <- sum(grid[i.min:(i.min+1)])/2
  predicted.models <- validation_models[predicted.lambda>min.log.lambda & predicted.lambda<max.log.lambda,]
  acc <- (1 - sum(predicted.models$total.errors) / sum(predicted.models$labels)) * 100
  data.table(samples = i, accuracy = acc)
}

constant_lambda_constant_phi <- function(validation_chunks, training_models, i){
  setkey(all.modelSelection,chunk.name)
  validation_models <- all.modelSelection[validation_chunks]
  training_models <- training_models[!is.na(min.log.lambda),]
  validation_models <- validation_models[!is.na(min.log.lambda),]
  training_models <- split(training_models, by="phi")
  dt <- do.call(rbind, map(training_models, function(dt_by_phi) {
    sorted.lambda <- sort(c(dt_by_phi$min.log.lambda, dt_by_phi$max.log.lambda))
    grid <- rle(sorted.lambda)$values
    errors <- vector(mode="integer", length= length(grid)-1)
    walk(1:nrow(dt_by_phi), function(i){
      a <- which(grid>=dt_by_phi$min.log.lambda[[i]] & grid<dt_by_phi$max.log.lambda[[i]])
      errors[a] <<- errors[a] + dt_by_phi$total.errors[[i]]
    })
    i.min <- which.min(errors)
    predicted.lambda <- sum(grid[i.min:(i.min+1)])/2
    data.table(lambda = predicted.lambda, phi = dt_by_phi$phi[[1]], errors = min(errors))
  }))
  predicted.lambda <- dt$lambda[which.min(dt$errors)]
  predicted.phi <- dt$phi[which.min(dt$errors)]
  predicted.models <- validation_models[predicted.phi == phi & predicted.lambda>min.log.lambda & predicted.lambda<max.log.lambda,]
  acc <- (1 - sum(predicted.models$total.errors) / sum(predicted.models$labels)) * 100
  data.table(samples= i, accuracy = acc)
}

methods <- list(
  "constant_lambda" = constant_lambda,
  "constant_lambda_constant_phi" = constant_lambda_constant_phi
)

get_model_name <- function(algo, lambda, phi)
{
  paste0(
    algo, "_",
    lambda,"_lambda",
    ifelse(is.na(phi), "", paste0("_", phi, "_phi"))
  )
}

get_method_name <- function(lambda, phi)
{
  paste0(
    lambda,"_lambda",
    ifelse(is.na(phi),"",paste0("_", phi, "_phi"))
  )
}

CV <- function(algo, targets, models, lambda, phi) {
  m <- get_method_name(lambda, phi)
  model_name <- get_model_name(algo, lambda, phi)
  if (!is.na(targets)) {
    load(paste0(path_to_input_data,"/",targets,".RData"), envir = .GlobalEnv)
  }
  load(paste0(path_to_input_data,"/",models,".RData"), envir = .GlobalEnv)
  res_cv <- do.call(rbind, imap(dp.peaks.sets, function(dataset, dataset.i){
    res_cv_dataset <- do.call(rbind, future_imap(dataset, function(fold, fold.i) {
      validation_chunks <- unlist(fold)
      training_chunks <- unlist(dataset[-fold.i])
      setkey(all.modelSelection,chunk.name)
      training_models <- all.modelSelection[training_chunks]
      setkey(training_models, chunk.name, sample.id)
      training_targets <- unique(training_models[,list(chunk.name, sample.id)])
      t <- sample.int(n=nrow(training_targets), replace=F)
      do.call(rbind,map(1:40, function(x){
        training_models_ <- training_models[training_targets[t[1:x],]]
        labels <- sum(unique(training_models_[,list(chunk.name, sample.id, labels)])$labels)
        res <- methods[[m]](validation_chunks, training_models_, x)
        res[,fold.i:=fold.i]
        res[,labels:=labels]
      }))
    }, .options = furrr_options(seed=NULL)))
    res_cv_dataset[,set.name := dataset.i]
    res_cv_dataset
  }))
  res_cv[,model := model_name]
  cv <- length(dp.peaks.sets[[1]])
  save(res_cv, file = paste0("data/test_error_vs_nb_labels_res_cv_",model_name,"_",cv,".RData"))
}

params <- tibble::tribble(
  ~ algo,               ~ targets,                         ~ models,                               ~ lambda,     ~ phi,       
  "PDPA_gaussian",      "all_targets_PDPA",                "all_modelSelection_PDPA",              "constant",   NA,        
  "PDPA_poisson",       "all_targets_PDPA_poisson",        "all_modelSelection_PDPA_poisson",      "constant",   NA,       
  "updown_gaussian",    "all_targets_updown",              "all_modelSelection_updown",            "constant",   NA,        
  "updown_poisson",     "all_targets_updown_poisson",      "all_modelSelection_updown_poisson",    "constant",   NA,        
  "updown_negbin",      "all_targets_updown_negbin",       "all_modelSelection_updown_negbin",     "constant",   "constant",
  "PDPA_negbin",        "all_targets_PDPA_negbin",         "all_modelSelection_PDPA_negbin",       "constant",   "constant",
)

pwalk(params, CV)