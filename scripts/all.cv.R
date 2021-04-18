library(data.table)
library(penaltyLearning)
library(furrr)
library(purrr)
library(Rcpp)

plan(multiprocess, workers=20)

path_to_input_data <- "data"
path_to_output_data <- "data"
load(paste0(path_to_input_data, "/dp.peaks.sets.RData"))
load(paste0(path_to_input_data, "/all.features.RData"))

linear_lambda <- function(validation_chunks, training_chunks, i){
  setkey(all.targets, chunk.name)
  setkey(all.features, chunk.name)
  training_targets <- all.targets[training_chunks]
  training_features <- all.features[training_chunks]
  validation_targets <- all.targets[validation_chunks]
  validation_features <- all.features[validation_chunks]
  training_targets <- training_targets[!is.na(min.log.lambda),]
  validation_targets <- validation_targets[!is.na(min.log.lambda),]
  training_targets <- training_targets[is.finite(min.log.lambda) | is.finite(max.log.lambda),]
  setkey(training_targets, chunk.name, sample.id)
  setkey(training_features, chunk.name, sample.id)
  training_features_targets <- training_features[training_targets]
  setkey(validation_targets, chunk.name, sample.id)
  setkey(validation_features, chunk.name, sample.id)
  validation_features_targets <- validation_features[validation_targets]
  training_targets_mat <- as.matrix(training_features_targets[,list(min.log.lambda, max.log.lambda)])
  training_features_mat <- as.matrix(training_features_targets[,!c("chunk.name", "sample.id", "phi", "min.log.lambda", "max.log.lambda", "errors")])
  validation_features_mat <- as.matrix(validation_features_targets[,!c("chunk.name", "sample.id", "phi", "min.log.lambda", "max.log.lambda", "errors")])
  fit <- IntervalRegressionCV(
    target.mat = training_targets_mat,
    feature.mat = training_features_mat,
  )
  predictions <- validation_features_targets[, list(chunk.name, sample.id, fit$predict(validation_features_mat))]
  names(predictions) <- c("chunk.name", "sample.id", "prediction")
  setkey(predictions, chunk.name, sample.id)
  setkey(all.modelSelection, chunk.name, sample.id)
  predictions <- merge(predictions, all.modelSelection)
  errors <- predictions[prediction>min.log.lambda & prediction<max.log.lambda,]
  acc <- (1 - sum(errors$total.errors) / sum(errors$labels)) * 100
  list(
    acc = data.table(fold.i = i, accuracy = acc),
    predicted_models = errors
  )
}

constant_lambda <- function(validation_chunks, training_chunks, i){
  setkey(all.modelSelection,chunk.name)
  training_models <- all.modelSelection[training_chunks]
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
  list(
    acc = data.table(fold.i = i, accuracy = acc),
    predicted_models = predicted.models
  )
}


constant_lambda_constant_phi <- function(validation_chunks, training_chunks, i){
  setkey(all.modelSelection,chunk.name)
  training_models <- all.modelSelection[training_chunks]
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
    list(
    acc = data.table(fold.i = i, accuracy = acc),
    predicted_models = predicted.models
  )
}


linear_lambda_constant_phi <- function(validation_chunks, training_chunks, i){
  # constant phi
  setkey(all.modelSelection,chunk.name)
  training_models <- all.modelSelection[training_chunks]
  training_models <- training_models[!is.na(min.log.lambda),]
  training_models <- split(training_models, by="phi")
  dt <- do.call(rbind, map(training_models, function(dt_by_phi) {
    sorted.lambda <- sort(c(dt_by_phi$min.log.lambda, dt_by_phi$max.log.lambda))
    grid <- c(sorted.lambda[sorted.lambda[2:length(sorted.lambda)]!=sorted.lambda[1:(length(sorted.lambda)-1)]],sorted.lambda[length(sorted.lambda)])
    errors <- vector(mode="integer", length= length(grid)-1)
    walk(1:nrow(dt_by_phi), function(i){
      a <- which(grid>=dt_by_phi$min.log.lambda[[i]] & grid<dt_by_phi$max.log.lambda[[i]])
      errors[a] <<- errors[a] + dt_by_phi$total.errors[[i]]
    })
    i.min <- which.min(errors)
    predicted.lambda <- sum(grid[i.min:(i.min+1)])/2
    data.table(phi = dt_by_phi$phi[[1]], errors = min(errors))
  }))
  predicted.phi <- dt$phi[which.min(dt$errors)]
  # linear lambda
  setkey(all.targets, chunk.name)
  setkey(all.features, chunk.name)
  training_targets <- all.targets[training_chunks]
  training_features <- all.features[training_chunks]
  validation_targets <- all.targets[validation_chunks]
  validation_features <- all.features[validation_chunks]
  training_targets <- training_targets[phi==predicted.phi,]
  validation_targets <- validation_targets[phi==predicted.phi,]
  training_targets <- training_targets[!is.na(min.log.lambda),]
  validation_targets <- validation_targets[!is.na(min.log.lambda),]
  training_targets <- training_targets[is.finite(min.log.lambda) | is.finite(max.log.lambda),]
  setkey(training_targets, chunk.name, sample.id)
  setkey(training_features, chunk.name, sample.id)
  training_features_targets <- training_features[training_targets]
  setkey(validation_targets, chunk.name, sample.id)
  setkey(validation_features, chunk.name, sample.id)
  validation_features_targets <- validation_features[validation_targets]
  training_targets_mat <- as.matrix(training_features_targets[,list(min.log.lambda, max.log.lambda)])
  training_features_mat <- as.matrix(training_features_targets[,!c("chunk.name", "sample.id", "phi", "min.log.lambda", "max.log.lambda", "errors")])
  validation_features_mat <- as.matrix(validation_features_targets[,!c("chunk.name", "sample.id", "phi", "min.log.lambda", "max.log.lambda", "errors")])
  fit <- IntervalRegressionCV(
    target.mat = training_targets_mat,
    feature.mat = training_features_mat
  )
  predictions <- validation_features_targets[, list(chunk.name, sample.id, fit$predict(validation_features_mat))]
  names(predictions) <- c("chunk.name", "sample.id", "prediction")
  setkey(predictions, chunk.name, sample.id)
  setkey(all.modelSelection, chunk.name, sample.id)
  predictions <- merge(predictions, all.modelSelection[phi == predicted.phi])
  errors <- predictions[prediction>min.log.lambda & prediction<max.log.lambda,]
  acc <- (1 - sum(errors$total.errors) / sum(errors$labels)) * 100
  list(
    acc = data.table(fold.i = i, accuracy = acc),
    predicted_models = errors
  )
}

methods <- list(
  "linear_lambda" = linear_lambda,
  "constant_lambda" = constant_lambda,
  "constant_lambda_constant_phi" = constant_lambda_constant_phi,
  "linear_lambda_constant_phi" = linear_lambda_constant_phi
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
  res_cv_list <- imap(dp.peaks.sets, function(dataset, dataset.i){
    res_cv_dataset_list <- future_imap(dataset, function(fold, fold.i) {
      validation_chunks <- unlist(fold)
      training_chunks <- unlist(dataset[-fold.i])
      methods[[m]](validation_chunks, training_chunks, fold.i)
    }, .options = furrr_options(seed=NULL))
    res_cv_dataset <- rbindlist(map(res_cv_dataset_list, ~.x$acc))
    res_cv_dataset[,set.name := dataset.i]
    predicted_models <- rbindlist(map(res_cv_dataset_list, ~.x$predicted_models))
    predicted_models[,set.name := dataset.i]
    list(
      acc = res_cv_dataset,
      predicted_models = predicted_models
    )
  })
  res_cv <- rbindlist(map(res_cv_list, ~.x$acc))
  res_cv[,model := model_name]
  predicted_models <- rbindlist(map(res_cv_list, ~.x$predicted_models))
  predicted_models[,model := model_name]
  cv <- length(dp.peaks.sets[[1]])
  save(res_cv, file = paste0(path_to_output_data, "/res_cv_",model_name,"_",cv,".RData"))
  save(predicted_models, file = paste0(path_to_output_data, "/predicted_models_",model_name,"_",cv,".RData"))
}

params <- tibble::tribble(
  ~ algo,                              ~ targets,                                        ~ models,                                                ~ lambda,     ~ phi,        
  "PDPA_poisson",                      "all_targets_PDPA_poisson",                       "all_modelSelection_PDPA_poisson",                       "linear",     NA,           
  "PDPA_poisson_largest_peak_rule",    "all_targets_PDPA_poisson_largest_peak_rule",     "all_modelSelection_PDPA_poisson_largest_peak_rule",     "linear",     NA,           
  "PDPA_poisson_smallest_peak_rule",   "all_targets_PDPA_poisson_smallest_peak_rule",    "all_modelSelection_PDPA_poisson_smallest_peak_rule",    "linear",     NA,           
  "PDPA_poisson",                      "all_targets_PDPA_poisson",                       "all_modelSelection_PDPA_poisson",                       "constant",   NA,           
  "PDPA_poisson_largest_peak_rule",    "all_targets_PDPA_poisson_largest_peak_rule",     "all_modelSelection_PDPA_poisson_largest_peak_rule",     "constant",   NA,           
  "PDPA_poisson_smallest_peak_rule",   "all_targets_PDPA_poisson_smallest_peak_rule",    "all_modelSelection_PDPA_poisson_smallest_peak_rule",    "constant",   NA,           
  "PDPA_gaussian",                     "all_targets_PDPA",                               "all_modelSelection_PDPA",                               "linear",     NA,           
  "PDPA_gaussian_largest_peak_rule",   "all_targets_PDPA_largest_peak_rule",             "all_modelSelection_PDPA_largest_peak_rule",             "linear",     NA,           
  "PDPA_gaussian_smallest_peak_rule",  "all_targets_PDPA_smallest_peak_rule",            "all_modelSelection_PDPA_smallest_peak_rule",            "linear",     NA,           
  "PDPA_gaussian",                     "all_targets_PDPA",                               "all_modelSelection_PDPA",                               "constant",   NA,           
  "PDPA_gaussian_largest_peak_rule",   "all_targets_PDPA_largest_peak_rule",             "all_modelSelection_PDPA_largest_peak_rule",             "constant",   NA,           
  "PDPA_gaussian_smallest_peak_rule",  "all_targets_PDPA_smallest_peak_rule",            "all_modelSelection_PDPA_smallest_peak_rule",            "constant",   NA,           
  "PDPA_negbin",                       "all_targets_PDPA_negbin",                        "all_modelSelection_PDPA_negbin",                        "linear",    "constant",   
  "PDPA_negbin_largest_peak_rule",     "all_targets_PDPA_negbin_largest_peak_rule",      "all_modelSelection_PDPA_negbin_largest_peak_rule",      "linear",    "constant",   
  "PDPA_negbin_smallest_peak_rule",    "all_targets_PDPA_negbin_smallest_peak_rule",     "all_modelSelection_PDPA_negbin_smallest_peak_rule",     "linear",    "constant",   
  "PDPA_negbin",                       "all_targets_PDPA_negbin",                        "all_modelSelection_PDPA_negbin",                        "constant",  "constant",   
  "PDPA_negbin_largest_peak_rule",     "all_targets_PDPA_negbin_largest_peak_rule",      "all_modelSelection_PDPA_negbin_largest_peak_rule",      "constant",  "constant",   
  "PDPA_negbin_smallest_peak_rule",    "all_targets_PDPA_negbin_smallest_peak_rule",     "all_modelSelection_PDPA_negbin_smallest_peak_rule",     "constant",  "constant",   
  "updown_poisson",                    "all_targets_updown_poisson",                     "all_modelSelection_updown_poisson",                     "constant",  NA,           
  "updown_poisson",                    "all_targets_updown_poisson",                     "all_modelSelection_updown_poisson",                     "linear",    NA,           
  "updown_gaussian",                   "all_targets_updown",                             "all_modelSelection_updown",                             "constant",  NA,           
  "updown_gaussian",                   "all_targets_updown",                             "all_modelSelection_updown",                             "linear",    NA,           
  "updown_negbin",                     "all_targets_updown_negbin",                      "all_modelSelection_updown_negbin",                      "constant",  "constant",   
  "updown_negbin",                     "all_targets_updown_negbin",                      "all_modelSelection_updown_negbin",                      "linear",    "constant",   
)

pwalk(params, CV)
