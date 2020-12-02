library(gfpop)
library(dequer)
library(data.table)

list_of_graphs <- list(
  PDPA = function(lambda){
    graph(
      Edge(
        state1 = "S1", 
        state2 = "S1", 
        type = "std", 
        penalty = lambda
      ),
      Edge(
        state1 = "S1", 
        state2 = "S1", 
        type = "null", 
        penalty = 0
      ) 
    )
  },
  updown = function(lambda){ 
    graph(
      Edge(state1 = "BG", 
        state2 = "UP", 
        type = "up", 
        penalty = lambda
      ),
      Edge(
        state1 = "UP", 
        state2 = "BG", 
        type = "down", 
        penalty = lambda
      ),
      Edge(
        state1 = "BG", 
        state2 = "BG", 
        type = "null", 
        penalty = 0
      ),
      Edge(
        state1 = "UP", 
        state2 = "UP", 
        type = "null", 
        penalty = 0
      ),
      StartEnd(
        start = "BG", 
        end = "BG"
      )
    )
  }
)

number_of_peaks <- function(theta){
  if (length(theta)>2) {
    sum(diff(theta[-length(theta)])>0 & diff(theta[-1])<0)
  } else {
    0
  }
}

print_model <- function(model, other="") {
  print(paste0(
    model$loss_revised, 
    " ", 
    model$nb_peaks,
    " ",
    model$nb_changepoints,
    " ",
    model$lambda,
    " ",
    other
  ))
}

list_of_losses <- list(
  mean = function(counts, theta, disp, weights) {
    sum((counts - theta)^2*weights)
  },
  poisson = function(counts, theta, disp, weights) {
    sum((counts - theta)^2*weights)
  }
)

list_of_solvers <- list(
  PDPA = function(counts, weights, lambda, loss_type, disp=0)
  {
    res_seg <- gfpop(mygraph = list_of_graphs[["PDPA"]](lambda), data=counts, weights=weights, type=loss_type)
    nb_peaks <- number_of_peaks(res_seg$parameters)
    theta <- rep(res_seg$parameters, times=diff(c(0,res_seg$changepoints)))
    loss <- list_of_losses[[loss_type]](counts, theta, 0, weights)
    list(
      loss = res_seg$globalCost, 
      loss_revised = loss, 
      nb_peaks = nb_peaks, 
      nb_changepoints = length(res_seg$changepoints)-1,
      lambda = lambda)
  }
)

sequential_search <-  function(counts, weights, target, solver, loss_type, disp=0)
{
  upper_bound_model <- list_of_solvers[[solver]](counts, weights, 0, loss_type, disp)
  lower_bound_model <- list_of_solvers[[solver]](counts, weights, Inf, loss_type, disp)
  while (target != lower_bound_model$nb_peaks & target != upper_bound_model$nb_peaks) {
    new_candidate_lambda <- (upper_bound_model$loss_revised - lower_bound_model$loss_revised) / (lower_bound_model$nb_changepoints - upper_bound_model$nb_changepoints)
    if (new_candidate_lambda > lower_bound_model$lambda | new_candidate_lambda < upper_bound_model$lambda) {
      warning(paste0("lambda out of bounds <", target,">"))
      new_candidate_lambda <- (lower_bound_model$lambda + upper_bound_model$lambda) / 2
    }
    new_candidate_model <- list_of_solvers[[solver]](counts, weights, new_candidate_lambda, loss_type, disp)
    if (new_candidate_model$nb_changepoints == lower_bound_model$nb_changepoints | new_candidate_model$nb_changepoints == upper_bound_model$nb_changepoints){
      break
    } else if (new_candidate_model$nb_peaks == target) {
      lower_bound_model <- new_candidate_model
    } else if (new_candidate_model$nb_peaks > target) {
      upper_bound_model <- new_candidate_model
    } else {
      lower_bound_model <- new_candidate_model
    } 
  }
  lower_bound_model
}

CROPS <- function(counts, weights, lower_bound_lambda, upper_bound_lambda, solver, loss_type, disp=0) {
  list_of_models <- list()
  last_index <- 1
  list_of_models[[last_index]] <- list_of_solvers[[solver]](counts, weights, lower_bound_lambda, loss_type, disp)
  last_index <- last_index + 1
  list_of_models[[last_index]] <- list_of_solvers[[solver]](counts, weights, upper_bound_lambda, loss_type, disp)
  q <- queue()
  pushback(q,c(1,2))
  while (length(q)!=0) {
    d <- pop(q)
    if (list_of_models[[d[[1]]]]$nb_changepoints > list_of_models[[d[[2]]]]$nb_changepoints + 1) {
      new_candidate_lambda <- (list_of_models[[d[[2]]]]$loss_revised - list_of_models[[d[[1]]]]$loss_revised) / (list_of_models[[d[[1]]]]$nb_changepoints - list_of_models[[d[[2]]]]$nb_changepoints)
      if (new_candidate_lambda < list_of_models[[d[[1]]]]$lambda| new_candidate_lambda > list_of_models[[d[[2]]]]$lambda) {
        warning(paste0("lambda out of bounds"))
        new_candidate_lambda <- (list_of_models[[d[[1]]]]$lambda + list_of_models[[d[[2]]]]$lambda) / 2
      }
      new_candidate_model <- list_of_solvers[[solver]](counts, weights, new_candidate_lambda, loss_type, disp)
      if (new_candidate_model$nb_changepoints < list_of_models[[d[[1]]]]$nb_changepoints & new_candidate_model$nb_changepoints > list_of_models[[d[[2]]]]$nb_changepoints) {
        last_index <- last_index + 1 
        list_of_models[[last_index]] <- new_candidate_model
        pushback(q,c(d[[1]],last_index))
        pushback(q,c(last_index, d[[2]]))
      }
    }
  }
  list_of_models
}

CROCS <- function(counts, weights, lower_bound_peak, upper_bound_peak, solver, loss_type, disp=0) {
  if (lower_bound_peak>0) {
    lower_bound_peak <- lower_bound_peak - 1
  }
  upper_bound_peak <- upper_bound_peak + 1
  lower_bound_lambda <- sequential_search(counts, weights, upper_bound_peak, solver, loss_type, disp)$lambda
  upper_bound_lambda <- sequential_search(counts, weights, lower_bound_peak, solver, loss_type, disp)$lambda
  CROPS(counts, weights, lower_bound_lambda, upper_bound_lambda, solver, loss_type, disp)
}


library(furrr)
furrr_options(seed=FALSE)
plan(multiprocess, workers=20)

load("data/chunks/H3K4me3_TDH_immune/1/counts.RData")
dt.counts <- data.table(counts)
res <- future_map(unique(dt.counts$sample.id), function(x) {
  dt.counts_ <- dt.counts[sample.id == x,]
  counts <- dt.counts_$coverage
  weights <- dt.counts_$chromEnd - dt.counts_$chromStart
  CROCS(counts, weights, 1, 9, "PDPA", "mean")
}, .options = furrr_options(seed=2020))

names(res) <- unique(dt.counts$sample.id)