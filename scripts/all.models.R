library(gfpop) # local version of gfpop which takes phi as parameter
library(dequer)
library(data.table)
library(purrr)
library(furrr)

graphs <- list(
  PDPA = function(lambda) {
    graph(
      Edge(state1 = "S1", state2 = "S1", type = "std", penalty = lambda),
      Edge(state1 = "S1", state2 = "S1", type = "null", penalty = 0)
    )
  },
  updown = function(lambda) {
    graph(
      Edge(state1 = "BG", state2 = "UP", type = "up", penalty = lambda),
      Edge(state1 = "UP", state2 = "BG", type = "down", penalty = lambda),
      Edge(state1 = "BG", state2 = "BG", type = "null", penalty = 0),
      Edge(state1 = "UP", state2 = "UP", type = "null", penalty = 0),
      StartEnd(start = "BG", end = "BG")
    )
  },
  updown_negbin = function(lambda) {
    graph(
      Edge(state1 = "BG", state2 = "UP", type = "down", penalty = lambda),
      Edge(state1 = "UP", state2 = "BG", type = "up", penalty = lambda),
      Edge(state1 = "BG", state2 = "BG", type = "null", penalty = 0),
      Edge(state1 = "UP", state2 = "UP", type = "null", penalty = 0),
      StartEnd(start = "UP", end = "UP")
    ) 
    # So that direction of changes are consistents with changes in the
    # mean, parametrization should be 1-theta instead of theta. While waiting for 
    # an update of the gfpop code, the trick is to reverse all edges of the up-down
    # graph.
  }
)

losses <- list(
  mean = function(counts, theta, phi, weights) {
    sum((counts-theta)^2*weights)
  },
  poisson = function(counts, theta, phi, weights) {
    theta[theta == 0] <- 1
    sum((theta-counts*log(theta))*weights)
  },
  negbin = function(counts, theta, phi, weights) {
    sum((-log(theta)-counts/phi*log(1-theta))*weights)
  }
)

transformations <- list(
  raw = function(x) x,
  anscombe_poisson = function(x) sqrt(x+3/8),
  logp1 = function(x) log(x+1)
)

peaks <- function(theta) {
  # 1,2,4,4,2 
  # 1,2,0,-2
  # 1,1,0,-1
  # 1,(1,-1)=peak
  if (length(theta) > 2) {
    sign_theta <- sign(diff(theta))
    sign_theta <- sign_theta[sign_theta!=0]
    if (length(sign_theta)>1){
      sum(sign_theta[1:(length(sign_theta)-1)] == 1 & 
        sign_theta[2:length(sign_theta)] == -1)
    } else {
      0
    }
  } else {
    0
  }
}

solver <- function(counts, weights, lambda, loss_type, graph, phi = 0) {
  fit <- gfpop(mygraph = graphs[[graph]](lambda), data = counts, weights = weights, type = loss_type, phi=phi)
  theta <- rep(fit$parameters, times = diff(c(0, fit$changepoints)))
  loss <- losses[[loss_type]](counts, theta, phi, weights)
  if (loss_type == "negbin") {
    fit$parameters <- 1-fit$parameters
  }
  peaks <- peaks(fit$parameters)
  list(
    peaks = peaks,
    changepoints = length(fit$changepoints) - 1,
    lambda = lambda,
    phi = phi,
    loss = loss,
    fit = list(fit)
  )
}

sequential_search <- function(counts, weights, target, loss_type, graph, phi = 0) {
  upper_bound_model <- solver(counts, weights, 0, loss_type, graph, phi)
  lower_bound_model <- solver(counts, weights, Inf, loss_type, graph, phi)
  while (target != lower_bound_model$peaks & target != upper_bound_model$peaks) {
    new_candidate_lambda <- (upper_bound_model$loss - lower_bound_model$loss) / (lower_bound_model$changepoints - upper_bound_model$changepoints)
    if (new_candidate_lambda > lower_bound_model$lambda | new_candidate_lambda < upper_bound_model$lambda) {
      warning("lambda out of bounds")
      new_candidate_lambda <- (lower_bound_model$lambda + upper_bound_model$lambda) / 2
    }
    new_candidate_model <- solver(counts, weights, new_candidate_lambda, loss_type, graph, phi)
    if (new_candidate_model$changepoints == lower_bound_model$changepoints) {
      if (new_candidate_model$loss < lower_bound_model$loss) { # gfpop misses sometimes the optimal model with k changepoints. It really improves the results with the negbin loss (the roots computation is definitely not straightforward).
        lower_bound_model <- new_candidate_model
      } else {
        break
      }
    } else if (new_candidate_model$changepoints == upper_bound_model$changepoints){ 
      if (new_candidate_model$loss < upper_bound_model$loss) { # gfpop misses sometimes the optimal model with k changepoints. It really improves the results with the negbin loss (the roots computation is definitely not straightforward).
        upper_bound_model <- new_candidate_model
      } else {
        break
      }
    } else if (new_candidate_model$peaks == target) {
      lower_bound_model <- new_candidate_model
    } else if (new_candidate_model$peaks > target) {
      upper_bound_model <- new_candidate_model
    } else {
      lower_bound_model <- new_candidate_model
    }
  }
  lower_bound_model
}

CROPS <- function(counts, weights, lower_bound_lambda, upper_bound_lambda, loss_type, graph, phi = 0) {
  list_of_models <- list()
  last_index <- 1
  list_of_models[[last_index]] <- solver(counts, weights, lower_bound_lambda, loss_type, graph, phi)
  last_index <- last_index + 1
  list_of_models[[last_index]] <- solver(counts, weights, upper_bound_lambda, loss_type, graph, phi)
  q <- queue()
  pushback(q, c(1, 2))
  while (length(q) != 0) {
    d <- pop(q)
    if (list_of_models[[d[[1]]]]$changepoints > list_of_models[[d[[2]]]]$changepoints + 1) {
      new_candidate_lambda <- (list_of_models[[d[[2]]]]$loss - list_of_models[[d[[1]]]]$loss) / (list_of_models[[d[[1]]]]$changepoints - list_of_models[[d[[2]]]]$changepoints)
      if (new_candidate_lambda < list_of_models[[d[[1]]]]$lambda | new_candidate_lambda > list_of_models[[d[[2]]]]$lambda) {
        new_candidate_lambda <- (list_of_models[[d[[1]]]]$lambda + list_of_models[[d[[2]]]]$lambda) / 2
      }
      new_candidate_model <- solver(counts, weights, new_candidate_lambda, loss_type, graph, phi)
      if (new_candidate_model$changepoints < list_of_models[[d[[1]]]]$changepoints & new_candidate_model$changepoints > list_of_models[[d[[2]]]]$changepoints) {
        last_index <- last_index + 1
        list_of_models[[last_index]] <- new_candidate_model
        pushback(q, c(d[[1]], last_index))
        pushback(q, c(last_index, d[[2]]))
      } else if (new_candidate_model$changepoints == list_of_models[[d[[1]]]]$changepoints & new_candidate_model$loss < list_of_models[[d[[1]]]]$loss) { # gfpop misses sometimes the optimal model with k changepoints. Really improves the results with the negbin loss (the roots computation is definitely not straightforward). 
        list_of_models[[d[[1]]]] <- new_candidate_model
        pushback(q, c(d[[1]], d[[2]]))
      } else if (new_candidate_model$changepoints == list_of_models[[d[[2]]]]$changepoints & new_candidate_model$loss < list_of_models[[d[[2]]]]$loss) { # gfpop misses sometimes the optimal model with k changepoints. Really improves the results with the negbin loss (the roots computation is definitely not straightforward).
        list_of_models[[d[[2]]]] <- new_candidate_model
        pushback(q, c(d[[1]], d[[2]]))
      }
    }
  }
  list_of_models
}

CROCS <- function(counts, weights, lower_bound_peak, upper_bound_peak, loss_type, graph, phi = 0, transformation="raw") {
  counts <- transformations[[transformation]](counts)
  if (lower_bound_peak > 0) {
    lower_bound_peak <- lower_bound_peak - 1
  }
  upper_bound_peak <- upper_bound_peak + 1
  lower_bound_lambda <- sequential_search(counts, weights, upper_bound_peak, loss_type, graph, phi)$lambda
  upper_bound_lambda <- sequential_search(counts, weights, lower_bound_peak, loss_type, graph, phi)$lambda
  CROPS(counts, weights, lower_bound_lambda, upper_bound_lambda, loss_type, graph, phi)
}

compute_models <- function(lower_bound_peak, upper_bound_peak, loss_type, graph, transformation, phi, directories, input, output) {
  future_walk(directories[[1]], function(directory) {
    load(file.path(directory, input))
    dt_counts <- setDT(counts)
    models <- map(phi, function(phi_p) {
      map(unique(dt_counts$sample.id)[1], function(sample_p) {
        dt_counts_c <- dt_counts[sample.id==sample_p,]
        counts <- dt_counts_c$coverage
        weights <- dt_counts_c$chromEnd - dt_counts_c$chromStart
        dt <- rbindlist(CROCS(
          counts = counts, 
          weights = weights, 
          lower_bound_peak = lower_bound_peak, 
          upper_bound_peak = upper_bound_peak, 
          loss_type = loss_type, 
          graph = graph, 
          transformation = transformation, 
          phi = phi_p
        ))
        dt[,sample:=sample_p]
        dt[order(changepoints),]
      })
    })
    dt.models <- rbindlist(map(models, rbindlist))
    save(dt.models, file = file.path(directory, output))
  })
}


d <- Sys.glob("data/chunks/H*/*")
m <- rbindlist(list(
  list(
    lower_bound_peak = 1, 
    upper_bound_peak = 9,
    loss_type = "poisson",
    graph = "PDPA",
    transformation = "raw",
    phi = list(0),
    directories = list(d),
    input = "counts.RData",
    output = "PDPA_poisson.RData"
  ),
  list(
    lower_bound_peak = 1, 
    upper_bound_peak = 9,
    loss_type = "mean",
    graph = "PDPA",
    transformation = "anscombe_poisson",
    phi = list(0),
    directories = list(d),
    input = "counts.RData",
    output = "PDPA.RData"
  ),
  list(
    lower_bound_peak = 1, 
    upper_bound_peak = 9,
    loss_type = "negbin",
    graph = "PDPA",
    transformation = "raw",
    phi = list(exp(seq(log(1),log(10000),length=16))),
    directories = list(d),
    input = "counts.RData",
    output = "PDPA_negbin.RData"
  ),
  list(
    lower_bound_peak = 1, 
    upper_bound_peak = 9,
    loss_type = "poisson",
    graph = "updown",
    transformation = "raw",
    phi = list(0),
    directories = list(d),
    input = "counts.RData",
    output = "updown_poisson.RData"
  ),
  list(
    lower_bound_peak = 1, 
    upper_bound_peak = 9,
    loss_type = "mean",
    graph = "updown",
    transformation = "anscombe_poisson",
    phi = list(0),
    directories = list(d),
    input = "counts.RData",
    output = "updown.RData"
  ),
  list(
    lower_bound_peak = 1, 
    upper_bound_peak = 9,
    loss_type = "negbin",
    graph = "updown",
    transformation = "raw",
    phi = list(exp(seq(log(1),log(10000),length=16))),
    directories = list(d),
    input = "counts.RData",
    output = "updown_negbin.RData"
  )
))

plan(multiprocess, workers = 20)
pwalk(m, compute_models)
plan(sequential)