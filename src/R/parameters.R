extract_param_or <- function(param_list, param_name, default) {
  v <- param_list[[param_name]]
  if (is.null(v)) {
    return(default)
  } else {
    return(v)
  }
}

cma_es_parameters <- R6::R6Class(
  "cma_es_parameters",
  public = list(
    dimension = NULL,
    lower = NULL,
    upper = NULL,
    log_level = NULL,
    fnscale = NULL,
    init_sigma = NULL,
    lambda = NULL,
    mu = NULL,
    stopfitness = NULL,
    budget = NULL,
    maxiter = NULL,
    weights = NULL,
    mueff = NULL,
    sc_tolx = NULL,
    keep.best = NULL,
    vectorized = NULL,
    cs = NULL,
    alpha_cov = NULL,
    cw = NULL,
    damps = NULL,
    initialize = function(dim, extra_params) {
      self$dimension <- dim
      self$lower <- extract_param_or(extra_params, "lower", -100)
      self$upper <- extract_param_or(extra_params, "upper", 100)
      self$log_level <- extract_param_or(extra_params, "log_level", logger::INFO)
      self$fnscale <- extract_param_or(extra_params, "fnscale", 1)
      self$stopfitness <-
        extract_param_or(extra_params, "stopfitness", 10^(-10))
      self$budget <-
        extract_param_or(extra_params, "budget", 10000 * self$dimension)
      self$init_sigma <-
        extract_param_or(extra_params, "init_sigma", 0.5)
      self$sc_tolx <-
        extract_param_or(extra_params, "stop.tolx", 1e-12 * self$init_sigma)
      self$keep.best <- extract_param_or(extra_params, "keep.best", TRUE)
      self$vectorized <- extract_param_or(extra_params, "vectorized", FALSE)
      self$lambda <-
        extract_param_or(
          extra_params,
          "lambda",
          4 + floor(3 * log(self$dimension))
        )
      self$mu <- extract_param_or(extra_params, "mu", floor(self$lambda / 2))
      self$maxiter <-
        extract_param_or(
          extra_params,
          "maxiter",
          round(self$budget / self$lambda)
        )
      self$weights <-
        extract_param_or(
          extra_params,
          "weights",
          log(self$mu + 1) - log(1:self$mu)
        )
      self$weights <- self$weights / sum(self$weights)
      self$mueff <- extract_param_or(
        extra_params,
        "mueff",
        1 / sum(self$weights^2)
      )
      self$cs <-
        extract_param_or(
          extra_params,
          "cs",
          (self$mueff + 2) / (self$dimension + self$mueff + 5)
        )
      self$alpha_cov <-
        extract_param_or(extra_params, "alpha_cov", 2)
      self$c1 <-
        extract_param_or(
          extra_params,
          "c1",
          self$alpha_cov / ((self$dimension + 1.3)^2 + self$mueff)
        )
      self$cw <-
        extract_param_or(
          extra_params,
          "cw",
          min(
            1 - self$c1,
            self$alpha_cov * ((self$mueff + self$mueff^(-1) - 2) / ((self$dimension + 2)^2 + .5 * self$alpha_cov * self$mueff))
          )
        )
      self$damps <-
        extract_param_or(
          extra_params,
          "damps",
          1 + 2 * max(
            0,
            sqrt((self$mueff - 1) / (self$dimension + 1)) - 1
          ) + self$cs
        )
    }
  )
)
