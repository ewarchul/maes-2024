cma_es_solutions <- R6::R6Class(
  "cma_es_solutions",
  public = list(
    initialize = function(dim, init_sigma = 1.0) {
      self$cov_evol_path <- rep(0.0, dim)
      self$sigma_evol_path <- rep(0.0, dim)
      self$cov_matrix <- diag(dim)
      self$sigma <- init_sigma
    },
    best_fitness = NULL,
    best_point = NULL,
    iter = 0,
    counteval = 0,
    cviol = 0,
    cov_evol_path = NULL,
    sigma_evol_path = NULL,
    cov_matrix = NULL,
    sigma = NULL,
    xmean = NULL,
    zmean = NULL
  )
)
