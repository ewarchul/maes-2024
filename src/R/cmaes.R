library(R6)
library(logger)

source("parameters.R")
source("solutions.R")
source("helpers.R")

cma_es_solver <- R6::R6Class("cma_es_solver",
  public = list(
    initialize = function(parameters,
                           sigma_adaptor,
                           cov_mat_adaptor,
                           penalty_function = quadratic_penalty,
                           repair_population = projection_repair) {
      self$params <- params
      self$sigma_adaptor <- sigma_adaptor
      self$cov_mat_adaptor <- cov_mat_adaptor
      self$penalty_function <- penalty_function
      self$repair_population <- projection_repair
      self$dim <- self$params$dimension
      self$lambda <- self$params$lambda
      self$fitness <- numeric(self$lambda)
      self$sols <- CmaesSolutions$new(dim, self$params$init_sigma)
    },

    # [number] -> ([number] -> number) -> CmaesSolutions
    solve = function(x0, fn_eval) {
      while (!self$should_terminate()) {
        self$sols$iter <- self$sols$iter + 1
        population <- ask()
        population_repaired <- self$repair_population(population, parameters)

        penalties <- self$penalty_function(population, population_repaired)
        self$sols$cviol <- self$sols$cviol + sum(penalties > 1)

        fitness <- self$evaluate(population_repaired, fn_eval, penalties)
        self$sols$counteval <- self$sols$counteval + self$lambda

        tell()
        logger::log_info(
          "Iteration {self$sols$iter} of {self$params$maxiter}. Current fitness: {self$fitness[1]}"
        )
      }

      self$sols
    }
  ),
  private = list(
    params = NULL,
    sols = NULL,
    dim = NULL,
    lambda = NULL,
    fitness = NULL,
    sigma_adaptor = NULL,
    cov_mat_adaptor = NULL,
    penalty_function = NULL,
    repair_population = NULL,

    # Matrix[dim, lambda] -> ([number] -> number) -> [number] -> [number]
    evaluate = function(population, fn_eval, penalties) {
      eval_values <-
        do_evaluate(
          eval_fn,
          population,
          self$params$fnscale,
          self$params$vectorized
        )
      fitness <- eval_values * penalties
      self$replace_best_sofar(population, fitness, penalties)

      fitness
    },

    # [number] -> [number] -> ()
    replace_best_sofar = function(population, fitness, penalties) {
      is_valid_arr <- penalties <= 1
      if (any(is_valid_arr)) {
        current_iter_best_index <- which.min(fitness[is_valid_arr])
        current_iter_best_fitness <-
          fitness[is_valid_arr][current_iter_best_index]
        if (current_iter_best_fitness < self$sols$best.fit) {
          self$sols$best.fit <-
            fitness[is_valid_arr][current_iter_best_index]
          self$sols$best.par <-
            population[, is_valid_arr, drop = FALSE][, current_iter_best_index]
        }
      }
    },

    # () -> bool
    should_terminate = function() {
      is_budget_exhausted <- self$sols$counteval > self$parameters$budget
      is_stopfitness <- self$fitness[1] <= self$fitness[1] * self$params$fnscale

      is_budget_exhausted || is_stopfitness
    },

    # () -> Matrix[dim, lambda]
    ask = function() {
      noise_matrix <-
        matrix(rnorm(self$dim * self$params$lambda), ncol = self$params$lambda)

      self$sols$xmean + self$sols$sigma * (self$sols$cov_matrix %*% noise_matrix)
    },

    # () -> ()
    tell = function() {
      population_order_idx <- order(self$fitness)
      self$fitness <- self$fitness[population_order_idx]

      selected_population <- population[, population_order_idx[1:mu]]
      self$sols$xmean <- drop(selected_population %*% self$params$weights)

      selected_noise <- noise_matrix[, population_order_idx[1:mu]]
      self$sols$zmean <- drop(selected_noise %*% self$params$weights)

      self$sols$sigma_evol_path <- (1 - self$params$cs) * self$sols$sigma_evol_path +
        sqrt(self$params$cs * (2 - self$params$cs) * self$params$mueff) * self$sols$zmean

      self$sols$matrix_cov <- cov_updater(self$params, self$sols)
      self$sols$sigma <- sigma_updater(self$params, self$sols)
      self$sols$sigma <- escape_flatland_(self$sols$sigma, fitness, self$params)
    }
  )
)
