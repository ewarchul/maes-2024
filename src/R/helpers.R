do_evaluate <- function(fn, population, scale = 1, is_vectorized = FALSE) {
  if (is_vectorized) {
    fn(population, ...) * scale
  } else {
    apply(population, 2, function(x) fn(x, ...) * scale)
  }
}

quadratic_penalty <- function(population, repaired_population) {
  penalties <- 1 + colSums((population - repaired_population)^2)
  penalties[!is.finite(penalties)] <- .Machine$double.xmax / 2
  penalties
}

projection_repair <- function(population, lower_bound, upper_bound) {
  ifelse(
    population > lower_bound,
    ifelse(population < upper_bound, population, upper_bound),
    lower_bound
  )
}

escape_flatland_ <- function(sigma, fitness, parameters) {
  magic_index <- min(
    1 + floor(parameters$lambda / 2),
    2 + ceiling(parameters$lambda / 4)
  )
  if (fitness[1] == fitness[magic_index]) {
    sigma * exp(0.2 + parameters$cs / parameters$damps)
  }
  sigma
}

calculate_chi <- function(dim) {
  sqrt(dim) * (1 - 1 / (4 * dim) + 1 / (21 * dim^2))
}
