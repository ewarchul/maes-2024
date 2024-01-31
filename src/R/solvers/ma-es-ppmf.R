library(magrittr)
ma_es_ppmf <- function(par, fn, ..., lower, upper, control = list()) {
  label <- "ma_es_ppmf"
  norm <- function(x) {
    drop(sqrt(crossprod(x)))
  }

  controlParam <- function(name, default) {
    v <- control[[name]]
    if (is.null(v)) {
      return(default)
    } else {
      return(v)
    }
  }

  xmean <- par
  N <- length(xmean)
  if (missing(lower)) {
    lower <- rep(-Inf, N)
  } else if (length(lower) == 1) {
    lower <- rep(lower, N)
  }

  if (missing(upper)) {
    upper <- rep(Inf, N)
  } else if (length(upper) == 1) {
    upper <- rep(upper, N)
  }

  trace <- controlParam("trace", FALSE)
  fnscale <- controlParam("fnscale", 1)
  stopfitness <- controlParam("stopfitness", 10^(-10))
  budget      <- controlParam("budget", 10000*N ) 

  sigma <- controlParam("sigma", 0.5)
  sc_tolx <- controlParam("stop.tolx", 1e-12 * sigma)
  keep.best <- controlParam("keep.best", TRUE)
  vectorized <- controlParam("vectorized", FALSE)

  log.all <- controlParam("diag", TRUE)
  log.sigma <- controlParam("diag.sigma", log.all)
  log.eigen <- controlParam("diag.eigen", FALSE)
  log.value <- controlParam("diag.value", log.all)
  log.pop <- controlParam("diag.pop", log.all)

  lambda <- controlParam("lambda", 4 + floor(3 * log(N)))
  mu <- controlParam("mu", floor(lambda / 2))


  maxiter <- controlParam("maxit", round(budget/lambda))


  weights <- controlParam("weights", log(mu + 1) - log(1:mu))
  weights <- weights / sum(weights)

  mueff <- controlParam("mueff", 1 / sum(weights^2))

  cs <- controlParam("cs", (mueff + 2) / (N + mueff + 5))
  alpha_cov <- 2
  c1 <- alpha_cov / ((N + 1.3)^2 + mueff)
  cw <- min(
    1 - c1,
    alpha_cov * ((mueff + mueff^(-1) - 2)/((N + 2)^2 + .5*alpha_cov*mueff))
  )

  damps <- 2
  p_target = 0.1

  stopifnot(length(upper) == N)
  stopifnot(length(lower) == N)
  stopifnot(all(lower < upper))
  stopifnot(length(sigma) == 1)

  best.fit <- Inf
  best.par <- NULL

  if (log.sigma) {
    sigma.log <- numeric(maxiter)
  }
  if (log.eigen) {
    eigen.log <- matrix(0, nrow = maxiter, ncol = N)
  }
  if (log.value) {
    value.log <- matrix(0, nrow = maxiter, ncol = mu)
  }
  if (log.pop) {
    pop.log <- array(0, c(N, mu, maxiter))
  }

  bestVal.log <- matrix(0, nrow = 0, ncol = 1)

  pc <- rep(0.0, N)
  ps <- rep(0.0, N)

  I <- diag(N)
  M <- I

  chiN <- sqrt(N) * (1 - 1 / (4 * N) + 1 / (21 * N^2))

  iter <- 0L 
  counteval <- 0L 
  cviol <- 0L 
  msg <- NULL 
  nm <- names(par) 

  arx <- matrix(0.0, nrow = N, ncol = lambda)
  arfitness <- numeric(lambda)

  eval_mean = Inf
  eval_meanOld = Inf
  while (counteval < budget) {
    iter <- iter + 1L

    if (!keep.best) {
      best.fit <- Inf
      best.par <- NULL
    }
    if (log.sigma) {
      sigma.log[iter] <- sigma
    }

    arz <- matrix(rnorm(N * lambda), ncol = lambda)

    arx <- xmean + sigma * (M %*% arz)

    vx <- ifelse(arx > lower, ifelse(arx < upper, arx, upper), lower)
    if (!is.null(nm)) {
      rownames(vx) <- nm
    }
    pen <- 1 + colSums((arx - vx)^2)
    pen[!is.finite(pen)] <- .Machine$double.xmax / 2
    cviol <- cviol + sum(pen > 1)

    if (vectorized) {
      y <- fn(vx, ...) * fnscale
    } else {
      y <- apply(vx, 2, function(x) fn(x, ...) * fnscale)
    }
    counteval <- counteval + lambda

    
    arfitness <- y * pen
    valid <- pen <= 1
    if (any(valid)) {
      wb <- which.min(y[valid])
      if (y[valid][wb] < best.fit) {
        best.fit <- y[valid][wb]
        best.par <- arx[, valid, drop = FALSE][, wb]
      }
    }

    arindex <- order(arfitness)
    arfitness <- arfitness[arindex]

    bestVal.log <- rbind(bestVal.log, min(suppressWarnings(min(bestVal.log)), min(arfitness)))

    aripop <- arindex[1:mu]
    selx <- arx[, aripop]
    xmean <- drop(selx %*% weights)
    selz <- arz[, aripop]
    zmean <- drop(selz %*% weights)

    if (log.pop) pop.log[, , iter] <- selx
    if (log.value) value.log[iter, ] <- arfitness[aripop]

    ps <- (1 - cs) * ps +
      sqrt(cs * (2 - cs) * mueff) * zmean

    M <- M * (I + 0.5 * c1 * (ps %o% ps - I) + 0.5 * cw * (zmean %o% zmean - I))

    eval_meanOld <- eval_mean
    mean_point <- apply(vx, 1, mean) %>% t() %>% t()
    eval_mean <- apply(mean_point, 2, function(x) fn(x))
    counteval <- counteval + 1

     p_succ = 
      length(which(arfitness < eval_meanOld))/lambda
    sigma = 
      sigma * exp(damps * (p_succ - p_target) / (1 - p_target))


    if (arfitness[1] <= stopfitness * fnscale) {
      msg <- "Stop fitness reached."
      break
    }

    if (arfitness[1] == arfitness[min(1 + floor(lambda / 2), 2 + ceiling(lambda / 4))]) {
      sigma <- sigma * exp(0.2 + cs / damps)
      if (trace) {
        message("Flat fitness function. Increasing sigma.")
      }
    }
    if (trace) {
      message(sprintf(
        "Iteration %i of %i: current fitness %f",
        iter, maxiter, arfitness[1] * fnscale
      ))
    }
  }
  cnt <- c(`function` = as.integer(counteval), gradient = NA)

  log <- list()
  if (log.value) log$value <- value.log[1:iter, ]
  if (log.sigma) log$sigma <- sigma.log[1:iter]
  if (log.eigen) log$eigen <- eigen.log[1:iter, ]
  if (log.pop) log$pop <- pop.log[, , 1:iter]

  log$bestVal <- bestVal.log

  names(best.fit) <- NULL
  res <- list(
    par = best.par,
    value = best.fit / fnscale,
    label = label,
    n.evals = cnt,
    resets = 0,
    convergence = ifelse(iter >= maxiter, 1L, 0L),
    message = msg,
    constr.violations = cviol,
    diagnostic = log
  )
  class(res) <- "ma_es_csa.result"
  return(res)
}
