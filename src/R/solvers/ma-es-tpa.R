ma_es_tpa <- function(par, fn, ..., lower, upper, control = list()) {
  label <- "ma_es_tpa"
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
  stopfitness <- controlParam("stopfitness", -Inf)
  maxiter <- controlParam("maxit", 100 * N^2)
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

  alpha_tpa <- 0.7
  beta_tpa <- alpha_tpa^(-1)
  c_alpha_tpa <- 0.5


  damps <- controlParam(
    "damps",
    1 + 2 * max(0, sqrt((mueff - 1) / (N + 1)) - 1) + cs
  )

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

  pc <- rep(0.0, N)
  ps <- rep(0.0, N)

  I <- diag(N)
  M <- I


  iter <- 0L 
  counteval <- 0L 
  cviol <- 0L 
  msg <- NULL 
  nm <- names(par) 

  arx <- matrix(0.0, nrow = N, ncol = lambda)
  arfitness <- numeric(lambda)
  alpha_s <- 0L
  while (iter < maxiter) {
    iter <- iter + 1L

    if (!keep.best) {
      best.fit <- Inf
      best.par <- NULL
    }
    if (log.sigma) {
      sigma.log[iter] <- sigma
    }

    arz <- matrix(rnorm(N * lambda), ncol = lambda)
    ary <- M %*% arz
    arx <- xmean + sigma * ary

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

    aripop <- arindex[1:mu]
    selx <- arx[, aripop]
    xmean <- drop(selx %*% weights)
    selz <- arz[, aripop]
    zmean <- drop(selz %*% weights)

    sely <- ary[, aripop]
    ymean <- drop(sely %*% weights)

    if (log.pop) pop.log[, , iter] <- selx
    if (log.value) value.log[iter, ] <- arfitness[aripop]

    ps <- (1 - cs) * ps +
      sqrt(cs * (2 - cs) * mueff) * zmean

    M <- M * (I + 0.5 * c1 * (ps %o% ps - I) + 0.5 * cw * (zmean %o% zmean - I))

    f_alpha <- xmean + alpha_tpa * sigma * ymean
    f_beta <- xmean + beta_tpa * sigma * ymean

    f_alpha_eval <- fn(f_alpha)
    f_beta_eval <- fn(f_beta)
    counteval <- counteval + 2

    alpha_s <- (1 - c_alpha_tpa) * alpha_s + 
      c_alpha_tpa * log(alpha_tpa) * ifelse(f_alpha_eval < f_beta_eval, 1, 0) +
      c_alpha_tpa * log(beta_tpa) * ifelse(f_alpha_eval >= f_beta_eval, 1, 0)

    sigma <- sigma * exp(alpha_s)

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

  names(best.fit) <- NULL
  res <- list(
    par = best.par,
    value = best.fit / fnscale,
    label = label,
    counts = cnt,
    convergence = ifelse(iter >= maxiter, 1L, 0L),
    message = msg,
    constr.violations = cviol,
    diagnostic = log
  )
  class(res) <- "ma_es_csa.result"
  return(res)
}
