ma <- function(params, sols) {
  I <- diag(params$dim)
  sols$cov_mat *
    (I + 0.5 * params$c1 * (sols$sigma_evol_path %o% sols$sigma_evol_path - I) +
      0.5 * params$cw * (sols$zmean %o% sols$zmean - I))
}
