#' Blurb of what this function does
#' 
#' @param p description
#' @param xi description
#' @param mu description
#' @param phi description
#' @param power description
#' @return what it returns

qtweedie_obs_dd <- function(p, xi = NULL, mu, phi, power = NULL) {
  # Modified from tweedie::qtweedie(); rewritten for speed & robustness
  if (is.null(power) & is.null(xi)) stop("Either xi or power must be given\n")
  xi.notation <- TRUE
  if (is.null(power)) {
    power <- xi
  } else xi.notation <- FALSE
  if (is.null(xi)) {
    xi.notation <- FALSE
    xi <- power
  }
  if (xi != power) {
    cat("Different values for xi and power given; the value of xi used.\n")
    power <- xi
  }
  
  index.par.long <- ifelse(xi.notation, "xi", "power")
  if (any(power < 1)) stop(paste(index.par.long, "must be greater than 1.\n"))
  if (any(phi <= 0)) stop("phi must be positive.")
  if (any(p < 0) || any(p > 1)) stop("p must be between zero and one.\n")
  if (any(mu <= 0)) stop("mu must be positive.\n")
  
  # align lengths
  len <- length(p)
  if (length(mu) > 1 && length(mu) != len) stop("mu must be scalar, or the same length as p.\n")
  if (length(mu) == 1) mu <- rep(mu, len)
  
  if (length(phi) > 1 && length(phi) != len) stop("phi must be scalar, or the same length as p.\n")
  if (length(phi) == 1) phi <- rep(phi, len)
  
  if (length(power) == 1 && length(mu) > 1) power <- rep(power, length(mu))
  if (length(power) > 1 && length(power) != len) {
    # power may be scalar or same length as p
    stop("power must be scalar, or the same length as p.\n")
  }
  
  # prepare output and trivial positions
  out <- rep(NA_real_, len)
  if (any(p == 1)) out[p == 1] <- Inf
  if (any(p == 0)) out[p == 0] <- 0
  
  interior_idx <- which(p > 0 & p < 1)
  if (length(interior_idx) == 0) return(out)
  
  mu.vec <- mu[interior_idx]
  phi.vec <- phi[interior_idx]
  p.vec <- p[interior_idx]
  power.vec <- if (length(power) == 1) rep(power, length(p.vec)) else power[interior_idx]
  
  # local wrapper for ptweedie difference (avoid repeated env lookups)
  pt2_local <- function(q, mu, phi, pwr, p.given) {
    ptweedie(q = q, mu = mu, phi = phi, power = pwr) - p.given
  }
  
  # Precompute vectorized quantities (fast)
  # qgamma and qpois are cheap vectorized calls
  qg_all <- qgamma(p.vec, rate = 1 / (phi.vec * mu.vec), shape = 1 / phi.vec)
  qp_all <- qpois(p.vec, lambda = mu.vec / phi.vec)
  
  # dtweedie(0) only needed for 1 < pwr < 2
  need_dt_idx <- which((power.vec > 1) & (power.vec < 2))
  step_all <- rep(NA_real_, length(p.vec))
  if (length(need_dt_idx) > 0) {
    step_all[need_dt_idx] <- tweedie::dtweedie(y = rep(0,length(need_dt_idx)), mu = mu.vec[need_dt_idx], xi=NULL,
                                               phi = phi.vec[need_dt_idx],
                                               power = power.vec[need_dt_idx][1])
  }
  
  # loop (tight) over interior points
  for (k in seq_along(p.vec)) {
    idx_global <- interior_idx[k]
    mu.1 <- mu.vec[k]
    phi.1 <- phi.vec[k]
    pwr  <- power.vec[k]
    prob <- p.vec[k]
    qg   <- qg_all[k]
    qp   <- qp_all[k]
    
    # trivial power cases
    if (pwr == 1) {
      out[idx_global] <- qp
      next
    }
    if (pwr == 2) {
      out[idx_global] <- qg
      next
    }
    
    # compute initial bracket 'start'
    if ((pwr > 1) & (pwr < 2)) {
      start <- (qg - qp) * pwr + (2 * qp - qg)
    } else { # pwr > 2
      start <- qg
    }
    
    # if small probability covered by point mass at zero (tweedie mass at zero)
    if ((pwr > 1) & (pwr < 2)) {
      step_val <- step_all[k]
      if (!is.na(step_val) && prob <= step_val) {
        out[idx_global] <- 0
        next
      }
    }
    
    # evaluate pt at start
    pt_start <- pt2_local(start, mu = mu.1, phi = phi.1, pwr = pwr, p.given = prob)
    if (pt_start == 0) {
      out[idx_global] <- start
      next
    }
    
    # find bracket start.2 such that pt(start) and pt(start.2) have opposite sign
    max_iter <- 120L
    success_bracket <- FALSE
    
    if (pt_start > 0) {
      # reduce start.2 by repeated halving until sign flips or max iterations
      start.2 <- start
      iter <- 0L
      while (iter < max_iter) {
        start.2 <- 0.5 * start.2
        pt_tmp <- pt2_local(start.2, mu = mu.1, phi = phi.1, pwr = pwr, p.given = prob)
        iter <- iter + 1L
        if (pt_tmp < 0) { success_bracket <- TRUE; break }
      }
      # if we didn't find a bracket by halving, attempt a small positive epsilon near zero
      if (!success_bracket) {
        start.2 <- start * (0.5 ^ max_iter)
        pt_tmp <- pt2_local(start.2, mu = mu.1, phi = phi.1, pwr = pwr, p.given = prob)
        if (pt_tmp < 0) success_bracket <- TRUE
      }
    } else { # pt_start < 0
      # expand start.2 upward until sign flips or max iterations
      start.2 <- start
      iter <- 0L
      while (iter < max_iter) {
        start.2 <- 1.5 * (start.2 + 2)
        pt_tmp <- pt2_local(start.2, mu = mu.1, phi = phi.1, pwr = pwr, p.given = prob)
        iter <- iter + 1L
        if (pt_tmp > 0) { success_bracket <- TRUE; break }
      }
      # final fallback: if still not bracketed, try a very large number
      if (!success_bracket) {
        start.2 <- start + 1e6
        pt_tmp <- pt2_local(start.2, mu = mu.1, phi = phi.1, pwr = pwr, p.given = prob)
        if (pt_tmp > 0) success_bracket <- TRUE
      }
    }
    
    # If bracket found, call uniroot; else fallback to nearest estimate
    if (success_bracket) {
      # Ensure endpoints have opposite signs to satisfy uniroot
      f.left <- pt2_local(start, mu = mu.1, phi = phi.1, pwr = pwr, p.given = prob)
      f.right <- pt2_local(start.2, mu = mu.1, phi = phi.1, pwr = pwr, p.given = prob)
      if (f.left * f.right < 0) {
        root_val <- tryCatch(
          uniroot(pt2_local, c(start, start.2), mu = mu.1, phi = phi.1, pwr = pwr, p.given = prob,
                  tol = 1e-12, maxiter = 1000L)$root,
          error = function(e) NA_real_
        )
        out[idx_global] <- root_val
      } else {
        # numerical pathology: fallback to start if uniroot bracket requirement not met
        out[idx_global] <- start
      }
    } else {
      # failed to bracket; return start as fallback (conservative)
      out[idx_global] <- start
    }
  } # end loop over interior pts
  
  out
}


