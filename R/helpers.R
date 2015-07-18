
Rhat1 <- function(mat) {
  # mat must be a iteration X chain matrix
  m <- ncol(mat) # chains
  n <- nrow(mat) # iterations per chain
  b <- apply(mat,2,mean)
  B <- sum((b-mean(mat))^2)*n/(m-1)
  w <- apply(mat,2,var)
  W <- mean(w)
  s2hat <- (n-1)/n*W + B/n
  Vhat <- s2hat + B/m/n 
  covWB <- n /m * (cov(w,b^2)-2*mean(b)*cov(w,b))
  varV <- (n-1)^2 / n^2 * var(w)/m +
    (m+1)^2 / m^2 / n^2 * 2*B^2/(m-1) +
    2 * (m-1)*(n-1)/m/n^2 * covWB
  df <- 2 * Vhat^2 / varV
  R <- sqrt((df+3) * Vhat / (df+1) / W)
  return(R)
}

# Little function to calculate posterior variances from simulation
colVars <- function(a){
  n <- dim(a)[[1]]
  c <- dim(a)[[2]]
  return(.colMeans(((a - matrix(.colMeans(a, n, c),
                                nrow=n, ncol=c, 
                                byrow=TRUE)) ^ 2), 
                   n, c) * n / (n-1))
}

# The calculation of Waic!  Returns lppd, p_waic_1, p_waic_2, and waic, which we define
# as 2*(lppd - p_waic_2), as recommmended in BDA
waic <- function(stanfit){
  log_lik <- rstan::extract (stanfit, "log_lik")$log_lik
  dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
    c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2*elpd_waic
  loo_weights_raw <- 1/exp(log_lik-max(log_lik))
  loo_weights_normalized <- loo_weights_raw/
    matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
  loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/
                    colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
  total <- colSums(pointwise)
  se <- sqrt(n*colVars(pointwise))
  return(list(waic=total["waic"], elpd_waic=total["elpd_waic"],
              p_waic=total["p_waic"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
              pointwise=pointwise, total=total, se=se))
}
