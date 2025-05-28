#' Estimate treatment effect and the corresponding variance estimation on the treatment arm using different covariate adjustment methods.
#'
#' @description
#' Implements a unified framework for comparing covariate adjustment method for completely randomized experiments under randomization-based framework.
#'
#' @param X The n by p covariates matrix.
#' @param Y Vector of n dimensional observed response.
#' @param A Vector of n dimensional treatment assignment.
#' @param H The n by n hat projection matrix corresponding to X.
#'
#' @return A list with two named vectors:
#' \describe{
#'   \item{point_est}{Point estimates for all estimators:
#'     \itemize{
#'       \item{\code{unadj}:} Unadjusted estimator
#'       \item{\code{db}:} Debiased estimator (Lu et al., 2023)
#'       \item{\code{adj2c}:} HOIF-inspired debiased estimator (Zhao et al., 2024), the same as \code{db}
#'       \item{\code{adj2}:} HOIF-motivated adjusted estimator (Zhao et al., 2024)
#'       \item{\code{adj3}:} Bias-free adjusted estimator based on \code{adj2}
#'       \item{\code{lin}:} Covariate-adjusted estimator (Lin, 2013)
#'       \item{\code{lin_db}:} Debiased estimator with population leverage scores (Lei, 2020)
#'     }}
#'   \item{var_est}{Variance estimates corresponding to each estimator:
#'     \itemize{
#'       \item{\code{unadj}:} Variance estimate for unadjusted estimator
#'       \item{\code{db}:} Variance estimate for debiased estimator (Lu et al., 2023)
#'       \item{\code{adj2c}:} Variance for \code{adj2c}, using formulas given in (Lu et al., 2023)
#'       \item{\code{adj2c_v2}:} Conservative variance for \code{adj2c} (Zhao et al., 2024)
#'       \item{\code{adj2}:} Variance for \code{adj2}, with formulas  motivated by (Lu et al., 2023)
#'       \item{\code{adj2_v2}:} Conservative variance for \code{adj2} (Zhao et al., 2024)
#'       \item{\code{adj3}:} Variance for \code{adj3},  with formulas  motivated by (Lu et al., 2023)
#'       \item{\code{adj3_v2}:} Conservative variance for \code{adj3} (Zhao et al., 2024)
#'       \item{\code{lin}:} HC3-type variance for Lin's (2013) estimator
#'       \item{\code{lin_db}:} HC3-type variance for Lei's (2020) estimator
#'     }}
#' }
#' @export
#'
#' @references
#' Winston Lin. Agnostic notes on regression adjustments to experimental data: Reexamining Freedman's critique. The Annals of Statistics, 7(1):295–318, 2013.
#'
#' Lihua Lei, Peng Ding. Regression adjustment in completely randomized experiments with a diverging number of covariates. Biometrika, 815–828, 2020.
#'
#' Sihui Zhao, Xinbo Wang, Lin Liu, & Xin Zhang. Covariate adjustment in randomized experiments motivated by higher-order influence functions. arXiv preprint. https://arxiv.org/abs/2411.08491.
#'
#' Xin Lu, Fan Yang, and Yuhao Wang. Debiased regression adjustment in completely randomized experiments with moderately high-dimensional covariates. arXiv preprint arXiv:2309.02073, 2023.
#'
#' @examples
#' set.seed(100)
#' n <- 500
#' p <- n * 0.3
#' beta <- runif(p, -1 / sqrt(p), 1 / sqrt(p))
#'
#' X <- mvtnorm::rmvt(n, sigma = diag(1, p), df = 3)
#' Y1 <- as.numeric(X %*% beta)
#' Y0 <- rep(0, n)
#'
#' pi1 <- 2/3
#' n1 <- ceiling(n * pi1)
#' ind <- sample(n, size = n1)
#' A <- rep(0, n)
#' A[ind] <- 1
#' Y <- Y1 * A + Y0 * (1 - A)
#'
#' Xc_svd <- svd(X)
#' H <- Xc_svd$u %*% t(Xc_svd$u)
#'
#' result_ls <- esti_mean_treat(X, Y, A, H)
#' point_est <- result_ls$point_est
#' var_est <- result_ls$var_est
#' print(paste0('True mean treat:', round(mean(Y1), digits = 3), '.'))
#' print('Absolute bias:')
#' print(abs(point_est - mean(Y1)))
#' print('Estimate variance:')
#' print(var_est)
#'
esti_mean_treat <- function(X, Y, A, H = NULL) {

  n <- nrow(X)
  p <- ncol(X)
  pi1 <- mean(A)
  pi0 <- 1 - pi1

  n1 <- sum(A)
  n0 <- n - n1

  Xc <- scale(X, scale = FALSE)
  if (is.null(H)) {
    Xc_svd <- svd(Xc)
    H <- Xc_svd$u %*% t(Xc_svd$u)
  }

  tau_unadj <- mean(A * Y) / mean(A)


  if(n1 > p){
    fit1 <- lm(Y~Xc,subset = (A==1))
    tau_hat_lin <- fit1$coefficients[1]
    e1 <- resid(fit1)
    lscores <- hat(Xc)
    lscores1 <- lscores[A==1]
    Delta1 <- sum(lscores1 * e1) / n1
    tau_hat_lin_db <- tau_hat_lin + n0 / n1 * Delta1

    point_vec_lin <- c(tau_hat_lin,tau_hat_lin_db)
  }else{
    point_vec_lin <- c(NA,NA)
  }
  names(point_vec_lin) <- c('lin','lin_db')


  YcA <- as.numeric(A * (Y - tau_unadj))
  tau_adj <- tau_unadj - mean((A / pi1 - 1) * H %*% YcA * n / n1)


  tau_db <- tau_adj + (1 - pi1) / pi1 * mean(A / pi1 * diag(H) * (Y - tau_unadj))


  ps_resid <- A / pi1 - 1
  or_resid <- A * Y / pi1
  or_resid_c <- YcA / pi1
  IF22 <- (t(ps_resid) %*% H %*% or_resid - sum(ps_resid * diag(H) * or_resid)) / n
  tau_adj2 <- tau_unadj - IF22


  IF22_c <- - (t(ps_resid) %*% H %*% or_resid_c - sum(ps_resid * diag(H) * or_resid_c)) / n
  tau_adj2c <- tau_unadj + IF22_c


  tau_adj3 <- tau_adj2 + (1-pi1)/pi1/n/(n-1)*sum(diag(H)*or_resid)

  point_est <- c(tau_unadj,tau_db,tau_adj2c,tau_adj2,tau_adj3,point_vec_lin)
  names(point_est) <- c('unadj','db','adj2c','adj2','adj3','lin','lin_db')





  var_hat_unadj <- pi0/pi1/n*sum((Y[A==1] - tau_unadj)^2)/(n1-1)
  if(n1 > p){
    mod1 <- lm(Y[A==1] ~ Xc[A==1, ])
    e1_hat <- as.numeric(resid(mod1))
    lscores1 <- hat(Xc[A==1, ])
    s1_HC3 <- sqrt(mean(e1_hat^2 / (1 - pmin(lscores1, 0.99))^2))
    var_hat_lin <- (1/n1 - 1/n)*s1_HC3^2

  }else{
    var_hat_lin <- NA
  }

  h <- diag(H)
  Q <- H^2
  diag(Q) <- diag(H) - diag(H)^2
  ones <- rep(1,n)
  M <- diag(n) - 1/n * ones %*% t(ones) - H + (diag(n) - 1/n * ones %*% t(ones)) %*% diag(h)
  B <- t(M) %*% M

  get_vec_s2 <- function(q,x,z,z0=1){
    q_sub <- q[z==z0]
    x_sub <- x[z==z0]
    s2 <- mean(q_sub*(x_sub-mean(x_sub))**2)
    return(s2)
  }
  I1_hat <- pi1*pi0*(get_vec_s2(q=pi1*pi0*diag(Q),x=Y/pi1^2,z=A,z0=1) +
                       get_vec_s2(q=diag(B)/pi1^2,x=Y,z=A,z0=1))

  get_mat_s2 <- function(Q,x,z,z0=1){
    Q_sub <- Q[z==z0,z==z0]
    x_sub <- x[z==z0] - mean(x[z==z0])

    s2 <- 1/pi1/n1*(t(x_sub) %*% Q_sub %*% x_sub - sum(x_sub * diag(Q_sub) * x_sub))
    return(s2)
  }
  I2_hat <- pi1*pi0*(get_mat_s2(Q = pi1*pi0*Q,x=Y/pi1**2,z=A,z0=1) +
                       get_mat_s2(Q=1/pi1^2*B,x=Y,z=A,z0=1))
  I3_hat <- 0
  I4_hat <- 0
  var_hat_db <- (I1_hat + I2_hat + I3_hat + I4_hat)/n


  Y1_hat <- A*Y/pi1
  tr_half_square_hat <- sum(h^2*A*Y^2/pi1) + sum((h*Y1_hat) %*% t(h*Y1_hat)) - sum((h*Y1_hat)^2)
  tr_hat <- sum(h^2*A*Y^2/pi1) + t(Y1_hat) %*% H^2 %*% Y1_hat - sum(Y1_hat*h^2*Y1_hat)
  Y1_sub <- Y[A==1]; H1_sub <- H[A==1,A==1]
  P <- diag(n) - 1/n*rep(1,n) %*% t(rep(1,n))
  M <- P - H + P%*%diag(diag(H))
  B <- t(M) %*% M
  var_hat_adj2c_wang <- (pi0/pi1)/n/n*(
    sum(diag(B)*A*(Y - tau_unadj)^2/pi1) + t(A*(Y-tau_unadj)/pi1) %*% B %*% (A*(Y-tau_unadj)/pi1) - sum((A*(Y-tau_unadj)/pi1) * diag(B) * (A*(Y-tau_unadj)/pi1))
  ) + (pi0/pi1)^2/n*(
    1/n1*sum(A*h*(1-h)*(Y-tau_unadj)^2) + 1/n1/pi1*(t(Y1_sub-tau_unadj) %*% H1_sub^2 %*% (Y1_sub - tau_unadj) - sum(diag(H1_sub)^2 * (Y1_sub - tau_unadj)^2))
  )


  resi <- Y - tau_unadj - (H - diag(diag(H))) %*% (A*(Y - tau_unadj)/pi1) - mean((1+diag(H))*A*(Y-tau_unadj)/pi1)
  var_hat_adj2c_v2 <- (pi0/pi1)/n/n*sum(A*resi^2/pi1) +
    (pi0/pi1)^2/n*(
      1/n1*sum(A*h*(1-h)*(Y-tau_unadj)^2) + 1/n1/pi1*(t(Y1_sub-tau_unadj) %*% H1_sub^2 %*% (Y1_sub - tau_unadj) - sum(diag(H1_sub)^2 * (Y1_sub - tau_unadj)^2))
    )


  P <- diag(n) - 1/n*rep(1,n) %*% t(rep(1,n))
  M <- P - H + P%*%diag(diag(H))
  B <- t(M) %*% M
  Y1_sub <- Y[A==1]; H1_sub <- H[A==1,A==1]
  var_hat_adj2_wang <- pi0/pi1/n/n*(
    sum(diag(B)*A*Y^2/pi1) + t(A*Y/pi1) %*% B %*% (A*Y/pi1) - sum((A*Y/pi1) * diag(B) * (A*Y/pi1))
  ) + (pi0/pi1)^2/n*(
    1/n1*sum(A*h*(1-h)*Y^2) + 1/n1/pi1*(t(Y1_sub) %*% H1_sub^2 %*% Y1_sub - sum(diag(H1_sub)^2 * Y1_sub^2))
  )


  resi <- Y - (H - diag(diag(H))) %*% (A*(Y)/pi1) - mean((1+diag(H))*A*(Y)/pi1)
  var_hat_adj2_v2 <- (pi0/pi1)/n/n*sum(A*resi^2/pi1) +
    (pi0/pi1)^2/n*(
      1/n1*sum(A*h*(1-h)*(Y)^2) + 1/n1/pi1*(t(Y1_sub) %*% H1_sub^2 %*% (Y1_sub) - sum(diag(H1_sub)^2 * (Y1_sub)^2))
    )


  var_hat_vec <- c(var_hat_unadj,var_hat_db,var_hat_adj2c_wang,var_hat_adj2c_v2,
                   var_hat_adj2_wang,var_hat_adj2_v2,
                   var_hat_adj2_wang,var_hat_adj2_v2,
                   var_hat_lin,var_hat_lin)
  names(var_hat_vec) <- c('unadj','db','adj2c','adj2c_v2',
                          'adj2','adj2_v2',
                          'adj3','adj3_v2',
                          'lin','lin_db')


  return(list(
    point_est = point_est,
    var_est = var_hat_vec
  ))

}



#' 	Covariate-Adjusted Treatment Effect Estimation under the Randomization-based Framework
#'
#' Implements the debiased estimators for average treatment effect (ATE) estimation
#' in completely randomized experiments with moderately high-dimensional covariates.
#'
#'
#' @param Y Numeric vector of length n containing observed responses.
#' @param X Numeric matrix (n x p) of covariates. Centering is required. May include intercept column.
#' @param A Binary vector of length n indicating treatment assignment (1 = treatment, 0 = control).
#'
#'
#' @return A list containing point estimates and variance estimates:
#' \describe{
#'    \item{tau_vec}{Point estimates for adj2 and adj2c:
#'        \itemize{
#'          \item{\code{adj2}:} Point estimation of the HOIF-inspired debiased estimator (Zhao et al., 2024).
#'          \item{\code{adj2c}:} Point estimation of the debiased estimator given by Lu et al. (2023), which is also the HOIF-inspired debiased estimator (Zhao et al., 2024)
#'        }
#'    }
#'    \item{var_vec}{Variance estimates for adj2 and adj2c:
#'      \itemize{
#'         \item{\code{adj2}:} Variance estimation of adj2, with formulas inspired by Lu et al. (2023).
#'         \item{\code{adj2c}:} Variance estimation of adj2c, with formulas given by Lu et al. (2023).
#'         \item{\code{adj2_v2}:} Variance estimation of adj2, with formulas given in Zhao et al. (2024), which is more conservative  compared to \code{adj2}.
#'         \item{\code{adj2c_v2}:} Variance estimation of adj2c, with formulas given in Zhao et al. (2024), which is more conservative  compared to \code{adj2c}.
#'      }
#'    }
#' }
#'
#' @export
#'
#' @references
#' Xin Lu, Fan Yang, and Yuhao Wang. Debiased regression adjustment in completely randomized experiments with moderately high-dimensional covariates. arXiv preprint arXiv:2309.02073, 2023.
#'
#' Sihui Zhao, Xinbo Wang, Lin Liu, & Xin Zhang. Covariate adjustment in randomized experiments motivated by higher-order influence functions. arXiv preprint. https://arxiv.org/abs/2411.08491.
#'
#' @examples
#' set.seed(100)
#' n <- 500
#' p <- n * 0.3
#' beta <- runif(p, -1 / sqrt(p), 1 / sqrt(p))
#'
#' X <- mvtnorm::rmvt(n, sigma = diag(1, p), df = 3)
#' Y1 <- as.numeric(X %*% beta)
#' Y0 <- as.numeric(X %*% beta - 1)
#'
#' pi1 <- 2/3
#' n1 <- ceiling(n * pi1)
#' ind <- sample(n, size = n1)
#' A <- rep(0, n)
#' A[ind] <- 1
#' Y <- Y1 * A + Y0 * (1 - A)
#'
#' Xc <- cbind(1, scale(X, scale = FALSE))
#' result.adj2.adj2c.random.ls <- fit.ate.adj2.adj2c.Random.CRE(Y, Xc, A)
#' point_est <- result.adj2.adj2c.random.ls$tau_vec
#' var_est <- result.adj2.adj2c.random.ls$var_vec
#' point_est
#' var_est

fit.ate.adj2.adj2c.Random.CRE <- function(Y, X, A) {

  # Calculate variance based on the randomization-based framework given in Wang et.al.
  # randomization scheme: completely randomized experiment
  # X is centered, can include constant or not

  n <- nrow(X)
  n1 <- sum(A == 1)
  n0 <- n - n1

  ### tau_unadj-
  tau_unadj <- mean(Y[A == 1]) - mean(Y[A == 0])
  var_unadj <- 1 / n * (var(Y[A == 1]) / mean(A) + var(Y[A == 0]) / (1 - mean(A)))

  ### tau_db (tau_adj2c)-
  r1 <- mean(A)
  r0 <- 1 - r1
  Y1_bar <- mean(Y[A == 1])
  Y0_bar <- mean(Y[A == 0])

  S_X2 <- t(X) %*% X / (n)
  S_X2_inv <- solve(S_X2)
  S_X1Y1 <- t(X[A == 1, ]) %*% (Y[A == 1] - Y1_bar) / (n1)
  S_X0Y0 <- t(X[A == 0, ]) %*% (Y[A == 0] - Y0_bar) / (n0)

  beta1 <- S_X2_inv %*% S_X1Y1
  beta0 <- S_X2_inv %*% S_X0Y0
  tau_adj <- sum(A * (Y - X %*% beta1)) / n1 - sum((1 - A) * (Y - X %*% beta0)) / n0

  H <- X %*% solve(t(X) %*% X)  %*% t(X)
  db <- r1 * r0 * (sum(A * diag(H) * (Y - Y1_bar)) / (r1 ^ 2 * n1)
                   - sum((1 - A) * diag(H) * (Y - Y0_bar)) / (r0 ^ 2 * n0))
  tau_adj2c <- tau_adj + db


  ### tau_db-
  mu1_adj2 <- X %*% (S_X2_inv %*% t(X[A == 1, ]) %*% (Y[A == 1]) / (n1))
  mu0_adj2 <- X %*% (S_X2_inv %*% t(X[A == 0, ]) %*% (Y[A == 0]) / (n0))
  tau_adj2 <- sum(A * (Y - mu1_adj2)) / n1 - sum((1 - A) * (Y - mu0_adj2)) / n0 +
    r1 * r0 * (sum(A * diag(H) * (Y)) / (r1 ^ 2 * n1) -  sum((1 - A) * diag(H) * (Y)) / (r0 ^
                                                                                           2 * n0))


  ### var_tau_db (tau_adj2c)
  Q <- H ^ 2
  diag(Q) <- diag(H) - diag(H) ^ 2
  M <- diag(n) - matrix(1 / n, n, n) - H + (diag(n) - matrix(1 / n, n, n)) %*% diag(diag(H))
  B <- t(M) %*% M
  I1 <- func_vec_s2(q = r1 * r0 * diag(Q), y = Y / r1 ** 2, Z = A, z = 1) + func_vec_s2(q = r1 * r0 * diag(Q), y = Y / r0 ** 2, Z = A, z = 0) +
    func_vec_s2(q=diag(B)/r1**2,y=Y,Z=A,z=1) + func_vec_s2(q=diag(B)/r0**2,y=Y,Z=A,z=0)
  I1 <- r1 * r0 * I1

  I2 <- func_mat_s2(Q = r1 * r0 * Q, Y = Y / r1 ** 2, Z = A, z = 1) + func_mat_s2(Q = r1 * r0 * Q, Y = Y / r0 ** 2, Z = A, z = 0) +
    func_mat_s2(Q=B/r1**2,Y=Y,Z=A,z=1) + func_mat_s2(Q=B/r0**2,Y=Y,Z=A,z=0)
  I2 <- r1 * r0 * I2


  I3_ub <- func_vec_s2(q = diag(B),y=Y,Z=A,z=1) - func_vec_s2(q=diag(Q),y=Y,Z=A,z=1) - func_mat_s2(Q=H,Y=Y,Z=A,z=1) +
    func_vec_s2(q = diag(B),y=Y,Z=A,z=0) - func_vec_s2(q=diag(Q),y=Y,Z=A,z=0) - func_mat_s2(Q=H,Y=Y,Z=A,z=0)
  I3_term <- func_mat_cross(H = H, Y = Y, Z = A)
  I3_ub <- I3_ub + 2 * I3_term

  I4_ub <- 2 * (func_mat_cross(H = B, Y, A) - func_mat_cross(H = Q, Y, A))

  I3_ub_new <- func_vec_s2(q = diag(B),y=Y,Z=A,z=1) - func_vec_s2(q=diag(Q),y=Y,Z=A,z=1) +
    func_vec_s2(q = diag(B),y=Y,Z=A,z=0) - func_vec_s2(q=diag(Q),y=Y,Z=A,z=0)


  var_adj2c <-  (I1 + I2 + min(I3_ub, I3_ub_new) + I4_ub) / n

  ##### Next consider another conservative variance estimator of HOIF, and plus I3, I4
  tau1_unadj <- mean(Y[A == 1])
  pi1_hat <- mean(A)
  pi0_hat <- 1 - pi1_hat
  Y1_sub <- Y[A == 1]
  H1_sub <- H[A == 1, A == 1]
  h <- diag(H)
  resi <- Y - tau1_unadj - (H - diag(diag(H))) %*% (A * (Y - tau1_unadj) / pi1_hat) - mean((1 + diag(H)) * A * (Y - tau1_unadj) / pi1_hat)
  var_tau1_adj2c_v2 <- (pi0_hat/pi1_hat)/n/n*sum(A*resi^2/pi1_hat) +
    (pi0_hat/pi1_hat)^2/n*(
      1/n1*sum(A*h*(1-h)*(Y-tau1_unadj)^2) +
        1/n1/pi1_hat*(t(Y1_sub-tau1_unadj) %*% H1_sub^2 %*% (Y1_sub - tau1_unadj) -
                        sum(diag(H1_sub)^2 * (Y1_sub - tau1_unadj)^2)))

  Y0_sub <- Y[A == 0]

  H0_sub <- H[A == 0, A == 0]
  tau0_unadj <- mean(Y[A == 0])
  resi <- Y - tau0_unadj - (H - diag(diag(H))) %*% ((1 - A) * (Y - tau0_unadj) / pi0_hat) - mean((1 + diag(H)) * (1 - A) * (Y - tau0_unadj) / pi0_hat)

  var_tau0_adj2c_v2 <- (pi1_hat/pi0_hat)/n/n*sum((1-A)*resi^2/pi0_hat) +
    (pi1_hat/pi0_hat)^2/n*(
      1/n0*sum((1-A)*h*(1-h)*(Y-tau0_unadj)^2) +
        1/n0/pi0_hat*(t(Y0_sub-tau0_unadj) %*% H0_sub^2 %*% (Y0_sub - tau0_unadj) -
                        sum(diag(H0_sub)^2 * (Y0_sub - tau0_unadj)^2)))


  var_adj2c_v2 <- as.numeric(var_tau1_adj2c_v2 + var_tau0_adj2c_v2 + (min(I3_ub, I3_ub_new) + I4_ub) / n)


  # var_tau_adj2
  I1 <- func_vec_s2(q=r1*r0*diag(Q),y=Y/r1**2,Z=A,z=1,is.center = F) + func_vec_s2(q=r1*r0*diag(Q),y=Y/r0**2,Z=A,z=0,is.center = F) +
    func_vec_s2(q=diag(B)/r1**2,y=Y,Z=A,z=1,is.center = F) + func_vec_s2(q=diag(B)/r0**2,y=Y,Z=A,z=0,is.center = F)
  I1 <- r1 * r0 * I1
  I2 <- func_mat_s2(Q=r1*r0*Q,Y=Y/r1**2,Z=A,z=1,is.center = F) + func_mat_s2(Q=r1*r0*Q,Y=Y/r0**2,Z=A,z=0,is.center = F) +
    func_mat_s2(Q=B/r1**2,Y=Y,Z=A,z=1,is.center = F) + func_mat_s2(Q=B/r0**2,Y=Y,Z=A,z=0,is.center = F)
  I2 <- r1 * r0 * I2
  (I1 + I2) / n

  I3_ub <- func_vec_s2(q = diag(B),y=Y,Z=A,z=1,is.center = F) - func_vec_s2(q=diag(Q),y=Y,Z=A,z=1,is.center = F) - func_mat_s2(Q=H,Y=Y,Z=A,z=1,is.center = F) +
    func_vec_s2(q = diag(B),y=Y,Z=A,z=0,is.center = F) - func_vec_s2(q=diag(Q),y=Y,Z=A,z=0,is.center = F) - func_mat_s2(Q=H,Y=Y,Z=A,z=0,is.center = F)
  I3_term <- func_mat_cross(H=H,Y=Y,Z=A,is.center = F)
  I3_ub <- I3_ub + 2*I3_term

  I3_ub_new <- func_vec_s2(q = diag(B),y=Y,Z=A,z=1,is.center = F) - func_vec_s2(q=diag(Q),y=Y,Z=A,z=1,is.center = F) +
    func_vec_s2(q = diag(B),y=Y,Z=A,z=0,is.center = F) - func_vec_s2(q=diag(Q),y=Y,Z=A,z=0,is.center = F)

  I4_ub <- 2*(func_mat_cross(H=B,Y,A,is.center = F) - func_mat_cross(H=Q,Y,A,is.center = F))

  var_adj2 <-  (I1 + I2 + min(I3_ub, I3_ub_new) + I4_ub) / n

  ##### Next consider another conservative variance estimator of HOIF, and plus I3, I4
  pi1_hat <- mean(A)
  pi0_hat <- 1 - pi1_hat
  Y1_sub <- Y[A == 1]
  H1_sub <- H[A == 1, A == 1]
  h <- diag(H)
  resi <- Y - (H - diag(diag(H))) %*% (A*(Y)/pi1_hat) - mean((1+diag(H))*A*(Y)/pi1_hat)
  var_tau1_adj2_v2 <- (pi0_hat/pi1_hat)/n/n*sum(A*resi^2/pi1_hat) +
    (pi0_hat/pi1_hat)^2/n*(
      1/n1*sum(A*h*(1-h)*(Y)^2) +
        1/n1/pi1_hat*(t(Y1_sub) %*% H1_sub^2 %*% (Y1_sub) -
                        sum(diag(H1_sub)^2 * (Y1_sub)^2)))

  Y0_sub <- Y[A == 0]

  H0_sub <- H[A == 0, A == 0]
  resi <- Y - (H - diag(diag(H))) %*% ((1 - A) * (Y) / pi0_hat) - mean((1 + diag(H)) * (1 - A) * (Y) / pi0_hat)

  var_tau0_adj2_v2 <- (pi1_hat/pi0_hat)/n/n*sum((1-A)*resi^2/pi0_hat) +
    (pi1_hat/pi0_hat)^2/n*(
      1/n0*sum((1-A)*h*(1-h)*(Y)^2) +
        1/n0/pi0_hat*(t(Y0_sub) %*% H0_sub^2 %*% (Y0_sub) -
                        sum(diag(H0_sub)^2 * (Y0_sub)^2)))

  var_adj2_v2 <- as.numeric(var_tau1_adj2_v2 + var_tau0_adj2_v2 + (min(I3_ub, I3_ub_new) + I4_ub) / n)


  tau_vec <- c(tau_adj2, tau_adj2c)
  names(tau_vec) <- c('adj2', 'adj2c')

  var_vec <- c(var_adj2, var_adj2c, var_adj2_v2, var_adj2c_v2)
  names(var_vec) <- c('adj2', 'adj2c', 'adj2_v2', 'adj2c_v2')

  return(list(tau_vec = tau_vec, var_vec = var_vec))

}


#' @keywords internal
func_vec_s2 <- function(q, y, Z, z, is.center = TRUE) {
  qz <- q[Z == z]
  if (is.center) {
    yz <- y[Z == z] - mean(y[Z == z])
  } else{
    yz <- y[Z == z]
  }
  return(mean(qz * yz ** 2))
}

#' @keywords internal
func_mat_s2 <- function(Q, Y, Z, z, is.center = TRUE) {
  Qz <- Q[Z == z, Z == z]

  if (is.center) {
    Yz <- Y[Z == z] - mean(Y[Z == z])
  } else{
    Yz <- Y[Z == z]
  }
  n <- length(Y)
  Qz <- Qz - diag(diag(Qz)) #对角线为0
  nz <- sum(Z == z)
  rz <- nz / n
  return(sum(tcrossprod(Yz, Yz) * Qz) / (rz * nz))
}

#' @keywords internal
func_mat_cross <- function(H, Y, Z, is.center = TRUE) {
  if (is.center) {
    Y1_tilde <- Y - mean(Y[Z == 1])
    Y0_tilde <- Y - mean(Y[Z == 0])
  } else{
    Y1_tilde <- Y
    Y0_tilde <- Y
  }
  n <- length(Y)
  r1 <- mean(Z)
  r0 <- 1 - r1
  tmp_result <- 0
  for (i in 1:n) {
    for (j in 1:n) {
      if (i != j & Z[i] == 1 & Z[j] == 0) {
        tmp_result <- tmp_result + H[i, j] * Y1_tilde[i] * Y0_tilde[j]
      }
    }
  }
  return(tmp_result / (n * r1 * r0))
}



#' Covariate-Adjusted Treatment Effect Estimation under the Super-Population Framework
#'
#' Implements HOIF-inspired debiased estimators for average treatment effect (ATE) with variance estimation
#' using influence function-based and asymptotic-variance. Designed for simple randomization experiments with moderately high-dimensional covariates.
#'
#'
#' @param Y Numeric vector of length n containing observed responses.
#' @param X Numeric matrix (n x p) of covariates. Centering is required. May include intercept column.
#' @param A Binary vector of length n indicating treatment assignment (1 = treatment, 0 = control).
#' @param intercept Logical. If TRUE (default), X already contains intercept. Set FALSE if X does not contain intercept.
#' @param pi1_hat Default is NULL. The assignment probability for the simple randomization.
#'
#' @return A list containing three named vectors:
#' \describe{
#'   \item{tau_vec}{Point estimates:
#'     \itemize{
#'       \item{\code{adj2}:} Point estimation of the HOIF-inspired debiased estimator (Zhao et al., 2024).
#'       \item{\code{adj2c}:} Point estimation of the the HOIF-inspired debiased estimator (Zhao et al., 2024), which is also the debiased estimator given by Lu et al. (2023).
#'     }}
#'   \item{var_infl_vec}{Influence function-based variance estimates:
#'     \itemize{
#'       \item{\code{adj2}:} Variance for \code{adj2} via the sample variance of its influence function formula.
#'       \item{\code{adj2c}:} Variance for \code{adj2c} via the sample variance of its influence function formula.
#'     }}
#'   \item{var_rb_vec}{Variance estimates inspired by RobinCar:
#'     \itemize{
#'       \item{\code{adj2}:} Variance  for \code{adj2} following the asymptotic variance given by Bannick et al. (2025).
#'       \item{\code{adj2c}:} Variance  for \code{adj2c} following the asymptotic variance given by Bannick et al. (2025).
#'     }}
#' }
#'
#' @export
#'
#' @references
#' Marlena S. Bannick, Jun Shao, Jingyi Liu, Yu Du, Yanyao Yi, Ting Ye (2025). A General Form of Covariate Adjustment in Randomized Clinical Trials. Biometrika.
#'
#' Xin Lu, Fan Yang, and Yuhao Wang. Debiased regression adjustment in completely randomized experiments with moderately high-dimensional covariates. arXiv preprint arXiv:2309.02073, 2023.
#'
#' Sihui Zhao, Xinbo Wang, Lin Liu, & Xin Zhang. Covariate adjustment in randomized experiments motivated by higher-order influence functions. arXiv preprint. https://arxiv.org/abs/2411.08491.
#'
#' @examples
#'
#' set.seed(120)
#' alpha0 <- 0.1;
#' n <- 400;
#'
#' p0 <- ceiling(n * alpha0)
#' beta0_full <- 1 / (1:p0) ^ (1 / 2) * (-1) ^ c(1:p0)
#' beta <- beta0_full / norm(beta0_full,type='2')
#'
#' Sigma_true <- matrix(0, nrow = p0, ncol = p0)
#' for (i in 1:p0) {
#'   for (j in 1:p0) {
#'     Sigma_true[i, j] <- 0.1 ** (abs(i - j))
#'   }
#' }
#'
#' X <- mvtnorm::rmvt(n, sigma = Sigma_true, df = 3)
#'
#' lp0 <- X %*% beta
#' delta_X <- 1  -  1/4 * X[, 2] -  1/8 * X[, 3]
#' lp1 <- lp0 + delta_X
#'
#' Y0 <- lp0 + rnorm(n)
#' Y1 <- lp1 + rnorm(n)
#'
#'
#' pi1 <- 1 / 2
#' A <- rbinom(n, size = 1, prob = pi1)
#' Y <- A * Y1 + (1 - A) * Y0
#'
#' Xc <- cbind(1, scale(X, scale = FALSE))
#' result.ate.adj2.adj2c.sp.ls <- fit.ate.adj2.adj2c.Super.SR(Y, Xc, A, intercept = TRUE)
#' result.ate.adj2.adj2c.sp.ls
#'

fit.ate.adj2.adj2c.Super.SR <- function(Y,
                                        X,
                                        A,
                                        intercept = TRUE,
                                        pi1_hat = NULL) {

  # Estimate ATE and its variance based on influence function and Robincar's variance estimator under the Super-Population based framework.
  # randomization scheme: simple randomization
  # X is centered, can include constant or not

  n <- nrow(X)
  n1 <- sum(A == 1)
  n0 <- n - n1

  if(is.null(pi1_hat)){
    pi1_hat <- mean(A)
  }


  ### tau_unadj-
  tau_unadj <- mean(Y[A == 1]) - mean(Y[A == 0])
  var_unadj <- 1 / n * (var(Y[A==1]) / mean(A) + var(Y[A==0]) / (1 - mean(A)))


  #### adj2
  H <- X %*% solve(t(X) %*% X) %*% t(X)

  yt_adj2_hat <- as.numeric((H - diag(diag(H))) %*% (A * Y / pi1_hat))
  yc_adj2_hat <- as.numeric((H - diag(diag(H))) %*% ((1 - A) * Y / (1 - pi1_hat)))

  psi1_adj2_vec <- A * Y / pi1_hat + (1 - A / pi1_hat) * yt_adj2_hat
  psi0_adj2_vec <- (1 - A) * Y / (1 - pi1_hat) + (1 - (1 - A) / (1 - pi1_hat)) * yc_adj2_hat
  tau1_adj2 <- mean(psi1_adj2_vec)
  tau0_adj2 <- mean(psi0_adj2_vec)

  if(intercept){
    infl1_adj2 <- psi1_adj2_vec - tau1_adj2
    infl0_adj2 <- psi0_adj2_vec - tau0_adj2
  }else{
    infl1_adj2 <- A * (Y - tau1_adj2) / pi1_hat + (1 - A / pi1_hat) * yt_adj2_hat
    infl0_adj2 <- (1 - A) * (Y - tau0_adj2) / (1 - pi1_hat) + (1 - (1 - A) / (1 - pi1_hat)) * yc_adj2_hat
  }

  tau_adj2 <- tau1_adj2 - tau0_adj2


  ## variance using IF
  var_infl_tau_adj2 <- 1 / n * mean((infl1_adj2 - infl0_adj2)**2)


  ### Codes from Bannick et al.
  # Get covariance between observed Y and predicted mu counterfactuals
  get.cov.Ya <- function(a, mumat){
    t_group <- A == a
    cv <- stats::cov(Y[t_group], mumat[t_group, ], use = 'complete.obs')
    return(cv)
  }
  get.cov.robincar <- function(mumat){
    # Covariance matrix between Y and mu
    cov_Ymu <- sapply(c(0, 1), get.cov.Ya, mumat)
    var_mu <- stats::var(mumat, na.rm = TRUE)
    Vhat <- (diag(c(var(Y[A == 0]), var(Y[A == 1]))) + diag(diag(var_mu)) - 2 * diag(diag(cov_Ymu))) %*% diag(1 / c(1 - pi1_hat, pi1_hat)) +
      cov_Ymu + t(cov_Ymu) - var_mu
    return(Vhat)
  }

  Vhat_adj2 <- get.cov.robincar(cbind(yc_adj2_hat, yt_adj2_hat))
  diff_mat <- matrix(c(-1,1),nrow=2)
  var_rb_tau_adj2 <- 1 / n * as.numeric(t(diff_mat) %*% Vhat_adj2 %*% diff_mat)


  # adj2c
  tau1_unadj <- mean(A * Y / pi1_hat); tau0_unadj <- mean((1 - A) * Y / (1 - pi1_hat))
  yt_adj2c_hat <- as.numeric((H - diag(diag(H))) %*% (A * (Y - tau1_unadj) / pi1_hat))
  yc_adj2c_hat <- as.numeric((H - diag(diag(H))) %*% ((1 - A) * (Y - tau0_unadj) / (1 - pi1_hat)))

  psi1_adj2c_vec <- A * Y / pi1_hat + (1 - A / pi1_hat) * yt_adj2c_hat
  psi0_adj2c_vec <- (1 - A) * Y / (1 - pi1_hat) + (1 - (1 - A) / (1 - pi1_hat)) * yc_adj2c_hat
  tau1_adj2c <- mean(psi1_adj2c_vec)
  tau0_adj2c <- mean(psi0_adj2c_vec)

  infl1_adj2c <- A * (Y - tau1_adj2c) / pi1_hat + (1 - A / pi1_hat) * yt_adj2c_hat
  infl0_adj2c <- (1 - A) * (Y - tau0_adj2c) / (1 - pi1_hat) + (1 - (1 - A) / (1 - pi1_hat)) * yc_adj2c_hat

  tau_adj2c <- tau1_adj2c - tau0_adj2c


  var_infl_tau_adj2c <- 1 / n * mean((infl1_adj2c - infl0_adj2c)**2)

  Vhat_adj2c <- get.cov.robincar(cbind(yc_adj2c_hat, yt_adj2c_hat))
  var_rb_tau_adj2c <- 1 / n * as.numeric(t(diff_mat) %*% Vhat_adj2c %*% diff_mat)


  tau_vec <- c(tau_adj2, tau_adj2c)
  names(tau_vec) <- c('adj2', 'adj2c')

  var_infl_vec <- c(var_infl_tau_adj2, var_infl_tau_adj2c)
  var_rb_vec <- c(var_rb_tau_adj2, var_rb_tau_adj2c)
  names(var_infl_vec) <- names(var_rb_vec) <- c('adj2','adj2c')

  return(list(
    tau_vec = tau_vec,
    var_infl_vec = var_infl_vec,
    var_rb_vec = var_rb_vec
  ))

}


#' Covariate-Adjusted Treatment Effect Estimation under the Super-Population Framework
#'
#' Implements HOIF-inspired debiased estimators for  treatment effect on the treatment/control arm, with variance estimation
#' using influence function-based and asymptotic-variance. Designed for simple randomization experiments with moderately high-dimensional covariates.
#'
#' @param Y Numeric vector of length n containing observed responses.
#' @param X Numeric matrix (n x p) of covariates. Centering is required. May include intercept column.
#' @param A Binary vector of length n indicating treatment assignment (1 = treatment, 0 = control).
#' @param intercept Logical. If TRUE (default), X already contains intercept. Set FALSE if X does not contain intercept.
#' @param Treated Logical. If TRUE (default), consider the treatment arm. Set FALSE if considering the control arm.
#' @param pi1_hat Default is NULL. The assignment probability for the simple randomization.
#'
#'
#' @return A list containing three named vectors:
#' \describe{
#'   \item{tau_vec}{Point estimates:
#'     \itemize{
#'       \item{\code{adj2}:} Point estimation of the HOIF-inspired debiased estimator (Zhao et al., 2024).
#'       \item{\code{adj2c}:} Point estimation of the the HOIF-inspired debiased estimator (Zhao et al., 2024), which is also the debiased estimator given by Lu et al. (2023).
#'     }}
#'   \item{var_infl_vec}{Influence function-based variance estimates:
#'     \itemize{
#'       \item{\code{adj2}:} Variance for \code{adj2} via the sample variance of its influence function formula.
#'       \item{\code{adj2c}:} Variance for \code{adj2c} via the sample variance of its influence function formula.
#'     }}
#'   \item{var_rb_vec}{Variance estimates inspired by RobinCar:
#'     \itemize{
#'       \item{\code{adj2}:} Variance  for \code{adj2} following the asymptotic variance given by Bannick et al. (2025).
#'       \item{\code{adj2c}:} Variance  for \code{adj2c} following the asymptotic variance given by Bannick et al. (2025).
#'     }}
#' }
#' @export
#'
#' @references
#' Marlena S. Bannick, Jun Shao, Jingyi Liu, Yu Du, Yanyao Yi, Ting Ye (2025). A General Form of Covariate Adjustment in Randomized Clinical Trials. Biometrika.
#'
#' Xin Lu, Fan Yang, and Yuhao Wang. Debiased regression adjustment in completely randomized experiments with moderately high-dimensional covariates. arXiv preprint arXiv:2309.02073, 2023.
#'
#' Sihui Zhao, Xinbo Wang, Lin Liu, & Xin Zhang. Covariate adjustment in randomized experiments motivated by higher-order influence functions. arXiv preprint. https://arxiv.org/abs/2411.08491.
#'
#' @examples
#'
#' set.seed(120)
#' alpha0 <- 0.1;
#' n <- 400;
#'
#' p0 <- ceiling(n * alpha0)
#' beta0_full <- 1 / (1:p0) ^ (1 / 2) * (-1) ^ c(1:p0)
#' beta <- beta0_full / norm(beta0_full,type='2')
#'
#' Sigma_true <- matrix(0, nrow = p0, ncol = p0)
#' for (i in 1:p0) {
#'   for (j in 1:p0) {
#'     Sigma_true[i, j] <- 0.1 ** (abs(i - j))
#'   }
#' }
#'
#' X <- mvtnorm::rmvt(n, sigma = Sigma_true, df = 3)
#'
#' lp0 <- X %*% beta
#' delta_X <- 1  -  1/4 * X[, 2] -  1/8 * X[, 3]
#' lp1 <- lp0 + delta_X
#'
#' Y0 <- lp0 + rnorm(n)
#' Y1 <- lp1 + rnorm(n)
#'
#'
#' pi1 <- 1 / 2
#' A <- rbinom(n, size = 1, prob = pi1)
#' Y <- A * Y1 + (1 - A) * Y0
#'
#' Xc <- cbind(1, scale(X, scale = FALSE))
#' result.treat.adj2.adj2c.sp.ls <- fit.mean.adj2.adj2c.Super.SR(Y, Xc, A, intercept = TRUE, Treated = TRUE)
#' result.control.adj2.adj2c.sp.ls <- fit.mean.adj2.adj2c.Super.SR(Y, Xc, A, intercept = TRUE, Treated = FALSE)
#' result.treat.adj2.adj2c.sp.ls
#' result.control.adj2.adj2c.sp.ls
#'
fit.mean.adj2.adj2c.Super.SR <- function(Y,
                                         X,
                                         A,
                                         intercept = FALSE,
                                         Treated = TRUE,
                                         pi1_hat = NULL) {

  # calculate variance based on influence function and Robincar's variance estimator under the super-population based framework
  # randomization scheme: simple randomization
  # X is centered, can include constant or not
  # Treated: if TRUE, estimate E[Y(1)] and its variance

  n <- nrow(X)

  A_new <- rep(0, n)
  if(Treated){
    A_new[A == 1] <- 1
  }else{
    A_new[A == 0] <- 1
  }
  A <- A_new
  Y[A == 0] <- 0

  result.ls <- fit.ate.adj2.adj2c.Super.SR(Y, X, A, intercept, pi1_hat)
  return(result.ls)
}
