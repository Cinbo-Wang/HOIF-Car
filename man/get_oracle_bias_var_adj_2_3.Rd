% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_oracle_infer.R
\name{get_oracle_bias_var_adj_2_3}
\alias{get_oracle_bias_var_adj_2_3}
\title{Estimate the oracle bias, the exact variance and approximated variance of the debiased estimator and the bias-free estimator motivated by HOIF (Zhao et al.(2024)).}
\usage{
get_oracle_bias_var_adj_2_3(X, Y1, n1 = NULL)
}
\arguments{
\item{X}{The n by p covariates matrix.}

\item{Y1}{Vector of n dimensional potential response Y(1).}

\item{n1}{The number of subjects in the treatment group.}
}
\value{
A list of oracle bias and variance of the adjusted estimator motivated by HOIF and the bias-free estimator.
\item{bias_adj2 }{The oracle bias of the estimator tau_adj2.}
\item{variance_exact_adj2}{The oracle exact variance of the estimator tau_adj2.}
\item{variance_approx_adj2}{The oracle approximated variance of the estimator tau_adj2 which omits the term of order o(1/n).}
\item{variance_exact_adj3}{The oracle exact variance of the bias-free estimator tau_adj3.}
\item{variance_unadj}{The oracle variance of the unadjusted estimator.}
}
\description{
Estimate the oracle bias, the exact variance and approximated variance of the debiased estimator and the bias-free estimator motivated by HOIF (Zhao et al.(2024)).
}
\examples{
# Linear setting
set.seed(100)
n <- 500
p <- 50
beta <- rt(p,3)

X <- mvtnorm::rmvt(n, sigma = diag(1, p), df = 3)
Y1 <- as.numeric(X \%*\% beta)
pi1 <- 0.50
n1 <- ceiling(n*pi1)

result_adj_db <- get_oracle_bias_var_adj_db(X = X,Y1=Y1,n1=n1)
result_adj2c <- get_oracle_bias_var_adj2c(X = X,Y1=Y1,n1=n1)
result_adj2_3 <- get_oracle_bias_var_adj_2_3(X = X,Y1=Y1,n1=n1)
unlist(result_adj_db)
unlist(result_adj2c)
unlist(result_adj2_3)



# Nonlinear setting
n <- 500;
alpha <- 0.2;
set.seed(1000)
p <- ceiling(n*alpha)
Sigma_true <- matrix(0,nrow=p,ncol=p)
for(i in 1:p){
  for(j in 1:p){
    Sigma_true[i,j] <- 0.1**(abs(i-j))
  }
}

X <- mvtnorm::rmvt(n, sigma = Sigma_true, df = 3)
beta <- rt(p,3)
or_baseline <- sign(X \%*\% beta) * abs(X \%*\% beta)^(1/2) + sin(X \%*\% beta)
epsilon1 <- epsilon0 <- rt(n,3)
Y1 <- 1 + as.numeric(or_baseline) + epsilon1


pi1 <- 0.50
n1 <- ceiling(n*pi1)

result_adj_db <- get_oracle_bias_var_adj_db(X = X,Y1=Y1,n1=n1) # from LYW paper
result_adj2c <- get_oracle_bias_var_adj2c(X = X,Y1=Y1,n1=n1)
result_adj2_3 <- get_oracle_bias_var_adj_2_3(X = X,Y1=Y1,n1=n1)
unlist(result_adj_db)
unlist(result_adj2c)
unlist(result_adj2_3)

}
\references{
Zhao, S., Wang, X., Liu, L., & Zhang, X. (2024). Covariate adjustment in randomized experiments motivated by higher-order influence functions. arXiv preprint. https://arxiv.org/abs/2411.08491
}
