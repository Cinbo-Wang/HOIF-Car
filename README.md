# HOIF-Car

The *HOIFCar* package facilitates the estimation of treatment effects through a range of covariate adjustment methods. It supports the calculation of point estimates and variance estimates for treatment effects. Additionally, the package provides oracle bias, oracle variance, and the corresponding approximated variance for the adjusted estimators derived from higher-order influence functions (HOIF).
   

Complete derivation for the exact variance of $\hat{\tau}_{\mathsf{adj}, 2}^{\dagger}$ under the randomization-based framework can be found [here](https://github.com/Cinbo-Wang/HOIF-Car/blob/main/var-db.pdf).

## Installation
```R
devtools::install_github("Cinbo-Wang/HOIF-Car")
```

## Example

We first consider the Completely Randomized Experiment (CRE) with moderately high-dimensional covariates under the randomization-based framework. 
```R
# Oracle Estimatation 
## Linear setting
set.seed(100)
n <- 500
p <- 50
beta <- rt(p,3)

X <- mvtnorm::rmvt(n, sigma = diag(1, p), df = 3)
Y1 <- as.numeric(X %*% beta)
pi1 <- 0.5
n1 <- ceiling(n*pi1)

require(HOIFCar)
result_adj_db <- get_oracle_bias_var_adj_db(X = X,Y1=Y1,n1=n1) # from LYW paper
result_adj2c <- get_oracle_bias_var_adj2c(X = X,Y1=Y1,n1=n1)
result_adj2_3 <- get_oracle_bias_var_adj_2_3(X=X, Y1=Y1,n1=n1)
unlist(result_adj_db)
unlist(result_adj2c)
unlist(result_adj2_3)



## Nonlinear setting
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
or_baseline <- sign(X %*% beta) * abs(X %*% beta)^(1/2) + sin(X %*% beta)
epsilon1 <- epsilon0 <- rt(n,3)
Y1 <- 1 + as.numeric(or_baseline) + epsilon1


pi1 <- 2/3
n1 <- ceiling(n*pi1)

result_adj_db <- get_oracle_bias_var_adj_db(X = X,Y1=Y1,n1=n1) # from LYW paper
result_adj2c <- get_oracle_bias_var_adj2c(X = X,Y1=Y1,n1=n1)
result_adj2_3 <- get_oracle_bias_var_adj_2_3(X=X, Y1=Y1,n1=n1)
unlist(result_adj_db)
unlist(result_adj2c)
unlist(result_adj2_3)



# Realistic estimation
set.seed(100)
n <- 500
p <- n * 0.3
beta <- runif(p, -1 / sqrt(p), 1 / sqrt(p))

X <- mvtnorm::rmvt(n, sigma = diag(1, p), df = 3)
Y1 <- as.numeric(X %*% beta)
Y0 <- rep(0, n)

pi1 <- 2/3
n1 <- ceiling(n * pi1)
ind <- sample(n, size = n1)
A <- rep(0, n)
A[ind] <- 1
Y <- Y1 * A + Y0 * (1 - A)

Xc_svd <- svd(X)
H <- Xc_svd$u %*% t(Xc_svd$u)

## Estimate mean treat on the treatment arm
result_ls <- esti_mean_treat(X, Y, A, H)
point_est <- result_ls$point_est
var_est <- result_ls$var_est
print(paste0('True mean treat:', round(mean(Y1), digits = 3), '.'))
print('Absolute bias:')
print(abs(point_est - mean(Y1)))
print('Estimate variance:')
print(var_est)

## Estimate ATE using adj2 and adj2c estimators
Xc <- cbind(1, scale(X, scale = FALSE))
result.adj2.adj2c.random.ls <- fit.ate.adj2.adj2c.Random.CRE(Y, Xc, A)
point_est <- result.adj2.adj2c.random.ls$tau_vec
var_est <- result.adj2.adj2c.random.ls$var_vec
point_est
var_est

```
We next consider the simple randomization experiments with moderately high-dimensional covariates under the Super-Population based framework. 
```R
set.seed(120)
alpha0 <- 0.1; 
n <- 400; 

p0 <- ceiling(n * alpha0)
beta0_full <- 1 / (1:p0) ^ (1 / 2) * (-1) ^ c(1:p0)
beta <- beta0_full / norm(beta0_full,type='2')

Sigma_true <- matrix(0, nrow = p0, ncol = p0)
for (i in 1:p0) {
  for (j in 1:p0) {
    Sigma_true[i, j] <- 0.1 ** (abs(i - j))
  }
}

X <- mvtnorm::rmvt(n, sigma = Sigma_true, df = 3)

lp0 <- X %*% beta
delta_X <- 1  -  1/4 * X[, 2] -  1/8 * X[, 3]
lp1 <- lp0 + delta_X

Y0 <- lp0 + rnorm(n)
Y1 <- lp1 + rnorm(n)


pi1 <- 1 / 2
A <- rbinom(n, size = 1, prob = pi1)
Y <- A * Y1 + (1 - A) * Y0


# Estimate ATE, EY1 and EY0 using adj2 and adj2c estimators
Xc <- cbind(1, scale(X, scale = FALSE))
result.ate.adj2.adj2c.sp.ls <- fit.ate.adj2.adj2c.Super.SR(Y, Xc, A, intercept = TRUE)
result.treat.adj2.adj2c.sp.ls <- fit.mean.adj2.adj2c.Super.SR(Y, Xc, A, intercept = TRUE, Treated = TRUE)
result.control.adj2.adj2c.sp.ls <- fit.mean.adj2.adj2c.Super.SR(Y, Xc, A, intercept = TRUE, Treated = FALSE)


result.ate.adj2.adj2c.sp.ls
result.treat.adj2.adj2c.sp.ls
```

## References
Zhao, S., Wang, X., Liu, L., & Zhang, X. (2024). Covariate adjustment in randomized experiments motivated by higher-order influence functions. arXiv preprint. https://arxiv.org/abs/2411.08491

Xin Lu, Fan Yang, and Yuhao Wang. Debiased regression adjustment in completely randomized experiments with moderately high-dimensional covariates. arXiv preprint arXiv:2309.02073, 2023.

Marlena S. Bannick, Jun Shao, Jingyi Liu, Yu Du, Yanyao Yi, Ting Ye (2025). A General Form of Covariate Adjustment in Randomized Clinical Trials. Biometrika.

