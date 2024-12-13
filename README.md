# HOIF-Car

The *HOIFCar* package facilitates the estimation of treatment effects through a range of covariate adjustment methods. It supports the calculation of point estimates and variance estimates for treatment effects within the treatment arm. Additionally, the package provides oracle bias, oracle variance, and the corresponding approximated variance for adjusted estimators derived from higher-order influence functions (HOIF).
   

Complete derivation for the exact variance of $\hat{\tau}_{\mathsf{adj}, 2}^{\dagger}$ under the randomization-based framework can be found [here](https://github.com/Cinbo-Wang/HOIF-Car/blob/main/var-db.pdf).

## Installation
```R
devtools::install_github("Cinbo-Wang/HOIF-Car")
```

## Example

```R
# Linear setting
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



# Realistic simulation
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

result_ls <- esti_mean_treat(X, Y, A, H)
point_est <- result_ls$point_est
var_est <- result_ls$var_est
print(paste0('True mean treat:', round(mean(Y1), digits = 3), '.'))
print('Absolute bias:')
print(abs(point_est - mean(Y1)))
print('Estimate variance:')
print(var_est)

```



## References
Zhao, S., Wang, X., Liu, L., & Zhang, X. (2024). Covariate adjustment in randomized experiments motivated by higher-order influence functions. arXiv preprint. https://arxiv.org/abs/2411.08491
