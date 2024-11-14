# HOIF-Car
Covariate Adjustment in Randomized Experiments Motivated by Higher-Order Influence Functions

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

X <- mvtnorm::rmvt(n, sigma = diag(1, p), df = 2)
Y1 <- as.numeric(X %*% beta)
pi1 <- 0.45

result_adj_db <- get_oracle_var_adj_db(X = X,Y1=Y1,pi1=pi1) # from LYW paper
result_adj2_dag <- get_oracle_bias_var_adj2_dag(X = X,Y1=Y1,pi1=pi1)
result_adj2_3 <- get_oracle_bias_var_adj_2_3(X=X, Y1=Y1,pi1=pi1)
unlist(result_adj_db)
unlist(result_adj2_dag)
unlist(result_adj2_3)



# Nonlinear setting
n <- 500;
alpha <- 0.2;
set.seed(1000)
p <- ceiling(n*alpha)
Sigma_true <- matrix(0,nrow=p,ncol=p)
for(i in 1:p){
  for(j in 1:p){
    Sigma_true[i,j] <- 0.3**(abs(i-j))
  }
}

X <- mvtnorm::rmvt(n, sigma = Sigma_true, df = 3)
beta <- rt(p,3)
or_baseline <- sign(X %*% beta) * abs(X %*% beta)^(1/2) + sin(X %*% beta)
epsilon1 <- epsilon0 <- rt(n,3)
Y1 <- 1 + as.numeric(or_baseline) + epsilon1


pi1 <- 0.45
result_adj_db <- get_oracle_var_adj_db(X = X,Y1=Y1,pi1=pi1)
result_adj2_dag <- get_oracle_bias_var_adj2_dag(X = X,Y1=Y1,pi1=pi1)
result_adj2_3 <- get_oracle_bias_var_adj_2_3(X=X, Y1=Y1,pi1=pi1)
unlist(result_adj_db)
unlist(result_adj2_dag)
unlist(result_adj2_3)


```



## References
Zhao, S., Wang, X., Liu, L., & Zhang, X. (2024). Covariate adjustment in randomized experiments motivated by higher-order influence functions. arXiv preprint. https://arxiv.org/abs/2411.08491
