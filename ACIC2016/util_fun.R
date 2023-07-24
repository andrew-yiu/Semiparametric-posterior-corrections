cov_freq <- function(mu,est, sd,sig = 0.05){
  return(as.numeric(mu > qnorm(sig/2, mean = est, sd = sd) && mu < qnorm(1- sig/2, mean = est, sd = sd))) 
}

len_freq <- function(est,sd, sig = 0.05){
  return(qnorm(1- sig/2, mean = est, sd = sd) - qnorm(sig/2, mean = est, sd = sd))
}

cov_bayes <- function(mu,mu_samp, sig =0.05){
  return(as.numeric((mu > quantile(mu_samp, probs = sig/2) && mu < quantile(mu_samp, probs = 1-sig/2))))
  }
    
len_bayes <- function(mu_samp, sig =0.05){  
  return(quantile(mu_samp, probs = 1-sig/2) - quantile(mu_samp, probs = sig/2))
}

rdirichlet <- function(n){
  wts <- rexp(n)
  return(wts/sum(wts))
}

expit <- function(x){
  return(1/(1+exp(-x)))
}

logit <- function(p){
  return(log(p/(1-p)))
}