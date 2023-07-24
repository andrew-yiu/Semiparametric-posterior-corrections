### Code to run ACIC 2016 simulations with BART+1step
library(aciccomp2016)
#load bartMachine
setwd("~/Downloads/")
library(BART)
# options(java.parameters = "-Xmx9000m")

# require(doSNOW)
library(parallel)
require(Kendall)
# library(doRNG)
num_cores <- 4
# cores <- 4
# cl <- makeSOCKcluster(cores)
# registerDoSNOW(cl)
source("baselines_ATT.R")
source("util_fun.R")
source("1step_ATT_fun.R")
param_num <- 20
print(param_num)

mcsapply<-function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, mc.preschedule = TRUE,
                    mc.set.seed = TRUE, mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L), 
                    mc.cleanup = TRUE, mc.allow.recursive = TRUE, affinity.list = NULL )
{
  answer <- mclapply(X = X, FUN = FUN, ...,mc.preschedule = mc.preschedule, 
                     mc.set.seed = mc.set.seed, mc.silent = mc.silent, mc.cores = mc.cores, 
                     mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive, affinity.list = affinity.list)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (!isFALSE(simplify) && length(answer)) 
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}

#INITIALIZE
trials <- 4

results <- matrix(, nrow = 100, ncol = 30)

start_time <- Sys.time()

pb <- txtProgressBar(max = trials, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

for (i in 1:trials){
  progress(i)
  sim_data <- dgp_2016(input_2016, param_num, i)
  x <- input_2016
  r <- sim_data$z
  y <- sim_data$y
  n = length(y)
  
  sdy <- sd(y)
  
  SATT <- mean(sim_data$y.1[sim_data$z == 1] - sim_data$y.0[sim_data$z == 1])
  
  M <- 4000
  Mskip <- 2000
  
  
  w <- cbind(x,r)
  x1 <- cbind(x,rep(1,length(r)))
  x0 <- cbind(x,rep(0,length(r)))
  names(x1)[59] <- "r"
  names(x0)[59] <- "r"
  x_comb <- rbind(x1, x0)
  
  invisible(capture.output(post <- mc.wbart(data.frame(w), y, data.frame(x=x_comb), ndpost = M, mc.cores = num_cores, seed = 99, nskip = Mskip)))
  
  fdraws1 <- post$yhat.test[,1:n]
  fdraws0 <- post$yhat.test[,(n+1):(2*n)]
  
  #BART for propensity score model with probit link
  invisible(capture.output(post <- mc.pbart(x.train = data.frame(x), y.train = r, ndpost = M, mc.cores = num_cores, seed = 99, nskip = Mskip)))
  
  pidraws <- post$prob.train
  #Estimate ATT
  
  # compute oracle one-step SACTT estimator
  pitrue <- sim_data$e
  f1true <- sim_data$mu.1
  f0true <- sim_data$mu.0
  pr1_orac <- mean(pitrue)
  CATT <- mean(pitrue*(f1true- f0true))/pr1_orac
  
  eif_orac <- (r/pr1_orac - ((1-r)/pr1_orac)*(pitrue/(1-pitrue)))*(y- (r*f1true + (1-r)*f0true)) +
    ((r-pitrue)/pr1_orac)*(f1true - f0true - CATT)
  theta_orac <- CATT + mean(eif_orac)
  
  orac_var <- mean(eif_orac^2)/n
  cov_orac = cov_freq(mu = SATT, est = theta_orac, sd = sqrt(orac_var))
  catt_cov_orac = cov_freq(mu = CATT, est = theta_orac, sd = sqrt(orac_var))
  len_orac = len_freq(est = theta_orac, sd = sqrt(orac_var))/sdy
  bias_orac = (theta_orac - SATT)/sdy
  catt_bias_orac = (theta_orac - CATT)/sdy
  
  pimean <- colMeans(pidraws)
  f1mean <- colMeans(fdraws1)
  f0mean <- colMeans(fdraws0)
  
  # compute frequentist one-step sactt estimator
  pr1_fre <- mean(pimean)
  theta_plug_sactt <- mean(pimean*(f1mean- f0mean))/pr1_fre
  
  eif_fre <- (r/mean(r) - ((1-r)/mean(r))*(pimean/(1-pimean)))*(y- (r*f1mean + (1-r)*f0mean)) +
    ((r-pimean)/mean(r))*(f1mean - f0mean - theta_plug_sactt)
  theta_fre <- theta_plug_sactt + mean(eif_fre)
  
  fre_var <- mean(eif_fre^2)/n
  cov_fre = cov_freq(mu = SATT, est = theta_fre, sd = sqrt(fre_var))
  catt_cov_fre = cov_freq(mu = CATT, est = theta_fre, sd = sqrt(fre_var))
  len_fre = len_freq(est = theta_fre, sd = sqrt(fre_var))/sdy
  bias_fre = (theta_fre - SATT)/sdy
  catt_bias_fre = (theta_fre - CATT)/sdy
  
  output <- mcsapply(X = 1:M, FUN = function(i) onestep(fdraws0[i,],fdraws1[i,],pidraws[i,],r,y,n), mc.cores = 4)
  # output <- sapply(1:M, function(i) onestep(fdraws0[,i],fdraws1[,i],pidraws[,i],r,y,n))
  #output <- sapply(1:M, function(j) onestep(fdraws0[j,],fdraws1[j,],pidraws[j,],r,y,n))
  
  theta_sactt <- output[1,]
  theta_att <- output[2,]
  theta_plug_sactt <- output[3,]
  theta_plug_catt <- output[4,]
  
  bias_sactt = (mean(theta_sactt) - SATT)/sdy
  bias_att = (mean(theta_att) - SATT)/sdy
  bias_plug_sactt = (mean(theta_plug_sactt) - SATT)/sdy
  bias_plug_catt = (mean(theta_plug_catt) - SATT)/sdy
  
  cov_sactt = cov_bayes(SATT,theta_sactt)
  cov_att = cov_bayes(SATT,theta_att)
  cov_plug_sactt = cov_bayes(SATT,theta_plug_sactt)
  cov_plug_catt = cov_bayes(SATT,theta_plug_catt)
  
  catt_bias_sactt = (mean(theta_sactt) - CATT)/sdy
  catt_bias_att = (mean(theta_att) - CATT)/sdy
  catt_bias_plug_sactt = (mean(theta_plug_sactt) - CATT)/sdy
  catt_bias_plug_catt = (mean(theta_plug_catt) - CATT)/sdy
  
  catt_cov_sactt = cov_bayes(CATT,theta_sactt)
  catt_cov_att = cov_bayes(CATT,theta_att)
  catt_cov_plug_sactt = cov_bayes(CATT,theta_plug_sactt)
  catt_cov_plug_catt = cov_bayes(CATT,theta_plug_catt)
  
  len_sactt = len_bayes(theta_sactt)/sdy
  len_att = len_bayes(theta_att)/sdy
  len_plug_sactt = len_bayes(theta_plug_sactt)/sdy
  len_plug_catt = len_bayes(theta_plug_catt)/sdy
  
  results[i,] <- c(bias_sactt, bias_att, bias_orac, bias_plug_sactt, bias_plug_catt,bias_fre, cov_sactt, cov_att, cov_orac, cov_plug_sactt, cov_plug_catt,cov_fre, len_sactt, len_att, len_orac, len_plug_sactt, len_plug_catt,len_fre, catt_bias_sactt, catt_bias_att, catt_bias_orac, catt_bias_plug_sactt, catt_bias_plug_catt,catt_bias_fre, catt_cov_sactt, catt_cov_att, catt_cov_orac, catt_cov_plug_sactt, catt_cov_plug_catt, catt_cov_fre)
  print(results[i,])
}
end <- Sys.time()-start_time
print(end)
# stopCluster(cl)

