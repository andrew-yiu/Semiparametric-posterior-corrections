library(foreach)
require(doSNOW)
require(Kendall)
library(dirichletprocess)
library(ExtDist)
library(invgamma)

library(CVST)
library(MASS)
library(doRNG)
library(tictoc)

setwd("~/Downloads/")

source('DirichletProcessIndepVariance.R')

set.seed(1)

#Correction
correct_l2 <- function(l2_samp,pdf_n_samp){
  B <- dim(pdf_n_samp)[1]
  n <- dim(pdf_n_samp)[2]
  w <- matrix(rexp(n*B),nrow = B)
  w <- w/rowSums(w)
  
  eic = 2*( rowSums(w*pdf_n_samp) - l2_samp) 
  l2_samp_corrected <- l2_samp + eic
  return(l2_samp_corrected)
}

#Initialize y
n <- 2000
y_plot = seq(-4,4,length.out = 500)
dy = y_plot[2]- y_plot[1]
n_rep <- 1000

l2_truth <- sum(dLaplace(y_plot)^2*dy)

var_truth <- sum(4*dLaplace(y_plot)*((dLaplace(y_plot) - l2_truth)^2)*dy)
len_true <- sqrt(var_truth/1000)*1.96*2

cores <- 4
cl <- makeSOCKcluster(cores)
registerDoSNOW(cl)

trials <- 1000
pb <- txtProgressBar(max = trials, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

tic()
pred <- foreach(k=1:trials, .packages=c("dirichletprocess","ExtDist", "invgamma"),
                .combine='rbind',.options.RNG=1,.options.snow=opts) %dorng%{
                  
                  source('DirichletProcessIndepVariance.R')
                  #Generate data
                  y <- rLaplace(n)
                  
                  n_samp = 5000 
                  n_burn = 1000
                  #dp_fit <- DirichletProcessGaussian(y) 
                  #system.time(dp_fit <- Burn(Fit(dp_fit, n_samp),n_burn))
                  
                  #New DPMM object and fit function 
                  dp_fit <- DirichletProcessGaussianIndepVariance(y,g0Priors = c(0,1,1,1),alpha = 1,mhStepSize = 0.5) 
                  system.time(dp_fit <- Burn(Fit.gaussianIndepVariance(dp_fit, n_samp),n_burn))

                  #Compute posterior samples of pdf
                  ind = 1:(n_samp-n_burn)
                  #pdf_grid_samp <- t(sapply(ind, function(i) PosteriorFunction(dp_fit,i)(y_plot)))
                  #pdf_n_samp <- t(sapply(ind, function(i) PosteriorFunction(dp_fit,i)(y)))
                  
                  #New posterior function sampling
                  pdf_grid_samp <- t(sapply(ind, function(i) PosteriorFunction.normalIndepVariance(dp_fit,i)(y_plot)))
                  pdf_n_samp <- t(sapply(ind, function(i) PosteriorFunction.normalIndepVariance(dp_fit,i)(y)))
                  
                  #Compute posterior of integral squared
                  l2_samp <- sapply(1:dim(pdf_grid_samp)[1],function(i) sum(pdf_grid_samp[i,]^2*dy))
                  l2_samp_corrected <- correct_l2(l2_samp,pdf_n_samp)
                  
                  post_mean <- mean(l2_samp)
                  lower <-quantile(l2_samp,probs = 0.025)
                  upper <-quantile(l2_samp,probs = 0.975)
                  
                  bias <- post_mean - l2_truth
                  coverage <- (l2_truth >=lower) & (l2_truth <= upper)
                  length <- upper - lower
                  
                  post_mean_corrected <- mean(l2_samp_corrected)
                  lower_corrected <-quantile(l2_samp_corrected,probs = 0.025)
                  upper_corrected <-quantile(l2_samp_corrected,probs = 0.975)
                  
                  bias_corrected <- post_mean_corrected - l2_truth
                  coverage_corrected <- (l2_truth >=lower_corrected) & (l2_truth <= upper_corrected)
                  length_corrected <- upper_corrected - lower_corrected
                  
                  c(post_mean, bias, coverage, length, post_mean_corrected, bias_corrected, coverage_corrected, length_corrected)
                }
toc()
stopCluster(cl)

require(gridExtra)
# 
# #png("dpmm_both.png", units="in", width=11, height=6, res=300)
# 
dat_ate <-  data.frame(ate_post=pred[,1], ate_post2 = pred[,5])
library(ggplot2)
g <- ggplot()
g <- g + geom_histogram(data = dat_ate, aes(y=..density..,ate_post, fill= "DPMM"), alpha = 0.3, bins = 20)
g <- g + geom_density(data = dat_ate, aes(ate_post, colour= "DPMM"))
g <- g + geom_histogram(data = dat_ate, aes(y=..density..,ate_post2, fill= "DPMM+1step"), alpha = 0.3, bins = 20)
g <- g + geom_density(data = dat_ate, aes(ate_post2, colour= "DPMM+1step"))
g <- g + labs(fill = "Posterior")
g <- g+ geom_vline(xintercept=l2_truth, linetype="dashed", color = "black")
g <- g +guides(colour="none")
g <- g+ xlim(0.20,0.3)
g <- g+ ylim(0,34)
g <- g + xlab("Integrated squared density") + ylab("Posterior mean density")
g
# 
# plot1 <- g
# 
# dat_ate <-  data.frame(ate_post=res[,1], ate_post2 = res[,5])
# library(ggplot2)
# h <- ggplot()
# h <- h + geom_histogram(data = dat_ate, aes(y=..density..,ate_post, fill= "DPMM"), alpha = 0.3, bins = 20)
# h <- h + geom_density(data = dat_ate, aes(ate_post, colour= "DPMM"))
# h <- h + geom_histogram(data = dat_ate, aes(y=..density..,ate_post2, fill= "DPMM+1step"), alpha = 0.3, bins = 20)
# h <- h + geom_density(data = dat_ate, aes(ate_post2, colour= "DPMM+1step"))
# h <- h + labs(fill = "Posterior") 
# h <- h+ geom_vline(xintercept=l2_truth, linetype="dashed", color = "black")
# h <- h +guides(colour="none")
# h <- h+ xlim(0.2,0.3)
# h <- h+ ylim(0,34)
# h <- h + xlab("Integrated squared density") + ylab("Posterior mean density")
# h
# 
# plot2 <- h
# grid.arrange(plot1, plot2, ncol=2)

# dev.off()
# 
# 

pred <- readRDS("dplmm1000.rds")
res <- readRDS("dplmm_2000.rds")

png("dpmm_both.png", units="in", width=11, height=6, res=300)
dat_ate <-  data.frame(ate_post=pred[,1], ate_post2 = pred[,5], oracle = ora1000[,1])
library(ggplot2)
g <- ggplot()
g <- g + geom_histogram(data = dat_ate, aes(y=..density..,ate_post, fill= "DPMM"), alpha = 0.3, bins = 20)
g <- g + geom_density(data = dat_ate, aes(ate_post, colour= "DPMM"))
g <- g + geom_histogram(data = dat_ate, aes(y=..density..,ate_post2, fill= "DPMM+1step"), alpha = 0.3, bins = 20)
g <- g + geom_density(data = dat_ate, aes(ate_post2, colour= "DPMM+1step"))
g <- g + geom_density(data = dat_ate, aes(oracle, colour= "Oracle"), linetype = "dashed")
g <- g + labs(fill = "Posterior")
g <- g+ geom_vline(xintercept=l2_truth, linetype="dashed", color = "black")
g <- g +guides(colour="none")
g <- g+ xlim(0.22,0.27)
g <- g+ ylim(0,60)
g <- g + xlab("Integrated squared density") + ylab("Posterior mean density")
g

plot1 <- g

dat_ate <-  data.frame(ate_post=res[,1], ate_post2 = res[,5], oracle = ora_2000[,1])
library(ggplot2)
h <- ggplot()
h <- h + geom_histogram(data = dat_ate, aes(y=..density..,ate_post, fill= "DPMM"), alpha = 0.3, bins = 20)
h <- h + geom_density(data = dat_ate, aes(ate_post, colour= "DPMM"))
h <- h + geom_histogram(data = dat_ate, aes(y=..density..,ate_post2, fill= "DPMM+1step"), alpha = 0.3, bins = 20)
h <- h + geom_density(data = dat_ate, aes(ate_post2, colour= "DPMM+1step"))
h <- h + geom_density(data = dat_ate, aes(oracle, colour= "Oracle"), linetype = "dashed")
h <- h + labs(fill = "Posterior")
h <- h+ geom_vline(xintercept=l2_truth, linetype="dashed", color = "black")
h <- h +guides(colour="none")
h <- h+ xlim(0.22,0.27)
h <- h+ ylim(0,60)
h <- h + xlab("Integrated squared density") + ylab("Posterior mean density")
h

plot2 <- h
grid.arrange(plot1, plot2, ncol=2)

dev.off()

# dat_ate <-  data.frame(ate_post=pred[,1], ate_post2 = pred[,5])
# library(ggplot2)
# g <- ggplot()
# g <- g + geom_histogram(data = dat_ate, aes(y=..density..,ate_post, fill= "DPMM"), alpha = 0.3, bins = 20)
# g <- g + geom_density(data = dat_ate, aes(ate_post, colour= "DPMM"))
# g <- g + geom_histogram(data = dat_ate, aes(y=..density..,ate_post2, fill= "DPMM+1step"), alpha = 0.3, bins = 20)
# g <- g + geom_density(data = dat_ate, aes(ate_post2, colour= "DPMM+1step"))
# g <- g + labs(fill = "Posterior")
# g <- g+ geom_vline(xintercept=l2_truth, linetype="dashed", color = "black")
# g <- g +guides(colour="none")
# g <- g+ xlim(0.21,0.28)
# g <- g+ ylim(0,60)
# g <- g + xlab("Integrated squared density") + ylab("Posterior mean density")
# g

plot(y_plot, dLaplace(y_plot), type = "l", ylim = c(0,0.5), lty = "dashed")
lines(y_plot, pdf_grid_samp[100,], type = "l", col = "red")
lines(y_plot, pdf_grid_samp[500,], type = "l", col = "green")
lines(y_plot, pdf_grid_samp[1000,], type = "l", col = "blue")
lines(y_plot, pdf_grid_samp[1500,], type = "l", col = "purple")

plot(l2_samp, l2_samp_corrected, xlab = "uncorrected", ylab = "corrected")
scat_grid <- 0.2+c(1:10000)/100000
lines(scat_grid, scat_grid, col = "red")
abline(v = 0.25, lty = "dashed")
abline(h = 0.25, lty = "dashed")

png("scat.png", units="in", width=9, height=6, res=300)
library(ggplot2)
scat_dat <- data.frame(initial = as.numeric(l2_samp), corrected = as.numeric(l2_samp_corrected))
s <- ggplot(data = scat_dat) + geom_point(mapping = aes(x = as.numeric(initial), y = as.numeric(corrected)), size = 1, alpha = 0.4, colour = "violetred1")
s <- s + xlab("Initial posterior") + ylab("One-step posterior")
s <- s + geom_abline(intercept = 0, slope = 1, alpha = 0.3)
s <- s + geom_hline(yintercept = l2_truth, alpha = 0.5, linetype = "dashed")
s <- s + geom_vline(xintercept = l2_truth, alpha = 0.5, linetype = "dashed")
s <- s + coord_fixed()
s <- s + xlim(0.20,0.27)
s
dev.off()
