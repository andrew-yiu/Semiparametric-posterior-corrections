#hajek estimator
hajek <- function(r,y,x,n,pi_samp){
  est <- mean(r*y/pi_samp)/(mean(r/pi_samp))
  sd <- sqrt(mean((r*((y-est)/pi_samp))^(2))/n)
  sd <- sd/abs(mean(r/pi_samp)) #sandwich covariance
  return(data.frame(est,sd))
}

#oracle efficient estimator
oracle <- function(r,y,x,n,mfun,pr){
  mtrue <- mfun(x)
  est <- mean(r*(y-mtrue)/pr(x)+mtrue)
  sd <- sqrt(mean((r*(y-mtrue)/pr(x)+mtrue - est)^(2))/n) #sandwich
  return(data.frame(est,sd))
}

#BART; M is number of posterior samples
fit_bart <- function(r,y,x, M = 1000){
  xfull <- rbind(x[which(r==1),], x[which(r==0),]) #rearrange covariates by treatment and control
  rfull <- r[c(which(r==1), which(r!=1))]
  yfull <- y[c(which(r==1), which(r!=1))]  

  #BART for outcome regression
  #Treatment group
  (capture.output(post <- wbart(data.frame(x = x[which(r==1),]), y[which(r==1)], data.frame(x=x[which(r==0),]), ndpost = M)))
  fdraws1 <- cbind(post$yhat.train, post$yhat.test)
  
  #Control group
  (capture.output(post <- wbart(data.frame(x = x[which(r==0),]), y[which(r==0)], data.frame(x=x[which(r==1),]), ndpost = M)))
  fdraws0 <- cbind(post$yhat.test,post$yhat.train) #put treatment first
  
  #BART for propensity score model with logistic link
  invisible(capture.output(postpi <- lbart(x.train = data.frame(xfull), y.train = rfull, ndpost = M)))
  
  pidraws <- postpi$prob.train
  
  return(list("fdraws1" = fdraws1,"fdraws0" = fdraws0, "pidraws" = pidraws,"yfull" = yfull,"xfull"= xfull,"rfull"= rfull))
}

#Logistic BART; M is number of posterior samples
fit_bart_bin <- function(r,y,x, M = 1000){
  xfull <- c(x[which(r==1)], x[which(r!=1)]) #rearrange covariates by selected and missing
  rfull <- r[c(which(r==1), which(r!=1))]
  yfull <- y[c(which(r==1), which(r!=1))]  
  
  #BART for outcome regression
  invisible(capture.output(post <- lbart(data.frame(x = x[which(r==1)]), y[which(r==1)], data.frame(x=x[which(r!=1)]), ndpost = M)))
  
  fdraws <- plogis(cbind(post$yhat.train, post$yhat.test))
  
  #BART for propensity score model with logistic link
  invisible(capture.output(postpi <- lbart(x.train = data.frame(xfull), y.train = rfull, ndpost = M)))
  
  pidraws <- postpi$prob.train
  
  return(list("fdraws" = fdraws, "pidraws" = pidraws,"yfull" = yfull,"xfull"= xfull,"rfull"= rfull))
}