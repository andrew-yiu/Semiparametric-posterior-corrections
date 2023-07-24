bb <- function(fdraw0, fdraw1,r,n){
  wts <- rdirichlet(n)*n
  n_treated <- sum(wts*r)
  theta_boot <- sum(wts*r*(fdraw1- fdraw0))/n_treated
  return(theta_boot)
}

onestep <- function(fdraw0,fdraw1,pi,r,y,n){
  wts <- rdirichlet(n)
  # wts_x_att <- rdirichlet(n)
  
  #Compute probability of T=1 and plug-in
  pr1_sactt <- mean(pi)
  theta_plug_sactt <- mean(pi*(fdraw1- fdraw0))/pr1_sactt
  theta_plug_catt <- sum(r*(fdraw1- fdraw0))/sum(r)
  
  # pr1_att <- sum(wts_x_att*pi)
  
  # compute probability of treatment under tilde P
  ptil_1 <- sum(wts*r)
  
  # theta_plug_att <- sum(wts_x_att*pi*(fdraw1- fdraw0))/pr1_att
  
  # pr1_orac <- mean(pitrue)
  # theta_plug_orac <- mean(pitrue*(f1true- f0true))/pr1_orac
  
  #Compute bias and update
  # eif_sactt <- (r/pr1_sactt - ((1-r)/pr1_sactt)*(pi/(1-pi)))*(y- (r*fdraw1 + (1-r)*fdraw0)) +
  #   ((r-pi)/pr1_sactt)*(fdraw1 - fdraw0 - theta_plug_sactt)
  # theta_sactt = theta_plug_sactt + sum(wts*eif_sactt)
  
  eif_sactt <- (r/ptil_1 - ((1-r)/ptil_1)*(pi/(1-pi)))*(y- (r*fdraw1 + (1-r)*fdraw0)) +
    ((r-pi)/ptil_1)*(fdraw1 - fdraw0 - theta_plug_sactt)
  theta_sactt = theta_plug_sactt + sum(wts*eif_sactt)
  
  # eif_att <- (r/ptil_1 - ((1-r)/ptil_1)*(pi/(1-pi)))*(y- (r*fdraw1 + (1-r)*fdraw0)) +
  #   (r/ptil_1)*(fdraw1 - fdraw0 - theta_plug_att)
  # theta_att = theta_plug_att + sum(wts*eif_att)
  
  att_cor <- (r/ptil_1 - ((1-r)/ptil_1)*(pi/(1-pi)))*(y- (r*fdraw1 + (1-r)*fdraw0)) +
    (r/ptil_1)*(fdraw1 - fdraw0)
  theta_att = sum(wts*att_cor)
  # 
  # eif_orac <- (r/pr1_orac - ((1-r)/pr1_orac)*(pitrue/(1-pitrue)))*(y- (r*f1true + (1-r)*f0true)) +
  # ((r-pitrue)/pr1_orac)*(f1true - f0true - theta_plug_orac)
  # theta_orac <- theta_plug_orac + sum(wts*eif_orac)

  return(c(theta_sactt, theta_att, theta_plug_sactt, theta_plug_catt))
}

# plugin BART samples with empirical propensity score model
bart_emp <- function(fdraw0,fdraw1, r) {
  theta_emp <- sum(r*(fdraw1- fdraw0))/sum(r)
  return(theta_emp)
}

# 
# #weighted TMLE
# #fluctuated negative log likelihood
# nll_tmle <- function(t, wts,fdraw,pi,r,y){
#   m_ <- expit(logit(fdraw) + t/pi)
#   m_ <- pmax(1e-6, pmin(m_,1-1e-6)) #clip m_ to lie between 0 and 1
#   nll <- -sum((wts*r/pi)*(y*log(m_) + (1-y)*log(1-m_)))
#   return(nll)
# }
# 
# #gradient of the above
# gr_nll_tmle <- function(t,wts,fdraw,pi,r,y){
#   m_ <- expit(logit(fdraw) + t/pi)
#   gr_nll <- -sum((wts*r/pi)*(y-m_))
#   return(gr_nll)
# }
# 
# #compute TMLE estimate
# tmle <- function(fdraw,pi,r,y,n){
#   wts <- rdirichlet(n)
#   opt <- optim(par = 0, nll_tmle, gr = gr_nll_tmle,wts,fdraw,pi,r,y, method = "L-BFGS-B")
#   t <- opt$par
#   m_ <- expit(logit(fdraw) + t/pi)
#   wts_x <- rdirichlet(n)
#   theta_tilde <- sum(wts_x*m_) #posterior over x as well - do I need independent weights?
#   return(theta_tilde)
# }