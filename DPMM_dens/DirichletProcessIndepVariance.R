library(dirichletprocess)
library(invgamma)


# Mixing Distribution Functions ####

#' Create a Gaussian Mixing Distribution with fixed variance.
#'
#'
#' @param priorParameters The prior parameters for the base measure.
#' @param sigma The fixed variance of the model.
#' @return A mixing distribution object.
#' @export
GaussianIndepVarianceMixtureCreate <- function(priorParameters=c(0,1,1,1), mhStepSize = 1) #mu0,sigma0,a,b
{
  mdobj <- MixingDistribution("normalIndepVariance",
                              priorParameters,
                              "conjugate")
  mdobj$sigma <- rinvgamma(1, priorParameters[3], priorParameters[4]) #draw sigma randomly from prior
  #mdobj$sigma <- 1 #initialize sigma = 1
  mdobj$mhStepSize <- mhStepSize
  return(mdobj)
}

#' @export
Likelihood.normalIndepVariance <- function(mdObj, x, theta) {
  as.numeric(dnorm(x, theta[[1]],  mdObj$sigma))
}

#' @export
logLikelihood.normalIndepVariance <- function(dpobj, x, sigma) {
  theta_n <- dpobj$clusterParameters[[1]][dpobj$clusterLabels] #get theta_{1:n} allocated to each cluster
  as.numeric(dnorm(x, theta_n,  sigma,log = TRUE))
}

#' @export
PriorDraw.normalIndepVariance <- function(mdObj, n = 1) {
  
  priorParameters <- mdObj$priorParameters
  
  mu <- rnorm(n, priorParameters[1], priorParameters[2]) #draw from normal with mu0, sigma0
  theta <- list(array(mu, dim = c(1, 1, n)))
  return(theta)
}

#' @export
PosteriorParameters.normalIndepVariance <- function(mdObj, x) {
  
  priorParameters <- mdObj$priorParameters
  
  n.x <- length(x)
  sigma <- mdObj$sigma
  mu0 <- priorParameters[1]
  sigma0 <- priorParameters[2] #this is now different to sigma
  
  sigmaPosterior <- (1/sigma0^2 + n.x/sigma^2) ^ (-1)
  muPosterior <- sigmaPosterior * (mu0/sigma0^2 + sum(x)/sigma^2)
  posteriorParameters <- matrix(c(muPosterior, sigmaPosterior), ncol=2)
  return(posteriorParameters)
}

#' @export
PosteriorDraw.normalIndepVariance <- function(mdObj, x, n = 1, ...) {
  
  posteriorParameters <- PosteriorParameters(mdObj, x)
  
  mu <- rnorm(n,
              posteriorParameters[1],
              posteriorParameters[2])
  theta <- list(array(mu, dim = c(1, 1, n)))
  return(theta)
}

#' @export
Predictive.normalIndepVariance <- function(mdObj, x) {
  priorParameters <- mdObj$priorParameters
  sigma0 <- priorParameters[[2]] #base distribution variance
  sigma <- mdObj$sigma
  
  predictiveArray <- numeric(length(x))
  
  for (i in seq_along(x)) {
    
    posteriorParameters <- PosteriorParameters(mdObj, x[i])
    
    predictiveArray[i] <- dnorm(x[i],
                                posteriorParameters[1],
                                sigma0^2 + sigma^2)
  }
  return(predictiveArray)
}

#' @export #custom plot function
PosteriorFunction.normalIndepVariance <- function (dpobj, ind) 
{
  post_clusters <- PosteriorClusters(dpobj, ind)
  dpobj$mixingDistribution$sigma <- dpobj$sigmaChain[ind] #set variance to appropriate posterior sample
  base_function <- function(x, theta) Likelihood(dpobj$mixingDistribution, 
                                                 x, theta)
  post_func <- weighted_function_generator(base_function, post_clusters$weights, 
                                           post_clusters$params)
  return(post_func)
}


# Dirichlet Process Functions ####

#' Create a Dirichlet Mixture of the Gaussian Distribution with fixed variance.
#'
#'
#' @param y Data.
#' @param sigma The fixed variance
#' @param g0Priors Base Distribution Priors.
#' @param alphaPriors Prior parameter distributions for the alpha concentration parameter.
#' @return Dirichlet process object
#'
#' @export
DirichletProcessGaussianIndepVariance <- function(y,
                                                  g0Priors = c(0,1,1,1), #mu0,sigma0, a,b
                                                  alpha = 1,
                                                  mhStepSize = 1) { #fixed alpha
  
  mdobj <- GaussianIndepVarianceMixtureCreate(g0Priors, mhStepSize =mhStepSize)
  dpobj <- DirichletProcessCreate(y, mdobj, mhDraws = 1) #not used
  dpobj <- Initialise(dpobj)
  dpobj$alpha <- alpha #hard-code alpha, don't update
  return(dpobj)
}


# Fitting Functions

#Metropolis-Hastings
MhParameterProposal.sigma <- function(sigma_old,mhStepSize){
  sigma_prop <- abs(sigma_old + mhStepSize*rnorm(1))
  return(sigma_prop)
}

PriorDensity.sigma <- function(mdobj,sigma){
  priorParameters <- mdobj$priorParameters
  sigmaDensity <- dinvgamma(sigma, priorParameters[3], priorParameters[4])
  return(as.numeric(sigmaDensity))
}

UpdateSigma <- function(dpObj){
  #load
  mixingDistribution <- dpObj$mixingDistribution
  x <- dpObj$data
  
  #old value
  sigma_old <- mixingDistribution$sigma
  #simulate proposal (folded normal, symmetric proposal)
  sigma_prop <- MhParameterProposal.sigma(sigma_old, mixingDistribution$mhStepSize)
  
  #acceptance ratio
  old_log_prior <- log(PriorDensity.sigma(mixingDistribution, sigma_old))
  old_log_Likelihood <- sum((logLikelihood.normalIndepVariance(dpObj, x, sigma_old)))
  
  new_log_prior <- log(PriorDensity.sigma(mixingDistribution, sigma_prop))
  new_log_Likelihood <- sum((logLikelihood.normalIndepVariance(dpObj, x, sigma_prop)))
  
  #accept/reject
  accept_prob <- min(1, exp(new_log_prior + new_log_Likelihood - old_log_prior - old_log_Likelihood))
  
  if (runif(1) < accept_prob) {
    dpObj$mixingDistribution$sigma <- sigma_prop
  } else {
    dpObj$mixingDistribution$sigma <- sigma_old
  }
  dpObj
}

#Run MCMC
Fit.gaussianIndepVariance <- function(dpObj, 
                                      its, 
                                      updatePrior = FALSE, 
                                      progressBar = interactive()){
  
  if (progressBar){
    pb <- txtProgressBar(min=0, max=its, width=50, char="-", style=3)
  }
  
  alphaChain <- numeric(its)
  likelihoodChain <- numeric(its)
  weightsChain <- vector("list", length = its)
  clusterParametersChain <- vector("list", length = its)
  priorParametersChain <- vector("list", length = its)
  labelsChain <- vector("list", length = its)
  sigmaChain <- numeric(its)
  
  for (i in seq_len(its)) {
    
    alphaChain[i] <- dpObj$alpha
    weightsChain[[i]] <- dpObj$pointsPerCluster / dpObj$n
    clusterParametersChain[[i]] <- dpObj$clusterParameters
    priorParametersChain[[i]] <- dpObj$mixingDistribution$priorParameters
    labelsChain[[i]] <- dpObj$clusterLabels
    sigmaChain[[i]] <- dpObj$mixingDistribution$sigma
    
    likelihoodChain[i] <- sum(log(LikelihoodDP(dpObj)))
    
    dpObj <- ClusterComponentUpdate(dpObj)
    dpObj <- ClusterParameterUpdate(dpObj)
    #dpObj <- UpdateAlpha(dpObj) #keep alpha fixed! no sampler 
    dpObj <- UpdateSigma(dpObj)

    if (updatePrior) {
      dpObj$mixingDistribution <- PriorParametersUpdate(dpObj$mixingDistribution,
                                                        dpObj$clusterParameters)
    }
    if (progressBar){
      setTxtProgressBar(pb, i)
    }
  }
  
  dpObj$weights <- dpObj$pointsPerCluster / dpObj$n
  dpObj$alphaChain <- alphaChain
  dpObj$likelihoodChain <- likelihoodChain
  dpObj$weightsChain <- weightsChain
  dpObj$clusterParametersChain <- clusterParametersChain
  dpObj$priorParametersChain <- priorParametersChain
  dpObj$labelsChain <- labelsChain
  dpObj$sigmaChain <- sigmaChain
  
  if (progressBar) {
    close(pb)
  }
  return(dpObj)
}

# Testing the Implementation (commented as we want to import this) ####

# testData <- c(rnorm(200, -1, 0.5), rnorm(200,1,0.5))
# #testData <- testData/sd(testData)
# 
# n_samp <- 3000
# dp <- DirichletProcessGaussianIndepVariance(testData,g0Priors = c(0,1,1,1),alpha = 1,mhStepSize = 0.5)
# dp <- Fit.gaussianIndepVariance(dp, n_samp, progressBar = TRUE)
# 
# 
# #Plot 
# ind = 1:(n_samp)
# y_plot = seq(-4,4,length.out = 100)
# pdf_grid_samp <- t(sapply(ind, function(i) PosteriorFunction.normalIndepVariance(dp,i)(y_plot)))
# hist(testData,pro = TRUE, ylim = c(0,0.45),nclass = 20) 
# lines(y_plot,colMeans(pdf_grid_samp),type = 'l')

# #For reference, standard location-scale mixture DPMM ####
# dp_fit <- DirichletProcessGaussian(testData)
# dp_fit <- Fit(dp_fit, n_samp)
# 
# #Plot
# ind = 1:(n_samp)
# y_plot = seq(-4,4,length.out = 100)
# pdf_grid_samp <- t(sapply(ind, function(i) PosteriorFunction(dp_fit,i)(y_plot)))
# hist(testData,pro = TRUE, ylim = c(0,0.45),nclass = 20) 
# lines(y_plot,colMeans(pdf_grid_samp),type = 'l')



