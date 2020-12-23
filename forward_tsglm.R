#
# Function "forward_tsglm" implements (in R) the Block-forward (BF) approach with 
# INGARCH(p,q) models availabale in tscount R package.
#
#  Input:
#
#   * Y = list(response, model, link, distr) where
#
#         response is a vector with the time series of counts, 
#         model is a list with info on past_obs and past_men
#         link is a string with the linkage function ("identity" or "log")
#         distr is a string with the conditional distribution ("poisson" or "nbinom")
#
#         see help of tsglm function for further details
#
#    * covar is the list with the blocks structure for the BF approach
#   
#         In this application covar is a list with 6 blocks:
#         {Temp}, {PM25, PM10}, {NOx, NO2}, {O3}, {SO2}, {CO}. The blocks can
#         be customised, namely they can be changed, removed and/or reordered 
#         for other applications. The number of sub-blocks within a block can 
#         also be altered, e.g., {PM25, PM10, NewSubBlock}. 
#         See S2 File for practical examples on how to manipulate the blocks.
# 
# 
#
# Development
#   Ana Martins (a.r.martins@ua.pt)
#   University of Aveiro, Portugal, December 2020
#   R version 3.6.2; tscount version 1.4.2;
#
# License 
# GNU General Public License, either Version 2 or Version 3
#
# Martins A, Scotto M, Deus R, Monteiro A, Gouveia S (2021) 
# Association between respiratory hospital admissions and air pollution in Portugal: a count time series approach
#
#

forward_tsglm = function(Y = list(response, model, link, distr), 
                         covar = list(temp =  NULL, 
                                      PM = list(PM25 = NULL, PM10 = NULL), 
                                      NO = list(NOx = NULL, NO2 = NULL), 
                                      O3 =  NULL, 
                                      SO2 =  NULL, 
                                      CO =  NULL)){
  #needed for data manipulation
  library(data.table) 
  library(tscount)

  #helper function for covar preprocessing (remove empty sub(blocks))
  processCovarList <- function(covar){
    
    ##helper function to whether an object is either NULL or a list of NULLs
    is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))
    
    ## Recursively step down into list, removing NULL objets
    rmNullObs <- function(x) {
      x <- Filter(Negate(is.NullOb), x)
      lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
    }
    
    processCovar = rmNullObs(covar)
    
    library(plotrix)
    #if list depth > 1, unlist a level
    for(i in 1:length(processCovar)){
      if(listDepth(processCovar[[i]])> 1){
        processCovar[[i]] <- unlist(processCovar[[i]], recursive = FALSE)
      }
    }
    
    return(covar = processCovar)
  }
 
  #helper function to test for significance
  signif = function(ci){
    if(is.matrix(ci)==TRUE){
      signi = vector()
      for(i in 1:dim(ci)[1]){
        a = !0 %between% c(ci[i,1],ci[i,2])
        signi[length(signi)+1] = a  
      }
    }  
    else{
      signi = !0 %between% c(ci[1],ci[2])
    }
    return(signi)  
  }

  #helper function to fit tsglm models
  fit_model = function(Y, X = NULL){
    
    #if X is empty = baseline model
    if(is.null(X)){
      f = tsglm(ts = Y$response, model = Y$model , link = Y$link, distr = Y$distr)
      #compute CIs
      CIs = confint(f)
    }
    else{
      f = tsglm(ts = Y$response, model = Y$model , link = Y$link, distr = Y$distr, xreg = cbind(X))
      
      #compute CIs only for xreg
      n_coef = sum(1 + length(f$model$past_obs) + length(f$model$past_mean)) # 1 = interecept, total number of coefs
      n_xreg = length(f$model$external) #total number of x_reg coefs
      CIs = confint(f)
      CIs = CIs[(n_coef+1):(n_coef+n_xreg), ]
    } 
    
    #compute AIC
    AIC = -2*f$logLik+2*(length(f$coefficients)+1)
    #logical vector indicating significance, TRUE = significant, FALSE = not significant
    logi_signif = signif(CIs)
    
    model = list(f = f, AIC = AIC, CI = CIs, significance = logi_signif)
    
    return(model)
  }
  
  
  covar = processCovarList(covar)

  X = NULL
  AIC = vector()
  p_val = vector()
  L <- length(covar) #number of blocks
  
  current_model = fit_model(Y = Y, X = X) #baseline model, pure INGARCH
  current_AIC = current_model$AIC
  current_pval = current_model$significance
  for(i in 1:L){
    N <- ncol(as.data.frame(covar[[i]])) #number of variables within each block
    j_min <- 0 
    for(j in 1:N){
      if(N==1){ # if the block only has one covar
        Z <- as.data.frame(covar[[i]])
        
        names(Z) <- names(covar[[i]])
        if(is.null(names(Z))){ #ensures that name of covariate being tested passes through
          names(Z) <- names(covar)[i]
        }
        candidate_model <- fit_model(Y, X = cbind(X, as.matrix(Z)))
      }
      else if (N==0){
        print(paste(names(covar)[[i]], 'has dimension zero. Please check your input.'))
        break
      }
      else{
        Z <- as.matrix(as.data.frame(covar[[i]][j]))
        candidate_model <- fit_model(Y, X = cbind(X, Z)) 
      }
      candidate_AIC = candidate_model$AIC
      candidate_pval = candidate_model$significance
      if(current_pval[length(current_pval)]==TRUE){ #if candidate regressor is significant
        if(candidate_AIC < current_AIC & all(candidate_pval)==TRUE){ #if previous regressors are also significant and the AIC is lower
          current_model <- candidate_model
          current_AIC <- candidate_AIC
          j_min <- j
        
        }
      }
    }
    if(j_min > 0){
      if(N==1){
      X <- cbind(X, as.matrix(as.data.frame(covar[[i]])))
      colnames(X)[ncol(X)] = names(covar)[i]

      }
      else{
      X <- cbind(X, as.matrix(as.data.frame(covar[[i]][j_min])))
      }
    }

  }

  return(model = current_model$f)

}


