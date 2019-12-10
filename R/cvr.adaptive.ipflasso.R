### Adaptive IPF-Lasso
###
### Copyright (C) 2019 Gerhard Schulze
###
### 
###
###
### This file is part of the `ipflasso' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


cvr.adaptive.ipflasso <- function(X, Y, family, type.measure, standardize = TRUE,
                                  alpha, type.step1, blocks, nfolds, ncv)
{ 
  
  ulblocks <- as.numeric(unlist(blocks))
  if(!setequal(ulblocks, 1:ncol(X)))
  {
    stop("Each predictor should be included in exactly one block.")
  }
  
  if(alpha == 1 & missing(type.step1))
  {
    warning("type.step1 is set to sep")
    type.step1 <- "sep"
  }
  
  if(alpha == 0 & missing(type.step1))
  {
    warning("type.step1 is set to comb")
    type.step1 <- "comb"
  } 
  
  if(alpha > 0 & alpha < 1 & missing(type.step1))
  {
    stop("type.step1 should be sep or comb")
  } 
  
  if(family == "gaussian")
  {
    if(type.measure != "mse")
      warning("type.measure is set to mse.")
    type.measure <- "mse"
  }
  
  if(family == "cox")
  {
    if(type.measure != "deviance")
      warning("type.measure is set to partial likelihood.")
    type.measure <- "deviance"
  }
  
  if(family == "binomial" & !is.element(type.measure, c("auc","class")))
  {
    warning("type.measure is set to class")
    type.measure <- "class"
  }
  
  ###################################################################################
  #############################     Step 1      #####################################
  ###################################################################################
  
  ################## determine average coefficients for each modality ###############
  M <- length(blocks)  
  means <- c(rep(0, M))
  index <- numeric()
  for (i in 1:M)
  {
    index <-append(index, rep(i,length(blocks[[i]])))
  }        

  ################## separate step1 models ##########################################
  if (type.step1 == "sep")
  {
    for (i in 1:M)
    {
      assign(paste("model", i, sep = ""), cv.glmnet(x = X[,blocks[[i]]], y = Y, family = family,
                                                    type.measure = type.measure,
                                                    standardize = standardize,
                                                    nfolds = nfolds, alpha = alpha))
      if (family == "cox") 
        means[i] <- mean(abs(coef(get(paste("model", i, sep = "")), s = "lambda.min")[,1]))
      if (family == "gaussian" | family == "binomial") 
        means[i] <- mean(abs(coef(get(paste("model", i, sep = "")), s = "lambda.min")[-1,1]))
    }
  }
  ################## combined step1 model ############################################ 
  if (type.step1 == "comb")
  {
    model1 <- cv.glmnet(x = X, y = Y, family = family, type.measure = type.measure,
                        standardize = standardize, nfolds = nfolds, alpha = alpha)
    coffs <- as.matrix(abs(coef(model1, s = "lambda.min")))
    if (family == "cox") means <- as.numeric(tapply(coffs, index, mean)) 
    if (family == "gaussian" | family == "binomial") 
      means <- as.numeric(tapply(coffs[-1,], index, mean)) 
  }
  ################## check if there are average betas = 0 and determine #############
  ################## exclusion vector,second step blocks and pf #####################
  
  if (any(means != 0))          # at least one step 1 coefficient mean is > 0
  {
    means_check <- 1 / means  
    
    #################### (i) Exclusion vector #######################################
    
    exc <- NULL
    a = 1
    for (j in 1:M) 
    {
      if (is.infinite(means_check[j]))
      {
        exc <- append(exc, c(a : (a + length(blocks[[j]]) - 1)))
      }
      a = a + length(blocks[[j]])
    }
    
    ################### (ii) pf vector ##############################################
    
    pf <- means_check[is.finite(means_check)]
    
    ################### (iii) new blocks variable ###################################
    
    for (k in 1:M) 
    {
      if (is.infinite(means_check[k])) blocks[[paste("block", k, sep = "")]] <- NULL
    }
    a = 1
    for (l in 1:length(blocks)) 
    {
      blocks[[l]] <- a : (a + length(blocks[[l]]) - 1)
      a <- a + length(blocks[[l]])
    }
    
    
    #################################################################################
    #############################     Step 2      ###################################
    #################################################################################
    
    if (!is.null(exc))
    {
      res <- cvr.ipflasso(X = X[,-exc], Y = Y, family = family, type.measure = type.measure,
                          standardize = standardize, alpha = 1, blocks = blocks, pf = pf,
                          nfolds = nfolds, ncv = ncv)
    } else {
      
      res <- cvr.ipflasso(X = X, Y = Y, family = family, type.measure = type.measure,
                          standardize = standardize, alpha = 1, blocks = blocks, pf = pf,
                          nfolds = nfolds, ncv = ncv)
    }

    return(lapply(rapply(list(res, means.step1 = means, exc = as.integer(exc)), 
                         enquote, how = "unlist"), eval)) 
    
    
  } else {              # all step 1 coefficient means are = 0 
    
    if (type.step1 == "sep")
    {
      ind.bestlambda <- 1L 
      exc <- c(1:ncol(X))
      
      # calculate nzero and average cvm
      nzerolist <- list()
      cvmlist   <- list()
      for (i in 1:M)
      {
        nzerolist[[i]] <- get(paste("model", i, sep = ""))$nzero
        cvmlist[[i]]   <- get(paste("model", i, sep = ""))$cvm
      }
      nzeromat <- do.call(rbind, lapply(nzerolist, "[", 1:max(sapply(nzerolist, length))))
      colnames(nzeromat) <-  c(paste0("s", 0:(ncol(nzeromat)-1)))
      nzero  <- colSums(nzeromat, na.rm = TRUE)
      cvmmat <- do.call(rbind, lapply(cvmlist, "[", 1:max(sapply(cvmlist, length))))
      cvm    <- colMeans(cvmmat, na.rm = TRUE)
      
      # combine model coefficients
     if (family == "cox")
      {
        coeff <- as.numeric(rep(NA,length(nzero)))
        names(coeff) <- c(paste0("s", 0:(length(nzero)-1)))
       
         for (i in 1:M)
        {
          if (get(paste("model", i, sep = ""))$glmnet.fit$dim[2] < length(nzero))
          {  
          namat <- matrix(NA, nrow = get(paste("model", i, sep = ""))$glmnet.fit$dim[1], 
                          ncol = (length(nzero) - get(paste("model", i, sep = ""))$glmnet.fit$dim[2]))
          coeff <- rbind(coeff, cbind(as.matrix(get(paste("model", i, sep = ""))$glmnet.fit$beta), namat))
          
          } else {
            
           coeff <- rbind(coeff,as.matrix(get(paste("model", i, sep = ""))$glmnet.fit$beta))  
          }
        }
        rownames(coeff)[1] <- ""
      }
      
     if (family == "gaussian" | family == "binomial")
      {
        intlist <- list()
        for (i in 1:M)
        {
          intlist[[i]] <- get(paste("model", i, sep = ""))$glmnet.fit$a0
        }
        intmat <- do.call(rbind, lapply(intlist, "[", 1:max(sapply(intlist, length))))
        colnames(intmat) <-  c(paste0("s", 0:(ncol(intmat)-1)))
        coeff   <- colMeans(intmat, na.rm = TRUE)
        
        for (i in 1:M)
        {
          if (get(paste("model", i, sep = ""))$glmnet.fit$dim[2] < length(nzero))
          {  
          namat <- matrix(NA, nrow = get(paste("model", i, sep = ""))$glmnet.fit$dim[1], 
                         ncol = (length(nzero) - get(paste("model", i, sep = ""))$glmnet.fit$dim[2]))
          coeff <- rbind(coeff, cbind(as.matrix(get(paste("model", i, sep = ""))$glmnet.fit$beta), namat))
          
          } else {
            
           coeff <- rbind(coeff,as.matrix(get(paste("model", i, sep = ""))$glmnet.fit$beta))  
          }
        }
        rownames(coeff)[1] <- "intercept"
      }
      
      # determine lambda sequence
      lambda1 <- c(rep(0, M))
      for (i in 1:M)
      {
        lambda1[i] <- get(paste("model", i, sep = ""))$lambda[1]
      }
      indmaxl <- which.max(lambda1)
      lambda <- get(paste("model", indmaxl, sep = ""))$lambda

      warning("All step 1 model coefficients are zero! Step 2 not performed.")
      return(list(coeff = coeff, ind.bestlambda = ind.bestlambda, lambda = lambda, cvm = cvm, 
                  nzero = nzero, family = family, means.step1 = means, exc = exc))
    }
    
    if (type.step1 == "comb")
    {
      coeff <- rbind(as.numeric(model1$glmnet.fit$a0)[1:length(model1$lambda)],
                     as.matrix(model1$glmnet.fit$beta)[,1:length(model1$lambda)])
      exc <- c(1:ncol(X))
      if (family == "cox")
      {
        rownames(coeff)[1] <- " " 
      }
      if (family != "cox")
      {
        rownames(coeff)[1] <- "intercept" 
      }
      if (type.measure == "auc")
      {
        ind.bestlambda <- which.max(model1$cvm)
      }
      if (type.measure != "auc")
      {
        ind.bestlambda <- which.min(model1$cvm)
      }

      warning("All step 1 model coefficients are zero! Step 2 not performed.") 
      return(list(coeff = coeff, ind.bestlambda = ind.bestlambda, lambda = model1$lambda, cvm = model1$cvm, 
                  nzero = model1$nzero, family = family, means.step1 = means, exc = exc)) 
    }  
  }
}    #end function
