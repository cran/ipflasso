\name{cvr.adaptive.ipflasso}
\alias{cvr.adaptive.ipflasso}
\title{Cross-validated integrative lasso with adaptive penalty factors}
\description{
Runs cvr.ipflasso applying different data based penalty factors to predictors from different blocks.
}
\usage{
cvr.adaptive.ipflasso(X, Y, family, type.measure, standardize = TRUE,
                                  alpha, type.step1, blocks, nfolds, ncv)
}
\arguments{
  \item{X}{a (nxp) matrix of predictors with observations in rows and predictors in columns.}
  
  \item{Y}{n-vector giving the value of the response (either continuous, numeric-binary 0/1, or \code{Surv} object).}

  \item{family}{should be "gaussian" for continuous \code{Y}, "binomial" for binary \code{Y}, "cox" for \code{Y} of type \code{Surv}.}

  \item{type.measure}{the accuracy/error measure computed in cross-validation. If not specified, type.measure is "class" (classification error) if \code{family="binomial"}, "mse" (mean squared error) if \code{family="gaussian"} and partial likelihood if \code{family="cox"}. If \code{family="binomial"}, one may specify \code{type.measure="auc"} (area under the ROC curve).}

  \item{standardize}{whether the predictors should be standardized or not. Default is TRUE.}
 
  \item{alpha}{the elastic net mixing parameter for step 1: \code{alpha}=1 yields the L1 penalty (Lasso), \code{alpha}=0 yields the L2 penalty (Ridge).}

\item{type.step1}{whether the models of step 1 should be run on the whole data set \code{X} (\code{type.step1="comb"}) or separately for each block (\code{type.step1="sep"}).}
  
  \item{blocks}{a list of length M of the format \code{list(block1=...,block2=...,} where the dots should be replaced by the indices of the predictors included in this block. The blocks should form a partition of 1:p.} 

  \item{nfolds}{the number of folds of the CV procedure.}

  \item{ncv}{the number of repetitions of the CV. Not to be confused with \code{nfolds}. For example, if one repeats 50 times 5-fold-CV (i.e. considers 50 random partitions into 5 folds in turn and averages the results), \code{nfolds} equals 5 and \code{ncv} equals 50.}
}

\details{The penalty factors are the inverse arithmetic means of the absolute model coefficients per block, generated in a first step of the function. The user can choose to determine these coefficients by running a Lasso model (\code{alpha=1}) or a Ridge model (\code{alpha=0}) either on the whole data set (\code{type.step1="comb"}) or seperately for each block (\code{type.step1="sep"}). If \code{type.step1} is ommited, it will be set to \code{"sep"} for Lasso and to \code{"comb"} for Ridge. 
If a Lasso model in step 1 returns any zero coefficient mean, the corresponding block will be excluded from the input date set \code{X} and step 2 will be run with the remaining blocks. If all model coefficient means are zero, step 2 will not be performed.
}

\value{
A list with the following arguments:
  \item{coeff}{the matrix of coefficients with predictors corresponding to rows and lambda values corresponding to columns. The first row contains the intercept of the models (for all families other than \code{"cox"}). 
  
In the special case of separate step 1 Lasso models and all coefficient means equal to zero, the intercept is the average of the separate model intercepts per block.}
  \item{ind.bestlambda}{the index of the best lambda according to CV.}
  \item{lambda}{the lambda sequence.
  
In the special case of separate step 1 Lasso models and all coefficient means equal to zero, it is the lambda sequence with the highest lambda value among the lambda sequences of all blocks.}
  \item{cvm}{the CV estimate of the measure specified by \code{type.measure} for each candidate lambda value.
  
In the special case of separate step 1 Lasso models and all coefficient means equal to zero, cmv is the average of the separate model cvms per block.}
  \item{nzero}{the number of non-zero coefficients in the selected model.

In the special case of separate step 1 Lasso models and all coefficient means equal to zero, nzero is the sum of the non-zero coefficients of the separate models per block.}
  \item{family}{see arguments.}
  \item{means.step1}{the arithmetic means of the absolute model coefficients per block, returned by the first step of the function.}
  \item{exc}{the exclusion vector containing the indices of the block(s) to be excluded from \code{X}.}
}

\references{
 Schulze, Gerhard (2017): Clinical Outcome Prediction Based on Multi-Omics Data: Extension of IPF-LASSO. Masterarbeit, Ludwig-Maximilians-Universitaet Muenchen
(Department of Statistics: Technical Reports) https://doi.org/10.5282/ubm/epub.59092
}

\author{
Gerhard Schulze (g-schulze@t-online.de)
}

\examples{
# load ipflasso library
library(ipflasso)

# generate dummy data
X<-matrix(rnorm(50*200),50,200)
Y<-rbinom(50,1,0.5)

cvr.adaptive.ipflasso(X=X,Y=Y,family="binomial",type.measure="class",standardize=FALSE,
                      alpha = 1,blocks=list(block1=1:50,block2=51:200),nfolds=5,ncv=10)
}
