\name{cvr2.ipflasso}
\alias{cvr2.ipflasso}
\title{Cross-validated integrative lasso with cross-validated penalty factors}
\usage{
cvr2.ipflasso(X, Y, family, type.measure, standardize=TRUE, 
              alpha=1, blocks, pflist, nfolds, ncv, 
              nzeromax = +Inf, plot=FALSE)
}
\arguments{
\item{X}{a (nxp) matrix of predictors with observations in rows and predictors in columns}

\item{Y}{n-vector giving the value of the response (either continuous, numeric-binary 0/1, or \code{Surv} object)}

\item{family}{should be "gaussian" for continuous \code{Y}, "binomial" for binary \code{Y}, "cox" for \code{Y} of type \code{Surv}}

\item{type.measure}{The accuracy/error measure computed in cross-validation. If not specified, type.measure is "class" (classification error) if \code{family="binomial"}, "mse" (mean squared error) if \code{family="gaussian"} and partial likelihood if \code{family="cox"}. If \code{family="binomial"}, one may specify \code{type.measure="auc"} (area under the ROC curve).}


\item{standardize}{whether the predictors should be standardized or not. Default is TRUE.}

\item{alpha}{the elastic net mixing parameter: \code{alpha}=1 yields the L1 penalty (lasso), \code{alpha}=0 yields the L2 penalty. Default is \code{alpha}=1 (lasso).}

\item{blocks}{a list of length M the format \code{list(block1=...,block2=...,} where the dots should be replaced by the indices of the predictors included in this block. The blocks should form a partition of 1:p.}

\item{pflist}{a list of candidate penalty factors (see the argument \code{pf} of the function \code{\link{cvr.ipflasso}}) of the format weightslist=list(c(1,1),c(1,2),c(2,1),...).}


\item{nfolds}{the number of folds of CV procedure.}

\item{ncv}{the number of repetitions of CV. Not to be confused with \code{nfolds}. For example, if one repeats 50 times 5-fold-CV (i.e. considers 50 random partitions into 5 folds in turn and averages the results), \code{nfolds} equals 5 and \code{ncv} equals 50.}

\item{nzeromax}{the maximal number of predictors allowed in the final model. Default is +Inf, i.e. the best model is selected based on CV without restriction.}

\item{plot}{If \code{plot=TRUE}, the function outputs plots of CV errors and number of included predictors for each block.}
}


\value{
A list with the following arguments:
  \item{coeff}{the matrix of coefficients obtained with the best combination of penalty factors, with covariates corresponding to rows and lambda values corresponding to columns. The first row contains the intercept of the model.}
\item{ind.bestlambda}{the index of the best lambda as selected by CV for the best combination of penalty factors.}
\item{bestlambda}{the best lambda as selected by CV for the best combination of penalty factors.}
\item{ind.bestpf}{the index of the best penalty factor selected by CV from the list of candidates \code{pflist}.}
  \item{cvm}{the CV error for each candidate lambda value, averaged over the ncv runs of \code{cv.glmnet}.}
 \item{a}{a list of length \code{length(pflist)} containing the outputs of the function \code{\link{cvr.ipflasso}} for all candidate penalty factors from \code{pflist}.}
 \item{family}{See arguments.}
 }


\description{
Runs cvr.glmnet giving different penalty factors to predictors from different blocks and chooses the penalty factors by cross-validation from the list \code{pflist} of candidates.}


\references{
 Boulesteix AL, De Bin R, Jiang X, Fuchs M, 2015. IPF-lasso: integrative L1-penalized regression with penalty factors for prediction based on multi-omics data. Technical Report 187, Department of Statistics, University of Munich.
}

\author{
  Anne-Laure Boulesteix (\url{http://www.ibe.med.uni-muenchen.de/organisation/mitarbeiter/020_professuren/boulesteix/eng.html}),   

}


\examples{
# load ipflasso library
library(ipflasso)

# generate dummy data
X<-matrix(rnorm(50*200),50,200)
Y<-rbinom(50,1,0.5)

cvr2.ipflasso(X=X,Y=Y,family="binomial",type.measure="class",standardize=FALSE,
              blocks=list(block1=1:50,block2=51:200),
              pflist=list(c(1,1),c(1,2),c(2,1)),nfolds=5,ncv=10)
}
