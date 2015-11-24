\name{cvr.ipflasso}
\alias{cvr.ipflasso}
\title{Cross-validated integrative lasso with fixed penalty factors}
\usage{
cvr.ipflasso(X, Y, family, type.measure, standardize=TRUE, alpha=1, blocks, pf, nfolds,
  ncv)
}
\arguments{
\item{X}{a (nxp) matrix of predictors with observations in rows and predictors in columns}

\item{Y}{n-vector giving the value of the response (either continuous, numeric-binary 0/1, or \code{Surv} object)}

\item{family}{should be "gaussian" for continuous \code{Y}, "binomial" for binary \code{Y}, "cox" for \code{Y} of type \code{Surv}}

\item{type.measure}{The accuracy/error measure computed in cross-validation. If not specified, type.measure is "class" (classification error) if \code{family="binomial"}, "mse" (mean squared error) if \code{family="gaussian"} and partial likelihood if \code{family="cox"}. If \code{family="binomial"}, one may specify \code{type.measure="auc"} (area under the ROC curve).}


\item{standardize}{whether the predictors should be standardized or not. Default is TRUE.}

\item{alpha}{the elastic net mixing parameter: \code{alpha}=1 yields the L1 penalty (lasso), \code{alpha}=0 yields the L2 penalty. Default is \code{alpha}=1 (lasso).}

\item{blocks}{a list of length M the format \code{list(block1=...,block2=...,} where the dots should be replaced by the indices of the predictors included in this block. The blocks should form a partition of 1:p.}

\item{pf}{a vector of length equal to the number of blocks M. Each entry contains the penalty factor to be applied to the predictors of the corresponding block. Example: if \code{pf}=c(1,2), the penalty applied to the predictors of the 2nd block is twice as large as the penalty applied to the predictors of the first block.}


\item{nfolds}{the number of folds of CV procedure.}

\item{ncv}{the number of repetitions of CV. Not to be confused with \code{nfolds}. For example, if one repeats 50 times 5-fold-CV (i.e. considers 50 random partitions into 5 folds in turn and averages the results), \code{nfolds} equals 5 and \code{ncv} equals 50.}
}


\value{
A list with the following arguments:
  \item{coeff}{the matrix of coefficients with predictors corresponding to rows and lambda values corresponding to columns. The first rows contains the intercept of the model (for all families other than \code{"cox"}). }
  \item{ind.bestlambda}{the index of the best lambda according to CV.}
 \item{lambda}{the lambda sequence.}
  \item{cvm}{the CV estimate of the measure specified by \code{type.measure} for each candidate lambda value.}
 \item{nzero}{the number of non-zero coefficients in the selected model.}
 \item{family}{See arguments.}
 }


\description{
Runs cvr.glmnet giving different penalty factors to predictors from different blocks.}


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

cvr.ipflasso(X=X,Y=Y,family="binomial",standardize=FALSE,
            blocks=list(block1=1:50,block2=51:200), 
            pf=c(1,2),nfolds=5,ncv=10,type.measure="class")
}

