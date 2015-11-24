\name{cvr.glmnet}
\alias{cvr.glmnet}
\title{Repeating cv.glmnet}
\usage{
cvr.glmnet(X, Y, family, standardize=TRUE,alpha=1, nfolds, ncv, type.measure,...)
}
\arguments{
\item{X}{a (nxp) matrix of predictors with observations in rows and predictors in columns}

\item{Y}{n-vector giving the value of the response (either continuous, numeric-binary 0/1, or \code{Surv} object)}

\item{family}{should be "gaussian" for continuous \code{Y}, "binomial" for binary \code{Y}, "cox" for \code{Y} of type \code{Surv}}

\item{standardize}{whether the predictors should be standardized or not. Default is TRUE.}

\item{alpha}{the elastic net mixing parameter: \code{alpha}=1 yields the L1 penalty (lasso), \code{alpha}=0 yields the L2 penalty. Default is \code{alpha}=1 (lasso).}

\item{nfolds}{the number of folds of CV procedure.}

\item{ncv}{the number of repetitions of CV. Not to be confused with \code{nfolds}. For example, if one repeats 50 times 5-fold-CV (i.e. considers 50 random partitions into 5 folds in turn and averages the results), nfolds equals 5 and ncv equals 50.}

\item{type.measure}{The accuracy/error measure computed in cross-validation. If not specified, type.measure is "class" (classification error) if \code{family="binomial"}, "mse" (mean squared error) if \code{family="gaussian"} and partial likelihood if \code{family="cox"}. If \code{family="binomial"}, one may specify \code{type.measure="auc"} (area under the ROC curve).}

 \item{...}{Other arguments to be passed to the function \code{cv.glmnet}.}

}

\value{
A list with the following arguments:
  \item{coeff}{the matrix of coefficients with predictors corresponding to rows and lambda values corresponding to columns. The first rows contains the intercept of the model (for all families other than \code{"cox"}).}
  \item{lambda}{the lambda sequence}
 \item{cvm}{the CV estimate of the measure specified by \code{type.measure} for each candidate lambda value, averaged over the ncv runs of \code{cv.glmnet}}
 }


\description{
the same as cv.glmnet but with several ncv repetitions of CV: cross-validation repeated ncv times (i.e. for ncv different random partitions)
}


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

cvr.glmnet(X=X,Y=Y,family="binomial",standardize=FALSE,nfolds=5,ncv=10,type.measure="class")
}

