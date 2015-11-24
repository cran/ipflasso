\name{ipflasso.predict}
\alias{ipflasso.predict}
\title{Using an IPF-lasso model for prediction of new observations}
\usage{
ipflasso.predict(object, Xtest)
}

\arguments{
\item{object}{the output of either \code{\link{cvr.ipflasso}} (if the user chooses the penalty factor himself) or \code{\link{cvr2.ipflasso}} (if the user cross-validates the penalty factor).}
\item{Xtest}{a ntest x p matrix containing the values of the predictors for the test data. It should have the same number of columns as the matrix \code{X} used to obtain the model \code{result}.}
}


\value{
A list with the following arguments:
  \item{linpredtest}{a ntest-vector giving the value of the linear predictor for the test observations}
  \item{classtest}{a ntest-vector with values 0 or 1 giving the predicted class for the test observations (for binary Y).}
  \item{probabilitiestest}{a ntest-vector giving the predicted probability of Y=1 for the test observations (for binary Y).}
 }


\description{
Derives predictions for new observations from a model fitted by the functions \code{\link{cvr.ipflasso}} or \code{\link{cvr2.ipflasso}}.}


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
Xtest<-matrix(rnorm(20*200),20,200)
Y<-rbinom(50,1,0.5)

# fitting the IPF-lasso model
model1<-cvr.ipflasso(X=X,Y=Y,family="binomial",standardize=FALSE,
                    blocks=list(block1=1:50,block2=51:200),
                    pf=c(1,2),nfolds=5,ncv=10,type.measure="class")

# making predictions from Xtest
ipflasso.predict(object=model1,Xtest=Xtest)
}

