\name{my.auc}
\alias{my.auc}
\title{Area under the curve (AUC)}
\usage{
my.auc(linpred, Y)
}

\arguments{
\item{linpred}{n-vector giving the value of the marker.}

\item{Y}{n-vector giving the binary status, coded as 0/1.}
}

\value{
the area under the curve
}


\description{
computes the area under the ROC curve (AUC) for the marker 'linpred' and the binary status 'Y'.
}

