% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BSplineMI.R
\name{calcWeights}
\alias{calcWeights}
\title{Calculate bin weights using B-Spline}
\usage{
calcWeights(x, nBins, splineOrder)
}
\arguments{
\item{x}{Gene expression matrix with genes in rows and samples in columns}

\item{nBins}{Number of bins}

\item{splineOrder}{Spline order}
}
\value{
3-dimensional array of weights with dimensions: nBins x samples x genes
}
\description{
Takes an expression matrix and calculates for each value the weights for each
bin using B-Splines of a given order.
}
\examples{
x <- matrix(0:100,nrow=1)
w <- calcWeights(x,nBins = 5,splineOrder = 3)
plot(NULL,xlim=c(0,100),ylim=c(0,1),xlab="x",ylab="weight")
for(i in 1:dim(w)[1]){
  lines(x[1, ],w[i, ,1],col=rainbow(dim(w)[1])[i])
}
}
