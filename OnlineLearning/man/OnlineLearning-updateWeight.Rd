\docType{methods}
\name{updateWeight,OnlineLearning-method}
\alias{updateWeight}
\alias{updateWeight,OnlineLearning-method}
\title{The \code{updateWeight} method of \code{OnlineLearning} class}
\arguments{
  \item{Response}{integer number of a class label and
  FeatureVec numeric vector consist of feature values.}
}
\description{
  \code{updateWeight} method
}
\details{
  \code{updateWeight} method call
  \code{updatePassiveAggressive} or \code{updateAROW} using
  \code{Param}.
}
\examples{
library(Matrix)
library(kernlab)
library(iterators)
# create training data
data(spam)
InsY <- matrix(ifelse(as.integer(spam[, "type"]) == 1, 1, -1), ncol = 1)
InsX <- spam[, -ncol(spam)]
InstanceOL <- OnlineLearning$new(
               ncol(InsX),
               list("Mode" = "Class", "Algo" = "AROW", "Gamma" = 0.5, "IsFull" = FALSE)
              )
set.seed(10)
ShuffleIter <- iter(sample(length(InsY), length(InsY), replace = FALSE))
# train
while(all(tryCatch(i <- nextElem(ShuffleIter), error = function(Cond){ return (FALSE)}) != FALSE)){
  InstanceOL$updateWeight(Response = InsY[i], FeatureVec = as.numeric(InsX[i, ]))
}
# test
Predict <- sign(apply(InsX, 1, FUN = "\%*\%", InstanceOL$Weight))
table(Predict, as.numeric(InsY))
}

