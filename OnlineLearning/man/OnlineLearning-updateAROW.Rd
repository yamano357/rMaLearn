\docType{methods}
\name{updateAROW,OnlineLearning-method}
\alias{updateAROW}
\alias{updateAROW,OnlineLearning-method}
\title{The \code{updateAROW} method of \code{OnlineLearning} class}
\arguments{
  \item{Margin}{\code{Response} * \code{Weight} *
  \code{Feature}}

  \item{Norm}{sum(\code{Feature} ** 2)}
}
\description{
  \code{updateAROW} method
}
\details{
  \code{updateAROW} method updates \code{Weight} based on
  AROW(Adaptive Regularization Of Weight Vector).
}
\examples{
library(Matrix)
library(kernlab)
data(spam)
InsY <- matrix(ifelse(as.integer(spam[, "type"]) == 1, 1, -1), ncol = 1)
InsX <- spam[, -ncol(spam)]
InstanceOL <- OnlineLearning$new(
               ncol(InsX),
               list("Mode" = "Class", "Algo" = "AROW", "Gamma" = 0.5, "IsFull" = FALSE)
              )
print(InstanceOL$Weight) # identical numeric(ncol(InsX))
# set values
InstanceOL$Response <- InsY[1]
InstanceOL$Feature <- as.numeric(InsX[1, ])
Margin <- InstanceOL$Response * as.numeric(InstanceOL$Weight \%*\% InstanceOL$Feature)
Norm <- as.numeric(crossprod(InstanceOL$Feature))
# update
InstanceOL$updateAROW(Margin = Margin, Norm = Norm)
print(InstanceOL$Weight)
}

