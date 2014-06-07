\docType{methods}
\name{updatePassiveAggressive,OnlineLearning-method}
\alias{updatePassiveAggressive}
\alias{updatePassiveAggressive,OnlineLearning-method}
\title{The \code{updatePassiveAggressive} method of \code{OnlineLearning} class}
\arguments{
  \item{Margin}{\code{Response} * \code{Weight} *
  \code{Feature}}

  \item{Norm}{sum(\code{Feature} ** 2)}
}
\description{
  \code{updatePassiveAggressive} method
}
\details{
  \code{updatePassiveAggressive} method updates
  \code{Weight} based on Online Passive aggressive(this
  method is implement by PA-2).
}
\examples{
library(Matrix)
library(kernlab)
# create training data
data(spam)
InsY <- matrix(ifelse(as.integer(spam[, "type"]) == 1, 1, -1), ncol = 1)
InsX <- spam[, -ncol(spam)]
# make instance
InstanceOL <- OnlineLearning$new(
               ncol(InsX),
               list("Mode" = "Class", "Algo" = "PA", "C" = 1)
              )
print(InstanceOL$Weight) # identical numeric(ncol(InsX))
# set values
InstanceOL$Response <- InsY[1]
InstanceOL$Feature <- as.numeric(InsX[1, ])
Margin <- InstanceOL$Response * as.numeric(InstanceOL$Weight \%*\% InstanceOL$Feature)
Norm <- as.numeric(crossprod(InstanceOL$Feature))
# update
InstanceOL$updatePassiveAggressive(Margin = Margin, Norm = Norm)
print(InstanceOL$Weight)
}

