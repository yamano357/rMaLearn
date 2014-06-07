\docType{class}
\name{OnlineLearning-class}
\alias{OnlineLearning}
\alias{OnlineLearning-class}
\title{OnlineLearning (title)}
\description{
  Online learning reference class (description)
}
\details{
  \code{OnlineLearning} class is implement a online
  learning algorithm class that has methods for
  classification and Regression based on PA and AROW
  (details).
}
\section{Fields}{
  \describe{ \item{\code{Weight}: }{ the result of weight
  vector by training.} \item{\code{Param}: }{ setting
  parameter list according to online learning algorithm.}
  \item{\code{Feature}: }{ feature vector at last.}
  \item{\code{Response}: }{ response value at last.}
  \item{\code{LossFunc}: }{ online learning algorithm
  function apply to update the weight vector.} }
}

\section{Contains}{
  NULL
}

\section{Methods}{
  \describe{ \item{\code{updateWeight(Response,
  FeatureVec)}: }{ \code{updateWeight} method call
  \code{updatePassiveAggressive} or \code{updateAROW} by
  \code{Param}.  }
  \item{\code{updatePassiveAggressive(Margin = Margin, Norm
  = Norm)}: }{ \code{updatePassiveAggressive} method return
  \code{Weight} based on PA's algorithm.  }
  \item{\code{updateAROW(Margin = Margin, Norm = Norm)}: }{
  \code{updateAROW} method return \code{Weight} and
  \code{Param} based on AROW's algorithm.  } }
}
\keyword{documentation}

