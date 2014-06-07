#' OLRefclassPkg (title)
#'
#' A online learning package for reference class documents using 'roxygen2' (description)
#'
#' This package is implement a online learning algorithm (PA, AROW) by reference class methods. (details)
#'
#' @name OLRefclassPkg-package
#' @aliases OLRefclassPkg
#' @rdname OnlineLerning
#' @docType package
#' @keywords documentation
#'
#' @author
#' Author: Yoshiaki AMANO. \cr
#' Maintainer: Yoshiaki AMANO. \cr
#' E-Mail: \email{yamano357@@gmail.com}
#' @seealso \code{\link[kernlab:onlearn]{onlearn}}
#' @import compiler
#' @exportPattern '^[^\\.]'
NULL


#' OnlineLearning (title)
#'
#' Online learning reference class (description)
#'
#' \code{OnlineLearning} class is implement a online learning algorithm class that has methods for classification and Regression based on PA and AROW (details).
#'
#' @name OnlineLearning-class
#' @aliases OnlineLearning
#' @docType class
#'
#' @section Fields:
#' \describe{
#' \item{\code{Weight}: }{ the result of weight vector by training.} 
#' \item{\code{Param}: }{ setting parameter list according to online learning algorithm.} 
#' \item{\code{Feature}: }{ feature vector at last.} 
#' \item{\code{Response}: }{ response value at last.} 
#' \item{\code{LossFunc}: }{ online learning algorithm function apply to update the weight vector.} 
#' }
#' @section Contains:
#' NULL
#' @section Methods:
#' \describe{
#'  \item{\code{updateWeight(Response, FeatureVec)}: }{
#'    \code{updateWeight} method call \code{updatePassiveAggressive} or \code{updateAROW} by \code{Param}.
#'  }
#'  \item{\code{updatePassiveAggressive(Margin = Margin, Norm = Norm)}: }{
#'    \code{updatePassiveAggressive} method return \code{Weight} based on PA's algorithm.
#'  }
#'  \item{\code{updateAROW(Margin = Margin, Norm = Norm)}: }{
#'    \code{updateAROW} method return \code{Weight} and \code{Param} based on AROW's algorithm.
#'  }
#' }
#' @keywords documentation

library(compiler)
library(BBmisc)
library(iterators)
suppressPackageStartupMessages(library(Matrix))

# avoid cmpfun's message(Note: no visible ...)
suppressAll({
  OnlineLearning <- setRefClass(
    Class = "OnlineLearning", 
    fields = list(
      Weight = "vector",
      Param = "list",
      Response = "numeric",
      Feature = "numeric",
      LossFunc = "function"
    ),
    methods = list(
      initialize = cmpfun(function(FeatureNum, Param){
      
        if(sum(is.element(c("Class", "Regression"), Param$Mode)) != 1){
          warning()
        }
        
        Param <<- Param
        if(Param$Algo == "PA"){
          if(Param$Mode == "Class"){
            Param$E <<- 1
          } else if(Param$Mode == "Regression"){
            Param$E <<- Param$E
          }
          LossFunc <<- .self$updatePassiveAggressive
        } else if(Param$Algo == "AROW"){
          if(Param$Mode == "Class"){
            Param$E <<- 1
          } else if(Param$Mode == "Regression"){
            Param$E <<- Param$E
          }
          
          if(Param$IsFull){
            Param$Alpha <<- diag(FeatureNum)
          } else{
            Param$Alpha <<- sparseMatrix(i = seq(FeatureNum), j = seq(FeatureNum), x = 1)
          }
          LossFunc <<- .self$updateAROW
        }
        
        Weight <<- numeric(FeatureNum)
      }),
      finalize = cmpfun(function(){
        Weight <<- numeric(length(.self$Weight))
      }),
      
      updateWeight = cmpfun(function(Response, FeatureVec){
        if(.self$Param$Mode == "Class"){
          Margin <- Response * as.numeric(.self$Weight %*% FeatureVec)
        } else if(.self$Param$Mode == "Regression"){
  
        }
  
        if(Margin < Param$E){
          Response <<- Response
          Feature <<- FeatureVec
          .self$LossFunc(Margin = Margin, Norm = as.numeric(crossprod(FeatureVec)))
        }
      }),
      updatePassiveAggressive = cmpfun(function(Margin = Margin, Norm = Norm){
        if(.self$Param$Mode == "Class"){
          Res <- .self$Response
          Loss <- 1 - Margin
        } else if(.self$Param$Mode == "Regression"){
  #        Res <- as.numeric(sign(.self$Response - (.self$Weight %*% as.numeric(.self$Feature))))
  
        } else{
          warning()
        }
        
        # Update
        # Only PA-2
        # Param$Alpha <<- Loss / Norm # PA
        # Param$Alpha <<- min(.self$Param$C, (Loss / Norm)) # PA-1      
        Param$Alpha <<- Loss / (Norm + (1 / (2 * .self$Param$C))) # PA-2
        
        Weight <<- Weight + (Param$Alpha * Res * .self$Feature)
      }),
      updateAROW = cmpfun(function(Margin = Margin, Norm = Norm){
        Confidence <- as.numeric(t(Feature) %*% Param$Alpha %*% Feature)
        
        # compute parameter
        Beta <- as.numeric(1 / (Confidence + Param$Gamma))
        Alpha <- max(0, 1 - (Response * as.numeric(t(Feature) %*% Weight))) * Beta
        
        # Update
        Weight <<- Weight + (Alpha * as.numeric(Param$Alpha %*% (Response * Feature)))
        if(Param$IsFull){
          Param$Alpha <<- Param$Alpha - (Beta * (Param$Alpha %*% (Feature %*% t(Feature)) %*% Param$Alpha))
        } else{
          Param$Alpha <<- Param$Alpha - (Beta * (Param$Alpha %*% sparseMatrix(i = seq(length(Feature)), j = seq(length(Feature)), x =(Feature ** 2)) %*% Param$Alpha))
        }
      })
    )
  )
})


#' The \code{updateWeight} method of \code{OnlineLearning} class
#'
#' \code{updateWeight} method
#'
#' \code{updateWeight} method call \code{updatePassiveAggressive} or \code{updateAROW} using \code{Param}.
#'
#' @name updateWeight,OnlineLearning-method
#' @aliases updateWeight
#' @rdname OnlineLearning-updateWeight
#' @docType methods
#'
#' @param Response integer number of a class label and FeatureVec numeric vector consist of feature values. 
#' @usage \S4method{updateWeight}{OnlineLearning}(Response = 1, FeatureVec = c(0, 1, 1))
#' @examples
#' library(Matrix)
#' library(kernlab)
#' library(iterators)
#' # create training data
#' data(spam)
#' InsY <- matrix(ifelse(as.integer(spam[, "type"]) == 1, 1, -1), ncol = 1)
#' InsX <- spam[, -ncol(spam)]
#' InstanceOL <- OnlineLearning$new(
#'                ncol(InsX),
#'                list("Mode" = "Class", "Algo" = "AROW", "Gamma" = 0.5, "IsFull" = FALSE)
#'               )
#' set.seed(10)
#' ShuffleIter <- iter(sample(length(InsY), length(InsY), replace = FALSE))
#' # train
#' while(all(tryCatch(i <- nextElem(ShuffleIter), error = function(Cond){ return (FALSE)}) != FALSE)){
#'   InstanceOL$updateWeight(Response = InsY[i], FeatureVec = as.numeric(InsX[i, ]))
#' }
#' # test
#' Predict <- sign(apply(InsX, 1, FUN = "%*%", InstanceOL$Weight))
#' table(Predict, as.numeric(InsY))
NULL


#' The \code{updatePassiveAggressive} method of \code{OnlineLearning} class
#'
#' \code{updatePassiveAggressive} method
#'
#' \code{updatePassiveAggressive} method updates \code{Weight} based on Online Passive aggressive(this method is implement by PA-2).
#'
#' @name updatePassiveAggressive,OnlineLearning-method
#' @aliases updatePassiveAggressive
#' @rdname OnlineLearning-updatePassiveAggressive
#' @docType methods
#'
#' @param Margin \code{Response} * \code{Weight} * \code{Feature}
#' @param Norm sum(\code{Feature} ** 2)
#' @usage \S4method{updatePassiveAggressive}{OnlineLearning}(Margin = Margin, Norm = Norm)
#' @examples
#' library(Matrix)
#' library(kernlab)
#' # create training data
#' data(spam)
#' InsY <- matrix(ifelse(as.integer(spam[, "type"]) == 1, 1, -1), ncol = 1)  
#' InsX <- spam[, -ncol(spam)]
#' # make instance
#' InstanceOL <- OnlineLearning$new(
#'                ncol(InsX),
#'                list("Mode" = "Class", "Algo" = "PA", "C" = 1)
#'               )
#' print(InstanceOL$Weight) # identical numeric(ncol(InsX))
#' # set values
#' InstanceOL$Response <- InsY[1]  
#' InstanceOL$Feature <- as.numeric(InsX[1, ])  
#' Margin <- InstanceOL$Response * as.numeric(InstanceOL$Weight %*% InstanceOL$Feature)  
#' Norm <- as.numeric(crossprod(InstanceOL$Feature))  
#' # update
#' InstanceOL$updatePassiveAggressive(Margin = Margin, Norm = Norm)  
#' print(InstanceOL$Weight)  
NULL


#' The \code{updateAROW} method of \code{OnlineLearning} class
#'
#' \code{updateAROW} method
#'
#' \code{updateAROW} method updates \code{Weight} based on AROW(Adaptive Regularization Of Weight Vector).
#'
#' @name updateAROW,OnlineLearning-method
#' @aliases updateAROW
#' @rdname OnlineLearning-updateAROW
#' @docType methods
#'
#' @param Margin \code{Response} * \code{Weight} * \code{Feature}
#' @param Norm sum(\code{Feature} ** 2)
#' @usage \S4method{updateAROW}{OnlineLearning}(Margin = Margin, Norm = Norm)
#' @examples
#' library(Matrix)
#' library(kernlab)
#' data(spam)
#' InsY <- matrix(ifelse(as.integer(spam[, "type"]) == 1, 1, -1), ncol = 1)
#' InsX <- spam[, -ncol(spam)]
#' InstanceOL <- OnlineLearning$new(
#'                ncol(InsX), 
#'                list("Mode" = "Class", "Algo" = "AROW", "Gamma" = 0.5, "IsFull" = FALSE)
#'               )
#' print(InstanceOL$Weight) # identical numeric(ncol(InsX))
#' # set values
#' InstanceOL$Response <- InsY[1]
#' InstanceOL$Feature <- as.numeric(InsX[1, ])
#' Margin <- InstanceOL$Response * as.numeric(InstanceOL$Weight %*% InstanceOL$Feature)
#' Norm <- as.numeric(crossprod(InstanceOL$Feature))
#' # update
#' InstanceOL$updateAROW(Margin = Margin, Norm = Norm)
#' print(InstanceOL$Weight)
NULL
