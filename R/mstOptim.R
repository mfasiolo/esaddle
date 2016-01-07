#####
#' Solving multiple related problems using smart initialization
#' @description Solving multiple related problems using smart initialization
#'
#' @param y matrix where each row represents the parameters of a different problem
#' @param lamHat matrix where each row represents the approximate solution to each problem.
#'               This is used to calculate the distance between problems. For instance the distance
#'               between the i-th and the j-th problem is the Euclidean distance between lamHat[i, ] and lamHat[j, ].
#'               In some cases it is appropriate to choose lamHat = y.
#' @param objFun function that solves the problem. I has signature function(.y, .lambda, ...) and it returns a list
#'               that contains the vector "lambda", which is the solution to the problem. 
#' @param gradFun function that returns a matrix containing the derivatives of lambda wrt y. It has 
#'                signature function(.y, .lambda, ...).       
#' @param mst (optional) minimum spanning tree which can be provided by the user.
#' @param ... extra arguments to be passed to objFun() and gradFun()
#' @author Matteo Fasiolo <matteo.fasiolo@@gmail.com>.
#' @export mstOptim
#'

mstOptim  <- function(y, lamHat, objFun, gradFun, mst = NULL, ...)
{
  
  ny <- nrow(y)
  
  if( is.null(mst)  )
  {
    # Calculate Minimum Spanning Tree across the lamHat points
    mat <- .Call("mst", X_ = t(lamHat), PACKAGE = "esaddle");
    mat[1:2, ] <- mat[1:2, ] + 1
        
    # Convert mst to a list, which is easier to traverse
    mst <- rep(list(NULL), ny)
    names(mst) <- 1:ny
    
    for(ii in 1:ncol(mat))
    {
      mst[[ mat[1, ii] ]] <- c(mst[[ mat[1, ii] ]], mat[2, ii])
      mst[[ mat[2, ii] ]] <- c(mst[[ mat[2, ii] ]], mat[1, ii])
    }
  }
  
  # Storage for lambda, dLambda/dy and for gradient of saddlepoint
  DLamDyStore <- list()
  
  todoList <- list( c("child" = 1, "dad" = 1) )
  lambdaE <- NULL
  ii <- todoList[[1]]["child"]
  
  out <- list()
  
  # Traversing ny nodes of the mst and evaluating saddlepoint at each (NOTE kk != ii)
  for(kk in 1:ny){
    
    out[[ii]] <-  objFun(.y = y[ii, ], .lambda = lambdaE, ...)
        
    toAdd <- mst[[ii]][mst[[ii]] != todoList[[1]]["dad"]]
    todoList <- todoList[-1]
    todoList <- c(todoList, lapply(toAdd, function(x_) c("child" = x_, "dad" = unname(ii))))
    
    if(kk < ny)
    {
      ii <- todoList[[1]]["child"] 
      dad <- todoList[[1]]["dad"]
      
      # Compute gradient of dad, if not computed before
      if( is.null(DLamDyStore[[as.character(dad)]]) ){ 
        
        DLamDyStore[[as.character(dad)]] <- gradFun(.y = y[dad, ], .lambda = out[[dad]][["lambda"]], .extra = out[[dad]][["extra"]]) 
        
      }
      
      lambdaE <- out[[dad]][["lambda"]] + DLamDyStore[[as.character(dad)]] %*% drop( y[ii, ] - y[dad, ] )
    }
    
  }
  
  return( out )
  
}