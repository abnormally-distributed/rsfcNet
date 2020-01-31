#' Soft Thresholding Function for Signed Correlation Networks
#'
#' @description This function implements the soft thresholding function from WGCNA. This is used to downweight negative
#' edges while at the same time shifting the matrix entries to be >= 0. This is accomplished by adding 1 to the matrix
#' entries, then dividing by two. Then, the resulting values are raised to a power, denoted as lambda.
#' Recommended values for lambda are 4 to 8 for partial correlation networks, and 8 to 12 for correlation networks. Larger
#' values result in stronger thresholding. The formula used is given below. \cr
#'
#' \eqn{w_{ij}^{\text {signed }}=\left[\frac{\rho\left(i,j\right)+1}{2}\right]^{\lambda}}
#'
#' @param x a correlation or partial correlation matrix
#' @param lambda the power to which the matrix will be raised. defaults to 8.
#' setting lambda to "auto" will automatically choose a value based on the operator norm of the transformed matrix.
#' A small value is added to lambda at each iteration, and the algorithm stops when the change in operator norm is <= than 0.01.
#' @param zero.diag if TRUE (the default) the diagonal of the returned matrix will be set to zero. otherwise, it will be set to 1.
#'
#' @return a matrix
#' @export
#'
soft_thresh <- function(x, lambda = 8, zero.diag = T){

  stfun <- function(x, lambda = "auto"){
    if (lambda < 1) return(cat("lambda must be 1 or greater."))
    if (sum(diag(x))==0)  diag(x) <- 1
    x = as.matrix(Matrix::nearPD((zapsmall(((x+1)/2)^lambda, digits = 10)))$mat)
    return(x)
  }

  iter <- 0; normval <- norm(stfun(x,lambda=1),"O"); i = 2; change <- 10000000;
  limit.reached <- F
  if (lambda == "auto"){
    while(change > 0.01 && !limit.reached){
      iter <- iter+1
      new.normval <- norm(stfun(x, lambda = i), "O")
      limit.reached <- new.normval <= 1
      i <- i + 0.175
      change <- normval - new.normval
      normval <- new.normval
    }
    lambda <- i
    out <- stfun(x, lambda=lambda)
    if (zero.diag) diag(out) <- 0
    attr(out, "lambda") <- lambda; attr(out, "o.norm") <- norm(out, "O")
    return(out)
  }
  else{
    out <- stfun(x,lambda=lambda)
    if (zero.diag) diag(out) <- 0
    return(out)
  }
}
