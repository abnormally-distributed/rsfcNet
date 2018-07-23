#' Variable selection using LASSO (L1 penalty) regression.
#'
#' This function is a convenience wrapper for glmnet's lasso capability. It will automatically select the best lambda value (strength of the L1 penalty) via cross-validation and return a table of variables that are good predictors of the response variable(s).
#' @param x A matrix of predictor (independent) variables.
#' @param y A vector of response (dependent) variables. For gaussian models this will be a numeric vector. For poisson it will be a numeric vector of 0+ integers (count data). For binomial or multinomial regression this will be a vector containing factor labels. For mgaussian (multiple response linear regression; see ) y is a matrix of numeric response variables.
#' @param var.names The names of each node.
#' @param family Options include "gaussian", "poisson", "binomial" (logistic), "multinomial" (multiple outcome logistic), and "mgaussian"
#' @param output Options include "html" (the default), "markdown", "latex", or return to console.
#' @param return.glmnet If set to TRUE returns a copy of the cross-validated glmnet model to the global environment. Defaults to false. Set to TRUE if you want to test it on a holdout sample (only recommended if n is large).
#' @return A table of variables that survived the variable selection.
#' @export
#' @author Brandon Vaughan
#' @examples
#' lasso.results = LASSO(predictor_matrix, smoker_status, family="binomial")
#'
#' @references
#'
#' Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22. URL http://www.jstatsoft.org/v33/i01/.
#'
#' Tibshirani, Robert (1996). Regression Shrinkage and Selection via the lasso. Journal of the Royal Statistical Society. Series B (methodological). Wiley. 58 (1): 267–88. https://www.jstor.org/stable/2346178
#'
LASSO = function(x, y, var.names = NULL, family="gaussian", output="html", return.glmnet=FALSE) {

  if (is.null(var.names)==TRUE) {
    colnames(x) = paste(rep("Var", ncol(x)), as.character(seq(1:ncol(x))))
  } else if (is.null(var.names)==FALSE) {
    colnames(x) = var.names
  }
  x.train = as.matrix(x)
  y.train = as.matrix(y)
  set.seed(1)
  cv.glmnet = cv.glmnet(x.train, y.train, alpha=1, nlambda=200, family=family, lambda.min.ratio=0.0001)
  if (return.glmnet==TRUE) {
    glmnet.model <<- cv.glmnet
  }
  coefs.glmnet = as.matrix(coef(cv.glmnet, s="lambda.min"))
  rownames(coefs.glmnet) = c("(Intercept)",var.names)
  coefs.glmnet = coefs.glmnet[ rowSums(coefs.glmnet)!=0, ] # Remove the estimates equal to zero
  var.names.res = rownames(coefs.glmnet)
  coefs.glmnet = as.data.frame(coefs.glmnet)
  terms = rownames(coefs.glmnet)
  if (output == "console") {
    return(coefs.glmnet)
  } else if (output == "html"){
    coefs.glmnet  = cbind(terms, data.frame(coefs.glmnet, row.names=var.names.res))
    colnames(coefs.glmnet) =c("term", "estimate")
    coefs.glmnet[,1:2]  %>%
      kable("html",escape = FALSE, digits=3, align='r') %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      add_header_above(c("Glmnet LASSO Regression Coefficients" = 2))
  } else if (output == "markdown") {
    coefs.glmnet  = cbind(terms, data.frame(coefs.glmnet, row.names=var.names.res))
    colnames(coefs.glmnet) =c("term", "estimate")
    coefs.glmnet[,1:2]  %>%
      kable("markdown",escape = FALSE, digits=3, align='r')
  } else if (output=="latex") {
    coefs.glmnet  = cbind(terms, data.frame(coefs.glmnet, row.names=var.names.res))
    colnames(coefs.glmnet) =c("term", "estimate")
    coefs.glmnet[,1:2]  %>%
      kable("latex",escape = FALSE, digits=3, align='r') %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      add_header_above(c("Glmnet LASSO Regression Coefficients" = 2))
  }
}
