#' Fit a vector of variables to the best distribution
#'
#' Automatically detect which distributions are appropriate to the data using Maximum Likelihood Estimation, then find the log likelihood and return
#' the parameters of the best fitting distribution. For example, fit the degree distribution to see if a power law distribution best describes the data.
#'
#' @param variable A vector of values describing a graph theory metric of interest.
#' @param only.best Return only the parameters for the best fitting distribution (default is TRUE) or all distributions (FALSE).
#'
#' @details First the function detects whether or not any negative values are present. If so, distributions that require all
#' positive values are excluded. Next, the data are fit to the appropriate distributions and the log likelihood is used
#' to select the best fitting model. This function was primarily designed for assessing whether or not the degree
#' distribution or strength distribution of a graph follows the power law, but can also be used to fit distributions
#' to other metrics (ie, eigenvector centrality, clustering coefficient, etc).
#'
#' @export
#'
#' @examples
#' **## Not run:**
#' k = degree(graph)
#' distr = dist_fit(k)
#'
#' s = strength(graph)
#' distr = dist(s)
#' ## End(**Not run**)
#'
#'
#' @seealso
#' \code{\link[rsfcNet]{strength_distribution}}
#'
#'
dist_fit = function(variable, only.best=TRUE) {

  if (any(variable < 0)=="TRUE") {

    fits <- list(
      normal = MASS::fitdistr(variable,"normal"),
      t = MASS::fitdistr(variable, "t", lower=c(mean(variable), sd(variable), 1)),
      power.law = igraph::fit_power_law(variable),
      gamma = MASS::fitdistr(variable, "gamma", lower=c(.0001,.0001))

    )
    fits[[1]]$sd = NULL
    fits[[2]]$sd = NULL
    fits[[4]]$sd = NULL
    fits[[3]]$loglik = fits[[3]]$logLik
    fits[[3]]$continuous = NULL
    fits[[3]]$logLik  = NULL
    fits[[3]]$KS.stat = NULL
    fits[[3]]$KS.p = NULL
  } else {

    fits <-  list(
      normal = MASS::fitdistr(variable,"normal"),
      lognormal = MASS::fitdistr(variable,"lognormal"),
      t = MASS::fitdistr(variable, "t", lower=c(mean(variable), sd(variable), 1)),
      gamma = MASS::fitdistr(variable, "gamma", lower=c(.0001,.0001)),
      exponential = MASS::fitdistr(variable, "exponential"),
      power.law = igraph::fit_power_law(variable)
    )

    fits[[1]]$sd = NULL
    fits[[2]]$sd = NULL
    fits[[3]]$sd = NULL
    fits[[4]]$sd = NULL
    fits[[5]]$sd = NULL
    fits[[6]]$loglik = fits[[6]]$logLik
    fits[[6]]$continuous = NULL
    fits[[6]]$logLik  = NULL
    fits[[6]]$KS.stat = NULL
    fits[[6]]$KS.p = NULL
  }

  logLik = sapply(fits, function(i) i$loglik)
  logLik2 = logLik[-which(logLik==max(logLik))]

  print("Log Likelihoods for Fitted Models")
  print(sapply(fits, function(i) i$loglik))

  print("Log-Likelihood Ratio for Best Model against Second Best")
  max(logLik)/max(logLik2)

  if (only.best==TRUE) {
  fits[[6]]$loglik  = NULL
  fits[[6]] = unlist(fits[[6]])
  best_fit = fits[which(logLik==max(logLik))]
  return(best_fit)
  } else {
  fits[[6]]$loglik  = NULL
  fits[[6]] = unlist(fits[[6]])
  return(fits)
  }
}

