#' Use a shrinkage estimated t-test over all pairwise comparisons.
#'
#' This function returns the pairwise comparisons for all nodes using the shrinkage t-test from the st package. The shrinkage t-statistic is a regularized estimate of the t-statistic for large P situations. It accounts for correlations between the independent variables such that an over-abundance of significant results are not driven by the correlation of variables. This approach was originally applied to microarray data (see Opgen-Rhein and Strimmer, 2007). This function extends the original by providing Bayes Factor tests, p-values, and q-values (FDR adjusted p-values). Bayes Factors calculated with BayesFactor package and q-values with fdrtools package. If n is less than half of P the Bayes Factors may be better to use, as q-values will be too conservative. Bayes Factors quantify evidence for or against a hypothesis and even a relatively low Bayes Factor threshold (such as BF > 3) is more conservative than a p<=.05 criterion (see Wetzels et al, 2011; Jeong and Boeck, 2017)
#' @param variables a data frame or matrix containing all dependent variables.
#' @param demograhics a vector containing factor labels for each subject (the independent variable)
#' @param ngrp1 Sample size of group 1.
#' @param ngrp2 Sample size of group 2.
#' @param var.equal Should variances be assumed to be equal? Default is false, which maintains power when variances are equal anyway.
#' @param output Options equal print all evidence measures to the console (select this if you want to assign the results to a data frame), or "html.qvalue", "html.pvalue", and "html.bf" to print specific oututs to an html table. Replace html with latex for a latex table, with markdown for a markdown format table, or with console to print to the console (ie, "latex.bf", "markdown.pvalue", "console.qvalue")
#' @param alpha Significance level for p-values and q-values. Defaults to .05
#' @param bf.threshold The strength of evidence needed to threshold interesting results from uninteresting. Defaults to 5. 3 is minimum to be considered any evidence, 5 the minimum for substantial evidence, and 10 for very strong evidence.
#' @param rscale The prior distribution scale for expected effect size. Defaults to "wide". Set to "ultrawide" to set more stringent expectations of effect size, "medium" for a liberal threshold (perhaps not wise for big P scenarios). See documentation in BayesFactor package for more information.
#' @return A table of the significantly different nodes in a format of your choice with a definition of significance of your choice.
#' @export
#' @author Brandon Vaughan
#' @examples
#' shrinkage_ttest(eigen.centrality, demographics$gender, var.equal=FALSE, ngrp1=63, ngrp2=58, output="html.bf")
#' @references
#' Opgen-Rhein, R., and K. Strimmer. 2007. Accurate ranking of differentially expressed genes by a distribution-free shrinkage approach. Statist. Appl. Genet. Mol. Biol. 6:9. http://dx.doi.org/10.2202/1544-6115.1252
#' Wetzels, R., Matzke, D., Lee, M. D., Rouder, J. N., Iverson, G. J., & Wagenmakers, E. (2011). Statistical Evidence in Experimental Psychology. Perspectives on Psychological Science, 6(3), 291-298. doi:10.1177/1745691611406923
#' Jeon, M., & Boeck, P. D. (2017). Decision qualities of Bayes factor and p value-based hypothesis testing. Psychological Methods, 22(2), 340-360. doi:10.1037/met0000140

shrinkage_ttest = function(variables, demographics, ngrp1=NULL, ngrp2=NULL,var.equal = FALSE, output="full.console", bf.threshold=5, alpha=.05,rscale = "medium"){

  diff = colMeans(variables[1:ngrp1,])-colMeans(variables[(ngrp1+1):(ngrp1+ngrp2),])

  tstats = st::shrinkt.stat(as.matrix(variables), demographics, var.equal = var.equal)

  Bayes.Factor = sapply(tstats, function(t) BayesFactor::ttest.tstat(t, n1=ngrp1, n2 = ngrp2, nullInterval = NULL, rscale = rscale, complement = FALSE, simple = TRUE))

  p.value = 2*sapply(abs(tstats), function(q) pt(q, df=pmin(ngrp1,ngrp2)-1, lower.tail=FALSE))

  q.value = as.matrix(fdrtool::fdrtool(p.value, statistic="pvalue", plot=FALSE)$qval)
  results = data.frame(mean.difference = diff, shrinkage.statistic=tstats, Bayes.Factor=Bayes.Factor, p.value=p.value, q.value=q.value)

  if (output=="full.console") {
    return(results)
  } else if (output=="html.bf") {
    results = results[which(results$Bayes.Factor >= bf.threshold),1:3]
    results[,1:3]  %>%
      knitr::kable("html",escape = FALSE, digits=3, align='r') %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      add_header_above(c("Shrinkage-T Test Results" = 4))
  } else if (output=="html.pvalue") {
    results = results[which(results$p.value <= alpha),c(1,2,4)]
    results[,1:3]  %>%
      knitr::kable("html",escape = FALSE, digits=3, align='r') %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      add_header_above(c("Shrinkage-T Test Results" = 4))
  } else if (output=="html.qvalue") {
    results = results[which(results$q.value <= alpha),c(1,2,5)]
    results[,1:3]  %>%
      knitr::kable("html",escape = FALSE, digits=3, align='r') %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      add_header_above(c("Shrinkage-T Test Results" = 4))
  }  else if (output=="latex.bf") {
    results = results[which(results$Bayes.Factor >= bf.threshold),1:3]
    results[,1:3]  %>%
      knitr::kable("latex",escape = FALSE, digits=3, align='r') %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      add_header_above(c("Shrinkage-T Test Results" = 4))
  } else if (output=="latex.pvalue") {
    results = results[which(results$p.value <= alpha),c(1,2,4)]
    results[,1:3]  %>%
      knitr::kable("latex",escape = FALSE, digits=3, align='r') %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      add_header_above(c("Shrinkage-T Test Results" = 4))
  } else if (output=="latex.qvalue") {
    results = results[which(results$q.value <= alpha),c(1,2,5)]
    results[,1:3]  %>%
      knitr::kable("latex",escape = FALSE, digits=3, align='r') %>%
      kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      add_header_above(c("Shrinkage-T Test Results" = 4))
  } else if (output=="markdown.bf") {
    results = results[which(results$Bayes.Factor >= bf.threshold),1:3]
    results[,1:3]  %>%
      knitr::kable("markdown",escape = FALSE, digits=3, align='r')
  } else if (output=="markdown.pvalue") {
    results = results[which(results$p.value <= alpha),c(1,2,4)]
    results[,1:3]  %>%
      knitr::kable("markdown",escape = FALSE, digits=3, align='r')
  } else if (output=="markdown.qvalue") {
    results = results[which(results$q.value <= alpha),c(1,2,5)]
    results[,1:3]  %>%
      knitr::kable("markdown",escape = FALSE, digits=3, align='r')
  } else if (output=="console.bf") {
    results = results[which(results$Bayes.Factor >= bf.threshold),1:3]
    return(results)
  } else if (output=="console.pvalue") {
    results = results[which(results$p.value <= alpha),c(1,2,4)]
    return(results)
  } else if (output=="console.qvalue") {
    results = results[which(results$q.value <= alpha),c(1,2,5)]
    return(results)
  }

}

