#' Adjust the distance matrix by quantile normalization for data with batch effect
#'
#'
#' @description This function applies quantile normalization on the distance matrix (dissimilarity matrix) and return the corrected distance matrix.
#' @param dat The original p*n batch effect data with n subjects and p RNA-seq measurements.
#' @param batch The vector of length n indicating which batch the subjects belong to.
#' @param method Method for the quantile normalization. There are two options: "row/column" and "vectorize".
#' @param cor_method Method to calculate the correlation matrix, can be 'spearman'(default), 'pearson' or 'kendall'.
#' @param tol The tolerance for the iterative method "row/column", which is the Euclidean distance of the vectorized two dissimilarity matrices before and after each iteration.
#' @param max Maximum number of the iteration if the tolerance is not reached.
#' @param logdat Whether conducting log transformation to data or not.
#' @param standardize Whether conducting standardization [(dat - mean)/sqrt(var)] to data or not.
#' @return Returns the corrected 1-correlation matrix between subjects.
#' @author Teng Fei. Email: tfei@emory.edu
#' @references Fei et al (2018), Mitigating the adverse impact of batch effects in sample pattern detection, Bioinformatics, <https://doi.org/10.1093/bioinformatics/bty117>.
#' @export
#' @examples
#'
#' library(pheatmap) #drawing heatmap
#'
#' data("ENCODE") #load the ENCODE data
#'
#' #Before correction, the subjects are clustered by species
#' pheatmap(cor(ENCODE))
#'
#' #Assigning the batches based on species
#' batches <- c(rep(1,13),rep(2,13))
#'
#' #QuantNorm correction
#' corrected.distance.matrix <- QuantNorm(ENCODE,batches,method='row/column', cor_method='pearson',
#'                                        logdat=FALSE, standardize = TRUE, tol=1e-4)
#' pheatmap(1-corrected.distance.matrix)

QuantNorm <- function (dat, batch, method = "row/column", cor_method = 'spearman', tol = 1e-2, max = 50, logdat = TRUE,
                       standardize = FALSE)
{
  dist = 10
  iter = 0

  if (standardize == TRUE) {
    ccc <- standardization(dat, batch, method=cor_method)
  }
  else if (logdat == FALSE) {
    ccc <- 1 - stats::cor(dat, method = cor_method)
  }
  else {
    ccc <- 1 - stats::cor(log(dat + 1), method = cor_method)
  }

  if (method == "vectorize") {
    ccc <- qnorm2(ccc, batch)
  }

  else if (method == "row/column") {
    while (dist > tol && iter < max){
      ccc.0 <- ccc
      ccc <- qnorm1plus(ccc, batch)
      dist = sqrt(sum((as.vector(ccc)-as.vector(ccc.0))^2))
      iter = iter+1
    }

    if (dist <= tol){
      cat(paste("Algorithm converged after", iter, "iterations. \n"))
    }else{
      cat(paste("The algorithm did not converge after", max,"iteration. \n The difference between the last two iteration is", dist,".\n"))
    }

  }else{
    cat("method must be 'row/column' or 'vectorize'. \n")
  }
  return(ccc)
}
