#' Correct scRNA-seq count matrix subject to batch effects by sample distance matrix correction
#'
#'
#' @description As an extension of package 'QuantNorm', which corrects the distance matrix to obtain reasonable sample pattern, this package use gradient descent algorithm to correct count matrix by linear transformation.
#' @param count.mat The original p*n batch effect data with n subjects and p RNA-seq measurements.
#' @param dist.mat The n*n distance matrix obtained by QuantNorm.
#' @param weight.mat An initial n*n weight matrix to conduct linear transformation. Default to be identity matrix if not specified.
#' @param m Number of groups to be divided for coordinate gradient descent. 1 < m <= n. Default to be 0.1n if not specified.
#' @param max.iter Maximum number of the iteration if the tolerance is not reached.
#' @param step.size Step size of the gradient descent algorithm.
#' @param tol Stopping criteria of the algorithm. The algorithm stops if the step size is smaller than tol.
#' @return Returns the corrected count matrix.
#' @author Teng Fei. Email: tfei@emory.edu
#' @references Fei et al (2018), Mitigating the adverse impact of batch effects in sample pattern detection, Bioinformatics, <https://doi.org/10.1093/bioinformatics/bty117>.
#' @export


scBatch <- function(count.mat,dist.mat,weight.mat, m, max.iter=30,step.size=0.0001, tol=1e-10){
  pn = dim(count.mat)
  if (!exists("weight.mat")){
    weight.mat = diag(pn[2])
  }
  if (!exists("m")){
    m = floor(0.1*pn[2])
  }
  p = pn[1]
  core = t(count.mat)%*%t((diag(p) - matrix(1,p,p)/p))%*%(diag(p) - matrix(1,p,p)/p)%*%count.mat
  core = (core + t(core))/2

  start_time <- proc.time()
  for (i in 1:max.iter){
    group = sample(1:m,size=pn[2],replace=TRUE,prob=rep(1/m,m))

    for (k in 1:m){
      idx = which(group==k)
      fdf = derifcn(count.mat,weight.mat,dist.mat,core,idx)
      f = fdf[[1]]
      df = fdf[[2]]
      for (j in 1:5){
        update.mat = weight.mat - step.size*df
        update.mat = update.mat/(rep(1,pn[2]) %o% apply(abs(update.mat),2,max))
        new.count.mat = count.mat%*%update.mat;
        A = 1-cor(new.count.mat);
        fnew = sum((as.vector(A)-as.vector(dist.mat))^2)

        if(fnew >= f){
          step.size = 0.5*step.size
        }else{
          step.size = 1.5*step.size
          weight.mat = update.mat
          end_time <- proc.time()
          time <- (end_time-start_time)[3]
          cat(i,',',f,',',step.size,',',time,'\n')
          break
        }
      }
    }
    if (step.size < tol){
      break
    }
  }

  new.count = count.mat%*%weight.mat
  new.count = new.count - min(new.count)
  return(new.count)
}
