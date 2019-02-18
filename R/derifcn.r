derifcn <- function(count.mat,weight.mat,dist.mat,core,idx){
  pn = dim(count.mat)
  new.count.mat = count.mat%*%weight.mat
  A = 1-stats::cor(new.count.mat)
  D = dist.mat
  AmD = A-D
  B = weight.mat
  f = sum((as.vector(A)-as.vector(D))^2)

  pn = dim(count.mat)
  W = (t(B)%*%core%*%B);
  W[W<0]=0
  W = W^0.5

  M = core%*%B;
  A1 = diag(1/W) %x% (M/(rep(1,length(diag(W))) %o% diag(W)))
  B1 = as.vector(M/(rep(1,length(diag(W))) %o% diag(W^3)))
  B2 = as.vector(t(W^2/(rep(1,length(diag(W))) %o% diag(W^2))))

  df = matrix(0,pn[2],pn[2])

  for (k in 1:length(idx)){
    i = idx[k]
    D = A1[(1+pn[2]*(i-1)):(pn[2]*i),] - as.matrix(B1[(1+pn[2]*(i-1)):(pn[2]*i)]) %x% t(as.matrix(B2[(1+pn[2]*(i-1)):(pn[2]*i)]))
    df[,i] = t(D) %*% AmD[,i]
    df[i,i] = sum(diag(t(D) %*% AmD))
  }

  df = -df

  result = list(f,df)
  return(result)
}
