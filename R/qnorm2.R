qnorm2 <-
function(ccc,batches){

  #ccc.0<-ccc<-1-cor(newbatch,method="spearman")

  batch_num <- max(batches)
  max_dim <- 0

  for (i in 1:batch_num){
    dim <- dim(ccc[batches == i, batches == i])[1]
    if (dim > max_dim){
      max_dim = dim
      big = i
    }
  }

  qt.x <- as.vector(ccc[batches == big, batches == big])

  for(ba1 in unique(batches))
  {
    for(ba2 in unique(batches))
    {
      qt.y<-ccc[batches == ba1, batches == ba2]
      dy <- dim(qt.y)
      qt.y <- as.vector(qt.y)
      qt.y <- zp.quantile(qt.x,qt.y)
      qt.y <- as.matrix(qt.y,nrow = dy[1],ncol=dy[2])
      ccc[batches == ba1, batches == ba2] <- qt.y
    }
  }

  return(ccc)
}
