qnorm1 <-  function(ccc,batches){
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
  qt.x <- ccc[batches == big, batches == big]
  for(ba1 in unique(batches)){
    for(ba2 in unique(batches)){
      qt.y<-ccc[batches == ba1, batches == ba2]
      dy <- dim(qt.y)
      if (ba1 != big || ba2 != big){
        qt.y1 <- qt.y
        qt.y2 <- qt.y
        for(i in 1:nrow(qt.y)){
          qt.y1[i,] <- zp.quantile(as.vector(qt.x),qt.y1[i,])
        }
        for(i in 1:ncol(qt.y)){
          qt.y2[,i] <- zp.quantile(as.vector(qt.x),qt.y2[,i])
        }
        qt.y <- (qt.y1 + qt.y2)/2
      }
      ccc[batches == ba1, batches == ba2] <- qt.y[1:dy[1],1:dy[2]]
    }
  }
  return(ccc)
}
