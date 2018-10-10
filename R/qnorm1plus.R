qnorm1plus <-  function(ccc,batches){
  batch_num <- max(batches)
  max_dim <- 0
  combnbatch <- combn(unique(batches),2)
  for (i in 1:batch_num){
    dim <- dim(ccc[batches == i, batches == i])[1]
    if (dim > max_dim){
      max_dim = dim
      big = i
    }
    combnbatch = cbind(combnbatch,c(i,i))
  }
  qt.x <- ccc[batches == big, batches == big]


  for (j in 1:dim(combnbatch)[2]){
    ba1 = combnbatch[1,j]
    ba2 = combnbatch[2,j]
    qt.y<-ccc[batches == ba1, batches == ba2]
    dy <- dim(qt.y)
    maxdy <- max(dy)
    mindy <- min(dy)
    whichmax <- which.max(dy)

    if (ba1 != big || ba2 != big){
      qt.y1 <- qt.y
      qt.y2 <- qt.y

      for (i in 1:maxdy){
        if(whichmax == 1){
          if(i <= mindy){
            qt.y2[,i] <- zp.quantile(as.vector(qt.x),qt.y2[,i])
          }
          qt.y1[i,] <- zp.quantile(as.vector(qt.x),qt.y1[i,])
        }else{
          if(i <= mindy){
            qt.y1[i,] <- zp.quantile(as.vector(qt.x),qt.y1[i,])
          }
          qt.y2[,i] <- zp.quantile(as.vector(qt.x),qt.y2[,i])
        }
      }
      qt.y <- (qt.y1 + qt.y2)/2
    }
    ccc[batches == ba1, batches == ba2] <- qt.y[1:dy[1],1:dy[2]]
    if (ba1 != ba2){
      ccc[batches == ba2, batches == ba1] <- t(qt.y[1:dy[1],1:dy[2]])
    }
  }
  return(ccc)
}
