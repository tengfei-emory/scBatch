qnorm <-
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
  
  qt.x <- ccc[batches == big, batches == big]
  qt.xrow <- qt.x
  qt.xcol <- qt.x
  
  for(i in 1:nrow(qt.x)){
    qt.xrow[i,] <- zp.quantile(qt.x[1,],qt.x[i,])
  }
  for(i in 1:ncol(qt.x)){
    qt.xcol[,i] <- zp.quantile(qt.x[,1],qt.x[,i])
  }
  
  for(ba1 in unique(batches))
  {
    for(ba2 in unique(batches))
    {
      qt.y<-ccc[batches == ba1, batches == ba2]
      dim1 <- max(dim(ccc[batches == ba1, batches == ba2]))
      dim2 <- min(dim(ccc[batches == ba1, batches == ba2]))
      dy <- dim(qt.y)
      if (ba1 != big && ba2 != big){
        dim0 = dim(ccc[batches == ba1, batches == ba2])
        #qt.y <- rbind(qt.y,matrix(0,abs(max_dim-dim0[1]),min(max_dim,dim0[2])))
        #qt.y <- cbind(qt.y,matrix(0,max_dim,abs(max_dim-dim0[2])))
        qt.y1 <- qt.y
        qt.y2 <- qt.y
        for(i in 1:nrow(qt.y)){
          qt.y1[i,] <- zp.quantile(qt.xrow[i,],qt.y1[i,])
        }
        for(i in 1:ncol(qt.y)){
          qt.y2[,i] <- zp.quantile(qt.xcol[,i],qt.y2[,i])
        }
        qt.y <- (qt.y1 + qt.y2)/2
      }else if (ba1 != big){
        #qt.y <- rbind(qt.y,matrix(0,abs(max_dim-dim2),max(max_dim,dim2)))
        qt.y1 <- qt.y
        qt.y2 <- qt.y
        for(i in 1:nrow(qt.y)){
          qt.y1[i,] <- zp.quantile(qt.xrow[i,],qt.y1[i,])
        }
        for(i in 1:ncol(qt.y)){
          qt.y2[,i] <- zp.quantile(qt.xcol[,i],qt.y2[,i])
        }
        qt.y <- (qt.y1 + qt.y2)/2
      }else if(ba2 != big){
        #qt.y <- cbind(qt.y,matrix(0,max(max_dim,dim2),abs(max_dim-dim2)))
        qt.y1 <- qt.y
        qt.y2 <- qt.y
        for(i in 1:nrow(qt.y)){
          qt.y1[i,] <- zp.quantile(qt.xrow[i,],qt.y1[i,])
        }
        for(i in 1:ncol(qt.y)){
          qt.y2[,i] <- zp.quantile(qt.xcol[,i],qt.y2[,i])
        }
        qt.y <- (qt.y1 + qt.y2)/2
      }
      ccc[batches == ba1, batches == ba2] <- qt.y[1:dy[1],1:dy[2]]
    }
  }
  
  return(ccc)
}
