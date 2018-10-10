zp.quantile <-
function(x,y)
  # x and y are two vectors of equal length
  # x serves as the quantile normalization template. Only updated y is returned
{
  o.x<-order(x)
  #r.x<-rank(x, ties.method = "random")

  o.y<-order(y)
  r.y<-rank(y, ties.method = "average")

  x<-x[o.x]
  y<-y[o.y]

  x2<-x[x>0]
  y2<-y[y>0]

  z.x<-seq(0,1,length.out=length(x2))
  z.y<-seq(0,1,length.out=length(y2))

  new.y2<-stats::approx(x=z.x, y=x2, xout=z.y)$y
  y[y>0]<-new.y2
  y<-y[r.y]

  y
}
