detect <- function(x, threshold, n){
  a <- x*(x<threshold)
  b <- numeric(length(x))
  for (i in n:length(x)) {
    if(a[i]!=0){
      b[i:(i-n+1)] <- b[i:(i-n+1)]+1
    }
  }
  b <- c(1:length(x))*(b>0)
  d <- data.frame()
  
  status <-1 # in anomaly
  if(b[1]==0){
    status <-0 # no anomaly
  }
  
  for (i in 1:length(b)) {
    if(b[i]!=0){
      if(status ==0){
        start <-b[i]
        status <-1
      }
    }
      
    if(b[i]==0){
      if (status ==1){
        if(i>1){
          end <-b[i-1]
          d <- rbind(d, data.frame(start,end))
          status <-0
        }
      }
    }
  }
  d[,'start'] <- d[,'start']-5
  d[,'end'] <- d[,'end']+n-1
  return(d)
}
