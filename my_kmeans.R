my_kmeans<- function(x, clusters = 3, maxiter=20) {
  set.seed(1)
  # start with random cluster centers
  clus <- data.frame(name=1:clusters,x[floor(runif(n=clusters, min = 1,max = nrow(x))),])
  
  # data points and cluster assignment in "data"
  # cluster coordinates in "clus"
  data <- data.frame(x, clus = NA)
  
  finish <- FALSE
  clus_old <- data$clus
  max_iter <- 1
  
  while(finish == FALSE){
    # assign cluster with minimum distance to each data point
    d<-data.frame()
    for(i in 1:nrow(x)){
      y<-x[i,]
      dis_from_centers<-rowSums((y[rep(1,clusters),]-clus[,2:ncol(clus)])^2)
      data$clus[i]<-which.min(dis_from_centers)
      d<-rbind(d,sqrt(dis_from_centers))
    }
    colnames(d)<-sprintf("dist%d",seq(1:clusters))
    # calculate new cluster centers
    for(i in 1:clusters) {
      tmp <- colMeans(subset(data[,1:(ncol(data)-1)], data$clus == i))
      clus[i,2:ncol(clus)] <- (tmp-mean(tmp))/sd(tmp)
    }
    if(max_iter>maxiter) finish <- TRUE
    if(identical(data$clus,clus_old)) finish <- TRUE
    clus_old <- data$clus
    max_iter <- max_iter+1
  }
  wss <- 0
  for(i in 1:nrow(x)){
    wss <- wss + sum((clus[data$clus[i],2:ncol(clus)]-x[i,])^2)
  }
  centers <- t(clus)
  centers <- centers[2:nrow(centers),]
  colnames(centers) <-sprintf("center %d",seq(1:clusters))
  result <- list(labels=data$clus,centers=centers, wss=wss,d=d)
  return(result)
}