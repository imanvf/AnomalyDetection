# data pre-processing
setwd("YOUR PATH")
source('my_kmeans.r')
source('multiplot.r')

# For ECG dataset use this:
training <- c(500:1500) # clean data == No Anomaly
data <-read.csv("ECG.csv")
x <- data[training,2]

step <- 2 # step size
n <- 30 # sequence length
m <- 5 # m-step Markov path
k <- 6 # number of clusters/states

# initializing the subsequences matrix
st <- 1
en <- st+n-1
tmp <- x[st:en]
tmp<- (tmp-mean(tmp))/sd(tmp)
y<- tmp
# creating subsequences of x
while (en<=length(x)-n){
  st <- st+step
  en <- st+n-1
  tmp <- x[st:en]
  tmp<- (tmp-mean(tmp))/sd(tmp) # normalizing the sequence
  y<- cbind(y,tmp) # concatenation of the new sequence to the matrix
}
# creating dataframe out of the subsequences matrix
df <- as.data.frame(t(y))
row.names(df) <- 1:nrow(df)

# running the k-means clustering and finding the optimal K
clustering <- my_kmeans(df, clusters = k)

max_dist <- numeric(k)
for (i in 1:k) {
  max_dist[i] <- max(clustering$d[clustering$labels==i,i])
}

index <- data.frame()
for (i in 1:(nrow(data)-n+1)){
  newSub <- data[i:(i+n-1),2]
  newSub <- as.data.frame((newSub-mean(newSub))/sd(newSub))
  newSubDist <- as.data.frame(t(sqrt(colSums((newSub[,rep(1,k)]-clustering$centers)^2))))
  colnames(newSubDist) <- sprintf("dist%d",seq(1:k))
  label <- which.min(newSubDist)
  if (newSubDist[label]>1.2*max_dist[label]){
    index <- rbind(index,i)
  }
}

index2 <- numeric(nrow(data))
for (i in 1:nrow(index)) {
  if(index[i,1]>0){
    index2[index[i,1]:(index[i,1]+n-1)] <- index2[index[i,1]:(index[i,1]+n-1)]+1
  }
}

b <- c(1:length(index2))*(index2>0)
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

anomalies <- data.frame(d)

datetime = c(1:nrow(data))
df <- data.frame(datetime, signal = data[,2], group=1)
for(j in 1:nrow(anomalies)){
  for(jj in anomalies[j,'start']:anomalies[j,'end'])
    df[jj,'group'] = j+1
}
df[,'group'] <- factor(df[,'group'])
df2 <- df[df[,'group']!=1,]
library(ggplot2)
plots<-list()
layout <- matrix(c(1:3), nrow = 3, byrow = TRUE)
plots[[1]]<- ggplot(data=df,aes(x=datetime,y=signal)) + geom_line(color='gray25') +
  geom_line(data=df2,aes(x=datetime,y=signal,group=group),color='red') +
  xlab('index') + ylab('') +
  scale_x_continuous(limits=c(0, 2300)) + ggtitle("ECG dataset") +
  theme(plot.title = element_text(hjust = 0.5))

# For SAS dataset use this:
data <-read.csv("SASdataset.csv")
training <- c(2000:3000) # clean data == No Anomaly
x <- data[training,2]

step <- 2 # step size
n <- 30 # sequence length
m <- 5 # m-step Markov path
k <- 6 # number of clusters/states

# initializing the subsequences matrix
st <- 1
en <- st+n-1
tmp <- x[st:en]
tmp<- (tmp-mean(tmp))/sd(tmp)
y<- tmp
# creating subsequences of x
while (en<=length(x)-n){
  st <- st+step
  en <- st+n-1
  tmp <- x[st:en]
  tmp<- (tmp-mean(tmp))/sd(tmp) # normalizing the sequence
  y<- cbind(y,tmp) # concatenation of the new sequence to the matrix
}
# creating dataframe out of the subsequences matrix
df <- as.data.frame(t(y))
row.names(df) <- 1:nrow(df)

# running the k-means clustering and finding the optimal K
clustering <- my_kmeans(df, clusters = k)

max_dist <- numeric(k)
for (i in 1:k) {
  max_dist[i] <- max(clustering$d[clustering$labels==i,i])
}

index <- data.frame()
for (i in 1:(nrow(data)-n+1)){
  newSub <- data[i:(i+n-1),2]
  newSub <- as.data.frame((newSub-mean(newSub))/sd(newSub))
  newSubDist <- as.data.frame(t(sqrt(colSums((newSub[,rep(1,k)]-clustering$centers)^2))))
  colnames(newSubDist) <- sprintf("dist%d",seq(1:k))
  label <- which.min(newSubDist)
  if (newSubDist[label]>1.2*max_dist[label]){
    index <- rbind(index,i)
  }
}

index2 <- numeric(nrow(data))
for (i in 1:nrow(index)) {
  if(index[i,1]>0){
    index2[index[i,1]:(index[i,1]+n-1)] <- index2[index[i,1]:(index[i,1]+n-1)]+1
  }
}

b <- c(1:length(index2))*(index2>0)
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

anomalies <- data.frame(d)

datetime = c(1:nrow(data))
df <- data.frame(datetime, signal = data[,2], group=1)
for(j in 1:nrow(anomalies)){
  for(jj in anomalies[j,'start']:anomalies[j,'end'])
    df[jj,'group'] = j+1
}
df[,'group'] <- factor(df[,'group'])
df2 <- df[df[,'group']!=1,]

plots[[2]]<- ggplot(data=df,aes(x=datetime,y=signal)) + geom_line(color='gray25') +
  geom_line(data=df2,aes(x=datetime,y=signal,group=group),color='red') +
  xlab('index') + ylab('') +
  scale_x_continuous(limits=c(0, 3850)) + ggtitle("Utility usage dataset") +
  theme(plot.title = element_text(hjust = 0.5))

# For NYC taxi dataset use this:
data <-read.csv("nyc_taxi.csv")
days <- as.POSIXlt(data[,1])$wday
index <- days>=1&days<=5
data <- data[index,]
training <- c(1:2000) # clean data == No Anomaly
x <- data[training,2]

step <- 2 # step size
n <- 30 # sequence length
m <- 5 # m-step Markov path
k <- 6 # number of clusters/states

# initializing the subsequences matrix
st <- 1
en <- st+n-1
tmp <- x[st:en]
tmp<- (tmp-mean(tmp))/sd(tmp)
y<- tmp
# creating subsequences of x
while (en<=length(x)-n){
  st <- st+step
  en <- st+n-1
  tmp <- x[st:en]
  tmp<- (tmp-mean(tmp))/sd(tmp) # normalizing the sequence
  y<- cbind(y,tmp) # concatenation of the new sequence to the matrix
}
# creating dataframe out of the subsequences matrix
df <- as.data.frame(t(y))
row.names(df) <- 1:nrow(df)

# running the k-means clustering and finding the optimal K
clustering <- my_kmeans(df, clusters = k)

max_dist <- numeric(k)
for (i in 1:k) {
  max_dist[i] <- max(clustering$d[clustering$labels==i,i])
}

index <- data.frame()
for (i in 1:(nrow(data)-n+1)){
  newSub <- data[i:(i+n-1),2]
  newSub <- as.data.frame((newSub-mean(newSub))/sd(newSub))
  newSubDist <- as.data.frame(t(sqrt(colSums((newSub[,rep(1,k)]-clustering$centers)^2))))
  colnames(newSubDist) <- sprintf("dist%d",seq(1:k))
  label <- which.min(newSubDist)
  if (newSubDist[label]>1.2*max_dist[label]){
    index <- rbind(index,i)
  }
}

index2 <- numeric(nrow(data))
for (i in 1:nrow(index)) {
  if(index[i,1]>0){
    index2[index[i,1]:(index[i,1]+n-1)] <- index2[index[i,1]:(index[i,1]+n-1)]+1
  }
}

b <- c(1:length(index2))*(index2>0)
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

anomalies <- data.frame(d)

datetime = c(1:nrow(data))
df <- data.frame(datetime, signal = data[,2], group=1)
for(j in 1:nrow(anomalies)){
  for(jj in anomalies[j,'start']:anomalies[j,'end'])
    df[jj,'group'] = j+1
}
df[,'group'] <- factor(df[,'group'])
df2 <- df[df[,'group']!=1,]
plots[[3]]<- ggplot(data=df,aes(x=datetime,y=signal)) + geom_line(color='gray25') +
  geom_line(data=df2,aes(x=datetime,y=signal,group=group),color='red') +
  xlab('index') + ylab('') +
  scale_x_continuous(limits=c(0, 7400)) + ggtitle("New York city taxi demand dataset") +
  theme(plot.title = element_text(hjust = 0.5))

png("WrongAssumption.png", width = 9, height = 6, units = 'in', res = 300)
multiplot(plotlist = plots, layout = layout)
dev.off()