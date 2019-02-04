# data pre-processing
setwd("YOUR PATH")
source('multiplot.r')
source('detect.r')
source('my_kmeans.r')

# For SAS dataset use this:
# data <-read.csv("SASdataset.csv")
# training <- c(2000:3000) # clean data == No Anomaly
# x <- data[training,2]

# For ECG dataset use this:
# training <- c(500:1500) # clean data == No Anomaly
# data <-read.csv("ECG.csv")
# x <- data[training,2]

# For NYC taxi dataset use this:
# data <-read.csv("nyc_taxi.csv")
# days <- as.POSIXlt(data[,1])$wday
# index <- days>=1&days<=5
# data <- data[index,]
# training <- c(1:2000) # clean data == No Anomaly
# x <- data[training,2]

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

# running the k-means clustering with clusters = k
clustering <- my_kmeans(df, clusters = k)

# running the k-means clustering and finding the optimal K
# w <- rep(0,7)
# for(j in 1:7){
#  j <- k
#  clustering <- my_kmeans(df, clusters = j)
#  w[j] <- clustering$wss
#}

# after observing wss, optimal K is determined

# training the Markov chain
path <- clustering$labels
P <- matrix(0,nrow=k,ncol=k)
Pi <- rep(0, k)
for(i in 1:(length(path)-1)){
  P[path[i],path[i+1]] <- P[path[i],path[i+1]]+1
  Pi[path[i]] <- Pi[path[i]]+1
}
P <- P/rowSums(P)
Pi[path[length(path)]] <- Pi[path[length(path)]]+1
Pi <- Pi/length(path)

# calling the library caret and building the LDA classifier
library(caret)
feat <- data.frame(labels = factor(clustering$labels),clustering$d)
# downsampling to balance the class sizes
newfeat<-downSample(feat,feat$labels)
# dropping the Class Variable
newfeat$Class<-NULL

library(MASS)
lda.fit <- lda(labels~.,newfeat)

# scoring new subsequences based on lda.fit and Markov path probability
# outcomes is the data frame whose rows are all combinations of states
ProbS <- data.frame()
IndexS <- data.frame()
scoredf <- data.frame()
for (i in 1:m){
  newSub <- data[i:(i+n-1),2]
  newSub <- as.data.frame((newSub-mean(newSub))/sd(newSub))
  newSubDist <- as.data.frame(t(sqrt(colSums((newSub[,rep(1,k)]-clustering$centers)^2))))
  colnames(newSubDist) <- sprintf("dist%d",seq(1:k))
  newSubIndex <- sort(predict(lda.fit, newSubDist)$posterior,index.return=TRUE,decreasing = TRUE)$ix[1:2]
  newSubProb <- sort(predict(lda.fit, newSubDist)$posterior,index.return=TRUE,decreasing = TRUE)$x[1:2]
  IndexS <- rbind(IndexS,newSubIndex)
  ProbS <- rbind(ProbS,newSubProb)
}
colnames(IndexS) <- sprintf("top%d",seq(1:2))
colnames(ProbS) <- sprintf("top%d",seq(1:2))

for (i in (m+1):(nrow(data)-n+1)){
  newSub <- data[i:(i+n-1),2]
  newSub <- as.data.frame((newSub-mean(newSub))/sd(newSub))
  newSubDist <- as.data.frame(t(sqrt(colSums((newSub[,rep(1,k)]-clustering$centers)^2))))
  colnames(newSubDist) <- sprintf("dist%d",seq(1:k))
  newSubIndex <- sort(predict(lda.fit, newSubDist)$posterior,index.return=TRUE,decreasing = TRUE)$ix[1:2]
  newSubProb <- sort(predict(lda.fit, newSubDist)$posterior,index.return=TRUE,decreasing = TRUE)$x[1:2]
  for (j in 1:(m-1)){
    IndexS[j,] <- IndexS[(j+1),]
    ProbS[j,] <- ProbS[(j+1),]
  }
  IndexS[m,] <- newSubIndex
  ProbS[m,] <- newSubProb
  
  weightProb <- expand.grid(s1=as.numeric(ProbS[1,]), s2=as.numeric(ProbS[2,]), s3=as.numeric(ProbS[3,]), s4=as.numeric(ProbS[4,]), s5=as.numeric(ProbS[5,]))
  MarkovProb <- expand.grid(s1=as.numeric(IndexS[1,]), s2=as.numeric(IndexS[2,]), s3=as.numeric(IndexS[3,]), s4=as.numeric(IndexS[4,]), s5=as.numeric(IndexS[5,]))
  ProbPath <- data.frame(Pi[MarkovProb[,1]])
  for (j in 2:ncol(MarkovProb)){
    tmp <- diag(P[MarkovProb[,(j-1)],MarkovProb[,j]])
    ProbPath <- cbind(ProbPath,tmp)
  }
  ProbPath[,2:ncol(ProbPath)] <- ProbPath[,2:ncol(ProbPath)] +1e-4
  score <- 0
  for(j in 1:nrow(ProbPath)){
    score <- score + log10(prod(ProbPath[j,]))*(prod(weightProb[j,]))
  }
  scoredf <- rbind(scoredf,score)
}
colnames(scoredf) <- "score"
threshold_soft <- quantile(scoredf[training,], 0.001)
threshold_hard <- threshold_soft*1.2
anomalies_soft <- detect(scoredf[,1], threshold_soft, n)
anomalies <- detect(scoredf[,1], threshold_hard, n)

anomalies_soft <- data.frame(anomalies_soft, timeID = data[anomalies_soft[,'start'],1])
anomalies <- data.frame(anomalies, timeID = data[anomalies[,'start'],1])

datetime = c(1:nrow(data))
df <- data.frame(datetime, signal = data[,2], group=1)
for(j in 1:nrow(anomalies)){
  for(jj in anomalies[j,'start']:anomalies[j,'end'])
    df[jj,'group'] = j+1
}
df[,'group'] <- factor(df[,'group'])
df2 <- df[df[,'group']!=1,]
library(ggplot2)
scoredf <- data.frame(time = c(1:nrow(scoredf)), scoredf)
plots<-list()
# For SAS dataset use this:
plots[[1]]<-ggplot(data=scoredf,aes(x=time,y=score)) + geom_line(color='gray25') +
  geom_line(aes(x=time,y=threshold_soft),linetype = 'twodash',color='blue',size=0.5) +
  geom_line(aes(x=time,y=threshold_hard),linetype = 'twodash',color='red',size=0.5) +
  xlab('index') + ylab('') + ggtitle("Expected Log-likelihood") +
  scale_x_continuous(limits=c(0, 3850)) +
  theme(plot.title = element_text(hjust = 0.5))

plots[[2]]<-ggplot(data=df,aes(x=datetime,y=signal)) + geom_line(color='gray25') +
  geom_line(data=df2,aes(x=datetime,y=signal,group=group),color='red') +
  xlab('index') + ylab('') +
  scale_x_continuous(limits=c(0, 3850)) +
  scale_y_continuous(limits=c(0, 2500)) + ggtitle("Anomaly segments") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 674, y = 300, label = "JUN 20", size = 3) +
  annotate("text", x = 1095, y = 300, label = "JUL 4", size = 3) +
  annotate("text", x = 1815, y = 300, label = "JUL 25", size = 3) +
  annotate("text", x = 3229.5, y = 300, label = "SEP 4", size = 3) 

layout <- matrix(c(1,2), nrow = 2, byrow = TRUE)
png("SASdataset.png", width = 9, height = 6, units = 'in', res = 300)
multiplot(plotlist = plots, layout = layout)
dev.off()
# finish with SAS dataset
# For ECG dataset use this:
plots[[1]]<-ggplot(data=scoredf,aes(x=time,y=score)) + geom_line(color='gray25') +
  geom_line(aes(x=time,y=threshold_soft),linetype = 'twodash',color='blue',size=0.5) +
  geom_line(aes(x=time,y=threshold_hard),linetype = 'twodash',color='red',size=0.5) +
  xlab('index') + ylab('') + ggtitle("Expected Log-likelihood") +
  scale_x_continuous(limits=c(0, 2300)) +
  theme(plot.title = element_text(hjust = 0.5))

plots[[2]]<-ggplot(data=df,aes(x=datetime,y=signal)) + geom_line(color='gray25') +
  geom_line(data=df2,aes(x=datetime,y=signal,group=group),color='red') +
  xlab('index') + ylab('') +
  scale_x_continuous(limits=c(0, 2300)) +
  ggtitle("Anomaly segments") +
  theme(plot.title = element_text(hjust = 0.5))

layout <- matrix(c(1,2), nrow = 2, byrow = TRUE)
png("ECGdataset.png", width = 9, height = 6, units = 'in', res = 300)
multiplot(plotlist = plots, layout = layout)
dev.off()
# finish with ECG dataset
# For nyc taxi dataset use this:
df <- data.frame(datetime, signal = data[,2], group=1)
for(j in 1:nrow(anomalies_soft)){
  for(jj in anomalies_soft[j,'start']:anomalies_soft[j,'end'])
    df[jj,'group'] = j+1
}
df[,'group'] <- factor(df[,'group'])
df3 <- df[df[,'group']!=1,]

plots[[1]]<-ggplot(data=scoredf,aes(x=time,y=score)) + geom_line(color='gray25') +
  geom_line(aes(x=time,y=threshold_soft),linetype = 'twodash',color='blue',size=0.5) +
  geom_line(aes(x=time,y=threshold_hard),linetype = 'twodash',color='red',size=0.5) +
  xlab('index') + ylab('') + ggtitle("Expected Log-likelihood") +
  scale_x_continuous(limits=c(0, 7400)) +
  theme(plot.title = element_text(hjust = 0.5))

plots[[2]]<-ggplot(data=df,aes(x=datetime,y=signal)) + geom_line(color='gray25') +
  geom_line(data=df3,aes(x=datetime,y=signal,group=group),color='blue') +
  geom_line(data=df2,aes(x=datetime,y=signal,group=group),color='red') +
  xlab('index') + ylab('') +
  scale_x_continuous(limits=c(0, 7400)) +
  scale_y_continuous(limits=c(-10000, 32000)) +
  ggtitle("Anomaly segments") +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 160, y = -2000, label = "JUL 3-4", size = 3) +
  annotate("text", x = 0.5*(2092+2158), y = -2000, label = "AUG 29", size = 3) +
  annotate("text", x = 0.5*(5119+5183), y = -2000, label = "NOV 26-27", size = 3) +
  annotate("text", x = 0.5*(6025+6431), y = -2000, label = "DEC 23-26", size = 3) +
  annotate("text", x = 0.5*(6025+6431), y = -7000, label = "DEC 29-JAN 2", size = 3) +
  annotate("text", x = 0.5*(6890+7248), y = -2000, label = "JAN 16, 23-26", size = 3)
  
layout <- matrix(c(1,2), nrow = 2, byrow = TRUE)
png("NYCdataset.png", width = 9, height = 6, units = 'in', res = 300)
multiplot(plotlist = plots, layout = layout)
dev.off()
# finish with nyc taxi dataset