#====
#RPA data analysis with randomForest
#====
#====================================
library(randomForest)
library(Metrics)           #RMSE
#library(e1071)              #取1:n 全排列的所有组合
#====================================
setwd("C:\\Users\\Gustavo\\Documents\\Coding\\RStudioProjects\\RF-and-PC\\Rscript Code Format") #Enter the directory you want to work in here
RPAd <- read.csv("individual.data.csv") #Read in data
RPAd <- as.data.frame(RPAd) #Checks to make sure data is in a dataframe
colnames(RPAd) #pull column names of RPAd as a vector
#length(RPAd[,4])
savenames <- colnames(RPAd[4:13]) #save a subset of the column names as a new vector
m <- data.frame(0) # create empty dataframe
m[1:179,] <- 0 #create 179 empty rows in m dataframe 
pa <- data.frame(z=o)

#-------------------------------------------------
#ten-fold cross-validation model performance with overall dataset
#-------------------------------------------------
set.seed(12345);disorder <- sample(length(RPAd[,2]),replace=F) #extracting a random sample of rows from the second column of RPAD without replacement
for(k in 1:10){
  n <- data.frame(r2=0,rm=0) #create empty dataframe with columns r2 and rm
  o <- disorder[(2*(k-1)+1):(2*k)] #taking a subset of disorder with length 65 #set to length of data divided by 10
  rf.data <- RPAd[-o,] #removing the rows matching with the above subset from the Data and saving as new dataframe
#for loop above is conducting the tenfold cross validation, o is the testing data set and rf.data is the training set, currently is set to a 90/10 training/testing split. #set for loop to column numbers of proteins
  for (i in 24:201){
    RPAdatai <- cbind(rf.data[i],rf.data[,4:13]) #make a dataframe with 1 protein's relative abundance and nanoparticle details
    colnames(RPAdatai) <- c("RPA",savenames) #change column name from protein accession number to "RPA"
    set.seed(12345)
    rf <- randomForest(RPA~., data=RPAdatai, proximity=F, importance=F)
    p <- predict(rf,RPAd[o,4:13]);
    a <- RPAd[o,i];
    pa <- cbind(pa, cbind(p,a));
    r2 <- cor(p,a);
    rm <- rmse(p,a);
    n <- rbind(n,cbind(r2,rm))
  }
  m <- cbind(m,n)
}
write.csv(m,"individual.r2.rmse.10.csv")



#-------------------------------------------------
#factor selection
#ten-fold cross-validation model performance 
#-------------------------------------------------


set.seed(12345);disorder <- sample(length(RPAd[,2]),replace=F)
for(k in 1:10){
  n <- data.frame(r2=0,rm=0)
  o <- disorder[(56*(k-1)+1):(56*k)]
  rf.data <- RPAd[-o,]
  for (i in 1:73){
    RPAdatai <- cbind(rf.data[i+2],rf.data[,79:99])
    colnames(RPAdatai) <- c("RPA",savenames)
    set.seed(12345)
    rf <- randomForest(RPA~., data=RPAdatai, proximity=F, importance=F)
    p <- predict(rf,RPAd[o,79:99]);
    a <- RPAd[o,i+2];
    r2 <- cor(p,a);
    rm <- rmse(p,a);
    n <- rbind(n,cbind(r2,rm))
  }
  m <- cbind(m,n)
}
write.csv(m,"proteins.r2.rmse.10.csv")



#排列组合，递增
#=========================
#factor selection
#=========================
#====================================
setwd("C:\\Users\\Gustavo\\Documents\\Coding\\RStudioProjects\\RF-and-PC\\Rscript Code Format") #Enter the directory you want to work in here
RPAd <- read.csv("class.280.save.csv")
RPAd <- as.data.frame(RPAd)
colnames(RPAd)
#length(RPAd[,4])
savenames <- colnames(RPAd[79:99])
m <- data.frame(0)
m[,1:4] <- 0
n <- data.frame(j=0,l=0, k=0,i=0,r2=0,rm=0)

for (j in 1:3){    #先做到3个因素
  select.num <- t(combn(1:15, j))
  
  for (l in 1:length(select.num[, 1])){
    select.factor <- cbind(RPAd[79:84], RPAd[select.num[l, ]+84])
    colnames(select.factor) <- c(colnames(RPAd[79:84]), colnames(RPAd[select.num[l, ]+84]))
    select.savenames <- colnames(select.factor)
    
    set.seed(12345);disorder <- sample(length(select.factor[,2]),replace=F)
    
    for(k in 1:10){
      o <- disorder[(56*(k-1)+1):(56*k)]
      rf.data <- select.factor[-o,]
      
      for (i in 70:75){
        RPAdatai <- cbind(RPAd[-o,i],rf.data)
        colnames(RPAdatai) <- c("RPA",select.savenames)
        set.seed(12345)
        rf <- randomForest(RPA~., data=RPAdatai, proximity=F, importance=F)
        p <- predict(rf,select.factor[o,]);
        ob <- RPAd[o,i];
        r2 <- cor(p,ob);
        rm <- rmse(p,ob);
        n <- rbind(n,cbind(j,l,k,i,r2,rm))
      }
    }
  }
}
write.csv(n,"proteins.r2.rmse.10.csv")


#人工剔除不重要及重复因素
#============================
setwd("C:\\Users\\Gustavo\\Documents\\Coding\\RStudioProjects\\RF-and-PC\\Rscript Code Format") #Enter the directory you want to work in here
RPAd <- read.csv("class.280.save.csv")
RPAd <- as.data.frame(RPAd)
colnames(RPAd)
#length(RPAd[,4])
savenames <- colnames(RPAd[c(79:80,82:84,87,91,96:98)])
m <- data.frame(0)
m[1:74,] <- 0

#-------------------------------------------------
#ten-fold cross-validation model performance with overall dataset
#-------------------------------------------------
set.seed(12345);disorder <- sample(length(RPAd[,2]),replace=F)
for(k in 1:10){
  n <- data.frame(r2=0,rm=0)
  o <- disorder[(56*(k-1)+1):(56*k)]
  rf.data <- RPAd[-o,]
  for (i in 1:73){
    RPAdatai <- cbind(rf.data[, i+2],rf.data[,c(79:80,82:84,87,91,96:98)])
    colnames(RPAdatai) <- c("RPA",savenames)
    set.seed(12345)
    rf <- randomForest(RPA~., data=RPAdatai, proximity=F, importance=F)
    p <- predict(rf,RPAd[o,c(79:80,82:84,87,91,96:98)]);
    a <- RPAd[o,i+2];
    r2 <- cor(p,a);
    rm <- rmse(p,a);
    n <- rbind(n,cbind(r2,rm))
  }
  m <- cbind(m,n)
}
write.csv(m,"r2.rmse.10f.csv")

