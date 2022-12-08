library(stringr)
library(readr)

# Import Dataset ----------------------------------------------------------
## Import Time Series Dataset
df<-data_preprocessed  # Import Dataset
  
# Create List of Dataframes -----------------------------------------------
num.columns.original.df<-ncol(df)

# Select number of feedback points ----------------------------------------------------------
number.feedback.points<-9  # years

# create names
# Temperature
names.feedback.points<-rep("",number.feedback.points)

for (j in 1:number.feedback.points){
  #print(i)
  names.feedback.points[j]<-paste("carbon.m",j,sep="")
}

# create feedback
data<-df[(number.feedback.points+1):nrow(df),]
data[,names.feedback.points]<-NA

for (j in 1:nrow(data)){
  data[j,(num.columns.original.df+1):(num.columns.original.df+number.feedback.points)] <- 
    df$carbon[(j-1+number.feedback.points):j]
}

write.csv(data,"carbon.with.feedback.9.csv",row.names=FALSE)