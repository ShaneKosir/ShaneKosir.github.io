library(corrplot)
library(reshape2)
library(mice)
library(plyr)
library(caret)
library(readr)

# Load model dataframe ----
original<-read_csv('Swell Functional Groups.csv')

# Remove irrelevant columns ----
irrelevant<-c('Molecule','Molecular Group','Molecular Formula')
original.downselected<-original[ , -which(names(original) %in% irrelevant)]

# Eliminate columns with low variation ---------------------------------------------------
index.low.variance<-nearZeroVar(original.downselected,freqCut=95/5,uniqueCut=10,
                                saveMetrics=FALSE,names=FALSE)
head(original.downselected[,index.low.variance])
df.low.variance.removed<-subset(original.downselected,select=-c(index.low.variance))

# Pearson correlation ----
cor.results<-cor(df.low.variance.removed,method='pearson',use='complete.obs')
corrplot(cor.results,tl.cex=.5)

# Identify highly correlated variables ----
diag(cor.results)<-NA
cor.results[lower.tri(cor.results)]<-NA
max.cor.90<-subset(melt(cor.results),value>.9)

# Remove highly correlated variables ----
#df.highly.correlated.removed<-subset(df.low.variance.removed,select=-c(max.cor.90[1,2]))
df.highly.correlated.removed<-df.low.variance.removed[,-which(names(df.low.variance.removed) %in% c('# Aro Rings','# Cyclo Rings',
                                                                                       'Avg Cyclo Ring Size','Avg Subs per Aro',
                                                                                       'Avg Len Aro Sub'))]

# Check for missing data ----
md.pattern(df.low.variance.removed)

# Remove rows with missing data----
df.missing.removed<-na.omit(df.low.variance.removed)

# Impute missing data ----
imputed.data<-mice(df.low.variance.removed,m=1,maxit = 50,method='pmm',seed=500)
# m  – Refers to 5 imputed data sets
# maxit – Refers to no. of iterations taken to impute missing values
# method – Refers to method used in imputation. we used predictive mean matching
df.imputed<-complete(imputed.data,1)

# Scale data ----
# minValue <- sapply(df.missing.removed, min, na.rm = TRUE)
# maxValue <- sapply(df.missing.removed, max, na.rm = TRUE)
# df.scaled<-as.data.frame(scale(df.missing.removed,center=minValue,scale = (maxValue-minValue)))

# vif_func ----
vif_func<-function(in_frame,thresh=5,trace=T,...){
  
  require(fmsb)
  
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
  
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  var_names <- names(in_frame)
  for(val in var_names){
    regressors <- var_names[-which(var_names == val)]
    form <- paste(regressors, collapse = '+')
    form_in <- formula(paste(val, '~', form))
    vif_init<-rbind(vif_init, c(val, VIF(lm(form_in, data = in_frame, ...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]))
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(var_names)
  }
  else{
    
    in_dat<-in_frame
    
    #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
    while(vif_max >= thresh){
      
      vif_vals<-NULL
      var_names <- names(in_dat)
      
      for(val in var_names){
        regressors <- var_names[-which(var_names == val)]
        form <- paste(regressors, collapse = '+')
        form_in <- formula(paste(val, '~', form))
        vif_add<-VIF(lm(form_in, data = in_dat, ...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2])))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
    
  }
  
}

# Remove collinear variables and recreate dataframe ----
x<-subset(df.low.variance.removed)
terms.to.keep<-vif_func(in_frame=x,thresh=5,trace=T)
ttk.df<-data.frame(terms.to.keep)
df.collinearity.removed<-df.missing.removed[,terms.to.keep]

# Add original columns for output----
df.out<-cbind(original[,which(names(original) %in% irrelevant)],df.low.variance.removed)

# Remove outliers----
#df.freeze.outliers.removed<-subset(df.out,df.out['Freeze Point, K']<400)

# Save to csv ----
write.csv(file='Swell Functional Groups Low Variance Removed.csv',x=df.out,row.names=FALSE)
