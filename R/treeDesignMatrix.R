#Helper function for gregTree
#Written by Daniell Toth

#Creates a linear design matrix for the regression tree
#Useful for creating survey weights
treeDesignMatrix <- function(splits, data){
  x<-rep(0,length=nrow(data)) #first column = x satifies no splits
  
  if(length(splits)<=0) return(x) else
    for(i in 1:length(splits))
      x<-cbind(x, #if x satisfies split 1 in that column else 0 in that column
               eval(parse(text=paste("ifelse(", splits[i], ", 1, 0)")), data),
               deparse.level=0)
    
    #colnames(x)<-paste("E", seq(0:(ncol(x)-1)), sep="")
    return(as.matrix(x[,-1])) 
}