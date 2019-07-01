# Random Forest helper function from rpms package modified for mtry and out of bag error estimateion

rpmsForestt <- function(rp_equ, data, weights=~1, strata=~1, clusters=~1, 
                        e_equ=~1, mtry = NULL, bin_size = NULL, perm_reps=100, pval=.25, 
                        f_size=200, cores=1){
  
  
  if(is.null(bin_size)){
    bin_size <- ceiling(nrow(data)^(1/2))
  } 
  
  #============= format data set ================================================
  varlist <- unique(c(all.vars(rp_equ), all.vars(e_equ), all.vars(weights), 
                      all.vars(strata), all.vars(clusters)))  
  
  #------- check variables against data set ------
  if(length(all.vars(rp_equ))==0) 
    stop("no variable for recursive partition given in formula") else
      if(!all(all.vars(rp_equ) %in% names(data)))
        stop("e_equ contains variables that are not in the data set")
  
  if(length(all.vars(e_equ))==0)
    e_equ<-formula(paste(all.vars(rp_equ)[1], 1, sep="~")) else
      if(all.vars(e_equ)[1]!=all.vars(rp_equ)[1]) 
        stop("Dependent variable must be same for rp_equ as e_equ") else
          if(!all(all.vars(rp_equ) %in% names(data)))
            stop("e_equ contains variables that are not in the data set")
  
  #---------- check all variables are in data set --------------
  if(!all(varlist %in% names(data))){
    e1ind <- which(!(varlist %in% names(data)))[1]
    stop(paste(varlist[e1ind], " not in dataset"))
  }
  
  #-----check all variables are numeric or categorical
  if(length(which(sapply(data[,varlist], 
                         function(x) !(is.factor(x)|is.numeric(x)))))>0){
    stop("RPMS works only on numeric or factor data types.")
  }
  
  #------ recurisive partitioning variables ------
  vX=all.vars(rp_equ)[-1]
  
  
  #--------- reduce data to only y with observations and required variables 
  data <- na.omit(data[,varlist])
  
  #-----------------------------------------------------------------------
  
  #===================== Design Information =================================
  
  des_ind <- c(FALSE, FALSE, FALSE)
  des<-list(weights=~1, strata=~1, clusters=~1)
  
  #------ Sample Weights ----------------------------
  if(is.numeric(weights)) { 
    if(length(weights)!=nrow(data)) stop("Number of strata labels != rows of data")
  } else
    if(length(all.vars(weights))==0) {
      weights <- rep(1, nrow(data))
    } else 
      if(all.vars(weights)[1] %in% names(data))
      {
        des$weights=weights
        weights <- as.numeric(data[,all.vars(weights)[1]])
        if(var(weights)>0) {
          des_ind[1] <- TRUE 
        } 
      } else {stop(paste("Problem with weights:",
                         all.vars(weights)[1], "not in data set"))}
  
  #------ Strata Labels ----------------------------
  if(is.numeric(strata) | is.factor(strata)){
    strata<-as.integer(strata)
    if(length(strata)!=nrow(data)) 
      stop("Number of strata labels != rows of data")
  } else
    if(length(all.vars(strata))==0) {strata <- rep(1L, nrow(data))} else
      if(all.vars(strata)[1] %in% names(data)) {
        des$strata=strata
        des_ind[2] <- TRUE
        strata <- as.integer(data[,all.vars(strata)[1]])} else 
          stop(paste("Problem with strata:",
                     all.vars(strata)[1], "not in data set"))
  
  #------ Cluster Labels ---------------------------- 
  if(is.numeric(clusters) | is.factor(clusters)){
    clusters<-as.integer(clusters)
    if(length(clusters)!=nrow(data)) 
      stop("Number of cluster labels != rows of data")
  } else
    if(length(all.vars(clusters))==0) {clusters <- seq(1L:nrow(data))} else
      if(all.vars(clusters)[1] %in% names(data)) {
        des$clusters <- clusters
        des_ind[3] <- TRUE
        clusters <- as.integer(data[,all.vars(clusters)[1]])} else
          stop(paste("Problem with clusters:",
                     all.vars(clusters)[1], "not in data set")) 
  #===================================================================================  
  
  des_string <- if(sum(des_ind)==0) "Simple Random Sample"
  else paste(c("unequal probability of selection", "stratified", "cluster")[which(des_ind)], "sample design")
  
  y <- all.vars(rp_equ)[1] 
  N<-nrow(data)
  xp <- length(vX) # number of variables
  
  
  tree <- vector(mode="list", length=f_size) # forest
  
  get_trees<-function(){
    
    # ----- randomize tree growth ------------------------------
    
    #--randomize data set ----- 
    s <- sample(N, size= round((2/3)*N), replace = FALSE) #2/3 sample
    oob <- which(!(1:N %in% s)) #out of bag set
    #-------------------------------------------------------
    
    #---randomize variables used  
    fvars <- paste(vX[sample(xp, size=mtry)], collapse="+")
    f_equ <- as.formula(paste0(y, "~", fvars))
    #---------------------------------------------------------
    fvX=all.vars(f_equ)[-1] 
    ####################   handle categorical variables for C++  ########################
    
    # ------------------identify categorical variable ------
    # need to handle length 1 separately
    if(length(fvX)==1) cat_vec <- is.factor(data[,fvX]) 
    else
      cat_vec <- sapply(data[, fvX], FUN=is.factor)
    
    n_cats<-length(which(cat_vec))
    
    #---------- There are categorical variables ------------
    if(n_cats>0){
      
      # ----- function to turn NA into ? category
      fn_q <- function(x){
        #x<-as.integer(x)
        #x[which(is.na(x))]<- (min(x, na.rm=TRUE)-1)
        nas <- which(is.na(x))
        if(length(nas>0)){
          x <- factor(x, levels=c(levels(x), "?"))
          x[nas]<- factor("?")
        }
        return(as.factor(x))
      } # end internal function fn
      # 
      
      # ---------- turn each NA into ? category
      if(length(which(cat_vec))==1) {
        data[,fvX[which(cat_vec)]] <- fn_q(data[,fvX[which(cat_vec)]])
      }
      else{
        data[,fvX[which(cat_vec)]] <- lapply(data[,fvX[which(cat_vec)]], fn_q)
      }
      
      # ----- store original levels for each categorical variable ------------
      cat_table <- list()
      
      # ----- function to turn categories into integers for C++ 
      for(i in 1:n_cats){
        cat_table[[i]] <- levels(data[,fvX[which(cat_vec)[i]]]) # --store original levels in cat_table
        #---- replace original levels with integer -------
        levels(data[,fvX[which(cat_vec)[i]]]) <- (1:length(levels(data[,fvX[which(cat_vec)[i]]])))
      }
      
    } # done with if categorical variables
    #------------------------------------------------------------------------------------
    
    #-------end randomize ------------------------------------
    ti <- rpms(rp_equ=f_equ, data[s,], weights[s], strata[s], clusters[s],
               e_equ=e_equ, perm_reps=perm_reps, pval=pval, bin_size = bin_size)
    ti$inclusion_ind <- s
    if(typeof(data$y) %in% c("numeric", "integer", "double")){
      #MSE
      ti$oob_error <- mean((predict(ti, newdata = data[oob,]) - data$y[oob])^2)
    }
    else{
      #Misclassification rate
      ti$oob_error <- mean(ifelse(predict(ti, newdata = data[oob,]) == data$y[oob], 0, 1))
    }
    return(ti)
    
  }#---end get_trees ---------------------------
  
  tree<-replicate(f_size, get_trees(), simplify = FALSE)
  #find average oob error
  oob_error <- sapply(1:f_size, FUN = function(x) tree[[x]]$oob_error) %>% mean()
  #variable importance
  var_imp <- suppressWarnings(as.expression( #warning from bind_rows
    lapply(1:f_size, function(x) tree[[x]]$frame %>%
             dplyr::select("var", "loss")) %>%
      bind_rows() %>% #rbind each tree's split info
      filter(var != "Root") %>% #ignore roots of each tree
      group_by(var) %>% #summarize decrease of loss function for each var.
      summarise(num_splits = n(), avg_loss = mean(loss))))
  #output
  f1<-list(tree = tree, oob_error = oob_error, var_imp = var_imp)
  class(f1)<-c("mase_rpms_forest")
  
  return(f1)
  
}

predict.mase_rpms_forest <- function(obj, newdata, oob = FALSE) {
  ntree <- length(obj$tree)
  p_matrix <- sapply(1:ntree, FUN = function(x)
    predict(obj$tree[[x]], newdata = newdata)) %>%
    matrix(nrow = ntree, byrow = TRUE)
  if(oob == TRUE){
    incl_mat <- sapply(1:ntree, FUN = function(x)
      1:nrow(newdata) %in% (obj$tree[[x]]$inclusion_ind)) %>%
      matrix(nrow = ntree, byrow = TRUE)
    p_matrix <- p_matrix * incl_mat
    bag_counts <- colSums(incl_mat)
    return(colSums(p_matrix)/bag_counts)
  }
  else{
    return(colSums(p_matrix)/ntree)
  }

}