#Helper function to compute linear GREG total for bootstrapping

logisticGregt <- function(data, x_pop_d, indices){
  #data: 1st column:y, 2nd column:weight, rest: x_sample_d
  d <- data[indices,]
  
  f <- paste(names(d)[1], "~", paste(names(d)[-c(1,2)], collapse=" + "))
  s_design<-survey::svydesign(ids=~1,weights=~weight,data=d)
  
  mod <- survey::svyglm(f, design=s_design,family=quasibinomial())
  
  y_hats_U <- as.matrix(predict(mod,newdata=data.frame(x_pop_d[,-1]),type="response",family=quasibinomial()))
  y_hats_s <- as.matrix(predict(mod,type="response",family=quasibinomial()))
  
  
  #Return estimator of total
  return(t(d$y-y_hats_s)%*%d$weight + sum(y_hats_U))
}

