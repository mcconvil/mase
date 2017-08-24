#Helper function to compute linear GREG total for bootstrapping

logisticGregt <- function(data, xpopd, indices){
  #data: 1st column:y, 2nd column:weight, rest: xsample_d
  d <- data[indices,]
  
  f <- paste(names(d)[1], "~", paste(names(d)[-c(1,2)], collapse=" + "))
  s.design<-survey::svydesign(ids=~1,weights=~weight,data=d)
  
  mod <- survey::svyglm(f, design=s.design,family=quasibinomial())
  
  y.hats.U <- as.matrix(predict(mod,newdata=data.frame(xpopd[,-1]),type="response",family=quasibinomial()))
  y.hats.s <- as.matrix(predict(mod,type="response",family=quasibinomial()))
  
  
  #Return estimator of total
  return(t(d$y-y.hats.s)%*%d$weight + sum(y.hats.U))
}

