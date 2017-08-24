#Helper function to compute postStrat total for bootstrapping

postStratt <- function(data, xpoptab, indices){

  #Define variables
  x= N_h = N_h_hats = ps_h = poptotal_h = strat_pop_total = NULL  
  
    #data: 3rd column:y, 2nd column:pis, 1st column: xsample
  d <- data[indices,]
  
  #y
  y <- d[,3]
  
  #pis 
  pis <- d[,2]
  
  #xsample
  xsample <- d[, 1]
  
  #Estimator
  #Compute estimator
  dat <- data.frame(x=xsample, pi=pis,y=y)
  tab <- dat %>%
    group_by(x) %>%
    summarize(poptotal_h = y%*%pi^(-1),N_h_hats = sum(pi^(-1)), var_h = var(y)) %>%
    inner_join(xpoptab, by=c("x")) %>%
    mutate(ps_h = N_h/N_h_hats, strat_pop_total = ps_h*poptotal_h, strat_pop_mean = strat_pop_total/N_h)
  
  
  #calculating the total estimate for y
  out <- sum(tab$strat_pop_total)
  
  return(out)
}

