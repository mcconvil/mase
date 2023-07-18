# Helper function to compute postStrat total for bootstrapping

postStratt <- function(data, xpoptab, indices){

  # Define variables
  x <- N_h <- N_h_hats <- ps_h <- poptotal_h <- strat_pop_total <- NULL  
  
  d <- data[indices,]
  
  #y
  y <- d[,3]
  
  #pis 
  pis <- d[,2]
  
  #xsample
  xsample <- d[, 1]
  
  #Estimator
  #Compute estimator
  dat <- data.frame(x = xsample, pi = pis, y = y)
  
  t <- aggregate(
    y ~ x,
    data = dat,
    FUN = function(t) data.frame(
      poptotal_h = sum(t * pi[1:length(t)]^(-1)),
      N_h_hats = sum(pi[1:length(t)]^(-1)),
      var_h = var(t)
    ),
    simplify = F)
  
  tab <- cbind(x = t[, "x"], do.call(rbind, t$y)) |> 
    merge(xpoptab, by = "x")
  
  tab$ps_h <- tab$N_h/tab$N_h_hats
  tab$strat_pop_total <- tab$ps_h * tab$poptotal_h
  tab$strat_pop_mean <- tab$strat_pop_total/tab$N_h
  
  # tab <- dat %>%
  #   group_by(x) %>%
  #   summarize(poptotal_h = y%*%pi^(-1),N_h_hats = sum(pi^(-1)), var_h = var(y)) %>%
  #   inner_join(xpoptab, by=c("x")) %>%
  #   mutate(ps_h = N_h/N_h_hats, strat_pop_total = ps_h*poptotal_h, strat_pop_mean = strat_pop_total/N_h)
  
  
  #calculating the total estimate for y
  out <- sum(tab$strat_pop_total)
  
  return(out)
}

