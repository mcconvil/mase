EBLUP <- function(formula, data, domains, pop_data) {
  # find response variable
  response_lab <- (formula)[2] %>% as.character()
  pred_lab <- all.vars(as.formula(formula))[-1L]
  r <- data %>% pull(response_lab) %>% as.vector()
  # fit LMM
  mod <- lme(formula, data, random = formula(paste0("~1|", domains)))
  # find estimates!
  out <- data %>% 
    dplyr::mutate(fixed = (mod$fitted %>% rowSums()),
           resid = (mod$residuals %>% rowSums()),
           response = r) %>%
    group_by_at(domains) %>%
    summarize(n_i = n(),
              srs = mean(response),
              srs_mse = var(response)/n_i,
              sre = mean(fixed),
              s = mean(resid),
              sd_nu = (mod$coefficients$random$DOMAIN) %>% sd()) %>%
    dplyr::mutate(gamma = (sd_nu^2 / (sd_nu^2 + mod$sigma^2/n_i)),
                  eblup = sre + s*gamma,
                  greg = sre + s) %>%
    dplyr::select(domains, n_i, srs, srs_mse, sre, greg, eblup, gamma) %>%
    arrange_at(domains)
  # EBLUP mse estimate things
  #make matrix of pop means
  X_bar <- cbind(rep(1,nrow(pop_data)),
                 pop_data %>% dplyr::select_at(pred_lab) %>%
                   as.matrix()) %>% unname()
  #make matrix of sample means
  x_bar <- cbind(rep(1,nrow(out)), data %>% dplyr::select_at(c(pred_lab, domains)) %>%
                   group_by_at(domains) %>%
                   summarize_all(funs(mean)) %>%
                   dplyr::select_at(pred_lab) %>%
                   as.matrix()) %>% unname()
  #check if population and sample domains are consistent
  if(nrow(X_bar) != nrow(x_bar)){
    stop("number of sample domains not equal to number of population domains, try filtering by domains %in% sample$domains, or vice versa.")
  }
  # make U matrix for each domain for MSE estimate
  makeU <- function(index, out){
    ni <- max(out[index,2], 1)
    g <- out[index,6]
    m1 <- diag(1, ni) %>% as.numeric()
    m2 <- matrix(rep(g/ni, ni^2), nrow = ni) %>% as.numeric()
    return(matrix((m1 - m2)/(mod$sigma^2), nrow = ni))
  }
  U_list <- lapply(1:nrow(out), FUN = makeU, out = out)
  #make list of design matricies
  design_list <- out %>% pull(domains) %>%
    lapply(function(dom) model.matrix(formula,
                            data = filter_at(data, .vars = vars(domains),
                                .vars_predicate = function(x) x == dom)) %>%
             as.matrix() )
  # Find the EBLUP MSE estimate!
  annoying_thing <- lapply(1:nrow(out), function(i) {
    t(design_list[[i]]) %*% U_list[[i]] %*% (design_list[[i]])
    }) %>% Reduce(f = "+", x = .) %>% solve()
  # Final
  out <- out %>%
    mutate(C1 = gamma*(mod$sigma^2/n_i),
           C2 = sapply(1:nrow(out), function(i) {
             t(X_bar[i,] - gamma[i]*x_bar[i,]) %*%
               annoying_thing %*% 
               (X_bar[i,] - gamma[i]*x_bar[i,])
           }),
           C3 = 1,
           eblup_mse = C1 + C2 + 2*C3) %>%
    dplyr::select(domains, n_i, srs, srs_mse, sre, greg, eblup, eblup_mse, gamma)
  return(out)
}
