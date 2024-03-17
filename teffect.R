teffect_ra <- function(formula, data, treatment.effect = "ATE", bootstrap_se = F, ipw = F, link = "probit", iter = 1000, alpha = 0.95, seed = 123) {
  # Dissect fml
  Y <- as.character(formula)[2]
  X <- unlist(strsplit(gsub(" ","", as.character(formula)[3]),"\\+"))
  
  D <- X[1]
  X <- X[-1]
  
  # Prep data
  dat <- data[c(Y, D, X)]
  dat <- na.omit(dat)
  
  # Prep treatment and outcome functions
  yf <- as.formula(paste0(Y, " ~ ", paste0(X, collapse = " + ")))
  xf <- as.formula(paste0(D, " ~ ", paste0(X, collapse = " + ")))
  
  # Estimate IPW weights
  if (ipw) {
    pscore <- predict(glm(xf, data = dat, family = binomial(link = link)))
    ipwts <- ifelse(dat[[D]] == 1, 1/pscore, 1/(1-pscore))
  }
  
  # Fit model
  fit_ra <- function(df) {
    # Fit model following https://clas.ucdenver.edu/marcelo-perraillon/sites/default/files/attached-files/matching_teffects_code_perraillon_0.do
    
    # Set weights (either ipw or constant weights)
    if (ipw) {
      df$wts <- ipwts
    } else {
      df$wts <- rep(1,nrow(df))
    }
    
    # Estimate potential outcome means
    if (treatment.effect == "ATE") {
      yhat_t <- predict(lm(yf, data = df[df[[D]] == 1,], weights = wts))
      yhat_c <- predict(lm(yf, data = df[df[[D]] == 0,], weights = wts))
      
      pom_t <- mean(yhat_t)
      pom_c <- mean(yhat_c)
      
    } else if (treatment.effect == "ATT") {
      yhat_atet <- predict(lm(yf, data = df[df[[D]] == 0,], weights = wts), newdata = df[df[[D]] == 1,])
      
      pom_t <- mean(df[df[[D]] == 1,][[Y]])
      pom_c <- mean(yhat_atet)
    } else {
      message("Error: {treatment.effect} unrecognised.")
      break
    }
    
    # E[\tau_i]  = E[Y_i(1)] - E[Y_i(0)] 
    return(pom_t - pom_c)
  }
  
  out <- list()
  out$estimate <- fit_ra(dat)
  
  # Bootstrap SEs
  if (bootstrap_se) {
    bootstrap_statistic <- function(data, indices) {
      resampled_data <- data[indices, ]
      return(fit_ra(resampled_data)) 
    }
    
    require(boot)
    set.seed(seed)
    bootstrap <- boot::boot(data = dat, statistic = bootstrap_statistic, R = iter)
    
    out$std.err <- as.numeric(summary(bootstrap)[4])
    out$conf.low <- boot::boot.ci(bootstrap, type = "perc", conf = alpha)$percent[4]
    out$conf.high <- boot::boot.ci(bootstrap, type = "perc", conf = alpha)$percent[5]
  }
  
  return(out)
}
