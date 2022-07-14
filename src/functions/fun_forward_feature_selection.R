# forward feature selection function

forward_feature_selection <- function(data,
                                      dep,            # the variable to be predicted
                                      vars,# the predictor variables
                                      random,# term for random effects     
                                      measure = "AIC"){      # either AIC or BIC  
  
  
  if (measure == "AIC") {
    
    # 1 - no variables are selected yet
    selected_vars <- NULL
    # therefore the model is calculated by using each variable individually
    
    fwd_fs <- lapply(seq(length(vars)), function(v){
      
      # calculate model
      formula <- paste(dep, " ~ ", vars[v],"+", random)
      lme <- lmer(formula, data = data)
      print(formula)
      
      # organize results and performance measure
      results <- data.frame(Variable   = vars[v],
                            #Adj_R_sqrd = round(summary(lme)$adj.r.squared, 4),
                            AIC        = AIC(lme))
      return(results)
    })
    
    fwd_fs <- do.call("rbind", fwd_fs)
    
    # based on the result I can now define the first best variable (which has the lowest AIC)
    
    AIC_current <- min(fwd_fs$AIC)
    
    best_var <- as.character(fwd_fs$Variable[which(fwd_fs$AIC == min(fwd_fs$AIC))])
    
    # this best variable is added to the selected variables and dropped from vars
    
    selected_vars <- best_var
    vars <- vars[-which(vars ==best_var)]
    
    # and I start the performance data.frame
    
    performance <- data.frame(n = 1,
                              predictor  = best_var,
                              AIC        = AIC_current)
    
    # 2 - now we can run a loop over the next variable selections, 
    # this will run until the new lowest AIC is higher than the previous minimum AIC
    
    finished = FALSE
    n = 1
    while (finished == FALSE) {
      
      
      fwd_fs <- lapply(seq(length(vars)), function(v){
        
        # calculate model
        formula <- paste(dep, " ~ ", 
                         paste(
                           paste(
                             selected_vars, collapse = " + "), 
                           vars[v],random, sep = " + "))
        
        lme <- lmer(formula, data = data)
        print(formula)
        
        # organize results and performance measure
        results <- data.frame(Variable   = vars[v],
                              #Adj_R_sqrd = round(summary(lmod)$adj.r.squared, 4),
                              AIC        = AIC(lme))
        return(results)
      })
      
      fwd_fs <- do.call("rbind", fwd_fs)
      
      # based on the result I can now define the next best variable (which has the lowest AIC)
      
      AIC_new <- min(fwd_fs$AIC)
      
      best_var <- as.character(fwd_fs$Variable[which(fwd_fs$AIC == min(fwd_fs$AIC))])
      
      # this best variable is added to the selected variables, if the new AIC is lower than the previous
      
      if (AIC_new < AIC_current) {
        selected_vars <- c(selected_vars, best_var)
        vars <- vars[-which(vars == best_var)]
        n <- n + 1
        
        tmp <- data.frame(n = n,
                          predictor = best_var,
                          AIC = AIC_new)
        
        performance <- rbind(performance, tmp)
        AIC_current <- AIC_new
        
      } else {
        
        print(paste("The optimum AIC has been reached. The selected variables are:", 
                    paste(selected_vars, collapse = ", ")))
        finished = TRUE
      }
      
      if (length(vars) == 0) {
        print("All variables have been included to the model")
        finished = TRUE
      }
      
    } # end of while statement
    
  } else if (measure == "BIC") {
    
    # 1 - no variables are selected yet
    selected_vars <- NULL
    # therefore the model is calculated by using each variable individually
    
    fwd_fs <- lapply(seq(length(vars)), function(v){
      
      # calculate model
      formula <- paste(dep, " ~ ", vars[v],"+", random)
      lme <- lmer(formula, data = data)
      print(formula)
      
      # organize results and performance measure
      results <- data.frame(Variable   = vars[v],
                            #Adj_R_sqrd = round(summary(lme)$adj.r.squared, 4),
                            BIC        = BIC(lme))
      return(results)
    })
    
    fwd_fs <- do.call("rbind", fwd_fs)
    
    # based on the result I can now define the first best variable (which has the lowest BIC)
    
    BIC_current <- min(fwd_fs$BIC)
    
    best_var <- as.character(fwd_fs$Variable[which(fwd_fs$BIC == min(fwd_fs$BIC))])
    
    # this best variable is added to the selected variables and dropped from vars
    
    selected_vars <- best_var
    vars <- vars[-which(vars ==best_var)]
    
    # and I start the performance data.frame
    
    performance <- data.frame(n = 1,
                              predictor  = best_var,
                              BIC        = BIC_current)
    
    # 2 - now we can run a loop over the next variable selections, 
    # this will run until the new lowest BIC is higher than the previous minimum BIC
    
    finished = FALSE
    n = 1
    while (finished == FALSE) {
      
      
      fwd_fs <- lapply(seq(length(vars)), function(v){
        
        # calculate model
        formula <- paste(dep, " ~ ", 
                         paste(
                           paste(
                             selected_vars, collapse = " + "), 
                           vars[v],random, sep = " + "))
        
        lme <- lmer(formula, data = data)
        print(formula)
        
        # organize results and performance measure
        results <- data.frame(Variable   = vars[v],
                              #Adj_R_sqrd = round(summary(lmod)$adj.r.squared, 4),
                              BIC        = BIC(lme))
        return(results)
      })
      
      fwd_fs <- do.call("rbind", fwd_fs)
      
      # based on the result I can now define the next best variable (which has the lowest BIC)
      
      BIC_new <- min(fwd_fs$BIC)
      
      best_var <- as.character(fwd_fs$Variable[which(fwd_fs$BIC == min(fwd_fs$BIC))])
      
      # this best variable is added to the selected variables, if the new BIC is lower than the previous
      
      if (BIC_new < BIC_current) {
        selected_vars <- c(selected_vars, best_var)
        vars <- vars[-which(vars == best_var)]
        n <- n + 1
        
        tmp <- data.frame(n = n,
                          predictor = best_var,
                          BIC = BIC_new)
        
        performance <- rbind(performance, tmp)
        BIC_current <- BIC_new
        
      } else {
        
        print(paste("The optimum BIC has been reached. The selected variables are:", 
                    paste(selected_vars, collapse = ", ")))
        finished = TRUE
      }
      
      if (length(vars) == 0) {
        print("All variables have been included to the model")
        finished = TRUE
      }
      
    } # end of while statement
    
    
  }
  
  return(performance)
}