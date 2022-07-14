# forward feature selection function with cross validation and RMSE


ffs_cv_rmse <- function(data,
                        dep,            # the variable to be predicted
                        vars,# the predictor variables
                        random,# term for random effects     
                        nfolds = 10){      # number of folds used for cv
  
  library(CAST)
  library(dplyr)
  
  folds <- CreateSpacetimeFolds(data,
                                spacevar = "cst_id", 
                                k = nfolds)
  
  
  # 1 - no variables are selected yet
  selected_vars <- NULL
  # therefore the model is calculated by using each variable individually
  
  
  cv <- lapply(seq(nfolds), function(i){
    
    # create training and test data set
    
    train <- data[folds[[1]][[i]],]
    test <- data[folds[[2]][[i]],]
    
    fwd_fs <- lapply(seq(length(vars)), function(v){
      
      # calculate model
      formula <- paste(dep, " ~ ", vars[v],"+", random)
      lme <- lmer(formula, data = train)
      #print(formula)
      
      pred <- predict(lme, newdata = test, allow.new.levels = TRUE)
      obsv <- test$temp
      resid <- obsv - pred
      
      # organize results and performance measure
      results <- data.frame(variable   = vars[v],
                            run = i,
                            pred = pred,
                            obsv = obsv,
                            resid = resid)
      return(results)
    })
    
    fwd_fs <- do.call("rbind", fwd_fs)
    
    })
    
    cv <- do.call("rbind", cv)
    
  rmse <- cv %>%  
    group_by(variable) %>% 
    summarize(rmse =  sqrt(mean((pred - obsv)^2, na.rm = TRUE)))
  
    
    # based on the result I can now define the first best variable (which has the lowest AIC)
    
    rmse_current <- min(rmse$rmse)
    
    best_var <- as.character(rmse$variable[which(rmse$rmse == min(rmse$rmse))])
    
    # this best variable is added to the selected variables and dropped from vars
    
    selected_vars <- best_var
    vars <- vars[-which(vars ==best_var)]
    
    # and I start the performance data.frame
    
    performance <- data.frame(n = 1,
                              predictor  = best_var,
                              RMSE        = rmse_current)
    
    # 2 - now we can run a loop over the next variable selections, 
    # this will run until the new lowest AIC is higher than the previous minimum AIC
    
    finished = FALSE
    n = 1
    while (finished == FALSE) {
      
      cv <- lapply(seq(nfolds), function(i){
        
        # create training and test data set
        
        train <- data[folds[[1]][[i]],]
        test <- data[folds[[2]][[i]],]
        
        fwd_fs <- lapply(seq(length(vars)), function(v){
          
          # calculate model
          formula <- paste(dep, " ~ ", 
                           paste(paste(selected_vars, collapse = " + "), 
                             vars[v],random, sep = " + "))
          lme <- lmer(formula, data = train)
          #print(formula)
          
          pred <- predict(lme, newdata = test, allow.new.levels = TRUE)
          obsv <- test$temp
          resid <- obsv - pred
          
          # organize results and performance measure
          results <- data.frame(variable   = vars[v],
                                run = i,
                                pred = pred,
                                obsv = obsv,
                                resid = resid)
          return(results)
        })
        
        fwd_fs <- do.call("rbind", fwd_fs)
        
      })
      
      cv <- do.call("rbind", cv)
      
      rmse <- cv %>%  
        group_by(variable) %>% 
        summarize(rmse =  sqrt(mean((pred - obsv)^2, na.rm = TRUE)))
    
      
      # based on the result I can now define the next best variable (which has the lowest AIC)
      
      rmse_new <- min(rmse$rmse)
      
      best_var <- as.character(rmse$variable[which(rmse$rmse == min(rmse$rmse))])
      
      # this best variable is added to the selected variables, if the new AIC is lower than the previous
      
      if (rmse_new < rmse_current) {
        selected_vars <- c(selected_vars, best_var)
        vars <- vars[-which(vars == best_var)]
        n <- n + 1
        
        tmp <- data.frame(n = n,
                          predictor = best_var,
                          RMSE = rmse_new)
        
        performance <- rbind(performance, tmp)
        rmse_current <- rmse_new
        
      } else {
        
        print(paste("The optimum RMSE has been reached. The selected variables are:", 
                    paste(selected_vars, collapse = ", ")))
        finished = TRUE
      }
      
      if (length(vars) == 0) {
        print("All variables have been included to the model")
        finished = TRUE
      }
      
    } # end of while statement
    
    return(performance)
}
