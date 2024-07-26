# load zoo library to apply rolling mean function as it already sums across columns(stocks)
library(zoo)
# load quad prog library since we want to solve a quadratic equation given in terms of weights, target return the covariance matrix of returns. 
library(quadprog)
# load tidyverse library since we want to read the stock data
library(tidyverse)    
# load the xts library for dates 
library(xts)

trade_strategy <- function(t0, D, I, x) {
  
  df <- read.csv("~/ST326/DATA/logprices.csv")            # reading in the .csv file
  dates <- as.Date(df[,1])                                # note that I made the first column dates when writing my csv file in Q1
  log_prices <- as.xts(df[,-1], order.by = dates)         # convert it back to an xts object 
  
  # Calculate daily log returns
  log_returns <- na.omit(diff(log_prices))               # remove first observation : NA value
  log_matrix <- as.matrix(log_returns)                   # convert into a matrix
  # Number of assets
  n_assets <- ncol(log_prices)
  # Number of days
  n_days <- ndays(log_returns)
  # Record daily portfolio returns, initialized as NA
  portfolio_daily_returns <- rep(NA, nrow(log_returns))
  
  # Record weights, initialized as NA
  weights_matrix <- matrix(NA, nrow = n_days, ncol = n_assets)
  
  
  # taking rolling means for each stock for given window D 
  rolling_mean_return <- rollapply(log_returns, width = D, FUN = mean , align = 'right', fill=NA)
  
    
    # Adjust the portfolio (update weights)
    # Note that our optimization problem is total ensure capital is fully allocated, to minimise variance and to ensure the portfolio returns satisfy the target return constraint
    # I reference section 4 of the lecture notes where our optimisation problem is equation 4.14 
  
   for (n in 0:(ceiling((nrow(log_returns) - t0) / I) - 1)) {           
      # Adjusting point
      rebalance_point <- t0 + n * I                                   # set day when the portfolio should be adjusted 
      end_point <- min(rebalance_point + I-1, nrow(weights_matrix))   # end of current adjustment period   
      if (rebalance_point <= nrow(log_returns)) {
        start_index <- max(1, rebalance_point - D + 1)                # start of window of historical data used in optimisation. 
        end_index <- rebalance_point                                  # end of window of historical data 
        
        # Ensure valid indexing
        if (start_index <= end_index) {
          current_window_returns <- log_returns[start_index:end_index, ]                 # log returns over a defined window i.e., sample returns  
          current_sigma <- cov(current_window_returns)                                   # captures unconditional sample covariances and variances of assets over window
          epsilon <- 1e-6                                                                # we apply a shrinkage technique to regularize the covariance matrix, 
                                                                                         # note that our shrinkage estimator needs to be small enough, so as to capture accurate covariance dynamics 
          current_sigma <- current_sigma + diag(epsilon, nrow(current_sigma))            # to ensure our covariance matrix is psd
          current_window_mean_vector <- as.vector(rolling_mean_return[rebalance_point,]) # mean return for each stock i on day t0 + nI, needed as term in optimization constraint 2) 
          
          target_mean <- sum(rolling_mean_return[rebalance_point, ], na.rm = TRUE)       # sum of mean return of all 10 stocks 
          target_return <- 0.01 + (x / 10) * target_mean                                 # optimization constraint for target return of the portfolio
          
          dmat <- 2 * current_sigma                                                      # due to how objective function is interpreted in the QP package, we multiply the covariance matrix by 2 
          dvec <-  rep(0, ncol(log_prices))
          amat <- rbind(rep(1, ncol(log_prices)), current_window_mean_vector)            # RHS constraints of 4.14  , 1) applying a coefficient of 1 to each weight , 2) mean returns of each stock (to give t(w) %*% mu )
          bvec <- c(1, target_return)                                                    # LHS constraints of 4.14  , 1) weights sum up to 1 *2) minimum target return to achieve 
          meq <- 1                                                                       # to ensure only the 1st constraint (1st row in amat)  is the equality constraint
          
          # Solve optimization problem
          result <- solve.QP(dmat, dvec, t(amat), bvec, meq)                # amat is transposed to enforce that dot product of this weights vector and the mean return vector must be equal to the target return 
          new_weights <- result$solution                                    # optimal weight solution Wopt
          
          # Update current weights
          for (day in rebalance_point:end_point) {
            weights_matrix[day, ] <- new_weights                            # adjust weight of the portfolio
          }
        }
      }
      
       for (day in rebalance_point:end_point) {                                           
      if (day <= nrow(log_returns)) {
        # Use new_weights if after adjusting point, otherwise use initial weights
        current_weights <- ifelse(day == rebalance_point, new_weights, weights_matrix[t0, ])
        daily_return <- sum(t(current_weights) %*% log_matrix[day,])                        # daily return t(w)*r
        portfolio_daily_returns[day] <- daily_return
      }
    }
  
  }
  # to calculate the sharpe ratio we need excess returns 
  excess_returns <- portfolio_daily_returns
  
  # Calculate average excess return and standard deviation
  average_excess_return <- mean(excess_returns, na.rm = TRUE)
  sd_excess_returns <- sd(excess_returns, na.rm = TRUE)
  
  # Annualize the Sharpe Ratio
  sharpe_ratio_annualized <- (average_excess_return / sd_excess_returns) * sqrt(252)       
  
  strategy <- as.data.frame(cbind(log_returns,weights_matrix,portfolio_daily_returns))
  # Return the log_prices, weights_matrix and daily portfolio returns
  return(list(sharpe_ratio= sharpe_ratio_annualized, strategy= strategy))
}

## Sample code for estimating relevant statistics in Q3

D_range <- seq(60,120, by= 10)
I_range <- seq(60,120, by= 10)
x_range <- seq(0.5,2.5, by=0.1)
t0_range <- seq(250, 2250, by=250)

for (D in D_range) { 
  for (t0 in t0_range) { 
    for (I in I_range) {
      for(x in x_range) {
        strategy_result <- trade_strategy(t0, D, I, x)
        dataframe_result <- strategy_result$strategy
        # Extract weights matrix from the strategy result
        weights_matrix <- dataframe_result[,c("weights_matrix","weights_matrix.1","weights_matrix.2", "weights_matrix.3", "weights_matrix.4", "weights_matrix.5", "weights_matrix.6","weights_matrix.7", "weights_matrix.8", "weights_matrix.9")]
        
        total_change = 0
        number_of_I_periods = 0
        
        # Calculate the total absolute change in weights
        for (n in 1:(ceiling((nrow(weights_matrix) - t0) / I) - 1)) {
          rebalance_point <- t0 + n * I
          if (rebalance_point + I <= nrow(weights_matrix)) {
            weight_change <- abs(weights_matrix[rebalance_point + I, ] - weights_matrix[rebalance_point, ])
            total_change <- total_change + sum(weight_change)
            number_of_I_periods <- number_of_I_periods + 1
          }
        }
        
        # Calculate the average absolute change in weights
        average_change <- if(number_of_I_periods > 0) total_change / number_of_I_periods else NA
        # Append the Sharpe Ratio result to the results data frame
        results <- rbind(results, data.frame(t0 = t0, D = D, I = I, x = x, Sharpe_Ratio = strategy_result$sharpe_ratio, AvsAchg= average_change ))
        print(paste("Processed t0:", t0, "D:", D, "I:", I, "x:", x, "Average Change:", average_change))
      }
    }
  }
}

write.csv(results, "~/ST326/DATA/results.final3.csv", row.names = FALSE)
