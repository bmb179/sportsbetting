# Install required packages if not already installed
if (!require('quadprog')) install.packages('quadprog', dependencies = TRUE)
options(scipen = 999)

library(quadprog)
library(plotly)
options(scipen = 999)

american.odds.conv <- function(american.odds) {
  ifelse(american.odds > 0, 
         (american.odds / 100) + 1, 
         (100 / abs(american.odds)) + 1)
}

vline <- function(x = 0, color = 'black') {
  list(
    type = 'line',
    y0 = 0,
    y1 = 1,
    yref = 'paper',
    x0 = x,
    x1 = x,
    line = list(color = color, dash = 'dot')
  )
}

interactive.plot <- function(team.weights, team.names, vline.num) {
  tooltip.text <- apply(team.weights[, team.names], 1, function(weights) {
    # Convert weights to percentages rounded to 2 decimal places
    percent.weights <- paste0(team.names,': ', round(weights * 100, 2), '%', collapse = '<br>')
    percent.weights
  })
  
  plot_ly(x = team.weights$risk, 
          y = team.weights$returns, 
          type = 'scatter', 
          mode = 'lines+markers',
          text = tooltip.text,
          hoverinfo = 'text+x+y',
          line = list(color = 'blue'),
          marker = list(color = 'blue')) %>%
    layout(
      title = 'Interactive Efficient Frontier',
      xaxis = list(title = 'Portfolio Risk (Standard Deviation)'),
      yaxis = list(title = 'Portfolio Return'),
      shapes = vline(vline.num)
    )
}

optimize.bets <- function(bet.data) {
  required.cols <- c('names', 'sharp.odds', 'american.odds')
  if (!all(required.cols %in% colnames(bet.data))) {
    stop('Input must contain "names", "sharp.odds", and "american.odds" columns')
  }
  
  bet.data$soft.odds <- american.odds.conv(bet.data$american.odds)
  bet.data$implied.probability <- 1 / bet.data$sharp.odds
  bet.data$expected.returns <- bet.data$soft.odds - 1
  
  bet.data$variances <- ((bet.data$sharp.odds - 1)^2 * 
                           bet.data$implied.probability * 
                           (1 - bet.data$implied.probability))
  
  n <- nrow(bet.data)
  covariance.matrix <- diag(bet.data$variances)
  
  portfolio.variance <- function(weights) {
    t(weights) %*% covariance.matrix %*% weights
  }
  
  target.returns <- seq(min(bet.data$expected.returns), 
                        max(bet.data$expected.returns), 
                        length.out = 50)
  
  efficient.weights <- vector('list', length(target.returns))
  efficient.variances <- numeric(length(target.returns))
  
  for (i in seq_along(target.returns)) {
    ones <- rep(1, n)
    Dmat <- covariance.matrix
    dvec <- rep(0, n)
    Amat <- cbind(ones, bet.data$expected.returns)
    bvec <- c(1, target.returns[i])
    
    result <- solve.QP(Dmat, dvec, Amat, bvec, meq = 2)
    
    efficient.weights[[i]] <- result$solution
    efficient.variances[i] <- portfolio.variance(result$solution)
  }
  
  team.weights <- data.frame(do.call(rbind, efficient.weights))
  colnames(team.weights) <- bet.data$names
  team.weights$risk <- sqrt(efficient.variances)
  team.weights$returns <- target.returns
  
  # Filter out rows with negative weights
  positive.weights <- team.weights[apply(team.weights[,bet.data$names] >= 0, 1, all),]
  xval <- positive.weights[which.max(positive.weights$returns / positive.weights$risk), ]$risk
  plotly.plot <- interactive.plot(positive.weights, bet.data$names, xval)
  
  return(list(
    bet.data = bet.data,
    weights = positive.weights,
    covariance_matrix = covariance.matrix,
    plot = plotly.plot
  ))
}

# Example usage
my.bets <- data.frame(
  names = c('WAS', 'LAR', 'BUF', 'NYJ'),
  sharp.odds = c(5.1, 3.5, 2.02, 5),
  american.odds = c(410, 250, 102, 400)
)

result <- optimize.bets(my.bets)
print(result$weights)
result$plot

# Max Risk Adjusted Returns
risk.adjusted <- result$weights[which.max(result$weights$returns / result$weights$risk), ]
risk.adjusted

