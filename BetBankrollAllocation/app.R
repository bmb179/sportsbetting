library(shiny)
library(quadprog)
library(plotly)
options(scipen = 999)

# Convert American odds to decimal European odds
american.odds.conv <- function(american.odds) {
  ifelse(american.odds >= 0, 
         (american.odds / 100) + 1, 
         (100 / abs(american.odds)) + 1)
}

# Add a dotted line to the plot to highlight optimal risk adjusted returns
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

# Creates an interactive plot showing the efficient frontier
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
  bet.data$implied.probability <- 1 / american.odds.conv(bet.data$sharp.odds)
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

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel('Sportsbetting Bankroll Manager'),
  div(
    strong('Assumptions:'),
    tags$ol(
      tags$li('At the moment a game starts, sharp sportsbooks (Pinnacle, BetFair, etc) provide the closest odds to the actual game probabilities.'),
      tags$li('The soft sportsbooks (DraftKings, FanDuel, etc) may create mispricings by giving better odds to acquire new customers.'),
      tags$li('Bettors have positive expected value if capital allocation optimizes risk-adjusted returns over a large number of mispriced bets.'),
      )
    ),
  div(
    strong('Instructions:'),
    tags$ol(
      tags$li('List 5 moneyline bets in the field below using the American odds (+/- 100) available within a half hour of the game start.'),
      tags$li('Leave unused columns blank.'),
      tags$li('The resulting graph will show the efficient allocation of capital on a curve showing the risk-return relationship.'),
      tags$li('The maximized potential profit relative to risk will be bisected by a dotted line.')
    )
  ),
  fluidRow(
    column(width = 2, textInput('bet1name', 'Bet 1: Name', value = 'WAS')),
    column(width = 2, numericInput('bet1sharp', 'Sharp Odds', value = -110)),
    column(width = 2, numericInput('bet1soft', 'Soft Odds', value = -105))
    ),
  fluidRow(
    column(width = 2, textInput('bet2name', 'Bet 2: Name', value = 'NYY')),
    column(width = 2, numericInput('bet2sharp', 'Sharp Odds', value = +200)),
    column(width = 2, numericInput('bet2soft', 'Soft Odds', value = +210))
    ),
  fluidRow(
    column(width = 2, textInput('bet3name', 'Bet 3: Name', value = 'COL')),
    column(width = 2, numericInput('bet3sharp', 'Sharp Odds', value = +500)),
    column(width = 2, numericInput('bet3soft', 'Soft Odds', value = +530))
    ),
  fluidRow(
    column(width = 2, textInput('bet4name', 'Bet 4: Name', value = NA)),
    column(width = 2, numericInput('bet4sharp', 'Sharp Odds', value = NA)),
    column(width = 2, numericInput('bet4soft', 'Soft Odds', value = NA))
    ),
  fluidRow(
    column(width = 2, textInput('bet5name', 'Bet 5: Name', value = NA)),
    column(width = 2, numericInput('bet5sharp', 'Sharp Odds', value = NA)),
    column(width = 2, numericInput('bet5soft', 'Soft Odds', value = NA))
    ),
  actionButton('calculate', 'Calculate'),
  tableOutput('table'),
  plotlyOutput(outputId = 'plot') 
)

server <- function(input, output) {
  observeEvent(input$calculate, {
    df <- data.frame(
      names = c(input$bet1name, input$bet2name, input$bet3name, input$bet4name, input$bet5name),
      sharp.odds = c(input$bet1sharp, input$bet2sharp, input$bet3sharp, input$bet4sharp, input$bet5sharp),
      american.odds = c(input$bet1soft, input$bet2soft, input$bet3soft, input$bet4soft, input$bet5soft)
    )
    df <- na.omit(df)
    result <- optimize.bets(df)
    risk.adjusted <- result$weights[which.max(result$weights$returns / result$weights$risk), ]
    output$table <- renderTable(risk.adjusted)
    output$plot <- renderPlotly(result$plot)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
