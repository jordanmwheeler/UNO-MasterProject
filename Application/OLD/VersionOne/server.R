### Loading packages
packages = c("shiny", "markdown", "shinythemes", "dplyr", "ggplot2", 
             "plotly", "tidyr", "kableExtra", "reshape2", "quantmod",
             "lubridate", "shinyWidgets", "Hmisc", "fGarch", "parallel")
invisible(lapply(packages, library, character.only = TRUE))

##################################################
##################################################
### Auto GARCH Code by Ivan Popivanov (https://www.r-bloggers.com/automatic-armagarch-selection-in-parallel/)
garchAutoTryFit = function(
  ll,
  data,
  trace=FALSE,
  forecast.length=1,
  with.forecast=TRUE,
  ic="AIC",
  garch.model="garch" )
{
  formula = as.formula( paste( sep="",
                               "~ arma(", ll$order[1], ",", ll$order[2], ")+",
                               garch.model,
                               "(", ll$order[3], ",", ll$order[4], ")" ) )
  fit = tryCatch( garchFit( formula=formula,
                            data=data,
                            trace=FALSE,
                            cond.dist=ll$dist ),
                  error=function( err ) TRUE,
                  warning=function( warn ) FALSE )
  pp = NULL
  if( !is.logical( fit ) ) {
    if( with.forecast ) {
      pp = tryCatch( predict( fit,
                              n.ahead=forecast.length,
                              doplot=FALSE ),
                     error=function( err ) FALSE,
                     warning=function( warn ) FALSE )
      if( is.logical( pp ) ) {
        fit = NULL
      }
    }
  } else {
    fit = NULL
  }
  if( trace ) {
    if( is.null( fit ) ) {
      cat( paste( sep="",
                  "   Analyzing (", ll$order[1], ",", ll$order[2],
                  ",", ll$order[3], ",", ll$order[4], ") with ",
                  ll$dist, " distribution done.",
                  "Bad model.\n" ) )
    } else {
      if( with.forecast ) {
        cat( paste( sep="",
                    "   Analyzing (", ll$order[1], ",", ll$order[2], ",",
                    ll$order[3], ",", ll$order[4], ") with ",
                    ll$dist, " distribution done.",
                    "Good model. ", ic, " = ", round(fit@fit$ics[[ic]],6),
                    ", forecast: ",
                    paste( collapse=",", round(pp[,1],4) ), "\n" ) )
      } else {
        cat( paste( sep="",
                    "   Analyzing (", ll[1], ",", ll[2], ",", ll[3], ",", ll[4], ") with ",
                    ll$dist, " distribution done.",
                    "Good model. ", ic, " = ", round(fit@fit$ics[[ic]],6), "\n" ) )
      }
    }
  }
  return( fit )
}
garchAuto = function(
  xx,
  min.order=c(0,0,1,1),
  max.order=c(5,5,1,1),
  trace=FALSE,
  cond.dists="sged",
  with.forecast=TRUE,
  forecast.length=1,
  arma.sum=c(0,1e9),
  cores=1,
  ic="AIC",
  garch.model="garch" )
{
  require( fGarch )
  require( parallel )
  len = NROW( xx )
  models = list( )
  for( dist in cond.dists )
    for( p in min.order[1]:max.order[1] )
      for( q in min.order[2]:max.order[2] )
        for( r in min.order[3]:max.order[3] )
          for( s in min.order[4]:max.order[4] )
          {
            pq.sum = p + q
            if( pq.sum <= arma.sum[2] && pq.sum >= arma.sum[1] )
            {
              models[[length( models ) + 1]] = list( order=c( p, q, r, s ), dist=dist )
            }
          }
  res = mclapply( models,
                  garchAutoTryFit,
                  data=xx,
                  trace=trace,
                  ic=ic,
                  garch.model=garch.model,
                  forecast.length=forecast.length,
                  with.forecast=TRUE,
                  mc.cores=cores )
  best.fit = NULL
  best.ic = 1e9
  for( rr in res )
  {
    if( !is.null( rr ) )
    {
      current.ic = rr@fit$ics[[ic]]
      if( current.ic < best.ic )
      {
        best.ic = current.ic
        best.fit = rr
      }
    }
  }
  if( best.ic < 1e9 ){
    return( best.fit )
  }
  return( NULL )
}
##################################################
##################################################

function(input, output, session) {
  output$overviewPlot <- renderPlotly({
    ### Collecting Stock Data from Yahoo! API
    stockSym = toupper(input$stockSymbol)
    stockSymDF = getSymbols(stockSym, src='yahoo', env = NULL)
    colnames(stockSymDF) = c("Open", "High", "Low", "Close", "Volume", "Adjusted")
    stockSymDF = as.data.frame(stockSymDF)
    stockSymDF$Date = as.Date(rownames(stockSymDF))
    stockSymDF = stockSymDF[stockSymDF$Date >= input$overviewDateRan[1] & stockSymDF$Date <= input$overviewDateRan[2],]
    stockSymDF$PercentChange = stockSymDF$Close / lag(stockSymDF$Close, 1) - 1
    stockSymDF$PercentChange = scales::percent(stockSymDF$PercentChange)
    stockSymDF = stockSymDF[,c(7,1,2,3,4,5,6,8)]
    rownames(stockSymDF) = NULL
    
    ### Creating specific graphics
    if (input$overviewGraphic == "cs"){
      stockSymDF %>%
        plot_ly(x = ~Date, type="candlestick",
                open = ~Open, close = ~Close,
                high = ~High, low = ~Low,
                increasing = list(line = list(color = "#377EB8")),
                decreasing = list(line = list(color = "#FF7F00"))) %>%
          layout(title = "Basic Candlestick Chart",
                 font = list(color = "#ffffff"),
                 yaxis = list(gridcolor = "#b3b3b3",
                              tickformat='$'),
                 xaxis = list(gridcolor = "#b3b3b3")) %>%
          layout(plot_bgcolor='#222222',
                 paper_bgcolor='#222222')
    }
    else if (input$overviewGraphic == "ts"){
      stockSymDF %>%
        plot_ly(x = ~Date, y = ~Close,mode = 'lines',
                line = list(color = '#377EB8'),
                hoverinfo = "text",
                text = paste0("Date: ", stockSymDF$Date, 
                              "\nOpen: ", stockSymDF$Open, 
                              "\nClose: ", stockSymDF$Close, 
                              "\nVolume: ", stockSymDF$Volume)) %>%
        layout(title = "Time Series Chart",
               font = list(color = "#ffffff"),
               yaxis = list(gridcolor = "#b3b3b3",
                            tickformat='$'),
               xaxis = list(gridcolor = "#b3b3b3")) %>%
        layout(plot_bgcolor='#222222',
               paper_bgcolor='#222222')
    }
    
    })
  
  output$overviewTable <- DT::renderDataTable({
    stockSym = toupper(input$stockSymbol)
    stockSymDF = getSymbols(stockSym, src='yahoo', env = NULL)
    colnames(stockSymDF) = c("Open", "High", "Low", "Close", "Volume", "Adjusted")
    stockSymDF = as.data.frame(stockSymDF)
    stockSymDF$Date = as.Date(rownames(stockSymDF))
    stockSymDF = stockSymDF[stockSymDF$Date >= input$overviewDateRan[1] & stockSymDF$Date <= input$overviewDateRan[2],]
    stockSymDF$PercentChange = stockSymDF$Close / lag(stockSymDF$Close, 1) - 1
    stockSymDF$PercentChange = scales::percent(stockSymDF$PercentChange)
    stockSymDF = stockSymDF[,c(7,1,2,3,4,5,6,8)]
    rownames(stockSymDF) = NULL
    DT::datatable(stockSymDF, rownames = FALSE)  %>% 
      DT::formatStyle(colnames(stockSymDF),color = '#377EB8', backgroundColor = '#222222', fontWeight = 'bold')
  })
  
  forecastPlotly = reactivePlot()
  observeEvent(input$forecastButton, {
    stockSym.For = toupper(input$stockSymbolForecast)
    stockSymDF.For = getSymbols(stockSym.For, src='yahoo', env = NULL)
    colnames(stockSymDF.For) = c("Open", "High", "Low", "Close", "Volume", "Adjusted")
    stockSymDF.For$Date = as.Date(rownames(stockSymDF.For))
    stockSymDF.For = stockSymDF[stockSymDF.For$Date >= input$forecastDateRan[1] & stockSymDF.For$Date <= input$forecastDateRan[2],]
    rownames(stockSymDF.For) = NULL
    stockSymDF.For = stockSymDF.For[,c(7,1:6)]
    
    stockSymDF.For$Return = log(stockSymDF.For$Close) - log(Lag(stockSymDF.For$Close))
    stockSymDF.For = stockSymDF.For[complete.cases(stockSymDF.For),]
    cpus = detectCores()
    
    numPredictions = input$numForecast
    numDay = numPredictions + (ceiling(numPredictions/7)*3)
    possibleDays = stockSymDF[nrow(stockSymDF),]$Date + 1:numDay
    possible.weekdays = possibleDays[!weekdays(possibleDays) %in% c("Saturday", "Sunday")]
    possible.weekdays = possible.weekdays[1:numPredictions]
    
    GARCHfit = garchAuto(stockSymDF.For$Return, cores = cpus, trace=TRUE)
    predictions = predict(GARCHfit, numPredictions)
    
    predictionDF = data.frame(Date = possible.weekdays, 
                              Close = NA, 
                              Return = as.numeric(unlist(predictions[1])),
                              Lower.Bound = NA,
                              Upper.Bound = NA,
                              StandardDeviation = as.numeric(unlist(predictions[3])))
    
    for(i in 1:numPredictions){
      if(i == 1){
        predictionDF[i,2] = exp(as.numeric(log(stockSymDF[nrow(stockSymDF),]$Close)) + as.numeric(predictionDF[i,3]))
        predictionDF[i,4] = exp(as.numeric(log(stockSymDF[nrow(stockSymDF),]$Close)) - as.numeric(predictionDF[i,3]+(1.96*predictionDF[i,6])))
        predictionDF[i,5] = exp(as.numeric(log(stockSymDF[nrow(stockSymDF),]$Close)) + as.numeric(predictionDF[i,3]+(1.96*predictionDF[i,6])))
      }
      else{
        predictionDF[i,2] = exp(as.numeric(log(predictionDF[i-1,2])) + as.numeric(predictionDF[i,3]))
        predictionDF[i,4] = exp(as.numeric(log(predictionDF[i-1,2])) - as.numeric(predictionDF[i,3]+(1.96*predictionDF[i,6])))
        predictionDF[i,5] = exp(as.numeric(log(predictionDF[i-1,2])) + as.numeric(predictionDF[i,3]+(1.96*predictionDF[i,6])))
      }
    }
    
    predictionDF$Forecast = "Forecast"
    GARCHDF = data.frame(Date = stockSymDF[,c(1)], Close = stockSymDF[,c(5)], Lower.Bound = NA, Upper.Bound = NA, Forecast = "Actual")
    GARCHDF = rbind(GARCHDF, predictionDF[,c(1,2,4,5,7)])
    GARCHDF$Line.Color = ifelse(GARCHDF$Forecast=="Actual", "#0033cc", "#33cc33")
    GARCHDF$Forecast = factor(GARCHDF$Forecast)
    
    forecastPlotly(plot_ly(GARCHDF, x = ~Date, y = round(GARCHDF$Close, 2), type = 'scatter', mode = 'lines',
                           color = ~Forecast, colors = c("#377EB8", "#39ac39"),
                           hoverinfo = "text",
                           hovertext = paste0(GARCHDF$Forecast, "\n", 
                                              "Date: ", GARCHDF$Date, "\n", 
                                              "Price: ", scales::dollar(GARCHDF$Close))) %>%
                     add_trace(y = round(GARCHDF$Upper.Bound, 2), type = 'scatter', mode = 'lines',
                               fillcolor='rgba(57,172,57,0.2)', line = list(color = 'rgba(57,172,57,0.5)', dash = 'dash'),
                               showlegend = FALSE, name = 'Upper Bound Error',
                               hoverinfo = "text",
                               hovertext = paste0("Upper Bound Error\n", 
                                                  "Date: ", GARCHDF$Date, "\n", 
                                                  "Price: ", scales::dollar(GARCHDF$Upper.Bound))) %>%
                     add_trace(y = round(GARCHDF$Lower.Bound, 2), type = 'scatter', mode = 'lines',
                               fill = 'tonexty', fillcolor='rgba(57,172,57,0.2)', line = list(color = 'rgba(57,172,57,0.5)', dash = 'dash'),
                               showlegend = FALSE, name = 'Lower Bound Error',
                               hoverinfo = "text",
                               hovertext = paste0("Lower Bound Error\n", 
                                                  "Date: ", GARCHDF$Date, "\n", 
                                                  "Price: ", scales::dollar(GARCHDF$Lower.Bound))) %>%
                     layout(xaxis = list(title = "Date"),
                            yaxis = list(title = "Stock Price",
                                         tickformat = "$.2f")))
  })
  
  output$forecastPlot <- renderPlotly({
    forecastPlotly()
  })
}
###  %>% DT::formatStyle(colnames(stockSymDF),color = '#377EB8', backgroundColor = '#222222', fontWeight = 'bold')
