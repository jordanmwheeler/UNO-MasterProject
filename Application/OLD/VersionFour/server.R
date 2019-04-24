### Loading packages
packages = c("shiny", "markdown", "shinythemes", "dplyr", "ggplot2", 
             "plotly", "tidyr", "kableExtra", "reshape2", "quantmod",
             "lubridate", "shinyWidgets", "Hmisc", "fGarch", "parallel",
             "shinycssloaders")
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

predict_swift=function(model, n.ahead){
  fit=model@fit
  u = fit$series$order[1]
  v = fit$series$order[2]
  p = fit$series$order[3]
  q = fit$series$order[4]		
  max.order = max(u, v, p, q)
  h.start = fit$series$h.start
  llh.start = fit$series$llh.start
  index = fit$params$index
  params = fit$params$params
  par = fit$par
  Names = names(index)
  for (Name in Names) params[Name] = par[Name]
  Names = names(params)
  cond.dist = fit$params$cond.dist
  leverage = fit$params$leverage
  mu = params["mu"]
  if (u > 0) {
    ar = params[substr(Names, 1, 2) == "ar"]
  }
  else {
    ar = c(ar1 = 0)
  }
  if (v > 0) {
    ma = params[substr(Names, 1, 2) == "ma"]
  }
  else {
    ma = c(ma1 = 0)
  }
  omega = params["omega"]
  if (p > 0) {
    alpha = params[substr(Names, 1, 5) == "alpha"]
  }
  else {
    alpha = c(alpha1 = 0)
  }
  if (p > 0 & leverage) {
    gamma = params[substr(Names, 1, 5) == "gamma"]
  }
  else {
    gamma = c(gamma1 = 0)
  }
  if (q > 0) {
    beta = params[substr(Names, 1, 4) == "beta"]
  }
  else {
    beta = c(beta1 = 0)
  }
  delta = params["delta"]
  skew = params["skew"]
  shape = params["shape"]
  
  M = n.ahead
  N = length(model@data)
  x = c(model@data, rep(mu, M))
  h = c(model@h.t, rep(0, M))
  z = c(model@fit$series$z, rep(mu, M))
  var.model = model@fit$series$model[2]
  if (var.model == "garch") {
    
    for (i in 1:M) {
      h[N + i] = omega + sum(beta * h[N + i - (1:q)])
      for (j in 1:p) {
        if (i - j > 0) {
          s = h[N + i - j]
        }
        else {
          s = z[N + i - j]^2
        }
        h[N + i] = h[N + i] + alpha[j] * s
      }
    }
  }
  mu <- mu/(1 - sum(ar))
  ARMA <- arima(x = model@data, order = c(max(u, 1), 0, 
                                          max(v, 1)), init = c(ar, ma, mu), transform.pars = FALSE, 
                optim.control = list(maxit = 0))
  prediction = predict(ARMA, n.ahead)
  meanForecast = as.vector(prediction$pred)
  
  
  hhat <- h[-(1:N)]^(2/delta[[1]])
  u2 <- length(ar)
  meanError <- hhat[1]
  
  if (n.ahead > 1) {
    for (i in 2:n.ahead) {
      temp_error=sum(hhat[1:i]*rev(c(1,(ARMAtoMA(ar=ar,ma=ma,lag.max=(i-1))^2))))
      meanError=c(meanError, temp_error)
    }
  }
  meanError <- sqrt(meanError)
  
  standardDeviation = h^(1/delta)
  
  forecast = data.frame(meanForecast = meanForecast, 
                        meanError = meanError, standardDeviation = standardDeviation[-(1:N)])
  
  forecast
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
                            tickformat='$.2f'),
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
                            tickformat='$.2f'),
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
    stockSymDF$Date = as.Date(rownames(stockSymDF), format = "%Y-%m-%d")
    stockSymDF = stockSymDF[stockSymDF$Date >= input$overviewDateRan[1] & stockSymDF$Date <= input$overviewDateRan[2],]
    stockSymDF$PercentChange = stockSymDF$Close / Lag(stockSymDF$Close, 1) - 1
    stockSymDF$PercentChange = scales::percent(stockSymDF$PercentChange)
    stockSymDF$Open = scales::dollar(stockSymDF$Open)
    stockSymDF$High = scales::dollar(stockSymDF$High)
    stockSymDF$Low = scales::dollar(stockSymDF$Low)
    stockSymDF$Close = scales::dollar(stockSymDF$Close)
    stockSymDF$Adjusted = scales::dollar(stockSymDF$Adjusted)
    stockSymDF$Volume = scales::comma(stockSymDF$Volume)
    stockSymDF = stockSymDF[,c(7,1,2,3,4,5,6,8)]
    rownames(stockSymDF) = NULL
    DT::datatable(stockSymDF, rownames = FALSE)  %>% 
      DT::formatStyle(colnames(stockSymDF),color = '#377EB8', backgroundColor = '#222222', fontWeight = 'bold')
  })
  
  man.forecast.Modal <- function(failed = FALSE) {
    modalDialog(
      title = "Specify ARMA and GARCH Parameters",
      column(6,
             h2("ARMA:"),
             numericInput("man.AR", label = "Choose AR Parameter:", value = 1),
             numericInput("man.MA", label = "Choose MA Parameter:", value = 1)),
      column(6,
             h2("GARCH:"),
             numericInput("man.GARCH", label = "Choose GARCH Parameter:", value = 1),
             numericInput("man.ARCH", label = "Choose ARCH Parameter:", value = 1)),
      if (failed)
        div(tags$b("Invalid Model Parameters", style = "color: red;")),
      easyClose = TRUE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("man.forecast", "Forecast")
      )
    )
  }
  
  observeEvent(input$man.forecast.button, {
    showModal(man.forecast.Modal())
  })
  
  GARCHDF <- eventReactive(c(input$forecastButton,
                             input$man.forecast), 
                           {if (input$forecastButton)
                           {withProgress(message = 'Running GRACH Forecasting', 
                                         detail = "This may take a few minutes...", 
                                         value = 0,
                                         {
                                           stockSym.For = toupper(input$stockSymbolForecast)
                                           stockSymDF.For = getSymbols(stockSym.For, src='yahoo', env = NULL)
                                           date_vals = index(stockSymDF.For)
                                           stockSymDF.For = data.frame(stockSymDF.For)
                                           colnames(stockSymDF.For) = c("Open", "High", "Low", "Close", "Volume", "Adjusted")
                                           stockSymDF.For$Date = date_vals
                                           stockSymDF.For = stockSymDF.For[stockSymDF.For$Date >= input$forecastDateRan[1] & stockSymDF.For$Date <= input$forecastDateRan[2],]
                                           rownames(stockSymDF.For) = NULL
                                           stockSymDF.For = stockSymDF.For[,c(7,1:6)]
                                           
                                           stockSymDF.For$Return = log(stockSymDF.For$Close) - log(Lag(stockSymDF.For$Close))
                                           stockSymDF.For = stockSymDF.For[complete.cases(stockSymDF.For),]
                                           cpus = detectCores()
                                           
                                           numPredictions = input$numForecast
                                           numDay = numPredictions + (ceiling(numPredictions/7)*3)
                                           possibleDays = stockSymDF.For[nrow(stockSymDF.For),]$Date + 1:numDay
                                           possible.weekdays = possibleDays[!weekdays(possibleDays) %in% c("Saturday", "Sunday")]
                                           possible.weekdays = possible.weekdays[1:numPredictions]
                                           
                                           incProgress(1/3)
                                           
                                           GARCHfit = garchAuto(stockSymDF.For$Return, trace=TRUE)
                                           
                                           incProgress(1/3)
                                           
                                           predictions = predict_swift(GARCHfit, numPredictions)
                                           
                                           predictionDF = data.frame(Date = possible.weekdays, 
                                                                     Close = NA, 
                                                                     Return = as.numeric(unlist(predictions[1])),
                                                                     Lower.Bound = NA,
                                                                     Upper.Bound = NA,
                                                                     StandardDeviation = as.numeric(unlist(predictions[3])))
                                           
                                           for(i in 1:numPredictions){
                                             if(i == 1){
                                               predictionDF[i,2] = exp(as.numeric(log(stockSymDF.For[nrow(stockSymDF.For),]$Close)) + as.numeric(predictionDF[i,3]))
                                               predictionDF[i,4] = exp(as.numeric(log(stockSymDF.For[nrow(stockSymDF.For),]$Close)) - as.numeric(predictionDF[i,3]+(1.96*predictionDF[i,6])))
                                               predictionDF[i,5] = exp(as.numeric(log(stockSymDF.For[nrow(stockSymDF.For),]$Close)) + as.numeric(predictionDF[i,3]+(1.96*predictionDF[i,6])))
                                             }
                                             else{
                                               predictionDF[i,2] = exp(as.numeric(log(predictionDF[i-1,2])) + as.numeric(predictionDF[i,3]))
                                               predictionDF[i,4] = exp(as.numeric(log(predictionDF[i-1,2])) - as.numeric(predictionDF[i,3]+(1.96*predictionDF[i,6])))
                                               predictionDF[i,5] = exp(as.numeric(log(predictionDF[i-1,2])) + as.numeric(predictionDF[i,3]+(1.96*predictionDF[i,6])))
                                             }
                                           }
                                           
                                           predictionDF$Forecast = "Forecast"
                                           GARCHDF = data.frame(Date = stockSymDF.For[,c(1)], Close = stockSymDF.For[,c(5)], Lower.Bound = NA, Upper.Bound = NA, Forecast = "Actual")
                                           GARCHDF = rbind(GARCHDF, predictionDF[,c(1,2,4,5,7)])
                                           GARCHDF$Color.Cat = GARCHDF$Forecast
                                           temp.row = GARCHDF[nrow(dplyr::filter(GARCHDF, Color.Cat == "Actual")),]
                                           GARCHDF[nrow(dplyr::filter(GARCHDF, Color.Cat == "Actual")),c("Color.Cat")] = "Forecast"
                                           GARCHDF = rbind(GARCHDF, temp.row)
                                           GARCHDF$Color.Cat = factor(GARCHDF$Color.Cat)
                                           GARCHDF = GARCHDF[order(GARCHDF$Date, GARCHDF$Color.Cat),]
                                           
                                           incProgress(1/3)
                                         })
                             
                             return(GARCHDF)
                           }
                             if (input$man.forecast){
                               # Check that data object exists and is data frame.
                               if ((input$man.AR >= 0) & (input$man.MA >= 0) &
                                   (input$man.GARCH >= 0) & (input$man.ARCH >= 0)) {
                                 {removeModal()
                                   withProgress(message = 'Running GRACH Forecasting', 
                                                detail = "This may take a few minutes...", 
                                                value = 0,
                                                {
                                                  stockSym.For = toupper(input$stockSymbolForecast)
                                                  stockSymDF.For = getSymbols(stockSym.For, src='yahoo', env = NULL)
                                                  date_vals = index(stockSymDF.For)
                                                  stockSymDF.For = data.frame(stockSymDF.For)
                                                  colnames(stockSymDF.For) = c("Open", "High", "Low", "Close", "Volume", "Adjusted")
                                                  stockSymDF.For$Date = date_vals
                                                  stockSymDF.For = stockSymDF.For[stockSymDF.For$Date >= input$forecastDateRan[1] & stockSymDF.For$Date <= input$forecastDateRan[2],]
                                                  rownames(stockSymDF.For) = NULL
                                                  stockSymDF.For = stockSymDF.For[,c(7,1:6)]
                                                  
                                                  stockSymDF.For$Return = log(stockSymDF.For$Close) - log(Lag(stockSymDF.For$Close))
                                                  stockSymDF.For = stockSymDF.For[complete.cases(stockSymDF.For),]
                                                  cpus = detectCores()
                                                  
                                                  numPredictions = input$numForecast
                                                  numDay = numPredictions + (ceiling(numPredictions/7)*3)
                                                  possibleDays = stockSymDF.For[nrow(stockSymDF.For),]$Date + 1:numDay
                                                  possible.weekdays = possibleDays[!weekdays(possibleDays) %in% c("Saturday", "Sunday")]
                                                  possible.weekdays = possible.weekdays[1:numPredictions]
                                                  
                                                  incProgress(1/3)
                                                  
                                                  formula2 = as.formula(paste0(" ~arma(", input$man.AR, ", ", input$man.MA, ") + garch(", input$man.GARCH, ", ", input$man.ARCH, ")"))
                                                  GARCHfit = garchFit(formula = formula2, data = stockSymDF.For$Return, trace = FALSE, cond.dist = "sged")
                                                  
                                                  incProgress(1/3)
                                                  
                                                  predictions = predict_swift(GARCHfit, numPredictions)
                                                  
                                                  predictionDF = data.frame(Date = possible.weekdays, 
                                                                            Close = NA, 
                                                                            Return = as.numeric(unlist(predictions[1])),
                                                                            Lower.Bound = NA,
                                                                            Upper.Bound = NA,
                                                                            StandardDeviation = as.numeric(unlist(predictions[3])))
                                                  
                                                  for(i in 1:numPredictions){
                                                    if(i == 1){
                                                      predictionDF[i,2] = exp(as.numeric(log(stockSymDF.For[nrow(stockSymDF.For),]$Close)) + as.numeric(predictionDF[i,3]))
                                                      predictionDF[i,4] = exp(as.numeric(log(stockSymDF.For[nrow(stockSymDF.For),]$Close)) - as.numeric(predictionDF[i,3]+(1.96*predictionDF[i,6])))
                                                      predictionDF[i,5] = exp(as.numeric(log(stockSymDF.For[nrow(stockSymDF.For),]$Close)) + as.numeric(predictionDF[i,3]+(1.96*predictionDF[i,6])))
                                                    }
                                                    else{
                                                      predictionDF[i,2] = exp(as.numeric(log(predictionDF[i-1,2])) + as.numeric(predictionDF[i,3]))
                                                      predictionDF[i,4] = exp(as.numeric(log(predictionDF[i-1,2])) - as.numeric(predictionDF[i,3]+(1.96*predictionDF[i,6])))
                                                      predictionDF[i,5] = exp(as.numeric(log(predictionDF[i-1,2])) + as.numeric(predictionDF[i,3]+(1.96*predictionDF[i,6])))
                                                    }
                                                  }
                                                  
                                                  predictionDF$Forecast = "Forecast"
                                                  GARCHDF = data.frame(Date = stockSymDF.For[,c(1)], Close = stockSymDF.For[,c(5)], Lower.Bound = NA, Upper.Bound = NA, Forecast = "Actual")
                                                  GARCHDF = rbind(GARCHDF, predictionDF[,c(1,2,4,5,7)])
                                                  GARCHDF$Color.Cat = GARCHDF$Forecast
                                                  temp.row = GARCHDF[nrow(dplyr::filter(GARCHDF, Color.Cat == "Actual")),]
                                                  GARCHDF[nrow(dplyr::filter(GARCHDF, Color.Cat == "Actual")),c("Color.Cat")] = "Forecast"
                                                  GARCHDF = rbind(GARCHDF, temp.row)
                                                  GARCHDF$Color.Cat = factor(GARCHDF$Color.Cat)
                                                  GARCHDF = GARCHDF[order(GARCHDF$Date, GARCHDF$Color.Cat),]
                                                  
                                                  incProgress(1/3)
                                                })
                                   return(GARCHDF)}
                               } else {
                                 showModal(man.forecast.Modal(failed = TRUE))
                               }
                             }})
  
  output$forecastPlot <- renderPlotly({
    if (!exists("GARCHDF")) return()
    
    GARCHDF = GARCHDF()
    # GARCHDF$Date = as.Date(GARCHDF$Date.Char, format = "%Y-%m-%d")
    
     plot_ly(GARCHDF, x = ~Date, y = round(GARCHDF$Close, 2), type = 'scatter', mode = 'lines',
            color = ~Color.Cat, colors = c("#377EB8", "#39ac39"),
            hoverinfo = "text",
            hovertext = paste0(GARCHDF$Forecast, "\n", 
                               weekdays(GARCHDF$Date), ": ", GARCHDF$Date, "\n", 
                               "Price: ", scales::dollar(GARCHDF$Close))) %>%
      add_trace(y = round(GARCHDF$Upper.Bound, 2), type = 'scatter', mode = 'lines',
                fillcolor='rgba(57,172,57,0.2)', line = list(color = 'rgba(57,172,57,0.5)', dash = 'dash'),
                showlegend = FALSE, name = 'Upper Bound Error',
                hoverinfo = "text",
                hovertext = paste0("Upper Bound Error\n", 
                                   weekdays(GARCHDF$Date), ": ", GARCHDF$Date, "\n", 
                                   "Price: ", scales::dollar(GARCHDF$Upper.Bound))) %>%
      add_trace(y = round(GARCHDF$Lower.Bound, 2), type = 'scatter', mode = 'lines',
                fill = 'tonexty', fillcolor='rgba(57,172,57,0.2)', line = list(color = 'rgba(57,172,57,0.5)', dash = 'dash'),
                showlegend = FALSE, name = 'Lower Bound Error',
                hoverinfo = "text",
                hovertext = paste0("Lower Bound Error\n", 
                                   weekdays(GARCHDF$Date), ": ", GARCHDF$Date, "\n", 
                                   "Price: ", scales::dollar(GARCHDF$Lower.Bound))) %>%
      layout(title = "Time Series Chart",
             font = list(color = "#ffffff"),
             yaxis = list(title = "Stock Price",
                          gridcolor = "#b3b3b3",
                          tickformat='$.2f'),
             xaxis = list(title = "Date",
                          gridcolor = "#b3b3b3"),
             plot_bgcolor='#222222',
             paper_bgcolor='#222222')
  })
  
  output$forecastTable <- DT::renderDataTable({
    if (input$forecastButton == 0 & input$man.forecast == 0) return()
    
    GARCHDF = GARCHDF()
    GARCHDF.tab <<- GARCHDF
    
    if (input$forecastDisplay.button == "Actuals and Forecast"){
      GARCHDF.tab$WeekDay = weekdays(GARCHDF.tab$Date)
      GARCHDF.tab$Close = scales::dollar(GARCHDF.tab$Close)
      GARCHDF.tab$Lower.Bound = scales::dollar(GARCHDF.tab$Lower.Bound)
      GARCHDF.tab$Upper.Bound = scales::dollar(GARCHDF.tab$Upper.Bound)
      DT::datatable(GARCHDF.tab[,c(7,1,2,3,4,5)], rownames = FALSE, colnames = c("Weekday", "Date", "Close", "Lower Bound", "Upper Bound", "Actual/Forecast"))  %>% 
        DT::formatStyle(columns = c("WeekDay", "Date", "Close", "Lower.Bound", "Upper.Bound", "Forecast"), color = '#377EB8', backgroundColor = '#222222', fontWeight = 'bold')
    }
    
    else if (input$forecastDisplay.button == "Only Forecast"){
      GARCHDF.tab$WeekDay = weekdays(GARCHDF.tab$Date)
      GARCHDF.tab$Close = scales::dollar(GARCHDF.tab$Close)
      GARCHDF.tab$Lower.Bound = scales::dollar(GARCHDF.tab$Lower.Bound)
      GARCHDF.tab$Upper.Bound = scales::dollar(GARCHDF.tab$Upper.Bound)
      GARCHDF.tab <<- GARCHDF.tab[GARCHDF.tab$Forecast == "Forecast",]
      DT::datatable(GARCHDF.tab[GARCHDF.tab$Forecast == "Forecast",c(7,1,2,3,4,5)], rownames = FALSE, colnames = c("Weekday", "Date", "Close", "Lower Bound", "Upper Bound", "Actual/Forecast"))  %>% 
        DT::formatStyle(columns = c("WeekDay", "Date", "Close", "Lower.Bound", "Upper.Bound", "Forecast"), color = '#377EB8', backgroundColor = '#222222', fontWeight = 'bold')
    }
  })
  
  
  output$forecast.download <- downloadHandler(
    filename = function() {
      paste0(input$stockSymbolForecast, "_GARCHForecast", '.csv')
    },
    content = function(con) {
      write.csv(GARCHDF.tab, con, row.names = FALSE)
    }
  )
  
}
###  %>% DT::formatStyle(colnames(stockSymDF),color = '#377EB8', backgroundColor = '#222222', fontWeight = 'bold')
