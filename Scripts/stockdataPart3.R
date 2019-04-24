### Loading packages
packages = c("dplyr", "ggplot2", "plotly", "tidyr", "BBmisc", "kableExtra", "reshape2", "quantmod",
             "lubridate", "Hmisc", "fGarch", "parallel")
invisible(lapply(packages, library, character.only = TRUE))

##################################################
### API Key: ZGYWNP1KJ5A6FALM
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

x = function() {
stockSym.For = toupper("aapl")
stockSymDF.For = getSymbols(stockSym.For, src='yahoo', env = NULL)
date_vals = index(stockSymDF.For)
stockSymDF.For = data.frame(stockSymDF.For)
colnames(stockSymDF.For) = c("Open", "High", "Low", "Close", "Volume", "Adjusted")
stockSymDF.For$Date = date_vals
stockSymDF.For = stockSymDF.For[stockSymDF.For$Date >= "2019-02-22" & stockSymDF.For$Date <= "2019-03-22",]
rownames(stockSymDF.For) = NULL
stockSymDF.For = stockSymDF.For[,c(7,1:6)]

stockSymDF.For$Return = log(stockSymDF.For$Close) - log(Lag(stockSymDF.For$Close))
stockSymDF.For = stockSymDF.For[complete.cases(stockSymDF.For),]
cpus = detectCores()

numPredictions = 10
numDay = numPredictions + (ceiling(numPredictions/7)*3)
possibleDays = stockSymDF.For[nrow(stockSymDF.For),]$Date + 1:numDay
possible.weekdays = possibleDays[!weekdays(possibleDays) %in% c("Saturday", "Sunday")]
possible.weekdays = possible.weekdays[1:numPredictions]

GARCHfit = garchAuto(stockSymDF.For$Return, trace=TRUE, min.order = c(0,2,1,1), max.order = c(0,2,1,1))
GARCHfit
formula2 = as.formula(paste0(" ~arma(", 5, ", ", 5, ") + garch(", 1, ", ", 1, ")"))
x = garchFit(formula = formula2, data = stockSymDF.For$Return, trace = FALSE, cond.dist = "sged")

predictions = predict_swift(garchFit(formula = formula2, data = stockSymDF.For$Return, trace = FALSE, cond.dist = "sged"), 5)

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
GARCHDF[nrow(dplyr::filter(GARCHDF, Color.Cat == "Actual")),c("Lower.Bound")] = GARCHDF[nrow(dplyr::filter(GARCHDF, Color.Cat == "Actual")),c("Close")]
GARCHDF[nrow(dplyr::filter(GARCHDF, Color.Cat == "Actual")),c("Upper.Bound")] = GARCHDF[nrow(dplyr::filter(GARCHDF, Color.Cat == "Actual")),c("Close")]
GARCHDF[nrow(dplyr::filter(GARCHDF, Color.Cat == "Actual")),c("Color.Cat")] = "Forecast"
GARCHDF = rbind(GARCHDF, temp.row)
GARCHDF$Color.Cat = factor(GARCHDF$Color.Cat)
GARCHDF = GARCHDF[order(GARCHDF$Date, GARCHDF$Color.Cat),]


return(GARCHDF)
}

GARCHDF.test = x()

class(GARCHDF.test$Date)

GARCHDF
insertR

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
  layout(xaxis = list(title = "Date"),
         yaxis = list(title = "Stock Price",
                      tickformat = "$.2f"))




stockSym = getSymbols("AAPL",src='yahoo', env = NULL)
date_val = index(stockSym)
stockSymDF = data.frame(stockSym)
colnames(stockSymDF) = c("Open", "High", "Low", "Close", "Volume", "Adjusted")
index(stockSymDF)
stockSymDF$Date = as.Date(rownames(stockSymDF), format = "%Y-%m-%d")
rownames(stockSymDF) = NULL
stockSymDF = stockSymDF[,c(7,1:6)]

stockSymDF$Return = log(stockSymDF$Close) - log(Lag(stockSymDF$Close))
stockSymDF = stockSymDF[complete.cases(stockSymDF),]
cpus = detectCores()

GARCHfit = garchAuto(stockSymDF$Return, trace=TRUE)
numPredictions = 30

numDay = numPredictions + (ceiling(numPredictions/7)*3)
possibleDays = stockSymDF[nrow(stockSymDF),]$Date + 1:numDay
possible.weekdays = possibleDays[!weekdays(possibleDays) %in% c("Saturday", "Sunday")]
possible.weekdays = possible.weekdays[1:numPredictions]
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
GARCHDF$hover

# plot_ly(GARCHDF, x = ~Date, y = ~Upper.Bound, type = 'scatter', mode = 'lines',
#              line = list(color = 'transparent'),
#              showlegend = FALSE, name = 'Upper Bound Error') %>%
#   add_trace(y = ~Lower.Bound, type = 'scatter', mode = 'lines',
#             fill = 'tonexty', fillcolor='rgba(0,100,80,0.2)', line = list(color = 'transparent'),
#             showlegend = FALSE, name = 'Lower Bound Error') %>%
#   add_trace(data = stockSymDF, x = ~Date, y = ~Close, type = 'scatter', mode = 'lines',
#             line = list(color="#0033cc"),
#             name = 'Actuals', showlegend = TRUE) %>%
#   add_trace(data = predictionDF, x = ~Date , y = ~Close, type = 'scatter', mode = 'lines',
#             line = list(color="#33cc33"),
#             name = 'Forecasts', showlegend = TRUE) %>%
#   layout(xaxis = list(title = "Date"),
#          yaxis = list(title = "Stock Price",
#                       tickformat = "$.2f"))

GARCHDF = GARCHDF[GARCHDF$Date >= "2019-01-01",]

plot_ly(GARCHDF, x = ~Date, y = round(GARCHDF$Close, 2), type = 'scatter', mode = 'lines',
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
                      tickformat = "$.2f"))
stockSymDF


stockSym = toupper('hemp')
stockSymDF = getSymbols(stockSym, src='yahoo', env = NULL)
colnames(stockSymDF) = c("Open", "High", "Low", "Close", "Volume", "Adjusted")
stockSymDF = as.data.frame(stockSymDF)
stockSymDF$Date = as.Date(rownames(stockSymDF), format = "%Y-%m-%d")
stockSymDF = stockSymDF[stockSymDF$Date >= "2019-01-01" & stockSymDF$Date <= "2019-03-14",]
stockSymDF$PercentChange = stockSymDF$Close / Lag(stockSymDF$Close, 1) - 1
stockSymDF$PercentChange = scales::percent(stockSymDF$PercentChange)
stockSymDF$Open = scales::dollar(stockSymDF$Open)
stockSymDF$High = scales::dollar(stockSymDF$High)
stockSymDF$Low = scales::dollar(stockSymDF$Low)
stockSymDF$Close = scales::dollar(stockSymDF$Close)
stockSymDF$Adjusted = scales::dollar(stockSymDF$Adjusted)
stockSymDF$Adjusted = scales::comma(stockSymDF$Volume)
stockSymDF = stockSymDF[,c(7,1,2,3,4,5,6,8)] 
rownames(stockSymDF) = NULL
DT::datatable(stockSymDF, rownames = FALSE)  %>% 
  DT::formatStyle(colnames(stockSymDF),color = '#377EB8', backgroundColor = '#222222', fontWeight = 'bold')