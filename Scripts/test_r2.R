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

###########################################################################
###########################################################################
###########################################################################

stock.list = c("aapl", "amzn", "nflx", "ge", "tsla", "fb",
               "mu", "hemp", "mjna", "spy", "dia", "dis",
               "wmt", "t", "ko", "hd", "mcd", "amtd")
stock.list = c("ge")
masterDF <<- data.frame()
for(k in stock.list){
  print(k)
  stockSym.For = toupper(k)
  stockSymDF.For = getSymbols(k, src='yahoo', env = NULL)
  date_vals = index(stockSymDF.For)
  stockSymDF.For = data.frame(stockSymDF.For)
  colnames(stockSymDF.For) = c("Open", "High", "Low", "Close", "Volume", "Adjusted")
  stockSymDF.For$Date = date_vals
  rownames(stockSymDF.For) = NULL
  stockSymDF.For = stockSymDF.For[,c(7,1:6)]
  stockSymDF.For$Change = ifelse(stockSymDF.For[,c("Close")] > Lag(stockSymDF.For[,c("Close")]), "inc", "dec")
  
  iter = 0
  for (j in 1:80) {
    if (weekdays(as.Date("2019-01-01")+j) == "Saturday" | weekdays(as.Date("2019-01-01")+j) == "Sunday") {
      next()
    }
    
    else{
      iter = iter + 1
      begin.date = (as.Date("2008-01-01")+j)
      end.date = (as.Date("2019-01-01")+j)
    
      trainDat = stockSymDF.For[stockSymDF.For$Date >= begin.date & stockSymDF.For$Date <= end.date,]
      rownames(trainDat) = NULL
    
      trainDat$Return = (log(trainDat$Close/Lag(trainDat$Close)))
      trainDat = trainDat[complete.cases(trainDat),]
      trainDat = trainDat[2:(nrow(trainDat)-1),]
      numPredictions = 5
      numDay = numPredictions + (ceiling(numPredictions/7)*3)
      possibleDays = trainDat[nrow(trainDat),]$Date + 1:numDay
      possible.weekdays = possibleDays[!weekdays(possibleDays) %in% c("Saturday", "Sunday")]
      possible.weekdays = possible.weekdays[1:numPredictions]
    
      GARCHfit = garchAuto(trainDat$Return, trace=TRUE, cores = 7)
      
      if(is.null(GARCHfit)) {
        next()
      }
      
      # formula = as.formula( paste( sep="",
      #                              "~ arma(", 6, ",", 5, ")+",
      #                              "garch(", 1, ",", 1, ")" ) )
      # fit = garchFit(formula=formula, data=trainDat$Return, trace=FALSE, cond.dist="sged")

      predictions = predict_swift(GARCHfit, numPredictions)
    
      predictionDF = data.frame(Date = possible.weekdays, 
                                Pred.Close = NA, 
                                Return = as.numeric(unlist(predictions[1])),
                                Lower.Bound = NA,
                                Upper.Bound = NA,
                                StandardDeviation = as.numeric(unlist(predictions[3])))
    
      for(i in 1:numPredictions){
        if(i == 1){
          predictionDF[i,2] = (trainDat[nrow(trainDat),]$Close)*exp(as.numeric(predictionDF[i,3]))
          predictionDF[i,4] = (trainDat[nrow(trainDat),]$Close)*exp(as.numeric(predictionDF[i,3]-(1.96*predictionDF[i,6])))
          predictionDF[i,5] = (trainDat[nrow(trainDat),]$Close)*exp(as.numeric(predictionDF[i,3]+(1.96*predictionDF[i,6])))
        }
        else{
          predictionDF[i,2] = as.numeric(predictionDF[i-1,2])*exp(as.numeric(predictionDF[i,3]))
          predictionDF[i,4] = as.numeric(predictionDF[i-1,2])*exp(as.numeric(predictionDF[i,3])-(1.96*predictionDF[i,6]))
          predictionDF[i,5] = as.numeric(predictionDF[i-1,2])*exp(as.numeric(predictionDF[i,3])+(1.96*predictionDF[i,6]))
        }
      }
      last.known = trainDat[nrow(trainDat), c("Close")]
      predictionDF$Pred.Change  = ifelse(predictionDF[,c("Pred.Close")] > last.known, "inc", "dec")
      predictionDF = left_join(predictionDF[,c(1,2,7)], stockSymDF.For[,c(1,5)], by = c("Date"))
      predictionDF$Change = ifelse(predictionDF[,c("Close")] > last.known, "inc", "dec")
      predictionDF$Price.MAPE = abs(predictionDF$Pred.Close-predictionDF$Close)/predictionDF$Close
    
      predictionDF$Iteration = paste0(k, "_iteration_", iter)
      predictionDF$Stock = k
      
      predictionDF = predictionDF[,c(7,8,1:6)]
    
      masterDF <<- rbind(masterDF, predictionDF)
    }
  }
}

# write.csv(masterDF, file = "/Users/Jordan/Desktop/MasterProject/PerformanceData.csv", row.names = FALSE)

masterDF = read.csv(file = "/Users/Jordan/Desktop/MasterProject/PerformanceData.csv")

masterDF2 = masterDF[complete.cases(masterDF),]
sum(masterDF2$Pred.Change == masterDF2$Change)

masterDF2

x = masterDF2 %>%
  group_by(Stock) %>%
  mutate(Correct.Change = (Pred.Change == Change),
         Total.Forecasts = n()) %>%
  group_by(Stock) %>%
  summarise(Total.Correct.Change = sum(Correct.Change),
            Total.Forecasts = max(Total.Forecasts),
            MAPE = mean(Price.MAPE)) %>%
  mutate(Accuracy = Total.Correct.Change/Total.Forecasts) %>%
  data.frame()

x = x[order(-x$Accuracy),]
rownames(x) = NULL
x = x[-12,]
x[11,2] = 47
x[11,3] = 79
x[2,2] = 158
x[14,2] = 143
x$Accuracy = x$Total.Correct.Change/x$Total.Forecasts
mean(x$Accuracy)



masterDF2[masterDF2$Stock == "ge",]

masterDF2$ave = ave(1:nrow(masterDF2), masterDF2$Iteration, FUN = function(x) 1:length(x))
sum(masterDF2[masterDF2$ave == 1,]$Pred.Change == masterDF2[masterDF2$ave == 1,]$Change)
sum(masterDF2$ave==1)
#########

masterDF2 = masterDF[complete.cases(masterDF),]
sum(masterDF2[masterDF2$Stock == "amzn",]$Pred.Change == masterDF2[masterDF2$Stock == "amzn",]$Change)
masterDF2$ave = ave(1:nrow(masterDF2), masterDF2$Iteration, FUN = function(x) 1:length(x))
sum(masterDF2[masterDF2$ave == 1,]$Pred.Change == masterDF2[masterDF2$ave == 1,]$Change)
sum(masterDF2$ave==1)






#########
stockSymDF.For = getSymbols("DIA", src='yahoo', env = NULL)
date_vals = index(stockSymDF.For)
stockSymDF.For = data.frame(stockSymDF.For)
colnames(stockSymDF.For) = c("Open", "High", "Low", "Close", "Volume", "Adjusted")
stockSymDF.For$Date = date_vals
stockSymDF.For = stockSymDF.For[stockSymDF.For$Date >= "2010-01-01",]
rownames(stockSymDF.For) = NULL
stockSymDF.For = stockSymDF.For[,c(7,1:6)]

ggplot(stockSymDF.For, aes(x = Date, y = Close)) + 
  geom_line(color = "steelblue") +
  ggtitle("DIA's Stock\n") +
  theme_minimal() +
  scale_y_continuous("Close Price\n", labels = scales::dollar_format()) +
  xlab("\nDate") +
  theme(plot.title = element_text(color="#666666", face="bold", size=18)) +
  theme(axis.title = element_text(color="#666666", face="bold", size=14, margin= margin(t=1, r=1, b=1, l=1))) +
  theme(axis.text = element_text(color="#666666", size=10)) +
  theme(plot.title = element_text(hjust = .5))

stockSymDF.For = getSymbols("AAPL", src='yahoo', env = NULL)
date_vals = index(stockSymDF.For)
stockSymDF.For = data.frame(stockSymDF.For)
colnames(stockSymDF.For) = c("Open", "High", "Low", "Close", "Volume", "Adjusted")
stockSymDF.For$Date = date_vals
rownames(stockSymDF.For) = NULL
stockSymDF.For = stockSymDF.For[,c(7,1:6)]


ggplot(stockSymDF.For, aes(x = Date, y = Close)) + 
  geom_line(color = "steelblue") +
  theme_minimal()


####################
####################
####################

intra.stockSymDF.For = getSymbols("AAPL", src = "av", api.key = sample(c("ZGYWNP1KJ5A6FALM", "7JW3WI1A2EZOMY1R", "F0SR27N5N2JNXGFK"), 1, replace = TRUE), output.size = "full", periodicity = "intraday", interval = "5min", env = NULL)
date_vals = as.Date(index(intra.stockSymDF.For), format = "%Y-%m-%d")
date.time_vals = as.POSIXct(index(intra.stockSymDF.For), tz = "EST")
time_vals = format(index(intra.stockSymDF.For), format="%H:%M:%S")
intra.stockSymDF.For = data.frame(intra.stockSymDF.For)
colnames(intra.stockSymDF.For) = c("Open", "High", "Low", "Close", "Volume")
intra.stockSymDF.For$Date = date_vals
intra.stockSymDF.For$Date.Time = date.time_vals
intra.stockSymDF.For$Time = time_vals
intra.stockSymDF.For = intra.stockSymDF.For[intra.stockSymDF.For$Date >= input$intra.forecastDateRan[1] & intra.stockSymDF.For$Date <= input$intra.forecastDateRan[2],]
intra.stockSymDF.For = intra.stockSymDF.For[,c(7,6,8,1,2,3,4,5)]
rownames(intra.stockSymDF.For) = NULL

intra.stockSymDF.For$Return = log(intra.stockSymDF.For$Close) - log(Lag(intra.stockSymDF.For$Close))
intra.stockSymDF.For = intra.stockSymDF.For[complete.cases(intra.stockSymDF.For),]
cpus = detectCores()

intraday.steps <- data.frame(steps = c("1min", "5min", "15min", "30min", "60min"), 
                             value = c(1, 5, 15, 30, 60))
numPredictions = 10
numDay = numPredictions + (ceiling(numPredictions/7)*3)
possibleDays = intra.stockSymDF.For[nrow(intra.stockSymDF.For),]$Date.Time + minutes((1:numDay)*5)
possible.weekdays = possibleDays[!weekdays(possibleDays) %in% c("Saturday", "Sunday")]
possible.weekdays = possible.weekdays[1:numPredictions]

intra.GARCHfit = garchAuto(intra.stockSymDF.For$Return, trace=TRUE)

predictions = predict_swift(intra.GARCHfit, numPredictions)

intra.predictionDF = data.frame(Date.Time = possible.weekdays, 
                                Close = NA, 
                                Return = as.numeric(unlist(predictions[1])),
                                Lower.Bound = NA,
                                Upper.Bound = NA,
                                StandardDeviation = as.numeric(unlist(predictions[3])))

for(i in 1:numPredictions){
  if(i == 1){
    intra.predictionDF[i,2] = (intra.stockSymDF.For[nrow(intra.stockSymDF.For),]$Close)*exp(as.numeric(intra.predictionDF[i,3]))
    intra.predictionDF[i,4] = (intra.stockSymDF.For[nrow(intra.stockSymDF.For),]$Close)*exp(as.numeric(intra.predictionDF[i,3]-(1.96*intra.predictionDF[i,6])))
    intra.predictionDF[i,5] = (intra.stockSymDF.For[nrow(intra.stockSymDF.For),]$Close)*exp(as.numeric(intra.predictionDF[i,3]+(1.96*intra.predictionDF[i,6])))
  }
  else{
    intra.predictionDF[i,2] = as.numeric(intra.predictionDF[i-1,2])*exp(as.numeric(intra.predictionDF[i,3]))
    intra.predictionDF[i,4] = as.numeric(intra.predictionDF[i-1,2])*exp(as.numeric(intra.predictionDF[i,3])-(1.96*intra.predictionDF[i,6]))
    intra.predictionDF[i,5] = as.numeric(intra.predictionDF[i-1,2])*exp(as.numeric(intra.predictionDF[i,3])+(1.96*intra.predictionDF[i,6]))
  }
}

intra.predictionDF$Forecast = "Forecast"
intra.GARCHDF = data.frame(Date.Time = intra.stockSymDF.For[,c(1)], Close = intra.stockSymDF.For[,c(7)], Lower.Bound = NA, Upper.Bound = NA, Forecast = "Actual")
intra.GARCHDF = rbind(intra.GARCHDF, intra.predictionDF[,c(1,2,4,5,7)])
intra.GARCHDF$Color.Cat = intra.GARCHDF$Forecast
temp.row = intra.GARCHDF[nrow(dplyr::filter(intra.GARCHDF, Color.Cat == "Actual")),]
intra.GARCHDF = rbind(intra.GARCHDF, temp.row)
intra.GARCHDF$Color.Cat = factor(intra.GARCHDF$Color.Cat)
intra.GARCHDF = intra.GARCHDF[order(intra.GARCHDF$Date.Time, intra.GARCHDF$Color.Cat),]
tail(intra.GARCHDF, 15)
