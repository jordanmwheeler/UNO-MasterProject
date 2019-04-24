### Loading packages
packages = c("dplyr", "ggplot2", "plotly", "tidyr", "BBmisc", "kableExtra", "reshape2", "quantmod")
lapply(packages, library, character.only = TRUE)

##################################################
### API Key: ZGYWNP1KJ5A6FALM

getSymbols.av("AMZN", api.key = "ZGYWNP1KJ5A6FALM", interval = "15min")
getSymbols("AMZN", src = "av", api.key = "ZGYWNP1KJ5A6FALM", output.size = "full", periodicity = "intraday", interval = "15min")
AMZN
AMZN.Dat = as.data.frame(AMZN)
tail(AMZN.Dat)
plot(x = c(1:nrow(AMZN.Dat)), y = AMZN.Dat$AMZN.Open)

getSymbols("AAPL",src='yahoo')

# basic example of ohlc charts
stockDF <- data.frame(Date=index(AAPL),coredata(AAPL))
stockDF <- tail(stockDF, 30)

p = plot_ly(data = stockDF, x = ~Date, type="candlestick",
        open = ~AAPL.Open, close = ~AAPL.Close,
        high = ~AAPL.High, low = ~AAPL.Low,
        increasing = list(line = list(color = "steelblue")),
        decreasing = list(line = list(color = "#ff8000"))) %>%
  layout(title = "Basic Candlestick Chart",
         font = list(color = "white"),
         yaxis = list(gridcolor = "white")) %>%
  layout(plot_bgcolor='#666666',
         paper_bgcolor='#666666')
p


pltyobj = plotly_build(p)

pltyobj$x$data[[1]]$hoverinfo <- "text"

pltyobj

getSymbols("AAPL",src='yahoo')
data.frame(AAPL)

df <- data.frame(Date=index(AAPL),coredata(AAPL))
df = data.frame(AAPL)
df <- tail(df, 30)

p <- df %>%
  plot_ly(x = ~Date, type="candlestick",
          open = ~AAPL.Open, close = ~AAPL.Close,
          high = ~AAPL.High, low = ~AAPL.Low) %>%
  layout(title = "Basic Candlestick Chart",
         xaxis = list(rangeslider = list(visible = F)))
p

library(RColorBrewer)
brewer.pal(n = 8, name = "Set1")

library(plotly)

library(quantmod)



getSymbols("RAVE",src='yahoo')
df <- data.frame(Date=index(RAVE),coredata(RAVE))
df <- tail(df, 30)

class(stockDF$Date)
class(stockSymDF$Date)

x = "RAVE"
getSymbols(x,src='yahoo', env = NULL)

stockSym = getSymbols("AAPL",src='yahoo', env = NULL)
colnames(stockSym) = c("Open", "High", "Low", "Close", "Volume", "Adjusted")
stockSymDF = as.data.frame(stockSym)
stockSymDF$Date = as.Date(rownames(stockSymDF))
stockSymDF = stockSymDF[stockSymDF$Date >= input$Date[1] & stockSymDF$Date <= input$Date[2],]
stockSymDF$PercentChange = stockSymDF$Close / lag(stockSymDF$Close, 1) - 1
stockSymDF$PercentChange = scales::percent(stockSymDF$PercentChange)
stockSymDF = stockSymDF[,c(7,1,2,3,4,5,6,8)]
rownames(stockSymDF) = NULL
stockSymDF = stockSymDF[stockSymDF$Date >= as.Date("2019-01-20") & stockSymDF$Date <= as.Date("2019-02-20"),]
stockSymDF$PercentChange = stockSymDF$Close / lag(stockSymDF$Close, 1) - 1
stockSymDF$PercentChange = scales::percent(stockSymDF$PercentChange)
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
  

df %>%
  plot_ly(x = ~Date, type="candlestick",
          open = ~RAVE.Open, close = ~RAVE.Close,
          high = ~RAVE.High, low = ~RAVE.Low) %>%
  layout(title = "Basic Candlestick Chart",
         paper_bgcolor='#404040',
         plot_bgcolor='#404040',
         xaxis = list(color = '#ffffff'),
         yaxis = list(color = '#ffffff'),
         titlefont = list(size = 24, color = '#ffffff'))


library(lubridate)

Sys.Date() %m-% months(1)


library(fGarch)
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

spy = getSymbols("SPY", auto.assign=FALSE)
(log(Cl(spy)) - log(lag(Cl(spy))))

rets = ROC(Cl(spy), na.pad=FALSE)
fit = garchAuto(rets, cores=8, trace=TRUE)
Cl(spy)[length(Cl(spy))]
exp(as.numeric(predict(fit,1)[1]) + as.numeric(log(Cl(spy)[length(Cl(spy))])))
predict(fit,25)




library(fGarch)
sp5=read.table("http://faculty.chicagobooth.edu/ruey.tsay/teaching/fts/sp500.dat")#Load data
plot(sp5,type="l")
m1=garchFit(formula=~arma(3,0)+garch(1,1),data=sp5,trace=F)
summary(m1)
m2=garchFit(formula=~garch(1,1),data=sp5,trace=F,cond.dist="std")
summary(m2)
stresi=residuals(m2,standardize=T)
plot(stresi,type="l")
Box.test(stresi,10,type="Ljung")
predict(m2,5)
