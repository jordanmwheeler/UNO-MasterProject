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
 