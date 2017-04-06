library(forecast)
library(knitr)
library(tseries)
library(lawstat)
library(vars)

# checking residual diagnostics
residual_diagnostics <- function(model) {
  # test for zero mean
  par(mfrow=c(1,1))
  # get residuals
  modelresiduals <- model$residuals
  # informal
  plot(modelresiduals, main="Residuals vs t", ylab="")
  abline(h=0, col="red")
  # formal
  ttest <- t.test(modelresiduals)
  print(paste("T-Test Statistic:",round(ttest$statistic,4),", P-value:",round(ttest$p.value,4)))
  
  # test for homoskedasticity
  # formal
  group <- c(rep(1,72),rep(2,72),rep(3,72),rep(4,72))
  levenetest <- levene.test(modelresiduals,group)
  print(paste("Levene Test Statistic:",round(levenetest$statistic,4),
              ", P-value:",round(levenetest$p.value,4)))
  
  # test for correlatedness
  tsdiag(model)
  
  # test for normality
  # informal
  qqnorm(modelresiduals, main="QQ-plot of Residuals")
  qqline(modelresiduals, col="red")
  # formal
  shapirowilktest <- shapiro.test(modelresiduals)
  print(paste("Shapiro Wilk Test Statistic:",round(shapirowilktest$statistic,4),
              ", P-value:",round(shapirowilktest$p.value,4)))
}

HW_residual_diagnostics <- function(resids) {
  # residual diagnostics
  # test for zero mean
  par(mfrow=c(1,1))
  # get residuals
  modelresiduals <- resids - bankruptcy_rate_train[13:288]
  # informal
  plot(modelresiduals, main="Residuals vs t", ylab="")
  abline(h=0, col="red")
  # formal
  ttest <- t.test(modelresiduals)
  print(paste("T-Test Statistic:",round(ttest$statistic,4),", P-value:",round(ttest$p.value,4)))
  
  # test for homoskedasticity
  # formal
  group <- c(rep(1,69),rep(2,69),rep(3,69),rep(4,69))
  levenetest <- levene.test(modelresiduals,group)
  print(paste("Levene Test Statistic:",round(levenetest$statistic,4),
              ", P-value:",round(levenetest$p.value,4)))
  
  # test for correlatedness
  acf(modelresiduals)
}

# read in data
train <- read.csv('train.csv', header = TRUE)
test <- read.csv('test.csv', header = TRUE)
# create train time series
bankruptcy_rate_train <- ts(train$Bankruptcy_Rate, start = c(1987,1), end = c(2010,12), frequency = 12)
unemploy_rate_train <- ts(train$Unemployment_Rate, start = c(1987,1), end = c(2010,12), frequency = 12)
population_train <- ts(train$Population, start = c(1987,1), end = c(2010,12), frequency = 12)
housing_PI_train <- ts(train$House_Price_Index, start = c(1987,1), end = c(2010,12), frequency = 12)
# create test time series
unemploy_rate_test <- ts(test$Unemployment_Rate, start = c(2011,1), end = c(2012,12), frequency = 12)
population_test <- ts(test$Population, start = c(2011,1), end = c(2012,12), frequency = 12)
housing_PI_test <- ts(test$House_Price_Index, start = c(2011,1), end = c(2012,12), frequency = 12)
# split train time series into a train set and test set
bankruptcy_rate_train_train <- ts(bankruptcy_rate_train[1:228], start = c(1987,1), end = c(2005,12), frequency = 12)
unemploy_rate_train_train <- ts(unemploy_rate_train[1:228], start = c(1987,1), end = c(2005,12), frequency = 12)
population_train_train <- ts(population_train[1:288], start = c(1987,1), end = c(2005,12), frequency = 12)
housing_PI_train_train <- ts(housing_PI_train[1:288], start = c(1987,1), end = c(2005,12), frequency = 12)
bankruptcy_rate_train_test <- ts(bankruptcy_rate_train[229:288], start = c(2006,1), end = c(2010,12), frequency = 12)
unemploy_rate_train_test <- ts(unemploy_rate_train[229:288], start = c(2006,1), end = c(2010,12), frequency = 12)
population_train_test <- ts(population_train[229:288], start = c(2006,1), end = c(2010,12), frequency = 12)
housing_PI_train_test <- ts(housing_PI_train[229:288], start = c(2006,1), end = c(2010,12), frequency = 12)

# visualizing data
plot(bankruptcy_rate_train)

par(mfrow=c(4,1))
plot(bankruptcy_rate_train)
plot(unemploy_rate_train)
plot(population_train)
plot(housing_PI_train)

######################
#### Holt-Winters ####
######################

par(mfrow=c(1,1))
acf(bankruptcy_rate_train_train, lag.max = 72) # there is seasonality and trend

# create for loop to find the best values for the smoothing parameters
fit_model <- data.frame(character(0), character(0), character(0), stringsAsFactors=F)
for (a in seq(.1,.9,.1)) {
  for (b in seq(.1,.9,.1)) {
    for (g in seq(.1,.9,.1)) {
      m <- HoltWinters(x = bankruptcy_rate_train_train, alpha = a, beta = b, gamma = g, seasonal = "mult")
      m_forecast <- forecast(m, h = 60)
      m_RMSE <- sqrt(mean((bankruptcy_rate_train_test - m_forecast$mean)**2))
      m <- HoltWinters(x = bankruptcy_rate_train, alpha = a, beta = b, gamma = g, seasonal = "mult")
      m_residuals <- m$fitted[,1] - bankruptcy_rate_train[13:288]
      group <- c(rep(1,69),rep(2,69),rep(3,69),rep(4,69))
      levenetest <- levene.test(m_residuals,group)
      fit_model <- rbind.data.frame(fit_model, c(paste("alpha=",a,",beta=",b,",gamma=",g,sep=''), as.character(m_RMSE), as.character(levenetest$p.value)), stringsAsFactors=F)
    }
  }
}
colnames(fit_model) <- c("model", "RMSE", "levenetest_pval")
fit_model$RMSE <- as.numeric(fit_model$RMSE)
fit_model$levenetest_pval <- as.numeric(fit_model$levenetest_pval)
# to print out top 5 best models
head(fit_model[order(fit_model$RMSE),], 5)

fit_model_HW <- fit_model

# best model chosen by predicted RMSE - does not satisfy homoskedasticity, uncorrelatedness assumptions
m <- HoltWinters(x = bankruptcy_rate_train_train, alpha = .4, beta = .3, gamma = .5, seasonal = "mult")
m_forecast <- forecast(m, h = 60)
m_RMSE <- sqrt(mean((bankruptcy_rate_train_test - m_forecast$mean)**2)) # 0.003712317
# displaying forecasts vs true values for best model chosen by predicted RMSE
plot(m_forecast, main = "Bankruptcy Rate Forecasts Using H-W Model vs True Values", ylab = "Bankruptcy Rate", xlab = "Time")
points(bankruptcy_rate_train_test, type="l", col="red")

# best model chosen by predicted RMSE - does not satisfy homoskedasticity, uncorrelatedness assumptions
m <- HoltWinters(x = bankruptcy_rate_train, alpha = .4, beta = .3, gamma = .5, seasonal = "mult")
HW_residual_diagnostics(m$fitted[,1])

# second best model chosen by predicted RMSE - does not satisfy homoskedasticity, uncorrelatedness assumptions
m <- HoltWinters(x = bankruptcy_rate_train, alpha = .6, beta = .7, gamma = .8, seasonal = "mult")
HW_residual_diagnostics(m$fitted[,1])

# third best model chosen by predicted RMSE - does not satisfy homoskedasticity, uncorrelatedness assumptions
m <- HoltWinters(x = bankruptcy_rate_train, alpha = .1, beta = .5, gamma = .1, seasonal = "mult")
HW_residual_diagnostics(m$fitted[,1])

# fourth best model chosen by predicted RMSE - does not satisfy homoskedasticity, uncorrelatedness assumptions
m <- HoltWinters(x = bankruptcy_rate_train, alpha = .1, beta = .5, gamma = .2, seasonal = "mult")
HW_residual_diagnostics(m$fitted[,1])

# fifth best model chosen by predicted RMSE - does not satisfy homoskedasticity, uncorrelatedness assumptions
m <- HoltWinters(x = bankruptcy_rate_train, alpha = .3, beta = .5, gamma = .2, seasonal = "mult")
HW_residual_diagnostics(m$fitted[,1])

# forecast into the future
m <- HoltWinters(x = bankruptcy_rate_train, alpha = .5, beta = .1, gamma = .1, seasonal = "add")
m_forecast <- forecast(m, h = 24)
plot(m_forecast, main = "Bankruptcy Rate Forecasts Using Holt-Winters Model", ylab = "Bankruptcy Rate", xlab = "Time")
forecast_table_df <- data.frame(as.yearmon(time(housing_PI_test)), "Prediction" = m_forecast$mean, "95% CI Lower Bound" = m_forecast$lower[,1], "95% CI Upper Bound" = m_forecast$upper[,1])
kable(forecast_table_df, col.names = c("", "Prediction", "95% CI Lower Bound", "95% CI Upper Bound"))


################
#### SARIMA ####
################

#### 1 seasonal difference, 1 ordinary difference ####

# log transform the data
bankruptcy_rate_train_train_log <- log(bankruptcy_rate_train_train)
bankruptcy_rate_train_log <- log(bankruptcy_rate_train)

plot(bankruptcy_rate_train_log)

# check for trend and/or seasonality
acf(bankruptcy_rate_train_train_log, lag.max = 72) # there is trend and seasonality, period=12

# perform seasonal differencing to remove seasonality
bankruptcy_rate_train_train_log_sdiff1 <- diff(bankruptcy_rate_train_train_log, lag=12)
acf(bankruptcy_rate_train_train_log_sdiff1, lag.max = 72)

# perform ordinary differencing to remove trend
bankruptcy_rate_train_train_log_sdiff1_diff1 <- diff(bankruptcy_rate_train_train_log_sdiff1)
acf(bankruptcy_rate_train_train_log_sdiff1_diff1, lag.max = 72)

# 1 ordinary differece, 1 seasonal difference

par(mfrow=c(2,1))
acf(bankruptcy_rate_train_train_log_sdiff1_diff1, lag.max = 72) # q=1, Q=0
pacf(bankruptcy_rate_train_train_log_sdiff1_diff1, lag.max = 72) # p=2, P=0

m <- arima(bankruptcy_rate_train_train_log, order = c(1,1,2), seasonal = list(order = c(0,1,0), period = 24))
m$residuals

# fitting models in the same neighborhood as SARIMA model chosen by looking at ACF/PACF plots to find optimal one
fit_model <- data.frame(character(0), character(0), character(0),
                        character(0), character(0), stringsAsFactors=F)
for (p in c(0:3)) {
  for (q in c(0:3)) {
    for (P in c(0:2)) {
      for (Q in c(0:2)) {
        for (methodname in c("CSS", "ML")) {
          try({m <- arima(bankruptcy_rate_train_train_log, order=c(p,1,q), seasonal=list(order=c(P,1,Q), period=12), method=methodname)
          m_forecast <- forecast(m, h=60, level=0.95)
          m_predictions <- exp(m_forecast$mean)
          m_RMSE <- sqrt(mean((bankruptcy_rate_train_test - m_predictions)**2))
          fit_model <- rbind.data.frame(fit_model,
                                        c(paste("method=",methodname,",p=",p,",q=",q,",P=",P,",Q=",Q,sep=''),
                                          as.character(m$loglik), as.character(m$sigma2),
                                          as.character(m$aic), as.character(m_RMSE)), stringsAsFactors=F)})
        }
      }
    }
  }
}
colnames(fit_model) = c("model", "loglik", "sigma2", "AIC", "predictRMSE")
fit_model$loglik <- as.numeric(fit_model$loglik)
fit_model$sigma2 <- as.numeric(fit_model$sigma2)
fit_model$AIC <- as.numeric(fit_model$AIC)
fit_model$predictRMSE <- as.numeric(fit_model$predictRMSE)

fit_model_SARIMA <- fit_model

# best by predicted RMSE - doesn't satisfy uncorrelatedness assumption
m010_211 <- arima(bankruptcy_rate_train_log, order=c(0,1,0), seasonal=list(order=c(2,1,1), period=12), method="ML")
residual_diagnostics(m010_211)

# second best by predicted RMSE - doesn't satisfy uncorrelatedness assumption
m110_211 <- arima(bankruptcy_rate_train_log, order=c(1,1,0), seasonal=list(order=c(2,1,1), period=12), method="ML")
residual_diagnostics(m110_211)

# third best by predicted RMSE - doesn't satisfy uncorrelatedness assumption
m010_011 <- arima(bankruptcy_rate_train_log, order=c(0,1,0), seasonal=list(order=c(0,1,1), period=12), method="ML")
residual_diagnostics(m010_011)

# fourth best by predicted RMSE - doesn't satisfy uncorrelatedness assumption
m310_012 <- arima(bankruptcy_rate_train_log, order=c(3,1,0), seasonal=list(order=c(0,1,2), period=12), method="ML")
residual_diagnostics(m310_012)

# fifth best by predicted RMSE - doesn't satisfy uncorrelatedness assumption
m211_012 <- arima(bankruptcy_rate_train_log, order=c(2,1,1), seasonal=list(order=c(0,1,2), period=12), method="ML")
residual_diagnostics(m211_012)

#### 1 ordinary difference, 0 seasonal difference ####

# perform ordinary difference to remove trend
bankruptcy_rate_train_train_log_diff1 <- diff(bankruptcy_rate_train_train_log)
acf(bankruptcy_rate_train_train_log_diff1, lag.max = 72)

par(mfrow=c(2,1))
acf(bankruptcy_rate_train_train_log_diff1, lag.max = 72) # q=0
pacf(bankruptcy_rate_train_train_log_diff1, lag.max = 72) # p=1

# fitting models in the same neighborhood as SARIMA model chosen by looking at ACF/PACF plots to find optimal one
fit_model <- data.frame(character(0), character(0), character(0),
                        character(0), character(0), stringsAsFactors=F)
for (p in c(0:3)) {
  for (q in c(0:3)) {
    for (methodname in c("CSS", "ML")) {
      try({m <- arima(bankruptcy_rate_train_train_log, order=c(p,1,q), method=methodname)
      m_forecast <- forecast(m, h=60, level=0.95)
      m_predictions <- exp(m_forecast$mean)
      m_RMSE <- sqrt(mean((bankruptcy_rate_train_test - m_predictions)**2))
      fit_model <- rbind.data.frame(fit_model,
                                    c(paste("method=",methodname,",p=",p,",q=",q,sep=''),
                                      as.character(m$loglik), as.character(m$sigma2),
                                      as.character(m$aic), as.character(m_RMSE)), stringsAsFactors=F)})
    }
  }
}
colnames(fit_model) = c("model", "loglik", "sigma2", "AIC", "predictRMSE")
fit_model$loglik <- as.numeric(fit_model$loglik)
fit_model$sigma2 <- as.numeric(fit_model$sigma2)
fit_model$AIC <- as.numeric(fit_model$AIC)
fit_model$predictRMSE <- as.numeric(fit_model$predictRMSE)

fit_model_SARIMA_1diff <- fit_model

# best model chosen by predicted RMSE - doesn't satisfy uncorrelatedness assumption
m011 <- arima(bankruptcy_rate_train_log, order=c(0,1,1), method="CSS")
residual_diagnostics(m011)


######################################
#### SARIMA + exogenous variables ####
######################################

#### 1 ordinary difference, 1 seasonal difference ####

unemploy_rate_train_train_log <- log(unemploy_rate_train_train)
unemploy_rate_train_test_log <- log(unemploy_rate_train_test)
unemploy_rate_train_log <- log(unemploy_rate_train)
unemploy_rate_test_log <- log(unemploy_rate_test)
population_train_train_log <- log(population_train_train)
population_train_test_log <- log(population_train_test)
population_train_log <- log(population_train)
population_test_log <- log(population_test)
housing_PI_train_train_log <- log(housing_PI_train_train)
housing_PI_train_test_log <- log(housing_PI_train_test)
housing_PI_train_log <- log(housing_PI_train)
housing_PI_test_log <- log(housing_PI_test)

# to determine if we should log transform the exogenous variables
cor(bankruptcy_rate_train_log, unemploy_rate_train)
cor(bankruptcy_rate_train_log, unemploy_rate_train_log) # log of unemployment rate is better
cor(bankruptcy_rate_train_log, population_train)
cor(bankruptcy_rate_train_log, population_train_log) # log of population is better
cor(bankruptcy_rate_train_log, housing_PI_train)
cor(bankruptcy_rate_train_log, housing_PI_train_log) # log of housing PI is better

# fitting models in the same neighborhood as SARIMA model chosen by looking at ACF/PACF plots to find optimal one
fit_model <- data.frame(character(0), character(0), character(0),
                        character(0), character(0), character(0), stringsAsFactors=F)
for (p in c(0:8)) {
  for (q in c(0:8)) {
    for (P in c(0:8)) {
      for (Q in c(0:8)) {
        for (methodname in c("CSS", "ML")) {
          for (df in list(data.frame(unemploy_rate_train_log), data.frame(population_train_log), 
                          data.frame(housing_PI_train_log), data.frame(unemploy_rate_train_log, population_train_log), 
                          data.frame(unemploy_rate_train_log, housing_PI_train_log), 
                          data.frame(population_train_log, housing_PI_train_log), 
                          data.frame(unemploy_rate_train_log, population_train_log, housing_PI_train_log))) {
            try({m <- arima(bankruptcy_rate_train_train_log, order=c(p,1,q), seasonal=list(order=c(P,1,Q), period=12), xreg = df[1:228,, drop=FALSE], method=methodname)
            m_forecast <- forecast(m, h=60, level=0.95, xreg = df[229:288,, drop=FALSE])
            m_predictions <- exp(m_forecast$mean)
            m_RMSE <- sqrt(mean((bankruptcy_rate_train_test - m_predictions)**2))
            fit_model <- rbind.data.frame(fit_model,
                                          c(paste("method=",methodname,",p=",p,",q=",q,",P=",P,",Q=",Q,sep=''),
                                            as.character(paste(colnames(df),collapse=',')),
                                            as.character(m$loglik), as.character(m$sigma2),
                                            as.character(m$aic), as.character(m_RMSE)), stringsAsFactors=F)})
          }
        }
      }
    }
  }
}
colnames(fit_model) = c("model", "exogenous variables", "loglik", "sigma2", "AIC", "predRMSE")
fit_model$loglik <- as.numeric(fit_model$loglik)
fit_model$sigma2 <- as.numeric(fit_model$sigma2)
fit_model$AIC <- as.numeric(fit_model$AIC)
fit_model$predRMSE <- as.numeric(fit_model$predRMSE)
head(fit_model[order(fit_model$predRMSE),], 5)

fit_model_ARIMAX <- fit_model

# best model chosen by predicted RMSE - doesn't satisfy uncorrelatedness assumption
m011_011 <- arima(bankruptcy_rate_train_log, order=c(0,1,1), seasonal=list(order=c(0,1,1), period=12), xreg=data.frame(unemploy_rate_train_log), method="CSS")
residual_diagnostics(m011_011)

# second best model chosen by predicted RMSE (0.002882778) - satisfies all assumptions! (shapiro-wilk fails, but QQ-plot looks good)
m213_212 <- arima(bankruptcy_rate_train_log, order=c(2,1,3), seasonal=list(order=c(2,1,2), period=12), xreg=data.frame(unemploy_rate_train_log, population_train_log), method="ML")
residual_diagnostics(m213_212)
# getting forecast into future
m213_212_forecast <- forecast(m213_212, h=24, level=0.95, xreg = data.frame(unemploy_rate_test_log, population_test_log))
m213_212_predictions <- exp(m213_212_forecast$mean)
m213_212_CIlower <- ts(exp(m213_212_forecast$lower), start = c(2011,1), frequency = 12)
m213_212_CIupper <- ts(exp(m213_212_forecast$upper), start = c(2011,1), frequency = 12)
par(mfrow=c(1,1))
plot(bankruptcy_rate_train, xlim=c(1987, 2013), ylim=c(0,0.06), main = "Forecasts from ARIMAX(2,1,3)(2,1,2)_12", ylab = "Bankruptcy Rate", xlab = "Time")
abline(v = 2011, lwd = 2, col = "black")
points(m213_212_predictions, type = "l", col = "blue")
points(m213_212_CIlower, type = "l", col = "red")
points(m213_212_CIupper, type = "l", col = "red")
legend("topleft", legend = c("Observed", "Predicted", "95% PI"), lty = 1, col = c("black", "blue", "red"), cex = 1)
# getting predicted RMSE for second best mdoel
m213_212 <- arima(bankruptcy_rate_train_train_log, order=c(2,1,3), seasonal=list(order=c(2,1,2), period=12), xreg=data.frame(unemploy_rate_train_train_log, population_train_train_log), method="ML")
m213_212_forecast <- forecast(m213_212, h=60, level=0.95, xreg = data.frame(unemploy_rate_train_test_log, population_train_test_log))
m213_212_predictions <- exp(m213_212_forecast$mean)
m213_212_RMSE <- sqrt(mean((bankruptcy_rate_train_test - m213_212_predictions)**2))

# third best model chosen by predicted RMSE - doesn't satisfy uncorrelatedness, normality assumptions
m210_111 <- arima(bankruptcy_rate_train_log, order=c(2,1,0), seasonal=list(order=c(1,1,1), period=12), xreg=data.frame(unemploy_rate_train_log, population_train_log), method="CSS")
residual_diagnostics(m210_111)

# fourth best model chosen by predicted RMSE - doesn't satisfy uncorrelatedness assumption
m111_011 <- arima(bankruptcy_rate_train_log, order=c(1,1,1), seasonal=list(order=c(0,1,1), period=12), xreg=data.frame(unemploy_rate_train_log), method="CSS")
residual_diagnostics(m111_011)

# seventh best model chosen by predicted RMSE (0.003045657) - satisfies all assumptions! (shapiro-wilk fails, but QQ-plot looks good)
m113_212 <- arima(bankruptcy_rate_train_log, order=c(1,1,3), seasonal=list(order=c(2,1,2), period=12), xreg=data.frame(unemploy_rate_train_log, population_train_log), method="ML")
residual_diagnostics(m113_212)

#### 1 ordinary difference, 0 seasonal difference ####

# fitting models in the same neighborhood as SARIMA model chosen by looking at ACF/PACF plots to find optimal one
fit_model <- data.frame(character(0), character(0), character(0),
                        character(0), character(0), character(0), stringsAsFactors=F)
for (p in c(0:3)) {
  for (q in c(0:3)) {
    for (methodname in c("CSS", "ML")) {
      for (df in list(data.frame(unemploy_rate_train_log), data.frame(population_train_log), 
                      data.frame(housing_PI_train_log), data.frame(unemploy_rate_train_log, population_train_log), 
                      data.frame(unemploy_rate_train_log, housing_PI_train_log), 
                      data.frame(population_train_log, housing_PI_train_log), 
                      data.frame(unemploy_rate_train_log, population_train_log, housing_PI_train_log))) {
        try({m <- arima(bankruptcy_rate_train_train_log, order=c(p,1,q), xreg = df[1:228,, drop=FALSE], method=methodname)
        m_forecast <- forecast(m, h=60, level=0.95, xreg = df[229:288,, drop=FALSE])
        m_predictions <- exp(m_forecast$mean)
        m_RMSE <- sqrt(mean((bankruptcy_rate_train_test - m_predictions)**2))
        fit_model <- rbind.data.frame(fit_model,
                                      c(paste("method=",methodname,",p=",p,",q=",q,sep=''),
                                        as.character(paste(colnames(df),collapse=',')),
                                        as.character(m$loglik), as.character(m$sigma2),
                                        as.character(m$aic), as.character(m_RMSE)), stringsAsFactors=F)})
      }
    }
  }
}
colnames(fit_model) = c("model", "exogenous variables", "loglik", "sigma2", "AIC", "predRMSE")
fit_model$loglik <- as.numeric(fit_model$loglik)
fit_model$sigma2 <- as.numeric(fit_model$sigma2)
fit_model$AIC <- as.numeric(fit_model$AIC)
fit_model$predRMSE <- as.numeric(fit_model$predRMSE)
head(fit_model[order(fit_model$predRMSE),], 5)

fit_model_ARIMAX_1diff <- fit_model

# best model chosen by predicted RMSE - does not satisfy uncorrelatedness assumption
m013 <- arima(bankruptcy_rate_train_log, order=c(0,1,3), xreg=data.frame(unemploy_rate_train_log, population_train_log), method = "ML")
residual_diagnostics(m013)

# second best model chosen by predicted RMSE - does not satisfy uncorrelatedness assumption
m31 <- arima(bankruptcy_rate_train_log, order=c(3,1,1), xreg=data.frame(unemploy_rate_train_log, population_train_log), method = "ML")
residual_diagnostics(m31)

# third best model chosen by predicted RMSE - does not satisfy uncorrelatedness assumption
m21 <- arima(bankruptcy_rate_train_log, order=c(2,1,1), xreg=data.frame(unemploy_rate_train_log, population_train_log), method = "CSS")
residual_diagnostics(m21)

# fourth best model chosen by predicted RMSE - does not satisfy uncorrelatedness assumption
m01 <- arima(bankruptcy_rate_train_log, order=c(0,1,1), xreg=data.frame(unemploy_rate_train_log, population_train_log), method = "ML")
residual_diagnostics(m01)


###################################
#### VAR + seasonal indicators ####
###################################

fit_model <- data.frame(character(0), character(0), character(0), stringsAsFactors=F)
for (pval in c(1:10)) {
  for (df in list(data.frame(bankruptcy_rate_train_train, unemploy_rate_train_train), data.frame(bankruptcy_rate_train_train, population_train_train), 
                  data.frame(bankruptcy_rate_train_train, housing_PI_train_train), data.frame(bankruptcy_rate_train_train, unemploy_rate_train_train, population_train_train), 
                  data.frame(bankruptcy_rate_train_train, unemploy_rate_train_train, housing_PI_train_train), 
                  data.frame(bankruptcy_rate_train_train, population_train_train, housing_PI_train_train), 
                  data.frame(bankruptcy_rate_train_train, unemploy_rate_train_train, population_train_train, housing_PI_train_train))) {
    try({m <- VAR(y = df, p = pval)
    m_forecast <- predict(m, n.ahead = 60, ci = 0.95)
    m_predictions <- m_forecast$fcst$bankruptcy_rate_train_train[,1]
    m_RMSE <- sqrt(mean((bankruptcy_rate_train_test - m_predictions)**2))
    fit_model <- rbind.data.frame(fit_model,
                                  c(paste("p=",pval),
                                    as.character(paste(colnames(df),collapse=',')),
                                    as.character(m_RMSE)), stringsAsFactors=F)})
  }
}
colnames(fit_model) = c("model", "endogenous variables", "predRMSE")
fit_model$predRMSE <- as.numeric(fit_model$predRMSE)
head(fit_model[order(fit_model$predRMSE),], 5)

fit_model_VAR <- fit_model

# best model chosen by predicted RMSE - doesn't satisfy uncorrelatedness, homoskedasticity assumptions for all endogenous variables
m8 <- VAR(y = data.frame(bankruptcy_rate_train, population_train, housing_PI_train), p = 8)
plot(m8)

# second best model chosen by predicted RMSE - doesn't satisfy uncorrelatedness
m8 <- VAR(y = data.frame(bankruptcy_rate_train, unemploy_rate_train, population_train, housing_PI_train), p = 8)
plot(m8)

# third best model chosen by predicted RMSE
m10 <- VAR(y = data.frame(bankruptcy_rate_train, housing_PI_train), p=10)
plot(m10)

# fourth best model chosen by predicted RMSE
m9 <- VAR(y = data.frame(bankruptcy_rate_train, population_train, housing_PI_train), p=9)
plot(m9)

# fifth best model chosen by predicted RMSE
m7 <- VAR(y = data.frame(bankruptcy_rate_train, population_train, housing_PI_train), p=7)
plot(m7)

