

##################### LIBRARY ###########################################################

library(tseries)
library(forecast)
library(Metrics)
library(ggplot2)
library(readr)
library(WaveletArima)
library(caret)
library(nnfor)
library(tsDyn)
library(fracdiff)
library(bsts)
library(forecastHybrid)
library(e1071)
library(tseriesChaos)
library(pracma)
library(Kendall)
library(nonlinearTseries)
library(GeneCycle)
library(fpp2)
library(ggfortify)
library(MTS)
library(WaveletArima)
library(smooth)

####################### SanJuan Data ########################

data = read.csv("Sanjuan_data_weekly.csv", stringsAsFactors = F)
df=data[,c("Cases","Rain")]
con_tr = df$Cases[37:1128]
xreg_tr = df$Rain[37:1128]
con_tst = df$Cases[1129:1180]

##################### Iquitos data ###########################

# data = read.csv("Iquitos_data_weekly.csv", stringsAsFactors = F)
# df=data[,c("Cases","Rain")]
# con_tr = df$Cases[1:520]
# xreg_tr = df$Rain[1:520]
# con_tst = df$Cases[521:572]

####################### Ahmedabad Data ########################

# data = read.csv("ahmedabad_dengue_data.csv", stringsAsFactors = F)
# df=data[,c("Cases","Rainfall")]
# con_tr = df$Cases[1:371]
# xreg_tr = df$Rain[1:371]
# con_tst = df$Cases[372:423]

######################## PLOTS #################################

autoplot(ts(con_tr))+ ggtitle('ahmedabad_dengue_data')+ylab('Cases')+xlab('week')+autolayer(ts(con_tst))
# diffset = diff(con_tr, differences = ndiffs(con_tr))
ggAcf(con_tr)+geom_point(color = 'navy blue') +ggtitle("ACF plot")
ggPacf(con_tr) + geom_point(color = 'navy blue') + ggtitle("PACF plot")

################# Statistical Tests ############################

fracdiff(con_tr)
hurstexp(con_tr)
Box.test(con_tr, lag = 1, type = c("Box-Pierce"))
skewness(con_tr)
kurtosis(con_tr)
max_lyapunov_expo <- lyap_k(con_tr, m=1, d=2, s=1, t=4, ref=length(con_tr), k=2, eps=4)
kpss.test(con_tr) 
nonlinearityTest(con_tr, verbose = TRUE)

################### ARIMA ###################################

fit_arima = auto.arima(as.ts(con_tr))
forecast_arima =  as.numeric(forecast(fit_arima, 52)$mean)
forecast::accuracy(forecast_arima, con_tst)
smape(con_tst, forecast_arima)*100
mase(con_tst, forecast_arima)

# plot(as.ts(c(con_tr)))
# lines(c(fitted(fit_arima)), type = "l", col = "red")

######################### Wavelet ARIMA ###########################

fit_warima = WaveletFittingarma(ts(con_tr), Waveletlevels = floor(log(length(con_tr))),boundary = "periodic", 
                                FastFlag = TRUE, MaxARParam = 10, MaxMAParam = 10, NForecast = 52)
fit_warima$FinalPrediction[fit_warima$FinalPrediction <0] <- 0
forecast_warima = as.data.frame(fit_warima$Finalforecast, h=52)
forecast::accuracy(fit_warima$Finalforecast, con_tst)
smape(con_tst, fit_warima$Finalforecast)*100
mase(con_tst, fit_warima$Finalforecast)

# plot(as.ts(c(xreg_tr)))
# lines(c(fit$FinalPrediction), type = "l", col = "red")

################## fitting ARIMAX ########################

fit_arimax =auto.arima(as.ts(con_tr), xreg=xreg_tr)
forecast_arimax = as.numeric(forecast(fit_arimax, 52,  xreg= xreg_tr[1:52])$mean)
forecast::accuracy(forecast_arimax, con_tst)
smape(con_tst, forecast_arimax)*100
mase(con_tst, forecast_arimax)

# plot(as.ts(c(con_tst)))
# lines(c(forecasts), type = "l", col = "red")

################ fitting ETS model ###########################

fit_ets = ets(as.ts(con_tr))
forecast_ets = forecast::forecast(fit_ets, h = 52)
forecast::accuracy(forecast_ets$mean, con_tst)
smape(con_tst, forecast_ets$mean)*100
mase(con_tst, forecast_ets$mean)

################## fitting ETSX ########################

fit_etsx = es(as.ts(con_tr), xreg=xreg_tr)
forecast_etsx = forecast.smooth(fit_etsx, h = 52)
forecast::accuracy(forecast_etsx$mean, con_tst)
smape(con_tst, forecast_etsx$mean)*100
mase(con_tst, forecast_etsx$mean)

################### SETAR ###########################

fit_SETAR = setar(as.ts(con_tr), m = 4)
forecast_SETAR = predict(fit_SETAR, n.ahead = 52)
forecast::accuracy(forecast_SETAR, con_tst)
smape(con_tst, forecast_SETAR)*100
mase(con_tst, forecast_SETAR)

################## TBATS ############################

fit_tbats = tbats(as.ts(con_tr))
forecast_TBATS=forecast::forecast(fit_tbats, h=52)
forecast::accuracy(forecast_TBATS, con_tst)
smape(con_tst, forecast_TBATS$mean)*100
mase(con_tst, forecast_TBATS$mean)

######################### Theta ################################

fit_theta=thetaf(as.ts(con_tr), h=52)
forecast::accuracy(fit_theta$mean, con_tst)
smape(con_tst, fit_theta$mean)*100
mase(con_tst, fit_theta$mean)

################### MLP/ANN ###################################

fit_ANN = mlp(as.ts(con_tr), hd = 5, reps = 1)
forecast_ANN = forecast::forecast(fit_ANN, h=52)
forecast::accuracy(forecast_ANN$mean, con_tst)
smape(con_tst, forecast_ANN$mean)*100
mase(con_tst, forecast_ANN$mean)

######################### BSTS ###############################

ss <- AddLocalLinearTrend(list(), as.ts(con_tr))
fit_bsts = bsts(as.ts(con_tr),state.specification = ss, niter = 1000)
forecast_BSTS <- predict(fit_bsts, horizon = 52)
burn <- SuggestBurn(0.1, fit_bsts)
fitted_bsts=as.numeric(-colMeans(fit_bsts$one.step.prediction.errors[-(1:burn),])+as.ts(con_tr))
forecast::accuracy(forecast_BSTS$mean, con_tst)
smape(con_tst, forecast_BSTS$mean)*100
mase(con_tst, forecast_BSTS$mean)

######################### ARFIMA ###############################

fit_ARFIMA=arfima( as.ts(con_tr))
forecast_ARFIMA = forecast::forecast(fit_ARFIMA, h=52)
forecast::accuracy(forecast_ARFIMA$mean, con_tst)
smape(con_tst, forecast_ARFIMA$mean)*100
mase(con_tst, forecast_ARFIMA$mean)

########################## ARNN ################################

fit_NNAR = nnetar(con_tr, repeats = 500, lambda="auto")
fc1 = forecast(fit_NNAR, h = 52)
fc1_df = data.frame(fc1)
forecast::accuracy(fc1_df$Point.Forecast, con_tst)
smape(con_tst, fc1_df$Point.Forecast)*100
mase(con_tst, fc1_df$Point.Forecast)

######################### NNARX model ##########################

fit_nnarx = nnetar(as.ts(con_tr), xreg = xreg_tr, repeats = 500, lambda="auto")
pred = as.numeric(forecast(fit_nnarx, h = 52, xreg = xreg_tr[1:52])$mean)
forecast::accuracy(pred, con_tst)
smape(con_tst, pred)*100
mase(con_tst, pred)

####################### Proposed WARNN ###########################

source("wnnar.R")
fit_wnnar = WaveletFittingnar(ts(con_tr), Waveletlevels = floor(log(length(con_tr))), boundary = "periodic", 
            FastFlag = TRUE, MaxARParam = 10, NForecast = 52)
fore_wnnar = as.data.frame(fit_wnnar$Finalforecast, h = 52)
forecast::accuracy(fore_wnnar$`fit_wnnar$Finalforecast`, con_tst)
smape(con_tst, fore_wnnar$`fit_wnnar$Finalforecast`)*100
mase(con_tst, fore_wnnar$`fit_wnnar$Finalforecast`)

####################### Proposed WARNNX ##########################

source("wavelet_nar.R")
fit_wnnarx = WaveletFittingnar(ts(con_tr), Waveletlevels = floor(log(length(con_tr))), 
             boundary = "periodic", FastFlag = TRUE, MaxARParam = 10, NForecast = 52)
fore_wnnarx = as.data.frame(fit_wnnarx$Finalforecast, h = 52)
forecast::accuracy(fore_wnnarx$`fit_wnnarx$Finalforecast`, con_tst)
smape(con_tst, fore_wnnarx$`fit_wnnarx$Finalforecast`)*100
mase(con_tst, fore_wnnarx$`fit_wnnarx$Finalforecast`)

######################### Hybrid Models ########################

### Hybrid ARIMA-ARNN

fit_res_ARNN = nnetar(fit_arima$residuals, repeats = 500, lambda="auto")
pred_res_ARNN = forecast::forecast(fit_res_ARNN, h = 52)
pred_arima_arnn = forecast(fit_arima, 52)$mean + pred_res_ARNN$mean
forecast::accuracy(pred_arima_arnn, con_tst)
smape(con_tst, pred_arima_arnn)*100
mase(con_tst, pred_arima_arnn)

### Hybrid ARIMA-WARIMA

fit_res_wbf = WaveletFittingarma(fit_arima$residuals, Waveletlevels = floor(log(length(con_tr))), 
boundary = 'periodic', FastFlag = TRUE, MaxARParam = 5, MaxMAParam = 5, NForecast = 52)
pred_arima_wbf = forecast(fit_arima, 52)$mean + fit_res_wbf$Finalforecast
forecast::accuracy(pred_arima_wbf, con_tst)
smape(con_tst, pred_arima_wbf)*100
mase(con_tst, pred_arima_wbf)

####################### Combination of Forecasts #####################

#### Ensemble ARIMA-ETS-Theta

avg_aef = hybridModel(con_tr, weights="equal",errorMethod = "MASE", models = "aef")
pred_aef = forecast::forecast(avg_aef, h = 52)
forecast::accuracy(pred_aef$mean, con_tst)
smape(con_tst, pred_aef$mean)*100
mase(con_tst, pred_aef$mean)

### Ensemble ETS-Theta-NNAR

avg_efn = hybridModel(con_tr, weights="equal",errorMethod = "MASE", models = "efn")
pred_efn = forecast::forecast(avg_efn, h = 52)
forecast::accuracy(pred_efn$mean, con_tst)
smape(con_tst, pred_efn$mean)*100
mase(con_tst, pred_efn$mean)

### Ensemble ARFIMA-Theta-ETS

pred_hybrid = 1/3*(forecast_ARFIMA$mean + fit_theta$mean 
                        + forecast_ets$mean)
forecast::accuracy(pred_hybrid, con_tst)
smape(con_tst, pred_hybrid)*100
mase(con_tst, pred_hybrid)

