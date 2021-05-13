WaveletFitting <- function(ts,Wvlevels,bndry,FFlag)
{
  mraout <- wavelets::modwt(ts, filter='haar', n.levels=Wvlevels,boundary=bndry, fast=FFlag)
  WaveletSeries <- cbind(do.call(cbind,mraout@W),mraout@V[[Wvlevels]])
  return(list(WaveletSeries=WaveletSeries,WVSeries=mraout))
}

WaveletFittingnar<- function(ts,Waveletlevels,boundary,FastFlag,MaxARParam,NForecast)
  
{
  WS <- WaveletFitting(ts=ts,Wvlevels=Waveletlevels,bndry=boundary,FFlag=FastFlag)$WaveletSeries
  AllWaveletForecast <- NULL;AllWaveletPrediction <- NULL
  
  for(WVLevel in 1:ncol(WS))
  {
    ts <- NULL
    ts <- WS[,WVLevel]
    # WaveletARMAFit <- forecast::auto.arima(x=as.ts(ts), d=NA, D=NA, max.p=MaxARParam, max.q=MaxMAParam,stationary=FALSE,
    #                                        seasonal=FALSE,ic=c("aic"), allowdrift=FALSE, allowmean=TRUE,stepwise = TRUE)
    # WaveletARIMAPredict <- WaveletARMAFit$fitted
    # WaveletARIMAForecast <- forecast::forecast(WaveletARMAFit,h=NForecast)
    # AllWaveletPrediction <- cbind(AllWaveletPrediction,WaveletARIMAPredict)
    # AllWaveletForecast <- cbind(AllWaveletForecast,as.matrix(WaveletARIMAForecast$mean))
    
    WaveletNARFit <- forecast::nnetar(y=as.ts(ts), xreg = xreg_tr, p = MaxARParam, repeats = 500)
    WaveletNARPredict <- WaveletNARFit$fitted
    WaveletNARForecast <- forecast::forecast(WaveletNARFit, xreg = xreg_tr[1:52], h=NForecast)
    AllWaveletPrediction <- cbind(AllWaveletPrediction,WaveletNARPredict)
    AllWaveletForecast <- cbind(AllWaveletForecast,as.matrix(WaveletNARForecast$mean))
    
    
  }
  Finalforecast <- rowSums(AllWaveletForecast,na.rm = T)
  FinalPrediction <- rowSums(AllWaveletPrediction,na.rm = T)
  return(list(Finalforecast=Finalforecast,FinalPrediction=FinalPrediction))
}