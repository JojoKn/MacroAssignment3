###Assignment 3###

###loading the required packages

library(BVAR)
library(vars)
library(fredr)
library(tseries)

###########################################
###Estimation and Identification of VARs###

###1###

###Open the API

fredr_set_key("3f87daa0298229fc49f2d28e0752dc58")
codes<-c("GDPC1",
         "GDPDEF",
         "FEDFUNDS")

###Setting the base variables

M <- length(codes)                   # how many variables (M)
t0 <- as.Date("1955-01-01")           # start date
t1 <- as.Date("2020-12-31")           # end date
t <- seq.Date(t0,t1,by="1 quarter")  # frequency, setting the increment of the sequence
bigT <- length(t)                       # length of time series (T)

###Generate an empty matrix for the time series to be put in

xraw<-matrix(NA, bigT, M)
colnames(xraw)<-codes
rownames(xraw)<-as.character(t)

###Programming a function to import the time series and make it part of the matrix

for (ii in 1:M){
  temp<-fredr(
  series_id<-codes[ii],
  observation_start=as.Date("1955-01-01"),
  observation_end=as.Date("2020-12-31"),
  frequency="q",
  aggregation_method="avg"
  )
  xraw[(t%in%temp$date),ii]<-(temp$value)
}
###Transform the matrix into a time series
xrawtseries<-xraw<-ts(xraw, start=c(1955,1), end=c(2020,4), frequency=4)

###Plot the Macro Time Series
plot.ts(xrawtseries, main="Macro Time Series")

###2###

#Real GDP
rgdpgr<-diff(log(xrawtseries[,1]), lag=4)*100
plot.ts(rgdpgr)
abline(h=0, col="green")

#GDP Deflator
gdpdefgr<-diff(log(xrawtseries[,2]), lag=4)*100
plot.ts(gdpdefgr)
abline(h=0, col="green")

###3###

summary(rgdpgr)
str(rgdpgr)
sd(rgdpgr)
adf.test(rgdpgr, alternative=c("stationary"))
#Explosiveness can be rejected for this time series, as p-value is smaller than 0.01. Furthermore,
#from eyeballing the time series, it seems to fluctuate around zero, with an only slightly decreasing variance.

summary(gdpdefgr)
str(gdpdefgr)
sd(gdpdefgr)
adf.test(gdpdefgr, alternative=c("stationary"))
#Explosiveness cannot be rejected for this time series. This also fits the plot of the time series, as the 
#variance in this time series shows strong outbursts, especially during the 70s and 80s.

#Therefore, take second differences for this time series:
gdpdefgr2<-diff(log(xrawtseries[,2]), lag=4, differences=2)*100
plot.ts(gdpdefgr2)
abline(h=0, col="green")
summary(gdpdefgr2)
str(gdpdefgr2)
sd(gdpdefgr2)
adf.test(gdpdefgr2, alternative=c("stationary"))
#Eyeballing the plot of the time series leads to the conclusion that it is stationary, as it is varying
#with almost constant variance around mean (almost) zero. The adf.test command can 
#reject the explosiveness of the process.

###4###

#Problem with the error terms in "normal" VARs: they are mutually correlated
#This does not allow for a clear economic interpretation. What we want to identify are structural shocks:
  #Orthogonal shocks (i.e. they are not mutually correlated)
  #with economic meaning
#To identify these shocks, we need the structural representation of the VAR:

#For this to be working, we need to find the matrix B0, which multiplies the whole VAR to obtain the
#structural representation. However, this leads to the problem that the model contains more parameters
#than we can estimate from the model, implying that we need to pose (Mx(M-1))/2 restrictions.
#The matrix B0 can be obtained using the Cholesky decomposition, which results in a lower triangular matrix
#that will therefore contain (Mx(M-1))/2 zeros above the main diagonal, which implies that there 
#are sufficient restrictions to estimate the model.

#Recursive Identification therefore means that the VAR will first be transformed using B0 to make 
#the error terms uncorrelated with each other before estimating the model. For this to work, we need to pin
#down B0, which is what we call identification. Recursive refers to the structure of B0 being a lower
#triangular matrix. The restrictions that are imposed by B0 not only have to include the minimum amount 
#of restrictions for the model to be identfied, but also need to have economically meaningful.

###5###

#How do you decide on the number of lags?

#What is the unit of measurement on each of the y-axes

#What is the dynamic impulse response of each of the variables to a monetary policy shock?

#Is this in accordance with what we would expect from theory?


###6###


##########################
###Bonus: Bayesian VARs###

###1###