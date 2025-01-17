###Assignment 3###

###loading the required packages

library(BVAR)
library(vars)
library(fredr)
library(tseries)
library(bvartools)

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

par(mfrow=c(2,1))

#Real GDP
rgdpgr<-diff(log(xrawtseries[,1]), lag=4)*100
plot.ts(rgdpgr)
abline(h=0, col="green")

#GDP Deflator
gdpdefgr<-diff(log(xrawtseries[,2]), lag=4)*100
plot.ts(gdpdefgr)
abline(h=0, col="green")

dev.off()

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

# Write down a VAR(p) process in its general form.
# y_t = c + \sum_{j=1}^{p}A_{ij}y_{t-j} + \epsilon,~~~ \epsilon \sim N(0,\Sigma)

# Transform the VAR(p) process to a SVAR(p) process with the according transformation.
# B_0y_t = B_0c + \sum_{j=1}^{p}B_0A_{ij}y_{t-j} + e_t,~~~ e_t \sim N(0,I)~

#Problem with the error terms in "normal" VARs: they are mutually correlated.
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
#of restrictions for the model to be identified, but also need to have economically meaningful.

# The recursive identification process is sensitive to the ordering of variables. In the Cholesky decomposition, the restrictions imply that the first variable is not contemporaneously correlated with any other variables, the second variable is only contemporaneously correlated with the first variable, and so on.
# In the context of Monetary Policy identification, a typical implementation is ordering the variables as output, inflation, then interest rate. This implies that output is only dependent its own structural shocks contemporaneously, inflation depends on itself and output contemporaneously and interest rate is endogenously dependent on output, inflation and interest shocks contemporaneously. This ordering has to be defended using economic intuition and cannot be verified empirically in the model.


###5###

#Create a new matrix which contains the growth rates of real GDP and the GDP deflator

y<-cbind(rgdpgr, gdpdefgr2, xraw[,3])
#create a subset as the first eight observations cannot be used due to differencing

yest<-y[9:nrow(y),]
colnames(yest)<-c("Outputgap", "Inflation", "FFR")
yest

t_<-t[9:length(t)]
t_

rownames(yest)<-as.character(t_)

#How do you decide on the number of lags?
#This can be done by using the Information criteria. In R, a function is implemented that returns the
#calculated selection criteria. Select minimum BIC from there (results in smallest number of lags).

lagselection<-VARselect(yest)
lag<-lagselection$selection[3]

#Estimation procedure of the VAR(2)

var<-VAR(yest, p=lag, type=c("const"))
summary(var)

#Storing the coefficients of the VAR(2)
A <- as.matrix(Bcoef(var))

#Notational Issues: Looking for the A representation of the SVAR; therefore 
#define a 3x3 identity matrix for the A matrix with NA elements at the spots that 
#one wants to have estimates for; set bmat=NULL in the estimation
#Compared to the notation from task 4, A corresponds to matrix B0.

amat<-diag(3)
amat[1,1]=amat[2,2]=amat[3,3]=NA
amat[2,1]=amat[3,1]=amat[3,2]=NA
amat


svar<-SVAR(var, Amat=amat, Bmat=NULL)

svar$type

#calculate the impulse response functions

impulseOutputgap<-vars::irf(svar, response="Outputgap", impulse="FFR", n.ahead=50)
impulseInflation<-vars::irf(svar, response="Inflation", impulse="FFR", n.ahead=50)
impulseFFR<-vars::irf(svar, response="FFR", impulse="FFR", n.ahead=50)

plot(impulseOutputgap)
plot(impulseInflation)
plot(impulseFFR)

###What is the unit of measurement on each of the y-axes

#FFR: interest rate in percentage points
#Output gap: GDP growth rate in percent
#Inflation rate: Change in inflation rate year on year in percent (?)


###What is the dynamic impulse response of each of the variables to a monetary policy shock?

#The shock is a contractionary monetary policy shock, i.e. the interest rate is increased
#This can be seen when looking at the impulse responses for the federal funds rate

#For the inflation, this would imply a decrease; however, this can only be observed in the
#long run; in the short run, inflation increases due to the monetary policy shock.

#The same holds true for the output gap (i.e. GDP growth), as this also increases in the
#short run before decreasing in the long run.

# Increase in interest rate is persistently high for 24 periods and returns to zero over time. There is an initial increase in the interest rate after the shock, peaking after 2 periods.
# Initial increase in GDP growth rate for approximately 5 periods, subsequently becoming negative, then tenting back towards zero over time.
# Similar to the GDP growth rate, the GDP deflector initially increases, peaking at period 3, then reduces to below zero. This stays persistently below zero for all 24 periods, tending back towards zero over time.

###Is this in accordance with what we would expect from theory?

# The responses are somewhat consistent with economic theory. Consider the basic New Keynesian model. With an contractionary monetary policy shock, interest rates stay persistently high rates over time, tending towards zero. This is consistent with our SVAR model.
# For GDP growth and the inflation responses, the basic New Keynesian model does not exhibit the initial hump-shaped increase that we observe in the IRFs, but the persistent negative GDP growth and inflation rates, tending back towards zero over time, is consistent with the theoretical model.

#No, for the output gap and the inflation it is not what we would expect, at least not in 
#the short run, as the impulse responses behave opposite to what we would expect.


###6###

#Instead of the output gap, one can look at the Industrial Production and its change
#over time. Industrial Production is a key figure in economics and should also react
#to monetary policy shocks, as a decrease in demand and a reduction in investment
#due to an increase in the interest rate should certainly affect industrial
#production as well.

Outputgap_Alt<-fredr(
  series_id="INDPRO",
  observation_start=as.Date("1955-01-01"),
  observation_end=as.Date("2020-12-31"),
  frequency="q",
  aggregation_method="avg"
)
Outputgap_Alt_<-ts(Outputgap_Alt$value, start=c(1955,1), end=c(2020,4), frequency=4)
plot.ts(Outputgap_Alt_)

Outputgap_Alt_gr<-diff(log(Outputgap_Alt_), lag=4)*100
plot(Outputgap_Alt_gr)
adf.test(Outputgap_Alt_gr, alternative=c("stationary"))
#Computed time series is stationary!
#Create an updated matrix with the exchanged data:
Outputgap_Alt_gr_<-Outputgap_Alt_gr[5:length(Outputgap_Alt_gr)]
yest1<-yest
yest1[,1]=Outputgap_Alt_gr_

#Selecting optimal lag length

lagselection1<-VARselect(yest1)
lag1<-lagselection1$selection[3]

#Estimation procedure of the VAR(2)

var1<-VAR(yest1, p=lag1, type=c("const")) 
summary(var1)   

#Storing the coefficients of the VAR(2)
A1 <- as.matrix(Bcoef(var1))

#Estimating the structural VAR:

svar1<-SVAR(var1, Amat=amat, Bmat=NULL)

#Impulse response functions:
library(vars)
impulseOutputgap1<-vars::irf(svar1, response="Outputgap", impulse="FFR", n.ahead=50)
impulseInflation1<-vars::irf(svar1, response="Inflation", impulse="FFR", n.ahead=50)
impulseFFR1<-vars::irf(svar1, response="FFR", impulse="FFR", n.ahead=50)

plot(impulseOutputgap1)
plot(impulseInflation1)
plot(impulseFFR1)
#Looks very similar to before

#Instead of the GDP Deflator as a measure for inflation, one can look at the CPI.
#Measures inflation not from the production side, but from the consumption side.
#While it can be expected that these inflation measures diverge, it seems to be a 
#reasonable assumption that both of the adequatly depict prices. 

Inflation_Alt<-GDPDEF_Alt<-fredr(
  series_id="CPIAUCSL",
  observation_start=as.Date("1955-01-01"),
  observation_end=as.Date("2020-12-31"),
  frequency="q",
  aggregation_method="avg"
)

Inflation_Alt_<-ts(Inflation_Alt$value, start=c(1955,1), end=c(2020,4), frequency=4)
plot.ts(Inflation_Alt_)

Inflation_Alt_gr<-diff(log(Inflation_Alt_), lag=4, diff=1)*100
plot(Inflation_Alt_gr)
adf.test(Inflation_Alt_gr, alternative=c("stationary"))
#Not stationary for p-value=0.05. Take second differences

Inflation_Alt_gr2<-diff(log(Inflation_Alt_), lag=4, diff=2)*100
plot(Inflation_Alt_gr2)
adf.test(Inflation_Alt_gr2, alternative=c("stationary")) #stationary

#Create an updated matrix with the exchanged data:

yest2<-yest
yest2[,2]=Inflation_Alt_gr2

#Selecting optimal lag length

lagselection2<-VARselect(yest2)
lag2<-lagselection2$selection[3]

#Estimation procedure of the VAR(6)

var2<-VAR(yest2, p=lag2, type=c("const")) 
summary(var2)   

#Storing the coefficients of the VAR(2)
A2 <- as.matrix(Bcoef(var2))

#Estimating the structural VAR:

svar2<-SVAR(var2, Amat=amat, Bmat=NULL)

#Impulse response functions:
library(vars)
impulseOutputgap2<-vars::irf(svar2, response="Outputgap", impulse="FFR", n.ahead=50)
impulseInflation2<-vars::irf(svar2, response="Inflation", impulse="FFR", n.ahead=50)
impulseFFR2<-vars::irf(svar2, response="FFR", impulse="FFR", n.ahead=50)

plot(impulseOutputgap2)
plot(impulseInflation2)
plot(impulseFFR2)
#Puzzle still there and seems to have gotten even stronger

#Instead of the Federal Funds Rate, one can look at 3-Months Treasury Bill Secondary
#Markets Rates. A rise in the interest rate should lead to a rise in the Treasury Bill
#Market rate, as investing becomes more attractive. 

FFR_Alt<-fredr(
  series_id="TB3MS",
  observation_start=as.Date("1955-01-01"),
  observation_end=as.Date("2020-12-31"),
  frequency="q",
  aggregation_method="avg"
)

FFR_Alt_<-ts(FFR_Alt$value, start=c(1955,1), end=c(2020,4), frequency=4)
plot(FFR_Alt_)
lines(xrawtseries[,3], col="green")
#Create an updated matrix with the exchanged data:

FFR_Alt__<-FFR_Alt_[9:length(FFR_Alt_)]
yest3<-yest
yest3[,3]=FFR_Alt__

#Selecting optimal lag length

lagselection3<-VARselect(yest3)
lag3<-lagselection3$selection[3]

#Estimation procedure of the VAR(6)

var3<-VAR(yest3, p=lag3, type=c("const")) 
summary(var3)   

#Storing the coefficients of the VAR(6)
A3 <- as.matrix(Bcoef(var3))

#Estimating the structural VAR:

svar3<-SVAR(var3, Amat=amat, Bmat=NULL)

#Impulse response functions:
library(vars)
impulseOutputgap3<-vars::irf(svar3, response="Outputgap", impulse="FFR", n.ahead=50)
impulseInflation3<-vars::irf(svar3, response="Inflation", impulse="FFR", n.ahead=50)
impulseFFR3<-vars::irf(svar3, response="FFR", impulse="FFR", n.ahead=50)

plot(impulseOutputgap3)
plot(impulseInflation3)
plot(impulseFFR3)


###Splitting the time series and only redo the analysis for the constrained sample:

#A suitable time point seems to be somewhere around the 1980s. During this time, the
#time series of the GDP deflator seems to be becoming stationary, which would enable 
#us to avoid taking second differences and solely work with first differences. This
#assumption will be tested in the following using data from 1985 Q1 onwards:

gdpdef<-xrawtseries[113:264,2]
gdpdefgr_<-diff(log(gdpdef), lag=4)*100
plot.ts(gdpdefgr_)
adf.test(gdpdefgr_, alternative=c("stationary"))
length(gdpdefgr_)
#ADF Test allows to reject explosiveness of the process; therefore, we can assume
#this time series to be stationary and do not need to take second differences. Allows
#for easier interpretation!

#Set up a new matrix for this subperiod: 1985 Q1 : 2020 Q4

yest4<-yest[109:256,]
yest4[,2]<-gdpdefgr_

#Estimation same as before:
lagselection4<-VARselect(yest4)
lag4<-lagselection4$selection[3]

#Estimation procedure of the VAR(2)

var4<-VAR(yest4, p=lag4, type=c("const")) 
summary(var4)   

#Storing the coefficients of the VAR(2)
A4 <- as.matrix(Bcoef(var4))

#Estimating the structural VAR:

svar4<-SVAR(var4, Amat=amat, Bmat=NULL)

#Impulse response functions:
library(vars)
impulseOutputgap4<-vars::irf(svar4, response="Outputgap", impulse="FFR", n.ahead=50)
impulseInflation4<-vars::irf(svar4, response="Inflation", impulse="FFR", n.ahead=50)
impulseFFR4<-vars::irf(svar4, response="FFR", impulse="FFR", n.ahead=50)

plot(impulseOutputgap4)
plot(impulseInflation4)
plot(impulseFFR4)

#Very poor performance of the whole system in terms of explaining the monetary
#policy shock; no negative impact on inflation, same holds true for the outputgap.

#Alternative 2: split just before the monetary policy shock

#Set up a new matrix for the subperiod 1957 Q1 : 2006 Q4
yest4a<-yest[9:200,]

#Estimation same as before:
lagselection4a<-VARselect(yest4a)
lag4a<-lagselection4a$selection[3]

#Estimation procedure of the VAR(2)

var4a<-VAR(yest4a, p=lag4a, type=c("const")) 
summary(var4a)   

#Storing the coefficients of the VAR(2)
A4a <- as.matrix(Bcoef(var4a))

#Estimating the structural VAR:

svar4a<-SVAR(var4a, Amat=amat, Bmat=NULL)

#Impulse response functions:
library(vars)
impulseOutputgap4a<-vars::irf(svar4a, response="Outputgap", impulse="FFR", n.ahead=50)
impulseInflation4a<-vars::irf(svar4a, response="Inflation", impulse="FFR", n.ahead=50)
impulseFFR4a<-vars::irf(svar4a, response="FFR", impulse="FFR", n.ahead=50)

plot(impulseOutputgap4a)
plot(impulseInflation4a)
plot(impulseFFR4a)


##########################
###Bonus: Bayesian VARs###

###1###
library(BVAR)
#For perfect reproducability, is is necessary to set the seed:
set.seed(123456789)

#Estimate the Bayesian VAR, lag selection from estimation of VAR/SVAR
x<-BVAR::bvar(yest, lags=lag)
print(x)
str(x)

#calculate the impulse response functions
setting<-bv_irf(horizon=24, identification = FALSE)
irf(x)<-BVAR::irf(x, setting, conf_bands=c(0.05, 0.1))
plot(BVAR::irf(x))

###Which prior did you use?
#Minnesota prior, but nothing adjusted, everything governed/optimised by the algorithm. Using this prior regularises the data. This helps to remove the underestimation of persistence common in frequentist approaches.

###What are the differences compared to the frequentist VAR?
#In estimation, we use the concept of belief updating; the prior is updated using the
#dataset to obtain the posterior. 
#Contentwise, the puzzle that was part of the frequentist VAR also appears in the 
#Impulse Response functions; even though there is a monetary policy shock which
#raises the interest rate, Inflation and Outputgap first increase before decreasing.
# Frequentist VAR (e.g. by estimating the VAR using OLS) has a small sample bias, where the persistence of parameter estimates is systematically underestimated. Frequentist VAR also faces the curse of dimensionality - the proliferation of parameters when new variables are added.
# Bayesian VAR approaches address these limitations through regularization and the choice of relevant priors, thus limiting the variance of parameter estimates.
# The interpretation of uncertainty bands is also different in both approaches - see below.

###How do you interpret the uncertainty bands?
#The uncertainty bands result from the estimation uncertainty introduced using the 
#Bayesian VAR. Between the bands, the credible set of impulse responses can be found.
# In Bayesian statistics, parameters are not considered to be fixed, but random variables. So, the uncertainty bands can be interpreted as the distribution of the parameter of interest. This reflects the fundamental epistemological uncertainty in the world.
# In frequentist statistics, we think of the parameter as a fixed and unknown value, and the uncertainty band is a set of values within which the 'true' value should lie with a given certainty.
