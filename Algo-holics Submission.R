library(tidyverse)
library(caret)
library(lubridate)

raw_dat <- read.csv(file="SPY Historical Data-3.csv", header=TRUE)

date = mdy(raw_dat$Date)
open = raw_dat$Open
high = raw_dat$High
low = raw_dat$Low
close = raw_dat$Close
volume = as.numeric(gsub('.{1}$', '',raw_dat$Volume))
growth = abs((close-open)/open)


L <- length(date[which(date=="2017-01-03"):which(date=="2020-09-30")])+1
S <- length(date[1:which(date=="2016-12-29")])


# parameters

# number of days ahead
n <- 4

# Short term EMA day ranges
N1 <- 14 

# Long term EMA day range
N2 <- 35

# How many days behind we fit the model
N3 <- 500

# Sell cut-off
c1 <- 0.45

# Buy cut-off
c2 <- 0.52

# Buy or sell exponentiator
h <- 45

# Set-up portfolio tracker, values are as per end of the day, with first entry relating to 2016-12-29
port <- data.frame(cash = rep(NA,L),units=rep(NA,L),carry=rep(NA,L))
port$cash[1] <- 1000000; port$units[1] <- 0; port$carry[1] <- 0;

# Set-up vector to store predictive scores
ps <- rep(NA,L)


for (i in 2:L) {
  
  # find the indices of the N3 days of biggest loss or gain since 2011
  indx <- tail(order(growth[1:(i+S)]), N3)
  indx <- indx[which(indx>N3)]
  len <- length(indx)
  
  # study the indicators on these days and 
  y <- ifelse(close[n+indx]-open[indx]>0,1,0)
  x1  <- OBV(close[indx],volume[indx])
  x2 <- x3 <- x4 <- rep(NA,len)
  for (j in 1:len) {
    x2[j] <- EMA(close[(indx[j]+i-N1):(indx[j]+i)],N1) - EMA(close[(indx[j]+i-N2):(indx[j]+i)],N2)  
    x3[j] <- RSI(close[(indx[j]+i-7):(indx[j]+i)],open[(indx[j]+i-7):(indx[j]+i)])
  }
  
  
  Mod <- glm(y ~ x1 + x2 + x3 + x1*x3, 
             family = binomial("logit"),na.action = na.omit) 
  
  newdat <- data.frame(x1 = RSI(close[(i+S-14):(i+S)],open[(i+S-14):(i+S)]),
                       x2 = OBV(close[(i+S-14):(i+S)],volume[(i+S-14):(i+S)])[14],
                       x3 = (EMA(close[(i+S-N1):(i+S)],N1) - EMA(close[(i+S-N2):(i+S)],N2)))                         
  
  # Compute our predictive score  
  ps[i] <- predict.glm(Mod, newdat, type = "response")

  # buy or sell units based on the score
  if (ps[i] < c1) {
    # Predictive score below cut-off, so we sell $trans
    trans <- h^(20*(c1-ps[i]))
    # Sell n units, based on opening price for the day, if we have n units in portfolio
    quan <- trans/open[i+S]
    if (port$units[i-1] > quan) {
      port$units[i] <- port$units[i-1] - quan
      port$cash[i] <- port$cash[i-1] + trans
    } else { 
      # If we do not have n units, sell all we can
      port$units[i] <- 0
      port$cash[i] <- port$cash[i-1] + port$units[i-1]*open[i+S]
    }
  } else if (ps[i] > c2) {
    # Predictive score above cut-off, so we buy $trans
    trans <- h^(20*(ps[i]-c2))
    quan <- trans/open[i+S]
    # Buy n units, based on opening price for the day, if we have enough cash
    if (port$cash[i-1] > trans) {
      port$units[i] <- port$units[i-1] + quan
      port$cash[i] <- port$cash[i-1] - trans
    } else { 
      # If we do not have $trans cash, buy as much as possible
      port$units[i] <- port$units[i-1] + port$cash[i-1]/open[i+S]
      port$cash[i] <- 0
    }
  } else {
    # Neither Buy nor sell
    port$units[i] <- port$units[i-1]
    port$cash[i] <- port$cash[i-1]
  }
  # calculate value at the end of the day
  port$carry[i] <- port$units[i]*close[i+S]
}  

# Visualise performance
v <- port$carry + port$cash
plot(date[(S+1):(L+S)],v[1:L]/v[1],type="l",xlab="Date",ylab="Growth in Value",col="red")
lines(date[(S+1):(L+S)],close[(S+1):(L+S)]/close[S+1])
legend(date[(S+1)],2,legend=c("Algo-holic's auto-trader", "SPDR S&P 500 ETF Trust"),
       col=c("red", "black"), lty=c(1,1), cex=0.8)

# Functions 

# RSI
RSI <- function(close, open) {
  mvmt <- close - open
  ave_up <- mean(mvmt[which(mvmt>0)]) 
  ave_down <- -mean(mvmt[which(mvmt<0)])
  res <- 100-100/(1+ave_up/ave_down)
  return(res)
}

# OBV
OBV <- function(close, volume) {
  N_OBV <- length(close)
  res <- rep(NA,N_OBV)
  res[1]<- volume[1]
  for (i in 2:N_OBV) {
    if (close[i] > close[i-1]) {
      res[i] <- res[i-1] + volume[i]
    } else if(close[i] < close[i-1]) {
      res[i] <- res[i-1] - volume[i]
    } else {res[i] <- res[i-1]}
  }
  return(res) 
}

# EMA
EMA <- function(close, N) {
  res <- rep(NA, N)
  k <- 2/(N+1)
  res[1] <- close[1]
  for (i in 2:N) {
    res[i]<- close[i]*k+res[i-1]*(1-k)
  }
  return(res[N])
}






