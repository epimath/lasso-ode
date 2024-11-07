
## Title: "SIWR to fit SIWR"
## Date: 10-25-2024
## Author: Jiale Tan

library(deSolve)
library(Matrix)
library(ggplot2)
library(reshape2)

#############  SIWR to SIWR

# 1. True model generator: SIWR
SIURode <- function(t, x, params){
  S = x[1]
  I = x[2]
  U = x[3]
  R = x[4]
  params <- abs(params)
  
  b_u = params[1]
  b_i = params[2]
  n_u = params[3]
  g = params[4]
  
  dS = -b_u*S*U-b_i*S*I
  dI = b_u*S*U + b_i*S*I - g*I
  dU = n_u*(I-U)
  dR = g*I
  
  list(c(dS, dI, dU, dR))
}

# 2. Generate data from the true model
Data_Generate <- function(par, initial, times){
  out <- ode(y = abs(initial), times = abs(times), func = SIURode, parms = abs(par))
  out <- as.data.frame(out)
  
  return(out)
}

# 3. Test model : SIWR
SIURode <- function(t, x, params){
  S = x[1]
  I = x[2]
  U = x[3]
  R = x[4]
  params <- abs(params)
  
  b_u = params[1]
  b_i = params[2]
  n_u = params[3]
  g = params[4]
  
  dS = -b_u*S*U-b_i*S*I
  dI = b_u*S*U + b_i*S*I - g*I
  dU = n_u*(I-U)
  dR = g*I
  
  list(c(dS, dI, dU, dR))
}

######### 4. Evaluation function
SIURML1=function(params,times,data,lamda,initial_value){
  params = abs(params)
  xcurr = ode(abs(initial_value), abs(times), SIURode, params, method='ode45')
  y = xcurr[,3]
  
  Lasso = sum((y - data)^2)  + lamda* sum(abs(params))
  
  return (Lasso)
}

# 5. Generate data from true model
b_u =  0.5
b_i = 0.25
n_u = 0.01
g = 0.25

S0 = 1- 1e-2
I0 = 1e-2/2
U0 = 1e-2/2
R0 = 0

times <- seq(0, 200, by = 1)
parameters <- c(b_u, b_i,n_u,g)
init       <- c(S0, I0 , U0, R0 )

DATA_full <- Data_Generate(parameters, init, times)
DATA <- DATA_full[,3]


# 6. Blocking CV to get MSE
Blocking_CV_1 <- function(Data_train,Data_test,lamda,initial_value,initial_par){
  # choose cutting point
  initial_par <- abs(initial_par)
  initial_value <- abs(initial_value)
  times = seq(1,length(Data_train),1)
  
  res = optim(par = initial_par,fn=SIURML1,times=times,data=Data_train, 
              lamda = lamda,initial_value = initial_value)
  
  parameters <- abs(res$par)
  #residule <- (res$value)/length(Data_train)
  #print(parameters)
  
  S0 = 1-Data_test[1]
  I0 = Data_test[1]/2
  U0 = Data_test[1]/2
  R0 = 0
  
  
  times <- seq(1,length(Data_test),1)
  
  init <- c(S0, I0 , U0, R0)
  
  out <- ode(y = abs(init), times = abs(times), func = SIURode, parms = parameters)
  out <- as.data.frame(out)
  
  MSE <- sum(abs(Data_test - out[,3])^2)/length(Data_test)
  
  
  return(list(MSE = MSE, parameters = parameters))
}

times <- seq(1, length(DATA), by = 1)
initial_value <- c(1- 1e-2, 1e-2/2 , 1e-2/2, 0 )
n <- length(DATA)


#################### Ramdonly split to test and training (Algorithm 1)
lamda = 0
c <- 1
mea_vector_SIWRsiwr<- c()
PARAMETERS_SIWRsiwr <- c()

while (c < 1000) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_SIWRsiwr <- c(mea_vector_SIWRsiwr,b$MSE)
  PARAMETERS_SIWRsiwr <- rbind(PARAMETERS_SIWRsiwr,b$parameters)
  
  c <- c+1
}


lamda = 0.1
c <- 1
mea_vector_1_SIWRsiwr <- c()
PARAMETERS_1_SIWRsiwr <- c()
while (c < 1000) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_1_SIWRsiwr <- c(mea_vector_1_SIWRsiwr,b$MSE)
  PARAMETERS_1_SIWRsiwr <- rbind(PARAMETERS_1_SIWRsiwr,b$parameters)
  
  c <- c+1
}


lamda = 0.2
c <- 1
mea_vector_2_SIWRsiwr <- c()
PARAMETERS_2_SIWRsiwr <- c()
while (c < 1000) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_2_SIWRsiwr <- c(mea_vector_2_SIWRsiwr,b$MSE)
  PARAMETERS_2_SIWRsiwr <- rbind(PARAMETERS_2_SIWRsiwr,b$parameters)
  
  c <- c+1
}



lamda = 0.3
c <- 1
mea_vector_3_SIWRsiwr <- c()
PARAMETERS_3_SIWRsiwr <- c()
while (c < 1000) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_3_SIWRsiwr <- c(mea_vector_3_SIWRsiwr,b$MSE)
  PARAMETERS_3_SIWRsiwr <- rbind(PARAMETERS_3_SIWRsiwr,b$parameters)
  
  c <- c+1
}


lamda = 0.4
c <- 1
mea_vector_4_SIWRsiwr <- c()
PARAMETERS_4_SIWRsiwr <- c()
while (c < 1000) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_4_SIWRsiwr <- c(mea_vector_4_SIWRsiwr,b$MSE)
  PARAMETERS_4_SIWRsiwr <- rbind(PARAMETERS_4_SIWRsiwr,b$parameters)
  
  c <- c+1
}


lamda = 0.5
c <- 1
mea_vector_5_SIWRsiwr <- c()
PARAMETERS_5_SIWRsiwr <- c()
while (c < 1000) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_5_SIWRsiwr <- c(mea_vector_5_SIWRsiwr,b$MSE)
  PARAMETERS_5_SIWRsiwr <- rbind(PARAMETERS_5_SIWRsiwr,b$parameters)
  
  c <- c+1
}




############################### Split around peak area (Algorithm 2)
times <- seq(1, length(DATA), by = 1)
initial_value <- c(1- 1e-2, 1e-2/2 , 1e-2/2, 0 )
n <- length(DATA)
tem <- which.max(DATA)


lamda = 0
c <- 1
peak_mea_vector_SIWRsiwr <- c()
peak_PARAMETERS_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_SIWRsiwr <- c(peak_mea_vector_SIWRsiwr,b$MSE)
  peak_PARAMETERS_SIWRsiwr <- rbind(peak_PARAMETERS_SIWRsiwr,b$parameters)
  
  c <- c+1
}


lamda = 0.1
c <- 1
peak_mea_vector_1_SIWRsiwr <- c()
peak_PARAMETERS_1_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_1_SIWRsiwr <- c(peak_mea_vector_1_SIWRsiwr,b$MSE)
  peak_PARAMETERS_1_SIWRsiwr <- rbind(peak_PARAMETERS_1_SIWRsiwr,b$parameters)
  
  c <- c+1
}

lamda = 0.2
c <- 1
peak_mea_vector_2_SIWRsiwr <- c()
peak_PARAMETERS_2_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_2_SIWRsiwr <- c(peak_mea_vector_2_SIWRsiwr,b$MSE)
  peak_PARAMETERS_2_SIWRsiwr <- rbind(peak_PARAMETERS_2_SIWRsiwr,b$parameters)
  
  c <- c+1
}

lamda = 0.3
c <- 1
peak_mea_vector_3_SIWRsiwr <- c()
peak_PARAMETERS_3_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_3_SIWRsiwr <- c(peak_mea_vector_3_SIWRsiwr,b$MSE)
  peak_PARAMETERS_3_SIWRsiwr <- rbind(peak_PARAMETERS_3_SIWRsiwr,b$parameters)
  
  c <- c+1
}

lamda = 0.4
c <- 1
peak_mea_vector_4_SIWRsiwr <- c()
peak_PARAMETERS_4_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_4_SIWRsiwr <- c(peak_mea_vector_4_SIWRsiwr,b$MSE)
  peak_PARAMETERS_4_SIWRsiwr <- rbind(peak_PARAMETERS_4_SIWRsiwr,b$parameters)
  
  c <- c+1
}

lamda = 0.5
c <- 1
peak_mea_vector_5_SIWRsiwr <- c()
peak_PARAMETERS_5_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_5_SIWRsiwr <- c(peak_mea_vector_5_SIWRsiwr,b$MSE)
  peak_PARAMETERS_5_SIWRsiwr <- rbind(peak_PARAMETERS_5_SIWRsiwr,b$parameters)
  
  c <- c+1
}





########################## Training and Test use different data (algorithms 3)

## Function to get residule in order to generate testing data
Blocking_CV_2 <- function(Data_train,lamda,initial_value,initial_par){
  
  initial_par <- abs(initial_par)
  initial_value <- abs(initial_value)
  times = seq(1,length(Data_train),1)
  
  res = optim(par = initial_par,fn=SIRML1,times=times,data=Data_train, 
              lamda = lamda,initial_value = initial_value)
  
  parameters <- abs(res$par)
  
  
  S0 = 1-Data_train[1]
  I0 = Data_train[1]/2
  U0 = Data_train[1]/2
  R0 = 0
  
  
  init <- c(S0, I0 , U0, R0)
  
  out <- ode(y = abs(init), times = abs(times), func = SIURode, parms = parameters)
  out <- as.data.frame(out)
  
  
  Residule <- Data_train - out[,3]
  
  return(Residule)
}


initial_value <- c(1- 1e-2, 1e-2/2 , 1e-2/2, 0 )
initial_par <- c(0.5, 0.25, 0.01, 0.25)
lamda <- 0
Residule <- Blocking_CV_2(DATA, lamda,initial_value,initial_par)

mean_residule <- mean(Residule)
sd_residule <- sd(Residule)

## Start running 
lamda = 0
c <- 1
residule_mea_vector_SIWRsiwr <- c()
residule_PARAMETERS_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_SIWRsiwr <- c(residule_mea_vector_SIWRsiwr,b$MSE)
  residule_PARAMETERS_SIWRsiwr <- rbind(residule_PARAMETERS_SIWRsiwr,b$parameters)
  
  c <- c+1
}



lamda = 0.1
c <- 1
residule_mea_vector_1_SIWRsiwr <- c()
residule_PARAMETERS_1_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_1_SIWRsiwr <- c(residule_mea_vector_1_SIWRsiwr,b$MSE)
  residule_PARAMETERS_1_SIWRsiwr <- rbind(residule_PARAMETERS_1_SIWRsiwr,b$parameters)
  
  c <- c+1
}


lamda = 0.2
c <- 1
residule_mea_vector_2_SIWRsiwr <- c()
residule_PARAMETERS_2_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_2_SIWRsiwr <- c(residule_mea_vector_2_SIWRsiwr,b$MSE)
  residule_PARAMETERS_2_SIWRsiwr <- rbind(residule_PARAMETERS_2_SIWRsiwr,b$parameters)
  
  c <- c+1
}


lamda = 0.3
c <- 1
residule_mea_vector_3_SIWRsiwr <- c()
residule_PARAMETERS_3_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_3_SIWRsiwr <- c(residule_mea_vector_3_SIWRsiwr,b$MSE)
  residule_PARAMETERS_3_SIWRsiwr <- rbind(residule_PARAMETERS_3_SIWRsiwr,b$parameters)
  
  c <- c+1
}


lamda = 0.4
c <- 1
residule_mea_vector_4_SIWRsiwr <- c()
residule_PARAMETERS_4_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_4_SIWRsiwr <- c(residule_mea_vector_4_SIWRsiwr,b$MSE)
  residule_PARAMETERS_4_SIWRsiwr <- rbind(residule_PARAMETERS_4_SIWRsiwr,b$parameters)
  
  c <- c+1
}


lamda = 0.5
c <- 1
residule_mea_vector_5_SIWRsiwr <- c()
residule_PARAMETERS_5_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_5_SIWRsiwr <- c(residule_mea_vector_5_SIWRsiwr,b$MSE)
  residule_PARAMETERS_5_SIWRsiwr <- rbind(residule_PARAMETERS_5_SIWRsiwr,b$parameters)
  
  c <- c+1
}





##################  Training and Test use different data (algorithms 4)

# Generate data for testing using different initial conditions 

DATA_full_test <- Data_Generate(parameters, init, times)
DATA_test <- DATA_full_test[,3]


b_u =  0.5
b_i = 0.25
n_u = 0.01
g = 0.25

S0 = 1- 1e-2
I0 = ((1e-2)*2)/3
U0 = (1e-2)/3
R0 = 0

times <- seq(0, 200, by = 1)
parameters <- c(b_u, b_i,n_u,g)
init       <- c(S0, I0 , U0, R0 )

DATA_full_test <- Data_Generate(parameters, init, times)
DATA_test <- DATA_full_test[,3]




lamda = 0
c <- 1
diff_mea_vector_SIWRsiwr <- c()
diff_PARAMETERS_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_SIWRsiwr <- c(diff_mea_vector_SIWRsiwr,b$MSE)
  diff_PARAMETERS_SIWRsiwr <- rbind(diff_PARAMETERS_SIWRsiwr,b$parameters)
  
  c <- c+1
}


lamda = 0.1
c <- 1
diff_mea_vector_1_SIWRsiwr <- c()
diff_PARAMETERS_1_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_1_SIWRsiwr <- c(diff_mea_vector_1_SIWRsiwr,b$MSE)
  diff_PARAMETERS_1_SIWRsiwr <- rbind(diff_PARAMETERS_1_SIWRsiwr,b$parameters)
  
  c <- c+1
}

lamda = 0.2
c <- 1
diff_mea_vector_2_SIWRsiwr <- c()
diff_PARAMETERS_2_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_2_SIWRsiwr <- c(diff_mea_vector_2_SIWRsiwr,b$MSE)
  diff_PARAMETERS_2_SIWRsiwr <- rbind(diff_PARAMETERS_2_SIWRsiwr,b$parameters)
  
  c <- c+1
}

lamda = 0.3
c <- 1
diff_mea_vector_3_SIWRsiwr <- c()
diff_PARAMETERS_3_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_3_SIWRsiwr <- c(diff_mea_vector_3_SIWRsiwr,b$MSE)
  diff_PARAMETERS_3_SIWRsiwr <- rbind(diff_PARAMETERS_3_SIWRsiwr,b$parameters)
  
  c <- c+1
}

lamda = 0.4
c <- 1
diff_mea_vector_4_SIWRsiwr <- c()
diff_PARAMETERS_4_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_4_SIWRsiwr <- c(diff_mea_vector_4_SIWRsiwr,b$MSE)
  diff_PARAMETERS_4_SIWRsiwr <- rbind(diff_PARAMETERS_4_SIWRsiwr,b$parameters)
  
  c <- c+1
}

lamda = 0.5
c <- 1
diff_mea_vector_5_SIWRsiwr <- c()
diff_PARAMETERS_5_SIWRsiwr <- c()
while (c < 500) {
  b_u <- rnorm(1, mean = 0.5, sd = 0.5*0.2)
  b_i <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  n_u <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.25, sd = 0.25*0.2)
  
  initial_par <- c('b_u '=b_u ,'b_i'=b_i, 'n_u'=n_u,'g'=g)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_5_SIWRsiwr <- c(diff_mea_vector_5_SIWRsiwr,b$MSE)
  diff_PARAMETERS_5_SIWRsiwr <- rbind(diff_PARAMETERS_5_SIWRsiwr,b$parameters)
  
  c <- c+1
}




######################  plots
data1 <- data.frame(
  beta1 = PARAMETERS_SIWRsiwr[,1],
  beta2 = PARAMETERS_SIWRsiwr[,2],
  group = "Without LASSO"
)

data2 <- data.frame(
  beta1 = PARAMETERS_3_SIWRsiwr[,1],
  beta2 = PARAMETERS_3_SIWRsiwr[,2],
  group = "With LASSO"
)

combined_data <- rbind(data1, data2)

extra_point <- data.frame(x = 0.5, y = 0.25)


ggplot(combined_data, aes(x = beta1, y = beta2, color = group)) +
  geom_point(alpha = 0.2) +
  labs(title = "b_u vs. b_i", x = "b_u", y = "b_i") +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 0.5)) +
  geom_point(data = extra_point, aes(x = x, y = y), color = "black", size = 3) +
  scale_color_manual(values = c("Without LASSO" = "blue", "With LASSO" = "red")) +
  theme_minimal()



######### MSE plot
vector1 <- mea_vector_SIWRsiwr
vector2 <- mea_vector_1_SIWRsiwr
vector3 <- mea_vector_2_SIWRsiwr
vector4 <- mea_vector_3_SIWRsiwr
vector5 <- mea_vector_4_SIWRsiwr
vector6 <- mea_vector_5_SIWRsiwr
df <- data.frame(`0` = vector1, `0.1` = vector2, `0.2` = vector3,`0.3` = vector4, `0.4` = vector5, `0.5` = vector6)

# Reshape the data to long format
df_melted <- melt(df, variable.name = "category", value.name = "value")

# Create the box plot
ggplot(df_melted, aes(x = as.numeric(category), y = value, group = category)) +
  geom_boxplot() +
  scale_x_continuous(breaks = c(1,2,3,4,5,6), labels = c("0", "0.1", "0.2","0.3","0.4","0.5")) + 
  ylim(0, 0.5) + 
  labs(title = "MSE in testing dataset (Algorithm 4)",
       x = "Lamda",
       y = "MSE") +
  theme_minimal()






#setwd("/Users/jialetan/Desktop/UM-class/Research paper/ML_Marisa/Final_Submission/Final_Results")
#load("PARAMETERS_3_SIWRsiwr.RData")
#load("PARAMETERS_SIWRsiwr.RData")

### plot model fitting 
Data_Generate <- function(par, initial, times){
  out <- ode(y = abs(initial), times = abs(times), func = SIURode, parms = abs(par))
  out <- as.data.frame(out)
  return(out)
}

DATA_PREDICTION_LASSO <- c()
init <- c(1- 1e-2, 1e-2/2 , 1e-2/2, 0 )

for (i in 1:dim(PARAMETERS_SIWRsiwr)[1]){
  
  parameters <- PARAMETERS_SIWRsiwr[i,]
  
  DATA_simulated <- Data_Generate(parameters, init, times)
  DATA_simulated <- DATA_simulated[,3]
  
  DATA_PREDICTION_LASSO <- rbind(DATA_PREDICTION_LASSO,DATA_simulated)
}


DATA_PREDICTION_LASSO <- as.matrix(DATA_PREDICTION_LASSO)
quantiles_low <- apply(DATA_PREDICTION_LASSO, 2, quantile, probs = 0.25)
quantiles_median <- apply(DATA_PREDICTION_LASSO, 2, quantile, probs = 0.5)
quantiles_high <- apply(DATA_PREDICTION_LASSO, 2, quantile, probs = 0.75)


data <- data.frame(
  x = times,
  y1 = DATA,
  y2 = quantiles_median,
  y2_lower = quantiles_low,
  y2_upper = quantiles_high
)


p <- ggplot(data, aes(x = x)) +
  geom_ribbon(aes(ymin = y2_lower, ymax = y2_upper, fill = "25%-75% quantile"), alpha = 0.5) +
  geom_line(aes(y = y2, color = "Prediction"), size = 0.5) +
  geom_point(aes(y = y1, color = "Data"), size = 0.5) +
  scale_color_manual(values = c("Data" = "blue", "Prediction" = "dark green")) +
  scale_fill_manual(values = c("25%-75% quantile" = "green")) +
  labs(color = " ", fill = " ") +
  labs(x = "Time (days)", y = "Infections", title = "SIWR to fit SIWR") +
  theme_minimal()

p + theme(
  plot.title = element_text(size = 15),          # Title font size
  axis.title.x = element_text(size = 12),        # X-axis title font size
  axis.title.y = element_text(size = 12),        # Y-axis title font size
  axis.text.x = element_text(size = 12),         # X-axis text (tick labels) font size
  axis.text.y = element_text(size = 12),         # Y-axis text (tick labels) font size
  legend.title = element_text(size = 12),        # Legend title font size
  legend.text = element_text(size = 10)          # Legend text font size
)








