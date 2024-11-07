
## Title: "Asymptomatic model to fit Exponential model"
## Date: 10-25-2024
## Author: Jiale Tan

library(deSolve)
library(Matrix)
library(ggplot2)
library(reshape2)

##### 1. True model generator: Exponential model
model1_exp <- function(t,x,params) {
  params <- abs(params)
  
  mu=0              
  betaI=params[1]           
  betaW=params[2]           
  gamma=0.25                 
  xi=params[3]              
  alpha=params[4]           
  kinv = params[5]
  
  S=x[1]
  I=x[2]
  W=x[3]
  R=x[4]
  
  dS = mu-betaI*S*I-betaW*S*W-mu*S+alpha*R
  dI = betaI*S*I+betaW*S*W-gamma*I-mu*I
  dW = xi*(I-W)
  dR = gamma*I-mu*R-alpha*R
  
  list(c(dS, dI, dW, dR))
}

##### 2.Generate data from the true model
Data_Generate <- function(par, initial, times){
  out <- ode(y = abs(initial), times = abs(times), func = model1_exp, parms = abs(par))
  out <- as.data.frame(out)
  
  return(out)
}


## 3. Generate dataset
betaI = 0.269
betaW = 1.62
xi = 0.00614
alpha = 0.001314
kinv = 0.0000116


S0 = 0.9986892
I0 = 1- S0
W0 = 0
R0 = 0

times <- seq(0, 200, by = 1)
parameters <- c(betaI, betaW, xi, alpha, kinv)
init       <- c(S0, I0 , W0 , R0 )

DATA <- Data_Generate(parameters, init, times)
DATA <- DATA[,3]/kinv   # CONVERT TO ACTUAL NUMBERS
n <- length(DATA)


##### 4. Test model generator: Asymptomatic model
model3_asymp_restrict <- function(t,x,params){
  
  params <- abs(params)
  
  mu=0 
  beta_IS = params[1]
  beta_IA = params[2]
  
  beta_W = params[3]
  alpha_S = params[4]
  alpha_A = params[5]
  gamma = 0.25
  ksi = params[6]
  kinv = params[7]
  q = params[8]
  
  S = x[1]
  I_S = x[2]
  I_A = x[3]
  R_S = x[4]
  R_A = x[5]
  W = x[6]
  
  dS = mu - beta_IS*S*I_S - beta_IA*S*I_A - beta_W*S*W - mu*S + alpha_S*R_S + alpha_A*R_A
  dI_S = min(q,1)*beta_IS*S*I_S + min(q,1)*beta_IA*S*I_A + min(q,1)*beta_W*S*W - gamma*I_S - mu*I_S
  dI_A = max((1-q),0)*beta_IS*S*I_S + max((1-q),0)*beta_IA*S*I_A + max((1-q),0)*beta_W*S*W - gamma*I_A - mu*I_A
  dR_S = gamma*I_S - alpha_S*R_S - mu*R_S
  dR_A = gamma*I_A - alpha_A*R_A - mu*R_A
  dW = ksi*(I_A + I_S - W)
  
  list(c(dS, dI_S, dI_A, dR_S, dR_A, dW))
}



######### 5. Evaluation function
SIURML1=function(params,times,data,lamda,initial_value){
  params = abs(params)
  xcurr = ode(abs(initial_value), abs(times), model3_asymp_restrict, params, method='ode45')
  y = xcurr[,3]/params['kinv']
  
  #Lasso = sum((y - data)^2)  + lamda* sum(abs(params))
  # Partially punish parameters
  Lasso = sum((y - data)^2)  + lamda*(abs(params[1])+abs(params[2])+abs(params[4])+abs(params[5])+abs(params[8]))
  
  return (Lasso)
}

##########  6. Blocking CV to get MSE
Blocking_CV_1 <- function(Data_train,Data_test,lamda,initial_value,initial_par){
  # choose cutting point
  initial_par <- abs(initial_par)
  initial_value <- abs(initial_value)
  times = seq(1,length(Data_train),1)
  
  res = optim(par = initial_par,fn=SIURML1,times=times,data=Data_train, 
              lamda = lamda,initial_value = initial_value)
  
  parameters <- abs(res$par)
  #print(parameters)
  
  S0 = 1 - Data_test[1]*parameters['kinv']
  I_S0 = Data_test[1]*parameters['kinv']
  I_A0 = 0
  R_S0 = 0
  R_A0 = 0
  W0 = 0
  
  times <- seq(1, length(Data_test), by = 1)
  
  init <- c(S0, I_S0 , I_A0, R_S0, R_A0, W0)
  
  out <- ode(y = abs(init), times = abs(times), func = model3_asymp_restrict, parms = parameters)
  out <- as.data.frame(out)
  
  MSE <- sum(abs(Data_test - out[,3]/parameters['kinv']))/length(Data_test)
  
  
  return(list(MSE = MSE, parameters = parameters))
}

times <- seq(1, length(DATA), by = 1)
initial_value = c(0.9986892, 0.0013108 , 0 , 0 , 0, 0 )
n <- length(DATA)


############################### Ramdonly split to test and training (Algorithm 1)
lamda = 0
c <- 1
mea_vector_ASYMPexp_367 <- c()
PARAMETERS_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q) 
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_ASYMPexp_367 <- c(mea_vector_ASYMPexp_367,b$MSE)
  PARAMETERS_ASYMPexp_367 <- rbind(PARAMETERS_ASYMPexp_367,b$parameters)
  c <- c+1
}



lamda = (1e+06)*1
c <- 1
mea_vector_1_ASYMPexp_367 <- c()
PARAMETERS_1_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q) 
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <-  Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_1_ASYMPexp_367 <- c(mea_vector_1_ASYMPexp_367,b$MSE)
  PARAMETERS_1_ASYMPexp_367 <- rbind(PARAMETERS_1_ASYMPexp_367,b$parameters)
  c <- c+1
}




lamda = (1e+06)*5
c <- 1
mea_vector_2_ASYMPexp_367 <- c()
PARAMETERS_2_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q) 
  
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <-  Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_2_ASYMPexp_367 <- c(mea_vector_2_ASYMPexp_367,b$MSE)
  PARAMETERS_2_ASYMPexp_367 <- rbind(PARAMETERS_2_ASYMPexp_367,b$parameters)
  c <- c+1
}



lamda = (1e+06)*10
c <- 1
mea_vector_3_ASYMPexp_367 <- c()
PARAMETERS_3_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q) 
  
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <-  Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_3_ASYMPexp_367 <- c(mea_vector_3_ASYMPexp_367,b$MSE)
  PARAMETERS_3_ASYMPexp_367 <- rbind(PARAMETERS_3_ASYMPexp_367,b$parameters)
  c <- c+1
}



lamda = (1e+06)*15
c <- 1
mea_vector_4_ASYMPexp_367 <- c()
PARAMETERS_4_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q) 
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <-  Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_4_ASYMPexp_367 <- c(mea_vector_4_ASYMPexp_367,b$MSE)
  PARAMETERS_4_ASYMPexp_367 <- rbind(PARAMETERS_4_ASYMPexp_367,b$parameters)
  c <- c+1
}




lamda = (1e+06)*20
c <- 1
mea_vector_5_ASYMPexp_367 <- c()
PARAMETERS_5_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q) 
  
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <-  Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_5_ASYMPexp_367 <- c(mea_vector_5_ASYMPexp_367,b$MSE)
  PARAMETERS_5_ASYMPexp_367 <- rbind(PARAMETERS_5_ASYMPexp_367,b$parameters)
  c <- c+1
}



lamda = (1e+07)*5
c <- 1
mea_vector_6_ASYMPexp_367 <- c()
PARAMETERS_6_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q) 
  
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <-  Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_6_ASYMPexp_367 <- c(mea_vector_6_ASYMPexp_367,b$MSE)
  PARAMETERS_6_ASYMPexp_367 <- rbind(PARAMETERS_6_ASYMPexp_367,b$parameters)
  c <- c+1
}







############################### Split around peak area (Algorithms 2)
times <- seq(1, length(DATA), by = 1)
initial_value = c(0.9986892, 0.0013108 , 0 , 0 , 0, 0 )
n <- length(DATA)
tem <- which.max(DATA)


lamda = 0
c <- 1
peak_mea_vector_ASYMPexp_367 <- c()
peak_PARAMETERS_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_ASYMPexp_367 <- c(peak_mea_vector_ASYMPexp_367,b$MSE)
  peak_PARAMETERS_ASYMPexp_367 <- rbind(peak_PARAMETERS_ASYMPexp_367,b$parameters)
  
  c <- c+1
}



lamda = (1e+06)*1
c <- 1
peak_mea_vector_1_ASYMPexp_367 <- c()
peak_PARAMETERS_1_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_1_ASYMPexp_367 <- c(peak_mea_vector_1_ASYMPexp_367,b$MSE)
  peak_PARAMETERS_1_ASYMPexp_367 <- rbind(peak_PARAMETERS_1_ASYMPexp_367,b$parameters)
  
  c <- c+1
}


lamda = (1e+06)*5
c <- 1
peak_mea_vector_2_ASYMPexp_367 <- c()
peak_PARAMETERS_2_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_2_ASYMPexp_367 <- c(peak_mea_vector_2_ASYMPexp_367,b$MSE)
  peak_PARAMETERS_2_ASYMPexp_367 <- rbind(peak_PARAMETERS_2_ASYMPexp_367,b$parameters)
  
  c <- c+1
}


lamda = (1e+06)*10
c <- 1
peak_mea_vector_3_ASYMPexp_367 <- c()
peak_PARAMETERS_3_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_3_ASYMPexp_367 <- c(peak_mea_vector_3_ASYMPexp_367,b$MSE)
  peak_PARAMETERS_3_ASYMPexp_367 <- rbind(peak_PARAMETERS_3_ASYMPexp_367,b$parameters)
  
  c <- c+1
}


lamda = (1e+06)*15
c <- 1
peak_mea_vector_4_ASYMPexp_367 <- c()
peak_PARAMETERS_4_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_4_ASYMPexp_367 <- c(peak_mea_vector_4_ASYMPexp_367,b$MSE)
  peak_PARAMETERS_4_ASYMPexp_367 <- rbind(peak_PARAMETERS_4_ASYMPexp_367,b$parameters)
  
  c <- c+1
}


lamda = (1e+06)*20
c <- 1
peak_mea_vector_5_ASYMPexp_367 <- c()
peak_PARAMETERS_5_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_5_ASYMPexp_367 <- c(peak_mea_vector_5_ASYMPexp_367,b$MSE)
  peak_PARAMETERS_5_ASYMPexp_367 <- rbind(peak_PARAMETERS_5_ASYMPexp_367,b$parameters)
  
  c <- c+1
}


lamda = (1e+07)*5
c <- 1
peak_mea_vector_6_ASYMPexp_367 <- c()
peak_PARAMETERS_6_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_6_ASYMPexp_367 <- c(peak_mea_vector_6_ASYMPexp_367,b$MSE)
  peak_PARAMETERS_6_ASYMPexp_367 <- rbind(peak_PARAMETERS_6_ASYMPexp_367,b$parameters)
  
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
  
  
  S0 = 1 - Data_train[1]*parameters['kinv']
  I_S0 = Data_train[1]*parameters['kinv']
  I_A0 = 0
  R_S0 = 0
  R_A0 = 0
  W0 = 0
  
  
  init <- c(S0, I_S0 , I_A0, R_S0, R_A0, W0)
  
  
  out <- ode(y = abs(init), times = abs(times), func = model3_asymp_restrict, parms = parameters)
  out <- as.data.frame(out)
  
  
  Residule <- Data_train - out[,3]/parameters['kinv']
  
  return(Residule)
}



initial_value = c(0.9986892, 0.0013108 , 0 , 0 , 0, 0 )
initial_par <- c(0.26,0.05,1.62,0.0013,0.0005,0.00614,0.0000116,1)
lamda <- 0
Residule <- Blocking_CV_2(DATA, lamda,initial_value,initial_par)

mean_residule <- mean(Residule)
sd_residule <- sd(Residule)


lamda = 0
c <- 1
residule_mea_vector_ASYMPexp_367 <- c()
residule_PARAMETERS_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_ASYMPexp_367 <- c(residule_mea_vector_ASYMPexp_367,b$MSE)
  residule_PARAMETERS_ASYMPexp_367 <- rbind(residule_PARAMETERS_ASYMPexp_367,b$parameters)
  c <- c+1
}



lamda = (1e+06)*1
c <- 1
residule_mea_vector_1_ASYMPexp_367 <- c()
residule_PARAMETERS_1_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_1_ASYMPexp_367 <- c(residule_mea_vector_1_ASYMPexp_367,b$MSE)
  residule_PARAMETERS_1_ASYMPexp_367 <- rbind(residule_PARAMETERS_1_ASYMPexp_367,b$parameters)
  c <- c+1
}


lamda = (1e+06)*5
c <- 1
residule_mea_vector_2_ASYMPexp_367 <- c()
residule_PARAMETERS_2_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_2_ASYMPexp_367 <- c(residule_mea_vector_2_ASYMPexp_367,b$MSE)
  residule_PARAMETERS_2_ASYMPexp_367 <- rbind(residule_PARAMETERS_2_ASYMPexp_367,b$parameters)
  c <- c+1
}


lamda = (1e+06)*10
c <- 1
residule_mea_vector_3_ASYMPexp_367 <- c()
residule_PARAMETERS_3_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_3_ASYMPexp_367 <- c(residule_mea_vector_3_ASYMPexp_367,b$MSE)
  residule_PARAMETERS_3_ASYMPexp_367 <- rbind(residule_PARAMETERS_3_ASYMPexp_367,b$parameters)
  c <- c+1
}


lamda = (1e+06)*15
c <- 1
residule_mea_vector_4_ASYMPexp_367 <- c()
residule_PARAMETERS_4_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_4_ASYMPexp_367 <- c(residule_mea_vector_4_ASYMPexp_367,b$MSE)
  residule_PARAMETERS_4_ASYMPexp_367 <- rbind(residule_PARAMETERS_4_ASYMPexp_367,b$parameters)
  c <- c+1
}


lamda = (1e+06)*20
c <- 1
residule_mea_vector_5_ASYMPexp_367 <- c()
residule_PARAMETERS_5_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_5_ASYMPexp_367 <- c(residule_mea_vector_5_ASYMPexp_367,b$MSE)
  residule_PARAMETERS_5_ASYMPexp_367 <- rbind(residule_PARAMETERS_5_ASYMPexp_367,b$parameters)
  c <- c+1
}


lamda = (1e+07)*5
c <- 1
residule_mea_vector_6_ASYMPexp_367 <- c()
residule_PARAMETERS_6_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_6_ASYMPexp_367 <- c(residule_mea_vector_6_ASYMPexp_367,b$MSE)
  residule_PARAMETERS_6_ASYMPexp_367 <- rbind(residule_PARAMETERS_6_ASYMPexp_367,b$parameters)
  c <- c+1
}



##################  Training and Test use different data (algorithms 4)
# Generate data for testing using different initial conditions

## Generate dataset
betaI = 0.269
betaW = 1.62
xi = 0.00614
alpha = 0.001314
kinv = 0.0000116


S0 = 0.97
I0 = 1- S0
W0 = 0
R0 = 0

times <- seq(0, 200, by = 1)
parameters <- c(betaI, betaW, xi, alpha, kinv)
init       <- c(S0, I0 , W0 , R0 )

DATA_full_test <- Data_Generate(parameters, init, times)
DATA_test <- DATA_full_test[,3]/kinv   # CONVERT TO ACTUAL NUMBERS


# Start running
lamda = 0
c <- 1
diff_mea_vector_ASYMPexp_367 <- c()
diff_PARAMETERS_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_ASYMPexp_367 <- c(diff_mea_vector_ASYMPexp_367,b$MSE)
  diff_PARAMETERS_ASYMPexp_367 <- rbind(diff_PARAMETERS_ASYMPexp_367,b$parameters)
  c <- c+1
}



lamda = (1e+06)*1
c <- 1
diff_mea_vector_1_ASYMPexp_367 <- c()
diff_PARAMETERS_1_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_1_ASYMPexp_367 <- c(diff_mea_vector_1_ASYMPexp_367,b$MSE)
  diff_PARAMETERS_1_ASYMPexp_367 <- rbind(diff_PARAMETERS_1_ASYMPexp_367,b$parameters)
  c <- c+1
}



lamda = (1e+06)*5
c <- 1
diff_mea_vector_2_ASYMPexp_367 <- c()
diff_PARAMETERS_2_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_2_ASYMPexp_367 <- c(diff_mea_vector_2_ASYMPexp_367,b$MSE)
  diff_PARAMETERS_2_ASYMPexp_367 <- rbind(diff_PARAMETERS_2_ASYMPexp_367,b$parameters)
  c <- c+1
}



lamda = (1e+06)*10
c <- 1
diff_mea_vector_3_ASYMPexp_367 <- c()
diff_PARAMETERS_3_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_3_ASYMPexp_367 <- c(diff_mea_vector_3_ASYMPexp_367,b$MSE)
  diff_PARAMETERS_3_ASYMPexp_367 <- rbind(diff_PARAMETERS_3_ASYMPexp_367,b$parameters)
  c <- c+1
}



lamda = (1e+06)*15
c <- 1
diff_mea_vector_4_ASYMPexp_367 <- c()
diff_PARAMETERS_4_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_4_ASYMPexp_367 <- c(diff_mea_vector_4_ASYMPexp_367,b$MSE)
  diff_PARAMETERS_4_ASYMPexp_367 <- rbind(diff_PARAMETERS_4_ASYMPexp_367,b$parameters)
  c <- c+1
}



lamda = (1e+06)*20
c <- 1
diff_mea_vector_5_ASYMPexp_367 <- c()
diff_PARAMETERS_5_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_5_ASYMPexp_367 <- c(diff_mea_vector_5_ASYMPexp_367,b$MSE)
  diff_PARAMETERS_5_ASYMPexp_367 <- rbind(diff_PARAMETERS_5_ASYMPexp_367,b$parameters)
  c <- c+1
}



lamda = (1e+07)*5
c <- 1
diff_mea_vector_6_ASYMPexp_367 <- c()
diff_PARAMETERS_6_ASYMPexp_367 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  beta_IA <- rnorm(1, mean = 0.04, sd = 0.04*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  q <- 1
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv, 'q' = q)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_6_ASYMPexp_367 <- c(diff_mea_vector_6_ASYMPexp_367,b$MSE)
  diff_PARAMETERS_6_ASYMPexp_367 <- rbind(diff_PARAMETERS_6_ASYMPexp_367,b$parameters)
  c <- c+1
}





######################  plots
data1 <- data.frame(
  beta1 = PARAMETERS_ASYMPexp_367[,1],
  beta2 = PARAMETERS_ASYMPexp_367[,2],
  group = "Without LASSO"
)

data2 <- data.frame(
  beta1 = PARAMETERS_3_ASYMPexp_367[,1],
  beta2 = PARAMETERS_3_ASYMPexp_367[,2],
  group = "With LASSO"
)

combined_data <- rbind(data1, data2)

extra_point <- data.frame(x = 0.269, y = 0)


ggplot(combined_data, aes(x = beta1, y = beta2, color = group)) +
  geom_point(alpha = 0.2) +
  labs(title = "beta_IS vs. beta_IA", x = "beta_IS", y = "beta_IA") +
  scale_x_continuous(limits = c(0, 0.5)) +
  scale_y_continuous(limits = c(0, 0.2)) +
  geom_point(data = extra_point, aes(x = x, y = y), color = "black", size = 3) +
  scale_color_manual(values = c("Without LASSO" = "blue", "With LASSO" = "red")) +
  theme_minimal()


########### MSE plots
vector1 <- mea_vector_ASYMPexp_367
vector2 <- mea_vector_1_ASYMPexp_367
vector3 <- mea_vector_2_ASYMPexp_367
vector4 <- mea_vector_3_ASYMPexp_367
vector5 <- mea_vector_4_ASYMPexp_367
vector6 <- mea_vector_5_ASYMPexp_367
df <- data.frame(`0` = vector1, `(1e+06)*1` = vector2, `(1e+06)*5` = vector3,`(1e+06)*10` = vector4, `(1e+06)*15` = vector5, `(1e+06)*20` = vector6)

# Reshape the data to long format
df_melted <- melt(df, variable.name = "category", value.name = "value")

# Create the box plot
ggplot(df_melted, aes(x = as.numeric(category), y = value, group = category)) +
  geom_boxplot() +
  scale_x_continuous(breaks = c(1,2,3,4,5,6), labels = c("0", "(1e+06)*1", "(1e+06)*5","(1e+06)*10","(1e+06)*15","(1e+06)*20")) + 
  ylim(0, 4000) + 
  labs(title = "MSE in testing dataset (Algorithm 4)",
       x = "Lamda",
       y = "MSE") +
  theme_minimal()





#setwd("/Users/jialetan/Desktop/UM-class/Research paper/ML_Marisa/Final_Submission/Final_Results")
#load("PARAMETERS_3_ASYMPexp_367.RData")
#load("PARAMETERS_ASYMPexp_367.RData")

### plot model fitting 
Data_Generate <- function(par, initial, times){
  out <- ode(y = abs(initial), times = abs(times), func = model3_asymp_restrict, parms = abs(par))
  out <- as.data.frame(out)
  return(out)
}

DATA_PREDICTION_LASSO <- c()
init <- c(0.9986892, 0.0013108 , 0 , 0 , 0, 0 )

for (i in 1:dim(PARAMETERS_3_ASYMPexp_367)[1]){
  
  parameters <- PARAMETERS_3_ASYMPexp_367[i,]
  
  DATA_simulated <- Data_Generate(parameters, init, times)
  DATA_simulated <- DATA_simulated[,3]/parameters['kinv']
  
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
  labs(x = "Time (days)", y = "Infections", title = "Asymptomatic to fit Exponential") +
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






