
## Title: "SI2R to fit SI2R"
## Date: 10-25-2024
## Author: Jiale Tan


library(deSolve)
library(Matrix)
library(ggplot2)
library(reshape2)

#############  SI2R to SIR

# 1. True model generator: SI2R
SIR1ode <- function(t, x, params){
  
  params <- abs(params)
  S = x[1]
  I1 = x[2]
  I2 = x[3]
  R = x[4]
  
  b1 = params[1]
  g1 = params[2]
  b2 = params[3]
  g2 = params[4]
  
  dS = -b1*S*I1-b2*S*I2
  dI1 = b1*S*I1 - g1*I1
  dI2 = b2*S*I2 - g2*I2
  dR = g1*I1 + g2*I2
  
  list(c(dS, dI1, dI2, dR))
}

# 2. Generate data from the true model
Data_Generate <- function(par, initial, times){
  out <- ode(y = abs(initial), times = abs(times), func = SIR1ode, parms = abs(par))
  out <- as.data.frame(out)
  
  return(out)
}

# 3. Test model : SI2R
SIR1ode <- function(t, x, params){
  
  params <- abs(params)
  S = x[1]
  I1 = x[2]
  I2 = x[3]
  R = x[4]
  
  b1 = params[1]
  g1 = params[2]
  b2 = params[3]
  g2 = params[4]
  
  dS = -b1*S*I1-b2*S*I2
  dI1 = b1*S*I1 - g1*I1
  dI2 = b2*S*I2 - g2*I2
  dR = g1*I1 + g2*I2
  
  list(c(dS, dI1, dI2, dR))
}


# 4. Evaluation Function
SIRML1=function(params,times,data,lamda,initial_value){
  params = abs(params)
  initial_value[2] = max(params[5],0)
  initial_value[3] = max((data[1] - params[5]), 0)
  initial_value[1] = min((1 - initial_value[2] - initial_value[3]), 1)
  # Simulate model
  xcurr = ode(abs(initial_value), abs(times), SIR1ode,params, method='ode45')
  
  # Measurement equation
  y = xcurr[,3] + xcurr[,4] 
  
  # Lasso objective function
  
  Lasso = sum((y - data)^2) + lamda* sum(abs(params))
  
  return (Lasso)
}


# 5. Generate data from true model
beta_1 =  1.4247
gamma_1 = 0.14286
beta_2 =  1.4247
gamma_2 = 0.14286

S0 = 1-1e-6
I10 = 1e-6/2
I20 = 1e-6/2 
R0 = 0

times <- seq(0, 70, by = 1)
parameters <- c(beta_1, gamma_1,beta_2,gamma_2, I10)
init       <- c(S0, I10 , I20, R0 )

DATA_full <- Data_Generate(parameters, init, times)
DATA <- DATA_full[,3]+DATA_full[,4]



# 6. Blocking CV to get MSE
Blocking_CV_1 <- function(Data_train,Data_test,lamda,initial_value,initial_par){
  # choose cutting point
  initial_par <- abs(initial_par)
  initial_value <- abs(initial_value)
  times_train = seq(1,length(Data_train),1)
  
  res = optim(par = initial_par,fn=SIRML1,times=times_train,data=Data_train, 
              lamda = lamda,initial_value = initial_value)
  
  parameters <- abs(res$par)
  #residule <- (res$value)/length(Data_train)
  #print(parameters)
  
  #Option 1 (less good?)
  #S0 = 1-Data_test[1]
  #I_10 = Data_test[1]*parameters[5]/(par[5] + par[6])
  #I_20 = Data_test[1]*parameters[6]/(par[5] + par[6])
  #R_0 = 0
  
  
  # option 2
  # run the model with the estimated parameters & Initial conditions, then take the final values for S, I1, I2, R and use as intiial
  I10_train = max(parameters[5],0)
  I20_train = max((Data_train[1] - parameters[5]), 0)
  S0_train = min((1-I10_train-I20_train),1)
  R0_train = 0
  init_train <- c(S0_train, I10_train, I20_train, R0_train)
  
  out_train <- ode(y = abs(init_train), times = abs(times_train), func = SIR1ode, parms = parameters)
  
  
  ## Testing part
  
  times_test <- seq(1,length(Data_test),1)
  
  S0 = out_train[length(times_train),2]
  I_10 = out_train[length(times_train),3]
  I_20 = out_train[length(times_train),4]
  R_0 = out_train[length(times_train),5]
    
  init <- c(S0, I_10 , I_20, R_0)
  
  out <- ode(y = abs(init), times = abs(times_test), func = SIR1ode, parms = parameters)
  out <- as.data.frame(out)
  
  MSE <- sum(abs(Data_test - out[,3]-out[,4])^2)/length(Data_test)
  
  
  return(list(MSE = MSE, parameters = parameters))
}


times <- seq(1, length(DATA), by = 1)
initial_value <- c(1-1e-6, (1e-6)/2 , (1e-6)/2, 0 )
n <- length(DATA)



#################### Ramdonly split to test and training (Algorithm 1)
lamda = 0
c <- 1
mea_vector_SInRsinr <- c()
PARAMETERS_SInRsinr <- c()
while (c < 500) {
  ## Uniformly 
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_SInRsinr <- c(mea_vector_SInRsinr,b$MSE)
  PARAMETERS_SInRsinr <- rbind(PARAMETERS_SInRsinr,b$parameters)
  
  
  c <- c+1
}


lamda = 0.4
c <- 1
mea_vector_1_SInRsinr <- c()
PARAMETERS_1_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_1_SInRsinr <- c(mea_vector_1_SInRsinr,b$MSE)
  PARAMETERS_1_SInRsinr <- rbind(PARAMETERS_1_SInRsinr,b$parameters)
  
  
  c <- c+1
}


lamda = 0.8
c <- 1
mea_vector_2_SInRsinr <- c()
PARAMETERS_2_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_2_SInRsinr <- c(mea_vector_2_SInRsinr,b$MSE)
  PARAMETERS_2_SInRsinr <- rbind(PARAMETERS_2_SInRsinr,b$parameters)
  
  
  c <- c+1
}

lamda = 1.2
c <- 1
mea_vector_3_SInRsinr <- c()
PARAMETERS_3_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_3_SInRsinr <- c(mea_vector_3_SInRsinr,b$MSE)
  PARAMETERS_3_SInRsinr <- rbind(PARAMETERS_3_SInRsinr,b$parameters)
  
  
  c <- c+1
}

lamda = 1.6
c <- 1
mea_vector_4_SInRsinr <- c()
PARAMETERS_4_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_4_SInRsinr <- c(mea_vector_4_SInRsinr,b$MSE)
  PARAMETERS_4_SInRsinr <- rbind(PARAMETERS_4_SInRsinr,b$parameters)
  
  
  c <- c+1
}

lamda = 2
c <- 1
mea_vector_5_SInRsinr <- c()
PARAMETERS_5_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_5_SInRsinr <- c(mea_vector_5_SInRsinr,b$MSE)
  PARAMETERS_5_SInRsinr <- rbind(PARAMETERS_5_SInRsinr,b$parameters)
  
  
  c <- c+1
}





############################### Split around peak area (Algorithm 2)
times <- seq(1, length(DATA), by = 1)
initial_value <- c(1-1e-6, (1e-6)/2 , (1e-6)/2, 0 )
n <- length(DATA)
tem <- which.max(DATA)


lamda = 0
c <- 1
peak_mea_vector_SInRsinr <- c()
peak_PARAMETERS_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_SInRsinr <- c(peak_mea_vector_SInRsinr,b$MSE)
  peak_PARAMETERS_SInRsinr <- rbind(peak_PARAMETERS_SInRsinr,b$parameters)
  
  c <- c+1
}


lamda = 0.4
c <- 1
peak_mea_vector_1_SInRsinr <- c()
peak_PARAMETERS_1_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_1_SInRsinr <- c(peak_mea_vector_1_SInRsinr,b$MSE)
  peak_PARAMETERS_1_SInRsinr <- rbind(peak_PARAMETERS_1_SInRsinr,b$parameters)
  
  c <- c+1
}

lamda = 0.8
c <- 1
peak_mea_vector_2_SInRsinr <- c()
peak_PARAMETERS_2_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_2_SInRsinr <- c(peak_mea_vector_2_SInRsinr,b$MSE)
  peak_PARAMETERS_2_SInRsinr <- rbind(peak_PARAMETERS_2_SInRsinr,b$parameters)
  
  c <- c+1
}

lamda = 1.2
c <- 1
peak_mea_vector_3_SInRsinr <- c()
peak_PARAMETERS_3_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_3_SInRsinr <- c(peak_mea_vector_3_SInRsinr,b$MSE)
  peak_PARAMETERS_3_SInRsinr <- rbind(peak_PARAMETERS_3_SInRsinr,b$parameters)
  
  c <- c+1
}

lamda = 1.6
c <- 1
peak_mea_vector_4_SInRsinr <- c()
peak_PARAMETERS_4_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_4_SInRsinr <- c(peak_mea_vector_4_SInRsinr,b$MSE)
  peak_PARAMETERS_4_SInRsinr <- rbind(peak_PARAMETERS_4_SInRsinr,b$parameters)
  
  c <- c+1
}

lamda = 2
c <- 1
peak_mea_vector_5_SInRsinr <- c()
peak_PARAMETERS_5_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  start <- 1
  end <- sample((tem-1):(tem+6), 1 , replace=F)
  
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  peak_mea_vector_5_SInRsinr <- c(peak_mea_vector_5_SInRsinr,b$MSE)
  peak_PARAMETERS_5_SInRsinr <- rbind(peak_PARAMETERS_5_SInRsinr,b$parameters)
  
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
  
  I_10 = max(parameters[5],0)
  I_20 = max((Data_train[1] - parameters[5]),0)
  S0 =  min((1 - I_10 - I_20),1)
  R_0 = 0
  init <- c(S0, I_10 , I_20, R_0)

  
  out <- ode(y = abs(init), times = abs(times), func = SIR1ode, parms = parameters)
  out <- as.data.frame(out)
  
  
  Residule <- Data_train - (out[,3] + out[,4])
  
  return(Residule)
}

initial_value <- c(1-1e-6, (1e-6)/2 , (1e-6)/2, 0 )
initial_par <- c(1.3, 0.12, 1.3, 0.12,  (1e-6)/2 )
lamda <- 0
Residule <- Blocking_CV_2(DATA, lamda,initial_value,initial_par)

mean_residule <- mean(Residule)
sd_residule <- sd(Residule)


lamda = 0
c <- 1
residule_mea_vector_SInRsinr <- c()
residule_PARAMETERS_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_SInRsinr <- c(residule_mea_vector_SInRsinr,b$MSE)
  residule_PARAMETERS_SInRsinr <- rbind(residule_PARAMETERS_SInRsinr,b$parameters)
  c <- c+1
}



lamda = 0.4
c <- 1
residule_mea_vector_1_SInRsinr <- c()
residule_PARAMETERS_1_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_1_SInRsinr <- c(residule_mea_vector_1_SInRsinr,b$MSE)
  residule_PARAMETERS_1_SInRsinr <- rbind(residule_PARAMETERS_1_SInRsinr,b$parameters)
  c <- c+1
}


lamda = 0.8
c <- 1
residule_mea_vector_2_SInRsinr <- c()
residule_PARAMETERS_2_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_2_SInRsinr <- c(residule_mea_vector_2_SInRsinr,b$MSE)
  residule_PARAMETERS_2_SInRsinr <- rbind(residule_PARAMETERS_2_SInRsinr,b$parameters)
  c <- c+1
}


lamda = 1.2
c <- 1
residule_mea_vector_3_SInRsinr <- c()
residule_PARAMETERS_3_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_3_SInRsinr <- c(residule_mea_vector_3_SInRsinr,b$MSE)
  residule_PARAMETERS_3_SInRsinr <- rbind(residule_PARAMETERS_3_SInRsinr,b$parameters)
  c <- c+1
}


lamda = 1.6
c <- 1
residule_mea_vector_4_SInRsinr <- c()
residule_PARAMETERS_4_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_4_SInRsinr <- c(residule_mea_vector_4_SInRsinr,b$MSE)
  residule_PARAMETERS_4_SInRsinr <- rbind(residule_PARAMETERS_4_SInRsinr,b$parameters)
  c <- c+1
}


lamda = 2
c <- 1
residule_mea_vector_5_SInRsinr <- c()
residule_PARAMETERS_5_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA+rnorm(length(DATA), mean = mean_residule, sd = sd_residule))  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  residule_mea_vector_5_SInRsinr <- c(residule_mea_vector_5_SInRsinr,b$MSE)
  residule_PARAMETERS_5_SInRsinr <- rbind(residule_PARAMETERS_5_SInRsinr,b$parameters)
  c <- c+1
}




##################  Training and Test use different data (algorithms 4)

#Generate data from true model using different initial conditions
beta_1 =  1.4247
gamma_1 = 0.14286
beta_2 =  1.4247
gamma_2 = 0.14286

S0 = 1-1e-6
I10 = 1e-6/3
I20 = (1e-6)*2/3
R0 = 0

times <- seq(0, 70, by = 1)
parameters <- c(beta_1, gamma_1,beta_2,gamma_2,I10)
init       <- c(S0, I10 , I20, R0 )

DATA_full_test <- Data_Generate(parameters, init, times)
DATA_test <- DATA_full_test[,3]+DATA_full_test[,4]



### Start running
lamda = 0
c <- 1
diff_mea_vector_SInRsinr <- c()
diff_PARAMETERS_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_SInRsinr <- c(diff_mea_vector_SInRsinr,b$MSE)
  diff_PARAMETERS_SInRsinr <- rbind(diff_PARAMETERS_SInRsinr,b$parameters)
  c <- c+1
}


lamda = 0.4
c <- 1
diff_mea_vector_1_SInRsinr <- c()
diff_PARAMETERS_1_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_1_SInRsinr <- c(diff_mea_vector_1_SInRsinr,b$MSE)
  diff_PARAMETERS_1_SInRsinr <- rbind(diff_PARAMETERS_1_SInRsinr,b$parameters)
  c <- c+1
}


lamda = 0.8
c <- 1
diff_mea_vector_2_SInRsinr <- c()
diff_PARAMETERS_2_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_2_SInRsinr <- c(diff_mea_vector_2_SInRsinr,b$MSE)
  diff_PARAMETERS_2_SInRsinr <- rbind(diff_PARAMETERS_2_SInRsinr,b$parameters)
  c <- c+1
}


lamda = 1.2
c <- 1
diff_mea_vector_3_SInRsinr <- c()
diff_PARAMETERS_3_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_3_SInRsinr <- c(diff_mea_vector_3_SInRsinr,b$MSE)
  diff_PARAMETERS_3_SInRsinr <- rbind(diff_PARAMETERS_3_SInRsinr,b$parameters)
  c <- c+1
}


lamda = 1.6
c <- 1
diff_mea_vector_4_SInRsinr <- c()
diff_PARAMETERS_4_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_4_SInRsinr <- c(diff_mea_vector_4_SInRsinr,b$MSE)
  diff_PARAMETERS_4_SInRsinr <- rbind(diff_PARAMETERS_4_SInRsinr,b$parameters)
  c <- c+1
}


lamda = 2
c <- 1
diff_mea_vector_5_SInRsinr <- c()
diff_PARAMETERS_5_SInRsinr <- c()
while (c < 500) {
  beta_1 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_1 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  beta_2 <- runif(1, min = 1.4247*0.8, max = 1.4247*1.2)
  gamma_2 <- runif(1, min = 0.14286*0.8, max = 0.14286*1.2)
  I1 <-  (1e-6)/2
  
  
  initial_par <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = I1)
  
  DATA_train <- DATA
  DATA_test <-  abs(DATA_test)  
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  diff_mea_vector_5_SInRsinr <- c(diff_mea_vector_5_SInRsinr,b$MSE)
  diff_PARAMETERS_5_SInRsinr <- rbind(diff_PARAMETERS_5_SInRsinr,b$parameters)
  c <- c+1
}




######### MSE plot
vector1 <- mea_vector_SInRsinr
vector2 <- mea_vector_1_SInRsinr
vector3 <- mea_vector_2_SInRsinr
vector4 <- mea_vector_3_SInRsinr
vector5 <- mea_vector_4_SInRsinr
vector6 <- mea_vector_5_SInRsinr
df <- data.frame(`0` = vector1, `0.4` = vector2, `0.8` = vector3,`1.2` = vector4, `1.6` = vector5, `2` = vector6)

# Reshape the data to long format
df_melted <- melt(df, variable.name = "category", value.name = "value")

# Create the box plot
ggplot(df_melted, aes(x = as.numeric(category), y = value, group = category)) +
  geom_boxplot() +
  scale_x_continuous(breaks = c(1,2,3,4,5,6), labels = c("0", "0.4", "0.8","1.2","1.6","2")) + 
  ylim(0 ,0.02) + 
  labs(title = "MSE in testing dataset (Algorithm 4)",
       x = "Lamda",
       y = "MSE") +
  theme_minimal()





############## plots
data1 <- data.frame(
  beta1 = PARAMETERS_SInRsinr[,2],
  beta2 = PARAMETERS_SInRsinr[,4],
  group = "Without LASSO"
)

data2 <- data.frame(
  beta1 = PARAMETERS_1_SInRsinr[,2],
  beta2 = PARAMETERS_1_SInRsinr[,4],
  group = "With LASSO (Optimal)"
)

combined_data <- rbind(data1, data2)

extra_point <- data.frame(x = 0.14286, y = 0.14286)

ggplot(combined_data, aes(x = beta1, y = beta2, color = group)) +
  geom_point(alpha = 0.15) +
  labs(title = "gamma1 vs. gamma2", x = "gamma1", y = "gamma2") +
  scale_x_continuous(limits = c(0, 0.25)) +
  scale_y_continuous(limits = c(0, 0.25)) +
  geom_point(data = extra_point, aes(x = x, y = y), color = "black", size = 2) +
  scale_color_manual(values = c("Without LASSO" = "blue", "With LASSO (Optimal)" = "red")) +
  theme_minimal()





#setwd("/Users/jialetan/Desktop/UM-class/Research paper/ML_Marisa/Final_Submission/Final_Results")
#load("PARAMETERS_1_SInRsinr.RData")
#load("PARAMETERS_SInRsinr.RData")

### plot model fitting 
Data_Generate <- function(par, initial, times){
  out <- ode(y = abs(initial), times = abs(times), func = SIR1ode, parms = abs(par))
  out <- as.data.frame(out)
  return(out)
}

DATA_PREDICTION_LASSO <- c()

for (i in 1:dim(PARAMETERS_1_SInRsinr)[1]){
  
  parameters <- PARAMETERS_1_SInRsinr[i,]
  
  I10 = max(parameters[5],0)
  I20 = max((DATA[1] - parameters[5]), 0)
  S0 =  min((1-I10-I20),1)
  R0 = 0
  init <- c(S0, I10, I20, R0)
  
  DATA_simulated <- Data_Generate(parameters, init, times)
  DATA_simulated <- DATA_simulated[,3]+DATA_simulated[,4]
  
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
  labs(x = "Time (days)", y = "Infections", title = "SI2R to fit SI2R") +
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




