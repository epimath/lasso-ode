
## Title: "SI2R to fit SIR"
## Date: 10-25-2024
## Author: Jiale Tan


library(deSolve)
library(Matrix)
library(ggplot2)
library(reshape2)


#############  SI2R to SIR

# 1. True model generator: SIR
SIRode <- function(t, x, params){
  
  params <- abs(params)
  
  S = x[1]
  I = x[2]
  R = x[3]
  
  b = params[1]
  g = params[2]
  
  dS = -b*S*I
  dI = b*S*I - g*I
  dR = g*I
  
  list(c(dS, dI, dR))
}


# 2. Generate data from the true model
Data_Generate <- function(par, initial, times){
  out <- ode(y = abs(initial), times = abs(times), func = SIRode, parms = abs(par))
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

  # Make initial conditions parameters 
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
beta =  1.4247
gamma = 0.14286

S0 = 1-1e-6
I0 = 1e-6
R0 = 0

times <- seq(0, 70, by = 1)
parameters <- c(beta, gamma)
init       <- c(S0, I0 , R0 )

DATA_full <- Data_Generate(parameters, init, times)
DATA <- DATA_full[,3]




# 6. Blocking CV to get MSE
Blocking_CV_1 <- function(Data_train,Data_test,lamda,initial_value,initial_par){
  
  initial_par <- abs(initial_par)
  initial_value <- abs(initial_value)
  times_train = seq(1,length(Data_train),1)
  
  res = optim(par = initial_par,fn=SIRML1,times=times_train,data=Data_train, 
              lamda = lamda,initial_value = initial_value)
  
  parameters <- abs(res$par)

  
  
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


####################################    Traditional Blocking CV: 2-fold randomly split

# Function to calculate the Euclidean distance between two points
distance <- function(p1, p2) {
  abs(p2-p1)
}

# The main sampling function
sample_three_points <- function(min_distance, lower_bound = 1, upper_bound = n-1) {
  repeat {
    # Sample three points 
    p <- sort(sample(lower_bound:upper_bound, 3 , replace=F))

    # Calculate distances between each pair of points
    d12 <- distance(lower_bound, p[1])
    d23 <- distance(p[1], p[2])
    d34 <- distance(p[2], p[3])
    d45 <- distance(p[3], upper_bound)
    
    # Check if all distances are greater than the specified minimum distance
    if (d12 > min_distance && d23 > min_distance && d34 > min_distance && d45 > min_distance) {
      return(p)
    }
  }
}



############ Start run the model 
lamda = 0
c <- 1
mea_vector_SInRsir_BCV <- c()
PARAMETERS_SInRsir_BCV <- c()

while (c < 500) {
  # initial parameters
  beta_1 <- rnorm(1, mean = 1.4247, sd = 1.4247*0.2)
  gamma_1 <- rnorm(1, mean = 0.14286, sd = 0.14286*0.2)
  beta_2 <- rnorm(1, mean = 0, sd = 1.4247*0.2)
  gamma_2 <- rnorm(1, mean = 0, sd = 0.14286*0.2)
  
  # random cuting off points
  Points <- sample_three_points(3)
  DATA_train_1 <- DATA[1:Points[1]]
  DATA_test_1 <-  DATA[(Points[1]+1):Points[2]]
  DATA_train_2 <- DATA[(Points[2]+1):Points[3]]
  DATA_test_2 <-  DATA[(Points[3]+1):n]
  
  
  initial_par_1 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_1[1]/2)
  initial_par_2 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_2[1]/2)
  

  b_1 <- Blocking_CV_1(DATA_train_1 , DATA_test_1, lamda,initial_value,initial_par_1)
  b_2 <- Blocking_CV_1(DATA_train_2 , DATA_test_2, lamda,initial_value,initial_par_2)
  
  q <- data.frame(b_1$parameters,b_2$parameters)
  b_parameters <- rowMeans(q)
  b_MSE <- mean(c(b_1$MSE,b_2$MSE))
  
  mea_vector_SInRsir_BCV <- c(mea_vector_SInRsir_BCV,b_MSE)
  PARAMETERS_SInRsir_BCV <- rbind(PARAMETERS_SInRsir_BCV,b_parameters)
  
  
  c <- c+1
}


c <- 1
while (c < 500) {
  
  beta_1 <- rnorm(1, mean = 0, sd = 1.4247*0.2)
  gamma_1 <- rnorm(1, mean = 0, sd = 0.14286*0.2)
  beta_2 <- rnorm(1, mean = 1.4247, sd = 1.4247*0.2)
  gamma_2 <- rnorm(1, mean = 0.14286, sd = 0.14286*0.2)

  # random cuting off points
  Points <- sample_three_points(3)
  DATA_train_1 <- DATA[1:Points[1]]
  DATA_test_1 <-  DATA[(Points[1]+1):Points[2]]
  DATA_train_2 <- DATA[(Points[2]+1):Points[3]]
  DATA_test_2 <-  DATA[(Points[3]+1):n]
  
  
  initial_par_1 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_1[1]/2)
  initial_par_2 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_2[1]/2)
  
  b_1 <- Blocking_CV_1(DATA_train_1 , DATA_test_1, lamda,initial_value,initial_par_1)
  b_2 <- Blocking_CV_1(DATA_train_2 , DATA_test_2, lamda,initial_value,initial_par_2)
  
  q <- data.frame(b_1$parameters,b_2$parameters)
  b_parameters <- rowMeans(q)
  b_MSE <- mean(c(b_1$MSE,b_2$MSE))
  
  mea_vector_SInRsir_BCV <- c(mea_vector_SInRsir_BCV,b_MSE)
  PARAMETERS_SInRsir_BCV <- rbind(PARAMETERS_SInRsir_BCV,b_parameters)
  
  
  c <- c+1
}




lamda = 0.4
c <- 1
mea_vector_1_SInRsir_BCV <- c()
PARAMETERS_1_SInRsir_BCV <- c()

while (c < 500) {
  # initial parameters
  beta_1 <- rnorm(1, mean = 1.4247, sd = 1.4247*0.2)
  gamma_1 <- rnorm(1, mean = 0.14286, sd = 0.14286*0.2)
  beta_2 <- rnorm(1, mean = 0, sd = 1.4247*0.2)
  gamma_2 <- rnorm(1, mean = 0, sd = 0.14286*0.2)
  
  # random cuting off points
  Points <- sample_three_points(3)
  DATA_train_1 <- DATA[1:Points[1]]
  DATA_test_1 <-  DATA[(Points[1]+1):Points[2]]
  DATA_train_2 <- DATA[(Points[2]+1):Points[3]]
  DATA_test_2 <-  DATA[(Points[3]+1):n]
  
  
  initial_par_1 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_1[1]/2)
  initial_par_2 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_2[1]/2)
  
  
  b_1 <- Blocking_CV_1(DATA_train_1 , DATA_test_1, lamda,initial_value,initial_par_1)
  b_2 <- Blocking_CV_1(DATA_train_2 , DATA_test_2, lamda,initial_value,initial_par_2)
  
  q <- data.frame(b_1$parameters,b_2$parameters)
  b_parameters <- rowMeans(q)
  b_MSE <- mean(c(b_1$MSE,b_2$MSE))
  
  mea_vector_1_SInRsir_BCV <- c(mea_vector_1_SInRsir_BCV,b_MSE)
  PARAMETERS_1_SInRsir_BCV <- rbind(PARAMETERS_1_SInRsir_BCV,b_parameters)
  
  
  c <- c+1
}


c <- 1
while (c < 500) {
  
  beta_1 <- rnorm(1, mean = 0, sd = 1.4247*0.2)
  gamma_1 <- rnorm(1, mean = 0, sd = 0.14286*0.2)
  beta_2 <- rnorm(1, mean = 1.4247, sd = 1.4247*0.2)
  gamma_2 <- rnorm(1, mean = 0.14286, sd = 0.14286*0.2)
  
  # random cuting off points
  Points <- sample_three_points(3)
  DATA_train_1 <- DATA[1:Points[1]]
  DATA_test_1 <-  DATA[(Points[1]+1):Points[2]]
  DATA_train_2 <- DATA[(Points[2]+1):Points[3]]
  DATA_test_2 <-  DATA[(Points[3]+1):n]
  
  
  initial_par_1 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_1[1]/2)
  initial_par_2 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_2[1]/2)
  
  b_1 <- Blocking_CV_1(DATA_train_1 , DATA_test_1, lamda,initial_value,initial_par_1)
  b_2 <- Blocking_CV_1(DATA_train_2 , DATA_test_2, lamda,initial_value,initial_par_2)
  
  q <- data.frame(b_1$parameters,b_2$parameters)
  b_parameters <- rowMeans(q)
  b_MSE <- mean(c(b_1$MSE,b_2$MSE))
  
  mea_vector_1_SInRsir_BCV <- c(mea_vector_1_SInRsir_BCV,b_MSE)
  PARAMETERS_1_SInRsir_BCV <- rbind(PARAMETERS_1_SInRsir_BCV,b_parameters)
  
  
  c <- c+1
}





lamda = 0.8
c <- 1
mea_vector_2_SInRsir_BCV <- c()
PARAMETERS_2_SInRsir_BCV <- c()

while (c < 500) {
  # initial parameters
  beta_1 <- rnorm(1, mean = 1.4247, sd = 1.4247*0.2)
  gamma_1 <- rnorm(1, mean = 0.14286, sd = 0.14286*0.2)
  beta_2 <- rnorm(1, mean = 0, sd = 1.4247*0.2)
  gamma_2 <- rnorm(1, mean = 0, sd = 0.14286*0.2)
  
  # random cuting off points
  Points <- sample_three_points(3)
  DATA_train_1 <- DATA[1:Points[1]]
  DATA_test_1 <-  DATA[(Points[1]+1):Points[2]]
  DATA_train_2 <- DATA[(Points[2]+1):Points[3]]
  DATA_test_2 <-  DATA[(Points[3]+1):n]
  
  
  initial_par_1 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_1[1]/2)
  initial_par_2 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_2[1]/2)
  
  
  b_1 <- Blocking_CV_1(DATA_train_1 , DATA_test_1, lamda,initial_value,initial_par_1)
  b_2 <- Blocking_CV_1(DATA_train_2 , DATA_test_2, lamda,initial_value,initial_par_2)
  
  q <- data.frame(b_1$parameters,b_2$parameters)
  b_parameters <- rowMeans(q)
  b_MSE <- mean(c(b_1$MSE,b_2$MSE))
  
  mea_vector_2_SInRsir_BCV <- c(mea_vector_2_SInRsir_BCV,b_MSE)
  PARAMETERS_2_SInRsir_BCV <- rbind(PARAMETERS_2_SInRsir_BCV,b_parameters)
  
  
  c <- c+1
}


c <- 1
while (c < 500) {
  
  beta_1 <- rnorm(1, mean = 0, sd = 1.4247*0.2)
  gamma_1 <- rnorm(1, mean = 0, sd = 0.14286*0.2)
  beta_2 <- rnorm(1, mean = 1.4247, sd = 1.4247*0.2)
  gamma_2 <- rnorm(1, mean = 0.14286, sd = 0.14286*0.2)
  
  # random cuting off points
  Points <- sample_three_points(3)
  DATA_train_1 <- DATA[1:Points[1]]
  DATA_test_1 <-  DATA[(Points[1]+1):Points[2]]
  DATA_train_2 <- DATA[(Points[2]+1):Points[3]]
  DATA_test_2 <-  DATA[(Points[3]+1):n]
  
  
  initial_par_1 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_1[1]/2)
  initial_par_2 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_2[1]/2)
  
  b_1 <- Blocking_CV_1(DATA_train_1 , DATA_test_1, lamda,initial_value,initial_par_1)
  b_2 <- Blocking_CV_1(DATA_train_2 , DATA_test_2, lamda,initial_value,initial_par_2)
  
  q <- data.frame(b_1$parameters,b_2$parameters)
  b_parameters <- rowMeans(q)
  b_MSE <- mean(c(b_1$MSE,b_2$MSE))
  
  mea_vector_2_SInRsir_BCV <- c(mea_vector_2_SInRsir_BCV,b_MSE)
  PARAMETERS_2_SInRsir_BCV <- rbind(PARAMETERS_2_SInRsir_BCV,b_parameters)
  
  
  c <- c+1
}




lamda = 1.2
c <- 1
mea_vector_3_SInRsir_BCV <- c()
PARAMETERS_3_SInRsir_BCV <- c()

while (c < 500) {
  # initial parameters
  beta_1 <- rnorm(1, mean = 1.4247, sd = 1.4247*0.2)
  gamma_1 <- rnorm(1, mean = 0.14286, sd = 0.14286*0.2)
  beta_2 <- rnorm(1, mean = 0, sd = 1.4247*0.2)
  gamma_2 <- rnorm(1, mean = 0, sd = 0.14286*0.2)
  
  # random cuting off points
  Points <- sample_three_points(3)
  DATA_train_1 <- DATA[1:Points[1]]
  DATA_test_1 <-  DATA[(Points[1]+1):Points[2]]
  DATA_train_2 <- DATA[(Points[2]+1):Points[3]]
  DATA_test_2 <-  DATA[(Points[3]+1):n]
  
  
  initial_par_1 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_1[1]/2)
  initial_par_2 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_2[1]/2)
  
  
  b_1 <- Blocking_CV_1(DATA_train_1 , DATA_test_1, lamda,initial_value,initial_par_1)
  b_2 <- Blocking_CV_1(DATA_train_2 , DATA_test_2, lamda,initial_value,initial_par_2)
  
  q <- data.frame(b_1$parameters,b_2$parameters)
  b_parameters <- rowMeans(q)
  b_MSE <- mean(c(b_1$MSE,b_2$MSE))
  
  mea_vector_3_SInRsir_BCV <- c(mea_vector_3_SInRsir_BCV,b_MSE)
  PARAMETERS_3_SInRsir_BCV <- rbind(PARAMETERS_3_SInRsir_BCV,b_parameters)
  
  
  c <- c+1
}


c <- 1
while (c < 500) {
  
  beta_1 <- rnorm(1, mean = 0, sd = 1.4247*0.2)
  gamma_1 <- rnorm(1, mean = 0, sd = 0.14286*0.2)
  beta_2 <- rnorm(1, mean = 1.4247, sd = 1.4247*0.2)
  gamma_2 <- rnorm(1, mean = 0.14286, sd = 0.14286*0.2)
  
  # random cuting off points
  Points <- sample_three_points(3)
  DATA_train_1 <- DATA[1:Points[1]]
  DATA_test_1 <-  DATA[(Points[1]+1):Points[2]]
  DATA_train_2 <- DATA[(Points[2]+1):Points[3]]
  DATA_test_2 <-  DATA[(Points[3]+1):n]
  
  
  initial_par_1 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_1[1]/2)
  initial_par_2 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_2[1]/2)
  
  b_1 <- Blocking_CV_1(DATA_train_1 , DATA_test_1, lamda,initial_value,initial_par_1)
  b_2 <- Blocking_CV_1(DATA_train_2 , DATA_test_2, lamda,initial_value,initial_par_2)
  
  q <- data.frame(b_1$parameters,b_2$parameters)
  b_parameters <- rowMeans(q)
  b_MSE <- mean(c(b_1$MSE,b_2$MSE))
  
  mea_vector_3_SInRsir_BCV <- c(mea_vector_3_SInRsir_BCV,b_MSE)
  PARAMETERS_3_SInRsir_BCV <- rbind(PARAMETERS_3_SInRsir_BCV,b_parameters)
  
  
  c <- c+1
}



lamda = 1.6
c <- 1
mea_vector_4_SInRsir_BCV <- c()
PARAMETERS_4_SInRsir_BCV <- c()

while (c < 500) {
  # initial parameters
  beta_1 <- rnorm(1, mean = 1.4247, sd = 1.4247*0.2)
  gamma_1 <- rnorm(1, mean = 0.14286, sd = 0.14286*0.2)
  beta_2 <- rnorm(1, mean = 0, sd = 1.4247*0.2)
  gamma_2 <- rnorm(1, mean = 0, sd = 0.14286*0.2)
  
  # random cuting off points
  Points <- sample_three_points(3)
  DATA_train_1 <- DATA[1:Points[1]]
  DATA_test_1 <-  DATA[(Points[1]+1):Points[2]]
  DATA_train_2 <- DATA[(Points[2]+1):Points[3]]
  DATA_test_2 <-  DATA[(Points[3]+1):n]
  
  
  initial_par_1 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_1[1]/2)
  initial_par_2 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_2[1]/2)
  
  
  b_1 <- Blocking_CV_1(DATA_train_1 , DATA_test_1, lamda,initial_value,initial_par_1)
  b_2 <- Blocking_CV_1(DATA_train_2 , DATA_test_2, lamda,initial_value,initial_par_2)
  
  q <- data.frame(b_1$parameters,b_2$parameters)
  b_parameters <- rowMeans(q)
  b_MSE <- mean(c(b_1$MSE,b_2$MSE))
  
  mea_vector_4_SInRsir_BCV <- c(mea_vector_4_SInRsir_BCV,b_MSE)
  PARAMETERS_4_SInRsir_BCV <- rbind(PARAMETERS_4_SInRsir_BCV,b_parameters)
  
  
  c <- c+1
}


c <- 1
while (c < 500) {
  
  beta_1 <- rnorm(1, mean = 0, sd = 1.4247*0.2)
  gamma_1 <- rnorm(1, mean = 0, sd = 0.14286*0.2)
  beta_2 <- rnorm(1, mean = 1.4247, sd = 1.4247*0.2)
  gamma_2 <- rnorm(1, mean = 0.14286, sd = 0.14286*0.2)
  
  # random cuting off points
  Points <- sample_three_points(3)
  DATA_train_1 <- DATA[1:Points[1]]
  DATA_test_1 <-  DATA[(Points[1]+1):Points[2]]
  DATA_train_2 <- DATA[(Points[2]+1):Points[3]]
  DATA_test_2 <-  DATA[(Points[3]+1):n]
  
  
  initial_par_1 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_1[1]/2)
  initial_par_2 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_2[1]/2)
  
  b_1 <- Blocking_CV_1(DATA_train_1 , DATA_test_1, lamda,initial_value,initial_par_1)
  b_2 <- Blocking_CV_1(DATA_train_2 , DATA_test_2, lamda,initial_value,initial_par_2)
  
  q <- data.frame(b_1$parameters,b_2$parameters)
  b_parameters <- rowMeans(q)
  b_MSE <- mean(c(b_1$MSE,b_2$MSE))
  
  mea_vector_4_SInRsir_BCV <- c(mea_vector_4_SInRsir_BCV,b_MSE)
  PARAMETERS_4_SInRsir_BCV <- rbind(PARAMETERS_4_SInRsir_BCV,b_parameters)
  
  
  c <- c+1
}





lamda = 2
c <- 1
mea_vector_5_SInRsir_BCV <- c()
PARAMETERS_5_SInRsir_BCV <- c()

while (c < 500) {
  # initial parameters
  beta_1 <- rnorm(1, mean = 1.4247, sd = 1.4247*0.2)
  gamma_1 <- rnorm(1, mean = 0.14286, sd = 0.14286*0.2)
  beta_2 <- rnorm(1, mean = 0, sd = 1.4247*0.2)
  gamma_2 <- rnorm(1, mean = 0, sd = 0.14286*0.2)
  
  # random cuting off points
  Points <- sample_three_points(3)
  DATA_train_1 <- DATA[1:Points[1]]
  DATA_test_1 <-  DATA[(Points[1]+1):Points[2]]
  DATA_train_2 <- DATA[(Points[2]+1):Points[3]]
  DATA_test_2 <-  DATA[(Points[3]+1):n]
  
  
  initial_par_1 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_1[1]/2)
  initial_par_2 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_2[1]/2)
  
  
  b_1 <- Blocking_CV_1(DATA_train_1 , DATA_test_1, lamda,initial_value,initial_par_1)
  b_2 <- Blocking_CV_1(DATA_train_2 , DATA_test_2, lamda,initial_value,initial_par_2)
  
  q <- data.frame(b_1$parameters,b_2$parameters)
  b_parameters <- rowMeans(q)
  b_MSE <- mean(c(b_1$MSE,b_2$MSE))
  
  mea_vector_5_SInRsir_BCV <- c(mea_vector_5_SInRsir_BCV,b_MSE)
  PARAMETERS_5_SInRsir_BCV <- rbind(PARAMETERS_5_SInRsir_BCV,b_parameters)
  
  
  c <- c+1
}


c <- 1
while (c < 500) {
  
  beta_1 <- rnorm(1, mean = 0, sd = 1.4247*0.2)
  gamma_1 <- rnorm(1, mean = 0, sd = 0.14286*0.2)
  beta_2 <- rnorm(1, mean = 1.4247, sd = 1.4247*0.2)
  gamma_2 <- rnorm(1, mean = 0.14286, sd = 0.14286*0.2)
  
  # random cuting off points
  Points <- sample_three_points(3)
  DATA_train_1 <- DATA[1:Points[1]]
  DATA_test_1 <-  DATA[(Points[1]+1):Points[2]]
  DATA_train_2 <- DATA[(Points[2]+1):Points[3]]
  DATA_test_2 <-  DATA[(Points[3]+1):n]
  
  
  initial_par_1 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_1[1]/2)
  initial_par_2 <- c('beta1'=beta_1,'gamma1'=gamma_1, 'beta2'=beta_2,'gamma2'=gamma_2, 'I1' = DATA_train_2[1]/2)
  
  b_1 <- Blocking_CV_1(DATA_train_1 , DATA_test_1, lamda,initial_value,initial_par_1)
  b_2 <- Blocking_CV_1(DATA_train_2 , DATA_test_2, lamda,initial_value,initial_par_2)
  
  q <- data.frame(b_1$parameters,b_2$parameters)
  b_parameters <- rowMeans(q)
  b_MSE <- mean(c(b_1$MSE,b_2$MSE))
  
  mea_vector_5_SInRsir_BCV <- c(mea_vector_5_SInRsir_BCV,b_MSE)
  PARAMETERS_5_SInRsir_BCV <- rbind(PARAMETERS_5_SInRsir_BCV,b_parameters)
  
  
  c <- c+1
}



######### MSE plot
vector1 <- mea_vector_SInRsir_BCV
vector2 <- mea_vector_1_SInRsir_BCV
vector3 <- mea_vector_2_SInRsir_BCV
vector4 <- mea_vector_3_SInRsir_BCV
vector5 <- mea_vector_4_SInRsir_BCV
vector6 <- mea_vector_5_SInRsir_BCV
df <- data.frame(`0` = vector1, `0.4` = vector2, `0.8` = vector3,`1.2` = vector4, `1.6` = vector5, `2` = vector6)

# Reshape the data to long format
df_melted <- melt(df, variable.name = "category", value.name = "value")

# Create the box plot
ggplot(df_melted, aes(x = as.numeric(category), y = value, group = category)) +
  geom_boxplot() +
  scale_x_continuous(breaks = c(1,2,3,4,5,6), labels = c("0", "0.4", "0.8","1.2","1.6","2")) + 
  ylim(0 ,0.5) + 
  labs(title = "MSE in testing dataset (Traditional Blocking CV Algorithm)",
       x = "Lamda",
       y = "MSE") +
  theme_minimal()






############## parameters plots
data1 <- data.frame(
  beta1 = PARAMETERS_SInRsir_BCV[,2],
  beta2 = PARAMETERS_SInRsir_BCV[,4],
  group = "Conventional CV"
)


combined_data <- data1

extra_point <- data.frame(x = 1.4247, y = 0)
extra_point1 <- data.frame(x = 0, y = 1.4247)

ggplot(combined_data, aes(x = beta1, y = beta2, color = group)) +
  geom_point(alpha = 0.15) +
  labs(title = "beta1 vs. beta2", x = "beta1", y = "beta2") +
  scale_x_continuous(limits = c(0, 4)) +
  scale_y_continuous(limits = c(0, 4)) +
  geom_point(data = extra_point, aes(x = x, y = y), color = "black", size = 2) +
  geom_point(data = extra_point1, aes(x = x, y = y), color = "black", size = 2) +
  theme_minimal()


