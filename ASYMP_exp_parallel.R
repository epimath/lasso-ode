## Title: "Model selection : Asymptomatic model VS. Exponential model"
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

##### 2. Generate data from the true model
Data_Generate <- function(par, initial, times){
  out <- ode(y = abs(initial), times = abs(times), func = model1_exp, parms = abs(par))
  out <- as.data.frame(out)
  
  return(out)
}


##########  Example

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



##### 4. Test model generator: Asymptomatic model + Exponential model
model3_asymp_restrict <- function(t,x,params){
  
  ##  Asymptomatic model
  params <- abs(params)
  
  mu=0 
  beta_IS = params[1]
  beta_IA = params[2]
  beta_W = params[3]
  alpha_S = params[4]
  alpha_A = params[5]
  gamma = 0.25
  ksi = params[6]
  kinv_0 = params[7]
  q = params[8]
  
  S_M3 = x[1]
  I_S = x[2]
  I_A = x[3]
  R_S = x[4]
  R_A = x[5]
  W_M3 = x[6]
  
  dS_M3 = mu - beta_IS*S_M3*I_S - beta_IA*S_M3*I_A - beta_W*S_M3*W_M3 - mu*S_M3 + alpha_S*R_S + alpha_A*R_A
  dI_S = min(q,1)*beta_IS*S_M3*I_S + min(q,1)*beta_IA*S_M3*I_A + min(q,1)*beta_W*S_M3*W_M3 - gamma*I_S - mu*I_S
  dI_A = max((1-q),0)*beta_IS*S_M3*I_S + max((1-q),0)*beta_IA*S_M3*I_A + max((1-q),0)*beta_W*S_M3*W_M3 - gamma*I_A - mu*I_A
  dR_S = gamma*I_S - alpha_S*R_S - mu*R_S
  dR_A = gamma*I_A - alpha_A*R_A - mu*R_A
  dW_M3 = ksi*(I_A + I_S - W_M3)
  
  
  
  ##  exponential model
  betaI=params[9]           
  betaW=params[10]           
  xi=params[11]              
  alpha=params[12]    
  kinv = params[13]
  wei = params[14]
  
  S=x[7]
  I=x[8]
  W=x[9]
  R=x[10]
  
  dS = mu-betaI*S*I-betaW*S*W-mu*S+alpha*R
  dI = betaI*S*I+betaW*S*W-gamma*I-mu*I
  dW = xi*(I-W)
  dR = gamma*I-mu*R-alpha*R
  
  list(c(dS_M3, dI_S, dI_A, dR_S, dR_A, dW_M3, dS, dI, dW, dR ))
}



######### 5. Evaluation function
SIURML1=function(params,times,data,lamda,initial_value){
  params = abs(params)
  xcurr = ode(abs(initial_value), abs(times), model3_asymp_restrict, params, method='ode45')
  y = min(params['wei'],1)*(xcurr[,3]/params['kinv_0']) + (1-min(params['wei'],1))*(xcurr[,9]/params['kinv'])
  
  Lasso = sum((y - data)^2)  + lamda*(abs(params[1])+abs(params[2])+abs(params[3])+abs(params[4])+abs(params[5])+abs(params[9])+abs(params[12]))
  #Lasso = sum((y - data)^2)  + lamda* sum(abs(params))
  
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
  
  S_0 = 1 - Data_test[1]*parameters['kinv_0']
  I_S0 = Data_test[1]*parameters['kinv_0']
  I_A0 = 0
  R_S0 = 0
  R_A0 = 0
  W_0 = 0
  
  S0 = 1 - Data_test[1]*parameters['kinv']
  I0 =Data_test[1]*parameters['kinv']
  W0 = 0
  R0 = 0
  
  
  times <- seq(1, length(Data_test), by = 1)
  
  init <- c(S_0, I_S0 , I_A0, R_S0, R_A0, W_0, S0, I0, W0, R0)
  
  out <- ode(y = abs(init), times = abs(times), func = model3_asymp_restrict, parms = parameters)
  out <- as.data.frame(out)
  
  p1 <- min(1,parameters['wei'])*(out[,3]/parameters['kinv_0'])
  p2 <- (1- min(1,parameters['wei']))*(out[,9]/parameters['kinv'])
  MSE <- sum(abs(Data_test - p1 - p2)^2)/length(Data_test)
  
  return(list(MSE = MSE, parameters = parameters))
}

times <- seq(1, length(DATA), by = 1)
initial_value = c((1 - DATA[1]*0.0000000000001), DATA[1]*0.0000000000001, 0 , 0 , 0, 0 ,
                  (1 - DATA[1]*0.0000116), DATA[1]*0.0000116, 0, 0)
n <- length(DATA)



############################### Ramdonly split to test and training (Algorithm 1)
lamda = 0
c <- 1
mea_vector_ASYMPexp_parallel <- c()
PARAMETERS_ASYMPexp_parallel <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.05, sd = 0.05*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  beta_W <- rnorm(1, mean = 0.3, sd = 0.3*0.1)
  alpha_S <- rnorm(1, mean = 0.0005, sd = 0.001*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.001*0.2)
  ksi <- 0
  kinv_0 <- 0.0000000000001
  q <- 1
  
  betaI <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.1)
  xi <- 0.00614
  alpha <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  kinv <- 0.0000116
  wei <- 0
  
  initial_par = c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=beta_W, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' = ksi, 'kinv_0' = kinv_0, 'q' = q,
                  'betaI' = betaI, 'betaW' = betaW, 'xi' = xi, 'alpha' = alpha,  'kinv' = kinv, 'wei' = wei)     
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_ASYMPexp_parallel <- c(mea_vector_ASYMPexp_parallel,b$MSE)
  PARAMETERS_ASYMPexp_parallel <- rbind(PARAMETERS_ASYMPexp_parallel,b$parameters)
  c <- c+1
}



lamda = (1e+06)*1
c <- 1
mea_vector_1_ASYMPexp_parallel <- c()
PARAMETERS_1_ASYMPexp_parallel <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.05, sd = 0.05*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  beta_W <- rnorm(1, mean = 0.3, sd = 0.3*0.1)
  alpha_S <- rnorm(1, mean = 0.0005, sd = 0.001*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.001*0.2)
  ksi <- 0
  kinv_0 <- 0.0000000000001
  q <- 1
  
  betaI <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.1)
  xi <- 0.00614
  alpha <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  kinv <- 0.0000116
  wei <- 0
  
  initial_par = c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=beta_W, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' = ksi, 'kinv_0' = kinv_0, 'q' = q,
                  'betaI' = betaI, 'betaW' = betaW, 'xi' = xi, 'alpha' = alpha,  'kinv' = kinv, 'wei' = wei)     
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <-  Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_1_ASYMPexp_parallel <- c(mea_vector_1_ASYMPexp_parallel,b$MSE)
  PARAMETERS_1_ASYMPexp_parallel <- rbind(PARAMETERS_1_ASYMPexp_parallel,b$parameters)
  c <- c+1
}





lamda = (1e+06)*5
c <- 1
mea_vector_2_ASYMPexp_parallel <- c()
PARAMETERS_2_ASYMPexp_parallel <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.05, sd = 0.05*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  beta_W <- rnorm(1, mean = 0.3, sd = 0.3*0.1)
  alpha_S <- rnorm(1, mean = 0.0005, sd = 0.001*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.001*0.2)
  ksi <- 0
  kinv_0 <- 0.0000000000001
  q <- 1
  
  betaI <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.1)
  xi <- 0.00614
  alpha <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  kinv <- 0.0000116
  wei <- 0
  
  initial_par = c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=beta_W, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' = ksi, 'kinv_0' = kinv_0, 'q' = q,
                  'betaI' = betaI, 'betaW' = betaW, 'xi' = xi, 'alpha' = alpha,  'kinv' = kinv, 'wei' = wei)     
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <-  Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_2_ASYMPexp_parallel <- c(mea_vector_2_ASYMPexp_parallel,b$MSE)
  PARAMETERS_2_ASYMPexp_parallel <- rbind(PARAMETERS_2_ASYMPexp_parallel,b$parameters)
  c <- c+1
}





lamda = (1e+06)*10
c <- 1
mea_vector_3_ASYMPexp_parallel <- c()
PARAMETERS_3_ASYMPexp_parallel <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.05, sd = 0.05*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  beta_W <- rnorm(1, mean = 0.3, sd = 0.3*0.1)
  alpha_S <- rnorm(1, mean = 0.0005, sd = 0.001*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.001*0.2)
  ksi <- 0
  kinv_0 <- 0.0000000000001
  q <- 1
  
  betaI <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.1)
  xi <- 0.00614
  alpha <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  kinv <- 0.0000116
  wei <- 0
  
  initial_par = c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=beta_W, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' = ksi, 'kinv_0' = kinv_0, 'q' = q,
                  'betaI' = betaI, 'betaW' = betaW, 'xi' = xi, 'alpha' = alpha,  'kinv' = kinv, 'wei' = wei)     
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <-  Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_3_ASYMPexp_parallel <- c(mea_vector_3_ASYMPexp_parallel,b$MSE)
  PARAMETERS_3_ASYMPexp_parallel <- rbind(PARAMETERS_3_ASYMPexp_parallel,b$parameters)
  c <- c+1
}




lamda = (1e+06)*15
c <- 1
mea_vector_4_ASYMPexp_parallel <- c()
PARAMETERS_4_ASYMPexp_parallel <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.05, sd = 0.05*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  beta_W <- rnorm(1, mean = 0.3, sd = 0.3*0.1)
  alpha_S <- rnorm(1, mean = 0.0005, sd = 0.001*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.001*0.2)
  ksi <- 0
  kinv_0 <- 0.0000000000001
  q <- 1
  
  betaI <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.1)
  xi <- 0.00614
  alpha <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  kinv <- 0.0000116
  wei <- 0
  
  initial_par = c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=beta_W, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' = ksi, 'kinv_0' = kinv_0, 'q' = q,
                  'betaI' = betaI, 'betaW' = betaW, 'xi' = xi, 'alpha' = alpha,  'kinv' = kinv, 'wei' = wei)     
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <-  Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_4_ASYMPexp_parallel <- c(mea_vector_4_ASYMPexp_parallel,b$MSE)
  PARAMETERS_4_ASYMPexp_parallel <- rbind(PARAMETERS_4_ASYMPexp_parallel,b$parameters)
  c <- c+1
}




lamda = (1e+06)*20
c <- 1
mea_vector_5_ASYMPexp_parallel <- c()
PARAMETERS_5_ASYMPexp_parallel <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.05, sd = 0.05*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  beta_W <- rnorm(1, mean = 0.3, sd = 0.3*0.1)
  alpha_S <- rnorm(1, mean = 0.0005, sd = 0.001*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.001*0.2)
  ksi <- 0
  kinv_0 <- 0.0000000000001
  q <- 1
  
  betaI <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.1)
  xi <- 0.00614
  alpha <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  kinv <- 0.0000116
  wei <- 0
  
  initial_par = c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=beta_W, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' = ksi, 'kinv_0' = kinv_0, 'q' = q,
                  'betaI' = betaI, 'betaW' = betaW, 'xi' = xi, 'alpha' = alpha,  'kinv' = kinv, 'wei' = wei)     
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <-  Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_5_ASYMPexp_parallel <- c(mea_vector_5_ASYMPexp_parallel,b$MSE)
  PARAMETERS_5_ASYMPexp_parallel <- rbind(PARAMETERS_5_ASYMPexp_parallel,b$parameters)
  c <- c+1
}



lamda = (1e+07)*5
c <- 1
mea_vector_6_ASYMPexp_parallel <- c()
PARAMETERS_6_ASYMPexp_parallel <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.05, sd = 0.05*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  beta_W <- rnorm(1, mean = 0.3, sd = 0.3*0.1)
  alpha_S <- rnorm(1, mean = 0.0005, sd = 0.001*0.2)
  alpha_A <- rnorm(1, mean = 0.0005, sd = 0.001*0.2)
  ksi <- 0
  kinv_0 <- 0.0000000000001
  q <- 1
  
  betaI <- rnorm(1, mean = 0.269, sd = 0.269*0.2)
  betaW <- rnorm(1, mean = 1.62, sd = 1.62*0.1)
  xi <- 0.00614
  alpha <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  kinv <- 0.0000116
  wei <- 0
  
  initial_par = c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=beta_W, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' = ksi, 'kinv_0' = kinv_0, 'q' = q,
                  'betaI' = betaI, 'betaW' = betaW, 'xi' = xi, 'alpha' = alpha,  'kinv' = kinv, 'wei' = wei)     
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <-  Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_6_ASYMPexp_parallel <- c(mea_vector_6_ASYMPexp_parallel,b$MSE)
  PARAMETERS_6_ASYMPexp_parallel <- rbind(PARAMETERS_6_ASYMPexp_parallel,b$parameters)
  c <- c+1
}


###### MSE plot
boxdata <- data.frame( A = mea_vector_ASYMPexp_parallel, 
                       B = mea_vector_1_ASYMPexp_parallel, 
                       C = mea_vector_2_ASYMPexp_parallel,
                       D = mea_vector_3_ASYMPexp_parallel,
                       E = mea_vector_4_ASYMPexp_parallel,
                       F = mea_vector_5_ASYMPexp_parallel,
                       G = mea_vector_6_ASYMPexp_parallel
)
boxplot(boxdata,ylim = c(0,(1e+08)*1))



############## plot the model fitting
Data_Generate <- function(par, initial, times){
  out <- ode(y = abs(initial), times = abs(times), func = model3_asymp_restrict, parms = abs(par))
  out <- as.data.frame(out)
  return(out)
}

DATA_PREDICTION_LASSO <- c()
init <- c((1 - DATA[1]*0.0000000000001), DATA[1]*0.0000000000001, 0 , 0 , 0, 0 ,
          (1 - DATA[1]*0.0000116), DATA[1]*0.0000116, 0, 0)

for (i in 1:dim(PARAMETERS_1_ASYMPexp_parallel)[1]){
  
  parameters <- PARAMETERS_1_ASYMPexp_parallel[i,]
  
  
  DATA_simulated <- Data_Generate(parameters, init, times)
  DATA_simulated <- min(parameters['wei'],1)*(DATA_simulated[,3]/parameters['kinv_0']) + (1-min(parameters['wei'],1))*(DATA_simulated[,9]/parameters['kinv'])
  
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
  labs(x = "Time (days)", y = "Infections", title = "Model selection (Asymptomatic vs Exponential)") +
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




######## Min-Max Normalize data and plot (make all parametyers in one scale)
# beta_IS
vector1 <- PARAMETERS_ASYMPexp_parallel[,1]   
vector2 <- PARAMETERS_6_ASYMPexp_parallel[,1]

# Get the overall min and max
combined_min <- min(c(vector1, vector2))
combined_max <- max(c(vector1, vector2))

# Normalize both vectors using the same scale
normalize_minmax <- function(x, min_val, max_val) {
  return ((x - min_val) / (max_val - min_val))
}

beta_IS_1_normalized <- normalize_minmax(vector1, combined_min, combined_max)
beta_IS_2_normalized <- normalize_minmax(vector2, combined_min, combined_max)
true_betaIS <- (0 - combined_min) / (combined_max - combined_min)

# beta_IA
vector1 <- PARAMETERS_ASYMPexp_parallel[,2]   
vector2 <- PARAMETERS_6_ASYMPexp_parallel[,2]
combined_min <- min(c(vector1, vector2))
combined_max <- max(c(vector1, vector2))
beta_IA_1_normalized <- normalize_minmax(vector1, combined_min, combined_max)
beta_IA_2_normalized <- normalize_minmax(vector2, combined_min, combined_max)
true_betaIA <- (0 - combined_min) / (combined_max - combined_min)

# alpha_S
vector1 <- PARAMETERS_ASYMPexp_parallel[,4]   
vector2 <- PARAMETERS_6_ASYMPexp_parallel[,4]
combined_min <- min(c(vector1, vector2))
combined_max <- max(c(vector1, vector2))
alpha_S_1_normalized <- normalize_minmax(vector1, combined_min, combined_max)
alpha_S_2_normalized <- normalize_minmax(vector2, combined_min, combined_max)
true_alphaS <- (0 - combined_min) / (combined_max - combined_min)

# alpha_A
vector1 <- PARAMETERS_ASYMPexp_parallel[,5]   
vector2 <- PARAMETERS_6_ASYMPexp_parallel[,5]
combined_min <- min(c(vector1, vector2))
combined_max <- max(c(vector1, vector2))
alpha_A_1_normalized <- normalize_minmax(vector1, combined_min, combined_max)
alpha_A_2_normalized <- normalize_minmax(vector2, combined_min, combined_max)
true_alphaA <- (0 - combined_min) / (combined_max - combined_min)

###### the first plot
##### combine the columns 
matrix1 <- cbind(beta_IS_1_normalized, beta_IA_1_normalized)
matrix1 <- cbind(matrix1, alpha_S_1_normalized)
matrix1 <- cbind(matrix1, alpha_A_1_normalized)
colnames(matrix1) <- c("beta_IS", "beta_IA", "alpha_S", "alpha_A")

matrix2 <- cbind(beta_IS_2_normalized, beta_IA_2_normalized)
matrix2 <- cbind(matrix2, alpha_S_2_normalized)
matrix2 <- cbind(matrix2, alpha_A_2_normalized)
colnames(matrix2) <- c("beta_IS", "beta_IA", "alpha_S", "alpha_A")

# Convert matrices into data frames in long format
df1 <- as.data.frame(matrix1)
df2 <- as.data.frame(matrix2)

# Add an ID column to differentiate the matrix origins
df1$group <- 'Without LASSO'
df2$group <- 'With LASSO (Optimal)'

# Convert the data frames to long format using the reshape2 or tidyr package
library(reshape2)
df1_long <- melt(df1, group.vars = 'group')
df2_long <- melt(df2, group.vars = 'group')

# Combine the two datasets
df_combined <- rbind(df1_long, df2_long)

ggplot(df_combined, aes(x = variable, y = value, fill = group)) + 
  geom_boxplot() +
  xlab("Parameters") +
  ylab("Values(Normalized)") +
  ggtitle("Parameters in asymptomatic model") +
  scale_y_continuous(limits = c(-0.0001, 0.2))+
  geom_segment(aes(x = 0.3, y = true_betaIS, xend = 1.4, yend = true_betaIS), 
               linetype="dashed", color="red", size=1)+
  geom_segment(aes(x = 1.6, y = true_betaIA, xend = 2.5, yend = true_betaIA), 
               linetype="dashed", color="red", size=1)+
  geom_segment(aes(x = 2.5, y = true_alphaS, xend = 3.5, yend = true_alphaS), 
               linetype="dashed", color="red", size=1)+
  geom_segment(aes(x = 3.5, y = true_alphaA, xend = 4.5, yend = true_alphaA), 
               linetype="dashed", color="red", size=1)+
  theme_minimal()








# betaI
vector1 <- PARAMETERS_ASYMPexp_parallel[,9]   
vector2 <- PARAMETERS_6_ASYMPexp_parallel[,9]
combined_min <- min(c(vector1, vector2))
combined_max <- max(c(vector1, vector2))
betaI_1_normalized <- normalize_minmax(vector1, combined_min, combined_max)
betaI_2_normalized <- normalize_minmax(vector2, combined_min, combined_max)
true_betaI <- (0.269 - combined_min) / (combined_max - combined_min)

# betaW
vector1 <- PARAMETERS_ASYMPexp_parallel[,10]   
vector2 <- PARAMETERS_6_ASYMPexp_parallel[,10]
combined_min <- min(c(vector1, vector2))
combined_max <- max(c(vector1, vector2))
betaW_1_normalized <- normalize_minmax(vector1, combined_min, combined_max)
betaW_2_normalized <- normalize_minmax(vector2, combined_min, combined_max)
true_betaW <- (1.62 - combined_min) / (combined_max - combined_min)

# alpha
vector1 <- PARAMETERS_ASYMPexp_parallel[,12]   
vector2 <- PARAMETERS_6_ASYMPexp_parallel[,12]
combined_min <- min(c(vector1, vector2))
combined_max <- max(c(vector1, vector2))
alpha_1_normalized <- normalize_minmax(vector1, combined_min, combined_max)
alpha_2_normalized <- normalize_minmax(vector2, combined_min, combined_max)
true_alpha <- (0.001314 - combined_min) / (combined_max - combined_min)


###### the second plot
##### combine the columns 
matrix1 <- cbind(betaI_1_normalized, betaW_1_normalized)
matrix1 <- cbind(matrix1, alpha_1_normalized)
colnames(matrix1) <- c("betaI", "betaW", "alpha")

matrix2 <- cbind(betaI_2_normalized, betaW_2_normalized)
matrix2 <- cbind(matrix2, alpha_2_normalized)
colnames(matrix2) <- c("betaI", "betaW", "alpha")

# Convert matrices into data frames in long format
df1 <- as.data.frame(matrix1)
df2 <- as.data.frame(matrix2)

# Add an ID column to differentiate the matrix origins
df1$group <- 'Without LASSO'
df2$group <- 'With LASSO (Optimal)'

# Convert the data frames to long format using the reshape2 or tidyr package
library(reshape2)
df1_long <- melt(df1, group.vars = 'group')
df2_long <- melt(df2, group.vars = 'group')

# Combine the two datasets
df_combined <- rbind(df1_long, df2_long)

ggplot(df_combined, aes(x = variable, y = value, fill = group)) + 
  geom_boxplot() +
  xlab("Parameters") +
  ylab("Values(Normalized)") +
  ggtitle("Parameters in exponential model") +
  scale_y_continuous(limits = c(0, 1))+
  geom_segment(aes(x = 0.3, y = true_betaI, xend = 1.4, yend = true_betaI), 
               linetype="dashed", color="red", size=1)+
  geom_segment(aes(x = 1.6, y = true_betaW, xend = 2.5, yend = true_betaW), 
               linetype="dashed", color="red", size=1)+
  geom_segment(aes(x = 2.5, y = true_alpha, xend = 3.5, yend = true_alpha), 
               linetype="dashed", color="red", size=1)+
  
  theme_minimal()













