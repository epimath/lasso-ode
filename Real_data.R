

## Title: "Real-world data application"
## Date: 10-25-2024
## Author: Jiale Tan


library(deSolve)
library(Matrix)
library(ggplot2)
library(reshape2)

#### 1. real data: 2006 cholera outbreak in Angola
DATA <- c(113, 60, 75, 148, 379, 2911, 4572, 5361, 5300, 6348, 5346,
          4412, 3558, 2271, 1931, 2251, 1692, 1184, 816, 748, 770, 522, 553, 379)
n <- length(DATA)


##### 2. Test model generator: Asymptomatic model 
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
  q = 0.2
  
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

######### 3. Evaluation function
SIURML1=function(params,times,data,lamda,initial_value){
  params = abs(params)
  xcurr = ode(abs(initial_value), abs(times), model3_asymp_restrict, params, method='ode45')
  y = xcurr[,3]/params['kinv']
  
  Lasso = sum((y - data)^2)  + lamda* sum(abs(params))
  #Lasso = sum((y - data)^2)  + lamda*(abs(params[1])+abs(params[2])+abs(params[4])+abs(params[5])+abs(params[8]))
  
  return (Lasso)
}


########## 4. Blocking CV to get MSE
Blocking_CV_1 <- function(Data_train,Data_test,lamda,initial_value,initial_par){
  # choose cutting point
  initial_par <- abs(initial_par)
  initial_value <- abs(initial_value)
  times = seq(1,length(Data_train),1)*7
  
  res = optim(par = initial_par,fn=SIURML1,times=times,data=Data_train, 
              lamda = lamda,initial_value = initial_value)
  
  parameters <- abs(res$par)
  #print(parameters)
  
  S0 = 1 - 5*Data_test[1]*parameters['kinv']
  I_S0 = Data_test[1]*parameters['kinv']
  I_A0 = 4*Data_test[1]*parameters['kinv']
  R_S0 = 0
  R_A0 = 0
  W0 = 0
  
  times <- seq(1, length(Data_test), by = 1)*7
  
  init <- c(S0, I_S0 , I_A0, R_S0, R_A0, W0)
  
  out <- ode(y = abs(init), times = abs(times), func = model3_asymp_restrict, parms = parameters)
  out <- as.data.frame(out)
  
  MSE <- sum(abs(Data_test - out[,3]/parameters['kinv']))/length(Data_test)
  
  
  return(list(MSE = MSE, parameters = parameters))
}


initial_value = c(1 - 5*DATA[1]*1.1212e-05, DATA[1]*1.1212e-05 , 4*DATA[1]*1.1212e-05 , 0 , 0, 0 )
n <- length(DATA)


############################### Ramdonly split to test and training (Algorithm 1)
## lamda == 0
lamda = 0
c <- 1
mea_vector_real2 <- c()
PARAMETERS_real2 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.3, sd = 0.3*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  betaW <- rnorm(1, mean = 3, sd = 3*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 2, sd = 2*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  #q <- 0.2
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv) 
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_real2 <- c(mea_vector_real2,b$MSE)
  PARAMETERS_real2 <- rbind(PARAMETERS_real2,b$parameters)
  c <- c+1
}


## lamda == (1e+06)*1
lamda = (1e+05)*5
c <- 1
mea_vector_1_real2 <- c()
PARAMETERS_1_real2 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.3, sd = 0.3*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  betaW <- rnorm(1, mean = 3, sd = 3*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 2, sd = 2*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  #q <- 0.2
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv) 
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_1_real2 <- c(mea_vector_1_real2,b$MSE)
  PARAMETERS_1_real2 <- rbind(PARAMETERS_1_real2,b$parameters)
  c <- c+1
}



## lamda == (1e+06)*5
lamda = (1e+05)*10
c <- 1
mea_vector_2_real2 <- c()
PARAMETERS_2_real2 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.3, sd = 0.3*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  betaW <- rnorm(1, mean = 3, sd = 3*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 2, sd = 2*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  #q <- 0.2
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv) 
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_2_real2 <- c(mea_vector_2_real2,b$MSE)
  PARAMETERS_2_real2 <- rbind(PARAMETERS_2_real2,b$parameters)
  c <- c+1
}



## lamda == (1e+06)*10
lamda = (1e+05)*15
c <- 1
mea_vector_3_real2 <- c()
PARAMETERS_3_real2 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.3, sd = 0.3*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  betaW <- rnorm(1, mean = 3, sd = 3*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 2, sd = 2*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  #q <- 0.2
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv) 
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_3_real2 <- c(mea_vector_3_real2,b$MSE)
  PARAMETERS_3_real2 <- rbind(PARAMETERS_3_real2,b$parameters)
  c <- c+1
}


## lamda == (1e+06)*15
lamda = (1e+05)*20
c <- 1
mea_vector_4_real2 <- c()
PARAMETERS_4_real2 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.3, sd = 0.3*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  betaW <- rnorm(1, mean = 3, sd = 3*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 2, sd = 2*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  #q <- 0.2
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv) 
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_4_real2 <- c(mea_vector_4_real2,b$MSE)
  PARAMETERS_4_real2 <- rbind(PARAMETERS_4_real2,b$parameters)
  c <- c+1
}




lamda = (1e+05)*25
c <- 1
mea_vector_5_real2 <- c()
PARAMETERS_5_real2 <- c()
while (c < 300) {
  beta_IS <- rnorm(1, mean = 0.3, sd = 0.3*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  betaW <- rnorm(1, mean = 3, sd = 3*0.01)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 2, sd = 2*0.2)
  ksi <- 0.00614
  kinv <- 0.0000116
  #q <- 0.2
  
  initial_par <- c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'betaW'=betaW, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' =  ksi, 'kinv' = kinv) 
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  b <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_5_real2 <- c(mea_vector_5_real2,b$MSE)
  PARAMETERS_5_real2 <- rbind(PARAMETERS_5_real2,b$parameters)
  c <- c+1
}




########### MSE plots
vector1 <- mea_vector_real2
vector2 <- mea_vector_1_real2
vector3 <- mea_vector_2_real2
vector4 <- mea_vector_3_real2
vector5 <- mea_vector_4_real2
vector6 <- mea_vector_5_real2
df <- data.frame(`0` = vector1, `(1e+05)*5` = vector2, `(1e+05)*10` = vector3,`(1e+05)*15` = vector4, `(1e+05)*20` = vector5, `(1e+05)*25` = vector6)

# Reshape the data to long format
df_melted <- melt(df, variable.name = "category", value.name = "value")

# Create the box plot
ggplot(df_melted, aes(x = as.numeric(category), y = value, group = category)) +
  geom_boxplot() +
  scale_x_continuous(breaks = c(1,2,3,4,5,6), labels = c("0", "(1e+05)*5", "(1e+05)*10","(1e+05)*15","(1e+05)*20","(1e+05)*25")) + 
  labs(title = "MSE in testing dataset",
       x = "Lamda",
       y = "MSE") +
  theme_minimal()




#setwd("/Users/jialetan/Desktop/UM-class/Research paper/ML_Marisa/Final_Submission/Final_Results")
#load("PARAMETERS_1_real2.RData")

####### plot model fitting
Data_Generate <- function(par, initial, times){
  out <- ode(y = abs(initial), times = abs(times), func = model3_asymp_restrict, parms = abs(par))
  out <- as.data.frame(out)
  
  return(out)
}

DATA <- c(113, 60, 75, 148, 379, 2911, 4572, 5361, 5300, 6348, 5346,
          4412, 3558, 2271, 1931, 2251, 1692, 1184, 816, 748, 770, 522, 553, 379)
n <- length(DATA)
times <-  seq(1,n,1)*7
init = c(1 - 5*DATA[1]*1.1212e-05, DATA[1]*1.1212e-05 , 4*DATA[1]*1.1212e-05 , 0 , 0, 0 )

DATA_PREDICTION_LASSO <- c()

for (i in 1:dim(PARAMETERS_1_real2)[1]){
  
  parameters <- PARAMETERS_1_real2[i,]

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
  labs(x = "Time (days)", y = "Infections", title = "Time series plot from real world data scenario") +
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






