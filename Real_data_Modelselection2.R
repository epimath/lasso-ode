## Title: "Model selection for real data: Asymptomatic model VS. Exponential VS. SIR"
## Date: 10-25-2024
## Author: Jiale Tan

library(deSolve)
library(Matrix)
library(ggplot2)
library(reshape2)

####      1. real data: 2006 cholera outbreak in Angola
DATA <- c(113, 60, 75, 148, 379, 2911, 4572, 5361, 5300, 6348, 5346,
          4412, 3558, 2271, 1931, 2251, 1692, 1184, 816, 748, 770, 522, 553, 379)
n <- length(DATA)


#####     2. Test model generator: Asymptomatic model + + Exponential + SIR model  
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
  # q = params[8]
  q = 0.2
  
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
  
  
  ##  Exponential model
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
  
  
  
  ## SIR model
  b = params[15]
  g = params[16]
  kinv_1 = params[17]
  wei_1 = params[18]
  
  S_1 = x[11]
  I_1 = x[12]
  R_1 = x[13]
  
  dS_1 = -b*S_1*I_1
  dI_1 = b*S_1*I_1 - g*I_1
  dR_1 = g*I_1
  
  
  list(c(dS_M3, dI_S, dI_A, dR_S, dR_A, dW_M3, dS, dI, dW, dR, dS_1, dI_1, dR_1))
}



#####    3. Evaluation function
SIURML1=function(params,times,data,lamda,initial_value){
  params = abs(params)
  xcurr = ode(abs(initial_value), abs(times), model3_asymp_restrict, params, method='ode45')
  y = min(params['wei'],1)*(xcurr[,3]/params['kinv_0']) + min(params['wei_1'],1)*(xcurr[,9]/params['kinv']) + (1- min(params['wei'],1) -min(params['wei_1'],1))*(xcurr[,13]/params['kinv_1'])
  
  Lasso = sum((y - data)^2)  + lamda*(abs(params[1])+abs(params[2])+abs(params[3])+abs(params[4])+abs(params[5])+abs(params[9])+abs(params[10])+abs(params[12])+abs(params[15])+abs(params[16]))
  
  return (Lasso)
}



#####    4. Blocking CV to get MSE
Blocking_CV_1 <- function(Data_train,Data_test,lamda,initial_value,initial_par){
  # choose cutting point
  initial_par <- abs(initial_par)
  initial_value <- abs(initial_value)
  times = seq(1,length(Data_train),1)*7
  
  res = optim(par = initial_par,fn=SIURML1,times=times,data=Data_train, 
              lamda = lamda,initial_value = initial_value)
  
  parameters <- abs(res$par)
  
  
  S_0 = 1 - 5*Data_test[1]*parameters['kinv_0']
  I_S0 = Data_test[1]*parameters['kinv_0']
  I_A0 = 4*Data_test[1]*parameters['kinv_0']
  R_S0 = 0
  R_A0 = 0
  W_0 = 0
  
  S0 = 1 - Data_test[1]*parameters['kinv']
  I0 = Data_test[1]*parameters['kinv']
  W0 = 0
  R0 = 0
  
  S_10 = 1 - Data_test[1]*parameters['kinv_1']
  I_10 =Data_test[1]*parameters['kinv_1']
  R_10 = 0
  
  
  times <- seq(1, length(Data_test), by = 1)*7
  
  init <- c(S_0, I_S0 , I_A0, R_S0, R_A0, W_0, S0, I0, W0, R0, S_10, I_10, R_10)
  
  out <- ode(y = abs(init), times = abs(times), func = model3_asymp_restrict, parms = parameters)
  out <- as.data.frame(out)
  
  p1 <- min(1,parameters['wei'])*(out[,3]/parameters['kinv_0'])
  p2 <- min(1,parameters['wei_1'])*(out[,9]/parameters['kinv'])
  p3 <- (1- min(parameters['wei'],1) -min(parameters['wei_1'],1))*(out[,13]/parameters['kinv_1'])
  MSE <- sum(abs(Data_test - p1 - p2 - p3)^2)/length(Data_test)
  
  return(list(MSE = MSE, parameters = parameters))
}


initial_value = c((1 - 5*DATA[1]*0.0000116), DATA[1]*0.0000116, 4*DATA[1]*0.0000116 , 0 , 0, 0 ,
                  (1 - DATA[1]*0.0000000000001), DATA[1]*0.0000000000001, 0, 0,
                  (1 - DATA[1]*0.0000000000001), DATA[1]*0.0000000000001, 0)
n <- length(DATA)



############################### Randomly split to test and training (Algorithm 1)
lamda = 0
c <- 1
mea_vector_ASYMPexpsir_parallel <- c()
PARAMETERS_ASYMPexpsir_parallel <- c()
while (c < 500) {
  beta_IS <- rnorm(1, mean = 0.3, sd = 0.3*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  beta_W <- rnorm(1, mean = 3, sd = 3*0.05)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 2, sd = 2*0.2)
  ksi <- 0.00614
  kinv_0 <- 0.0000116
  q <- 0.2
  
  betaI <- rnorm(1, mean = 0.01, sd = 0.01*0.2)         
  betaW <- rnorm(1, mean = 0.3, sd = 0.3*0.2)          
  xi <- 0.0000000000001       
  alpha <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  kinv <- 0.0000000000001
  wei <- 1
  
  b <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  kinv_1 <- 0.0000000000001
  wei_1 <- 0
  
  initial_par = c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'beta_W'=beta_W, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' = ksi, 'kinv_0' = kinv_0, 'q' = q,
                  'betaI' = betaI, 'betaW' = betaW, 'xi'= xi, 'alpha' = alpha, 'kinv' = kinv, 'wei' = wei,
                  'b' = b, 'g' = g, 'kinv_1' = kinv_1, 'wei_1' = wei_1)     
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  B <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_ASYMPexpsir_parallel <- c(mea_vector_ASYMPexpsir_parallel,B$MSE)
  PARAMETERS_ASYMPexpsir_parallel <- rbind(PARAMETERS_ASYMPexpsir_parallel,B$parameters)
  c <- c+1
}






lamda = (1e+05)*5
c <- 1
mea_vector_1_ASYMPexpsir_parallel <- c()
PARAMETERS_1_ASYMPexpsir_parallel <- c()
while (c < 500) {
  beta_IS <- rnorm(1, mean = 0.3, sd = 0.3*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  beta_W <- rnorm(1, mean = 3, sd = 3*0.05)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 2, sd = 2*0.2)
  ksi <- 0.00614
  kinv_0 <- 0.0000116
  q <- 0.2
  
  betaI <- rnorm(1, mean = 0.01, sd = 0.01*0.2)         
  betaW <- rnorm(1, mean = 0.3, sd = 0.3*0.2)          
  xi <- 0.0000000000001       
  alpha <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  kinv <- 0.0000000000001
  wei <- 1
  
  b <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  kinv_1 <- 0.0000000000001
  wei_1 <- 0
  
  initial_par = c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'beta_W'=beta_W, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' = ksi, 'kinv_0' = kinv_0, 'q' = q,
                  'betaI' = betaI, 'betaW' = betaW, 'xi'= xi, 'alpha' = alpha, 'kinv' = kinv, 'wei' = wei,
                  'b' = b, 'g' = g, 'kinv_1' = kinv_1, 'wei_1' = wei_1)     
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  B <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_1_ASYMPexpsir_parallel <- c(mea_vector_1_ASYMPexpsir_parallel,B$MSE)
  PARAMETERS_1_ASYMPexpsir_parallel <- rbind(PARAMETERS_1_ASYMPexpsir_parallel,B$parameters)
  c <- c+1
}





lamda = (1e+05)*10
c <- 1
mea_vector_2_ASYMPexpsir_parallel <- c()
PARAMETERS_2_ASYMPexpsir_parallel <- c()
while (c < 500) {
  beta_IS <- rnorm(1, mean = 0.3, sd = 0.3*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  beta_W <- rnorm(1, mean = 3, sd = 3*0.05)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 2, sd = 2*0.2)
  ksi <- 0.00614
  kinv_0 <- 0.0000116
  q <- 0.2
  
  betaI <- rnorm(1, mean = 0.01, sd = 0.01*0.2)         
  betaW <- rnorm(1, mean = 0.3, sd = 0.3*0.2)          
  xi <- 0.0000000000001       
  alpha <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  kinv <- 0.0000000000001
  wei <- 1
  
  b <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  kinv_1 <- 0.0000000000001
  wei_1 <- 0
  
  initial_par = c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'beta_W'=beta_W, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' = ksi, 'kinv_0' = kinv_0, 'q' = q,
                  'betaI' = betaI, 'betaW' = betaW, 'xi'= xi, 'alpha' = alpha, 'kinv' = kinv, 'wei' = wei,
                  'b' = b, 'g' = g, 'kinv_1' = kinv_1, 'wei_1' = wei_1)     
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  B <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_2_ASYMPexpsir_parallel <- c(mea_vector_2_ASYMPexpsir_parallel,B$MSE)
  PARAMETERS_2_ASYMPexpsir_parallel <- rbind(PARAMETERS_2_ASYMPexpsir_parallel,B$parameters)
  c <- c+1
}






lamda = (1e+05)*15
c <- 1
mea_vector_3_ASYMPexpsir_parallel <- c()
PARAMETERS_3_ASYMPexpsir_parallel <- c()
while (c < 500) {
  beta_IS <- rnorm(1, mean = 0.3, sd = 0.3*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  beta_W <- rnorm(1, mean = 3, sd = 3*0.05)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 2, sd = 2*0.2)
  ksi <- 0.00614
  kinv_0 <- 0.0000116
  q <- 0.2
  
  betaI <- rnorm(1, mean = 0.01, sd = 0.01*0.2)         
  betaW <- rnorm(1, mean = 0.3, sd = 0.3*0.2)          
  xi <- 0.0000000000001       
  alpha <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  kinv <- 0.0000000000001
  wei <- 1
  
  b <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  kinv_1 <- 0.0000000000001
  wei_1 <- 0
  
  initial_par = c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'beta_W'=beta_W, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' = ksi, 'kinv_0' = kinv_0, 'q' = q,
                  'betaI' = betaI, 'betaW' = betaW, 'xi'= xi, 'alpha' = alpha, 'kinv' = kinv, 'wei' = wei,
                  'b' = b, 'g' = g, 'kinv_1' = kinv_1, 'wei_1' = wei_1)     
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  B <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_3_ASYMPexpsir_parallel <- c(mea_vector_3_ASYMPexpsir_parallel,B$MSE)
  PARAMETERS_3_ASYMPexpsir_parallel <- rbind(PARAMETERS_3_ASYMPexpsir_parallel,B$parameters)
  c <- c+1
}






lamda = (1e+05)*20
c <- 1
mea_vector_4_ASYMPexpsir_parallel <- c()
PARAMETERS_4_ASYMPexpsir_parallel <- c()
while (c < 500) {
  beta_IS <- rnorm(1, mean = 0.3, sd = 0.3*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  beta_W <- rnorm(1, mean = 3, sd = 3*0.05)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 2, sd = 2*0.2)
  ksi <- 0.00614
  kinv_0 <- 0.0000116
  q <- 0.2
  
  betaI <- rnorm(1, mean = 0.01, sd = 0.01*0.2)         
  betaW <- rnorm(1, mean = 0.3, sd = 0.3*0.2)          
  xi <- 0.0000000000001       
  alpha <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  kinv <- 0.0000000000001
  wei <- 1
  
  b <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  kinv_1 <- 0.0000000000001
  wei_1 <- 0
  
  initial_par = c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'beta_W'=beta_W, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' = ksi, 'kinv_0' = kinv_0, 'q' = q,
                  'betaI' = betaI, 'betaW' = betaW, 'xi'= xi, 'alpha' = alpha, 'kinv' = kinv, 'wei' = wei,
                  'b' = b, 'g' = g, 'kinv_1' = kinv_1, 'wei_1' = wei_1)     
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  B <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_4_ASYMPexpsir_parallel <- c(mea_vector_4_ASYMPexpsir_parallel,B$MSE)
  PARAMETERS_4_ASYMPexpsir_parallel <- rbind(PARAMETERS_4_ASYMPexpsir_parallel,B$parameters)
  c <- c+1
}








lamda = (1e+05)*25
c <- 1
mea_vector_5_ASYMPexpsir_parallel <- c()
PARAMETERS_5_ASYMPexpsir_parallel <- c()
while (c < 500) {
  beta_IS <- rnorm(1, mean = 0.3, sd = 0.3*0.2)
  beta_IA <- rnorm(1, mean = 0.02, sd = 0.02*0.2)
  beta_W <- rnorm(1, mean = 3, sd = 3*0.05)
  alpha_S <- rnorm(1, mean = 0.001314, sd = 0.001314*0.2)
  alpha_A <- rnorm(1, mean = 2, sd = 2*0.2)
  ksi <- 0.00614
  kinv_0 <- 0.0000116
  q <- 0.2
  
  betaI <- rnorm(1, mean = 0.01, sd = 0.01*0.2)         
  betaW <- rnorm(1, mean = 0.3, sd = 0.3*0.2)          
  xi <- 0.0000000000001       
  alpha <- rnorm(1, mean = 0.0005, sd = 0.0005*0.2)
  kinv <- 0.0000000000001
  wei <- 1
  
  b <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  g <- rnorm(1, mean = 0.01, sd = 0.01*0.2)
  kinv_1 <- 0.0000000000001
  wei_1 <- 0
  
  initial_par = c('beta_IS'= beta_IS,'beta_IA'= beta_IA,'beta_W'=beta_W, 'alpha_S'=alpha_S,'alpha_A'=alpha_A, 'ksi' = ksi, 'kinv_0' = kinv_0, 'q' = q,
                  'betaI' = betaI, 'betaW' = betaW, 'xi'= xi, 'alpha' = alpha, 'kinv' = kinv, 'wei' = wei,
                  'b' = b, 'g' = g, 'kinv_1' = kinv_1, 'wei_1' = wei_1)     
  
  start <- 1
  end <- sample(5:(n-4), 1 , replace=F)
  DATA_train <- DATA[start:end]
  DATA_test <-  DATA[(end+1):n]
  B <- Blocking_CV_1(DATA_train , DATA_test, lamda,initial_value,initial_par)
  mea_vector_5_ASYMPexpsir_parallel <- c(mea_vector_5_ASYMPexpsir_parallel,B$MSE)
  PARAMETERS_5_ASYMPexpsir_parallel <- rbind(PARAMETERS_5_ASYMPexpsir_parallel,B$parameters)
  c <- c+1
}




#####  MSE plot
vector1 <- mea_vector_ASYMPexpsir_parallel
vector2 <- mea_vector_1_ASYMPexpsir_parallel
vector3 <- mea_vector_2_ASYMPexpsir_parallel
vector4 <- mea_vector_3_ASYMPexpsir_parallel
vector5 <- mea_vector_4_ASYMPexpsir_parallel
vector6 <- mea_vector_5_ASYMPexpsir_parallel
df <- data.frame(`0` = vector1, `(1e+05)*5` = vector2, `(1e+05)*10` = vector3,`(1e+05)*15` = vector4, `(1e+05)*20` = vector5, `(1e+05)*25` = vector6)

# Reshape the data to long format
df_melted <- melt(df, variable.name = "category", value.name = "value")

# Create the box plot
ggplot(df_melted, aes(x = as.numeric(category), y = value, group = category)) +
  geom_boxplot() +
  scale_x_continuous(breaks = c(1,2,3,4,5,6), labels = c("0", "(1e+05)*5", "(1e+05)*10","(1e+05)*15","(1e+05)*20","(1e+05)*25")) + 
  ylim(0 ,5e+07) +
  labs(title = "MSE in testing dataset",
       x = "Lamda",
       y = "MSE") +
  theme_minimal()




##### Parameters plot

df1 <- as.data.frame(PARAMETERS_5_ASYMPexpsir_parallel)
df1 <- df1[, c("beta_IS", "beta_IA","beta_W","alpha_S","alpha_A","b","g","betaI","betaW","alpha")]
df_long <- melt(df1, variable.name = "Variable", value.name = "Value")
df_long$ColorGroup <- ifelse(df_long$Variable %in% c("beta_IS", "beta_IA","beta_W","alpha_S","alpha_A"), "Asymptomatic model",
                 ifelse(df_long$Variable %in% c("b","g"), "SIR model",
                        "Exponential model"))
# Plot using ggplot2
ggplot(df_long, aes(x = Variable, y = Value, fill = ColorGroup)) +
  geom_boxplot() +
  scale_fill_manual(values = c("Asymptomatic model" = "blue", "SIR model"= "green","Exponential model" = "red")) +
  ylim(-0.1, 3.5)+
  labs(title = "Parameters estimation",
       x = "Parameters",
       y = "Values",
       fill = "Group") +  # Legend title
  theme_minimal()






##### calculate weights in each model
df1 <- as.data.frame(PARAMETERS_5_ASYMPexpsir_parallel)
df1 <- df1[, c("wei", "wei_1")]
df1[,1] <- ifelse(df1[,1] >= 1, 1, df1[,1])
df1[,2] <- ifelse(df1[,2] >= 1, 1, df1[,2])
df1$wei_2 <- apply(df1, 1, function(row) max(0, 1 - row[1] - row[2]))
df1 <- df1[, c("wei", "wei_2","wei_1")]

df_long <- melt(df1, variable.name = "Variable", value.name = "Value")
df_long$ColorGroup <- ifelse(df_long$Variable %in% c("wei"), "Asymptomatic model",
                             ifelse(df_long$Variable %in% c("wei_2"), "SIR model",
                                    "Exponential model"))




##### plot model fitting
Data_Generate <- function(par, initial, times){
  out <- ode(y = abs(initial), times = abs(times), func = model3_asymp_restrict, parms = abs(par))
  out <- as.data.frame(out)
  
  return(out)
}

DATA <- c(113, 60, 75, 148, 379, 2911, 4572, 5361, 5300, 6348, 5346,
          4412, 3558, 2271, 1931, 2251, 1692, 1184, 816, 748, 770, 522, 553, 379)
n <- length(DATA)
times <-  seq(1,n,1)*7
init = c((1 - 5*DATA[1]*0.0000116), DATA[1]*0.0000116, 4*DATA[1]*0.0000116 , 0 , 0, 0 ,
                  (1 - DATA[1]*0.0000000000001), DATA[1]*0.0000000000001, 0, 0,
                  (1 - DATA[1]*0.0000000000001), DATA[1]*0.0000000000001, 0)

DATA_PREDICTION_LASSO <- c()



for (i in 1:dim(PARAMETERS_5_ASYMPexpsir_parallel)[1]){
  
  parameters <- PARAMETERS_5_ASYMPexpsir_parallel[i,]
  
  DATA_simulated <- Data_Generate(parameters, init, times)
  DATA_simulated <- min(1,parameters['wei'])*DATA_simulated[,3]/parameters['kinv_0']+
                    min(1,parameters['wei_1'])*(DATA_simulated[,9]/parameters['kinv'])+
                    (1- min(parameters['wei'],1) -min(parameters['wei_1'],1))*(DATA_simulated[,13]/parameters['kinv_1'])
                    
  
  DATA_PREDICTION_LASSO <- rbind(DATA_PREDICTION_LASSO,DATA_simulated)
}


DATA_PREDICTION_LASSO <- as.matrix(DATA_PREDICTION_LASSO)
quantiles_low <- apply(DATA_PREDICTION_LASSO, 2, quantile, probs = 0.2)
quantiles_median <- apply(DATA_PREDICTION_LASSO, 2, quantile, probs = 0.5)
quantiles_high <- apply(DATA_PREDICTION_LASSO, 2, quantile, probs = 0.7)


data <- data.frame(
  x = times,
  y1 = DATA,
  y2 = quantiles_median,
  y2_lower = quantiles_low,
  y2_upper = quantiles_high
)



p <- ggplot(data, aes(x = x)) +
  geom_ribbon(aes(ymin = y2_lower, ymax = y2_upper, fill = "20%-70% quantile"), alpha = 0.5) +
  geom_line(aes(y = y2, color = "Prediction"), size = 0.5) +
  geom_point(aes(y = y1, color = "Data"), size = 0.5) +
  scale_color_manual(values = c("Data" = "blue", "Prediction" = "dark green")) +
  scale_fill_manual(values = c("20%-70% quantile" = "green")) +
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







