# delete everything
rm(list=ls(all=TRUE))
par(mfrow=c(1,1), mar = c(6.1, 7.1, 3.1, 2.1), lwd = 1, cex.axis = 2, cex.lab = 2.5, family = "Times New Roman")

library(deSolve) #solves differential equations
CBSD_CSS_behaviour <- function(time, y, parms)
{
  # extract state variables
  S_C <- y[1]  
  I_C <- y[2]
  S_N <- y[3]
  I_N <- y[4]
  # epidemiological parameters
  beta <- parms$beta
  gamma <- parms$gamma 
  p <- parms$p
  N <- parms$N
  
  # economic parameters
  Y <- parms$Y
  L <- parms$L
  phi <- parms$phi
  
  # behavioural parameter
  eta <- parms$eta
  modelType <- parms$modelType # this sets which model of behaviour will be used.
  #1 = "strategy vs. population" model
  #2 = "strategy vs. alternative" model
  #3 = "grower vs. population" model
  #4 = "grower vs. alternative" model
  
  # payoffs for each outcome
  P_SN <- Y
  P_IN <- Y - L
  P_SC <- Y - phi
  P_IC <- Y - phi - L
  
  # probabilities of infection
  pVert <- p*(I_C+I_N)/N
  
  qBeta <- parms$nu_beta * beta # parameter for mis-estimating beta
  nu_I <- parms$nu_I # parameter for mis-estimating the prevalence of infection 

  if(nu_I*(I_C+I_N) < N){ # need to ensure that estimate of infection does not exceed the total number of fields
    pVert <- p*nu_I*(I_C+I_N)/N
    q_C <- qBeta * nu_I*(I_C + I_N) / ( qBeta * nu_I*(I_C + I_N) + gamma) 
    q_N <- pVert + (1 - pVert) * q_C
  } else{
    pVert <- p*N/N
    q_C <- qBeta *N / ( qBeta * N + gamma) 
    q_N <- pVert + (1 - pVert) * q_C
  }
  
  
  q_C <- qBeta * (I_C + I_N) / ( qBeta * (I_C + I_N) + gamma) 
  q_N <- pVert + (1 - pVert) * q_C
  
  P_C <- q_C * P_IC + (1 - q_C) * P_SC
  P_N <- q_N * P_IN + (1 - q_N) * P_SN
  
  P <- P_C * (S_C + I_C) / N + P_N * (S_N + I_N) / N
  
  # switching terms
  z_SC <- 0
  z_IC <- 0
  z_SN <- 0
  z_IN <- 0
  
  if(modelType == 1) # Strategy vs. population
  {
    z_SC <- max(0, 1 - exp(-eta*(P - P_C)))
    z_IC <- z_SC
    z_SN <- max(0, 1 - exp(-eta*(P - P_N)))
    z_IN <- z_SN
  } else if(modelType == 2) # Strategy vs. alternative
  {
    z_SC <- max(0, 1 - exp(-eta*(P_N - P_C)))
    z_IC <- z_SC
    z_SN <- max(0, 1 - exp(-eta*(P_C - P_N)))
    z_IN <- z_SN
  } else if(modelType == 3) # Grower vs. population
  {
    z_SC <- max(0, 1 - exp(-eta*(P - P_SC)))
    z_IC <- 1 - exp(-eta*(P - P_IC))
    z_SN <- 0
    z_IN <- max(0, 1 - exp(-eta*(P - P_IN)))
  } else if(modelType == 4) # Grower vs. alternative
  {
    z_SC <- max(0, 1 - exp(-eta*(P_N - P_SC)))
    z_IC <- 1 - exp(-eta*(P_N - P_IC))
    z_SN <- 0
    z_IN <- max(0, 1 - exp(-eta*(P_C - P_IN)))
  }
  
  plantingRate_C <- gamma * (S_C * (1 - z_SC) + I_C * (1 - z_IC) + S_N * z_SN + I_N * z_IN)
  plantingRate_N <- gamma * (S_C * z_SC + I_C * z_IC + S_N * (1 - z_SN) + I_N * (1 - z_IN))
  
  dS_C <- plantingRate_C - beta * S_C * (I_C + I_N) - gamma * S_C
  dI_C <- beta * S_C * (I_C + I_N) - gamma * I_C
  dS_N <- plantingRate_N * (1 - pVert) - beta * S_N * (I_C + I_N) - gamma * S_N
  dI_N <- plantingRate_N * pVert + beta * S_N * (I_C + I_N) - gamma * I_N
  
  return(list(c(dS_C, dI_C, dS_N, dI_N)))
}


readParams <- function(N, beta, gamma, p, Y, phi, L, nu_beta, nu_I, eta, modelType) ## function to return parameters of the model
{
  retval <- list(N = N,
                 beta = beta, 
                 gamma = gamma, 
                 p = p,
                 Y = Y,
                 phi = phi,
                 L = L, 
                 nu_beta = nu_beta,
                 nu_I = nu_I,
                 eta = eta, 
                 modelType = modelType)
  return(retval)
}

#Below are default parameters without any mis-estimation using the "strategy vs. population" model. Leads to no control at equilibrium.
parameters <- readParams(beta = .0067, gamma = 1/300,p =.8, N = 1, Y = 1, phi = 0.25, L = 0.6, nu_beta = 1, nu_I = 1, eta = 10, modelType = 1) 
yini<- c(0.1, 0, 0.99*(1-0.1), (1-0.99)*(1-0.1)) #Initial conditions in format (S_C, I_C, S_N, I_N)
times <- seq(0,10*3300, .10) #Timespan to run model

out <- ode(func = CBSD_CSS_behaviour,  y = yini, time = times, parms = parameters)
plot(out[,1], out[,2], lwd = 3, ty = "l", ylim = c(0,1), col = "coral", main = "", ylab = "Proportion of growers", xlab = "Time (seasons)", xaxt = "n") #SC
axis(side = 1, at = seq(0, max(times ), max(times /5)),labels = seq(0, max(times /300), max(times /300)/5))
lines(out[,1], out[,3], lwd = 3, col = "firebrick") #IC
lines(out[,1], out[,4], lwd = 3, lty = 1, col = "dodgerblue3") #SN
lines(out[,1], out[,5], lwd = 3, lty = 1, col = "light blue") #IN
lines(out[,1], out[,2] + out[,3], col ="chartreuse3", lwd = 3) #All controllers
legend("right", ncol =1, cex = 0.9,legend = c( expression("S"[N]*""),expression("S"[C]*""), expression("I"[N]*""), expression("I"[C]*""), "Control"), col = c(" dodgerblue3","coral", "light blue", "firebrick", "chartreuse3"), lwd =2 )

######## Below is a scan over the rate of horizontal tranmsission ####
parameters <- readParams(beta = .0067, gamma = 1/300,p =.8, N = 1, Y = 1, phi = 0.25, L = 0.6, nu_beta = 1, nu_I = 1, eta = 10, modelType = 4) 
yini<- c(0.1, 0, 0.99*(1-0.1), (1-0.99)*(1-0.1)) #Initial conditions in format (S_C, I_C, S_N, I_N)
times <- seq(0,50*300, 10) #Timespan to run model

beta_sequence_det <- seq(0, 0.015, 0.0001)
control_equilibrium_det <- c()
infection_equilibrium_det <- c()
for(beta in beta_sequence_det){
  parameters$beta <- beta
  out <- ode(func =CBSD_CSS_behaviour,  y = yini, time = times, parms = parameters)
  S_CS <- out[length(out[,1]),2]
  I_CS <- out[length(out[,1]),3]
  S_NC <- out[length(out[,1]),4]
  I_NC <- out[length(out[,1]),5]
  control_equilibrium_det <- c(control_equilibrium_det, (S_CS + I_CS))
  infection_equilibrium_det <- c(infection_equilibrium_det, (I_NC + I_CS))
}
plot(beta_sequence_det, control_equilibrium_det, ty = "l", lwd = 3, col = "#80cc99", ylim = c(0,1), ylab = "", xlab = expression(paste("Rate of horizontal transmission (", beta, ")")))
lines(beta_sequence_det, infection_equilibrium_det, lwd = 3, col ="#9980cc")
title(main = expression("Deterministic model"), cex.main = 2.5)

