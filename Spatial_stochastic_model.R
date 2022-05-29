
#########Spatial stochastic 
################################ Libraries #############################################################################
library(colorspace)
library(fields)
library(matrixStats)
################################ Functions #############################################################################

## Parameters 
readParams <- function(N = 750,
                       beta = 110,
                       disp_i  = 150,
                       gamma = 1/300,
                       S_C = 0.1, 
                       I_C = 0, 
                       S_N = 0.9, 
                       I_N = 0.01,
                       yMax = 3160,
                       xMax = 3160, 
                       Y = 1,
                       phi = 0.1, 
                       L = 0.6, 
                       p = 0.9, 
                       eta = 10,
                       disp_i_est = 1, #parameter for mis-estimating the rate of whitefly dispersal
                       beta_s = 1){ 
  retval <- list(N = N,beta = beta, disp_i = disp_i, S_C = N*S_C, I_C = N*I_C, 
                 S_N = N*S_N*(1 - I_N), I_N = N*S_N*I_N, xMax = xMax, yMax =yMax,
                 Y = Y, phi =phi, L = L, p = p, gamma = gamma, eta= eta, disp_i_est = disp_i_est, beta_s = beta_s)
  return(retval)
  
}

## Generate location of hosts 
generateHosts <- function(params){
  x <- runif(params$N,0,params$xMax)  
  y <- runif(params$N,0,params$yMax)
  retVal <- data.frame(x,y)
  return(retVal)
}

## Initial state of the system 
Initial_States <- function(params){
  S_C <- rep("S_C", params$S_C)
  S_N <- rep("S_N", params$S_N)
  I_C <- rep("I_C", params$I_C)
  I_N <- rep("I_N", length.out = params$I_N + 1)
  Initial_states <- c(S_C, S_N, I_C, I_N)
  return(Initial_states)
}

## Dataframe to hold locations, state of each field, time of infection, previous yield etc. 
CBSD_dataframe <-  function(params, hosts, states){
  CBSD_dataframe   <- data.frame(matrix(0, nrow = params$N, ncol = 11))
  colnames(CBSD_dataframe  ) <- c("x_coord", "y_coord", "ID", "State", "Harvest_Rate", "Infection","tI", "Prev_State", "Y", "Profit", "Total_rate")    
  CBSD_dataframe  [,1] <- hosts$x  
  CBSD_dataframe  [,2] <- hosts$y
  CBSD_dataframe  [,3] <- 1:params$N  
  CBSD_dataframe  [,4] <- sample(states, params$N, replace = FALSE) 
  CBSD_dataframe  [,5] <- params$gamma
  return(CBSD_dataframe  )
}

## Distance between fields
Distance_Matrix <- function(params, dataframe, hosts){
  Distance_matrix <- matrix(0, nrow = params$N, ncol = params$N)      
  colnames(Distance_matrix) <- 1:params$N
  rownames(Distance_matrix) <- 1:params$N
  Distance_matrix <- as.matrix(dist(hosts))
  
  return(Distance_matrix)
}

## Probability of infection between fields - exponential kernel 
Kernel_Matrix_infection <- function(params, Dist_matrix){
  Kernel_matrix <- matrix(0, nrow = params$N, ncol = params$N)      
  colnames(Kernel_matrix) <- 1:params$N    
  rownames(Kernel_matrix) <- 1:params$N    
  for(i in 1:params$N){      
    for (j in 1:params$N){      
      Kernel_matrix[[i,j]] <- exp(-1*(Dist_matrix[[i,j]])/params$disp_i)
      Kernel_matrix[[i,j]] <- Kernel_matrix[i,j]*(1/(2*pi*as.numeric(params$disp_i)*as.numeric(params$disp_i))) 
    }      
  }      
  diag(Kernel_matrix) <- 0
  return(Kernel_matrix)
}

## Mis-estimated probability of infection between fields - exponential kernel 
Kernel_Matrix_Estimation<- function(params, Dist_matrix){
  Kernel_matrix <- matrix(0, nrow = params$N, ncol = params$N)      
  colnames(Kernel_matrix) <- 1:params$N    
  rownames(Kernel_matrix) <- 1:params$N    
  for(i in 1:params$N){      
    for (j in 1:params$N){      
      Kernel_matrix[[i,j]] <- exp(-1*(Dist_matrix[[i,j]])/(params$disp_i_est * params$disp_i))
      Kernel_matrix[[i,j]] <- Kernel_matrix[i,j]*(1/(2*pi*as.numeric((params$disp_i_est * params$disp_i))*as.numeric((params$disp_i_est * params$disp_i)))) 
    }      
  }      
  diag(Kernel_matrix) <- 0
  return(Kernel_matrix)
}



## Calculates the initial probability of infection for each field based on other infected fields 
Infection_pressure <- function(dataframe, params, infection_matrix){
  for(Field in 1:params$N){
    if(dataframe$State[Field] == "I_N" | dataframe$State[Field] == "I_C"){ ### All fields can be infected by all other fields 
      dataframe$Infection <-dataframe$Infection+params$beta*infection_matrix[,Field]*(dataframe$State=="S_C" | dataframe$State=="S_N") #Updates the probability of infection for all the susceptible fields
    }
  }
  
  return(dataframe)
}

## Calculates the initial probability of infection for each field based on other infected fields during the epidemic
epidemic_infection <- function(dataframe, params, kernel_infection, host){
  if(dataframe$State[host] == "S_C"){ #if a S_C field is infected, it becomes I_C
    dataframe$State[host] <- "I_C"
    dataframe$tI[host] <- t
    dataframe$Infection[host] <- 0
    inf_pressure <- kernel_infection[,host]
    dataframe$Infection <-dataframe$Infection+params$beta*inf_pressure*(dataframe$State=="S_N"|dataframe$State=="S_C" ) #only updates the infection pressure for susceptible fields
    
  }else if(dataframe$State[host] == "S_N"){ #if a S_C field is infected, it becomes I_N
    dataframe$State[host] <- "I_N"
    dataframe$tI[host] <- t
    dataframe$Infection[host] <- 0
    inf_pressure <- kernel_infection[,host]
    dataframe$Infection <-dataframe$Infection+params$beta*inf_pressure*(dataframe$State=="S_N"|dataframe$State=="S_C" ) 
  }
  return(dataframe)
}

## Calculate yields and profit 
Strategy_profits <- function(dataframe, params){
  dataframe$Profit <- params$Y - params$phi*(dataframe$State=="S_C"| dataframe$State=="I_C")
  dataframe$Profit <- dataframe$Profit - params$L*(dataframe$State=="I_N"| dataframe$State=="I_C")

  return(dataframe)
}


next_strategy_compare_grower_to_strategy <- function(dataframe, params,  host,prev_state_host, no_infected, which_I, kernel_matrix){
  profit_host <- dataframe$Profit[host]
  grower_states <-  dataframe$State
  controllers <- c( which(grower_states == "S_C" | grower_states == "I_C"), 0)
  non_controllers <- c(which(grower_states== "S_N"| grower_states == "I_N"),0)
  
  infected <- c(which(grower_states== "I_C"| grower_states == "I_N"),0)
  No_I <- length(infected)
  
  FOI <-  params$beta_s*params$beta*sum(rowSums(kernel_matrix[, infected, drop = FALSE ]))/params$N

  expected_PC <- ((FOI)/(FOI + params$gamma))*(P_IC) + ((params$gamma)/(FOI + params$gamma))*(P_SC)
  
  q <- (params$p*(No_I/params$N)) + (1 - params$p*(No_I/params$N))*(FOI)/(FOI+ params$gamma)
  
  expected_PN<- (q)*(P_IN) + (1-q)*(P_SN)
  
  
  if(prev_state_host == "S_N" | prev_state_host == "I_N"){
    ifelse(profit_host  > expected_PC, State_host <-  "S_N",State_host <- sample(c("S_C", "S_N"), 1, F, c(1-exp(-1*params$eta*(expected_PC - profit_host)), 1-(1-exp(-1*params$eta*(expected_PC - profit_host))))))
    
  }else if(prev_state_host == "S_C"){## maybe change if S_C
    ifelse(profit_host  > expected_PN,  State_host <- "S_C",dataframe$State[host]<-State_host <- sample(c("S_N", "S_C"), 1, F, c(1-exp(-1*params$eta*(expected_PN - profit_host)), 1-(1-exp(-1*params$eta*(expected_PN - profit_host))))))
    
  } else{
    State_host <- "S_N"
  }
  return(State_host)
}


####### Runs ##############

parameters <- readParams(beta =110, S_C = 0.1, I_N = 0.01, S_N = .9, eta = 10, phi = 0.25, L = 0.6, p= .8, N = 750, xMax =  3160, yMax =  3160, disp_i_est =  1, beta_s = 1)

host_locs <- generateHosts(parameters)
Initial_states<- Initial_States(parameters)
CBSD_CSS_Dataframe <- CBSD_dataframe(parameters, host_locs, Initial_states)
Distance_matrix <- Distance_Matrix(parameters, CBSD_CSS_Dataframe, host_locs)
Kernel_matrix <- Kernel_Matrix_infection(parameters, Distance_matrix) ## These seem very low, consider paramaterisation 
Kernel_Matrix_estimation <- Kernel_Matrix_Estimation(parameters, Distance_matrix) ## These seem very low, consider paramaterisation 
CBSD_CSS_Dataframe <- Infection_pressure(CBSD_CSS_Dataframe, parameters, Kernel_matrix)
No_Its <- 10
tmax <- 10*300
t_seq <- seq(0.1,tmax,10)


S_C_avg <-  data.frame(matrix(0, ncol = No_Its, nrow = length(t_seq))) ## Need to be able to store outcomes for each run. 
S_N_avg <- data.frame(matrix(0, ncol = No_Its, nrow = length(t_seq)))  ## Later used to calculate average
I_C_avg <-  data.frame(matrix(0, ncol = No_Its, nrow = length(t_seq)))
I_N_avg <- data.frame(matrix(0, ncol = No_Its, nrow = length(t_seq)))
average_profit <- data.frame(matrix(0, ncol = No_Its, nrow = length(t_seq)))

Extinction_events <- data.frame(matrix(0,nrow = 1, ncol = No_Its ))
P_SC <- parameters$Y - parameters$phi
P_IC <- parameters$Y - parameters$phi - parameters$L
P_SN <- parameters$Y 
P_IN <- parameters$Y - parameters$L
Total_harvest_rate <- parameters$N*parameters$gamma #This is constant so only need to calculate once
prev_infected <- 0
EXPECTED <- F



for(x in 1:No_Its){
  print(x)
  CBSD_CSS_Dataframe <- CBSD_dataframe(parameters, host_locs, Initial_states)
  CBSD_CSS_Dataframe <- Infection_pressure(CBSD_CSS_Dataframe, parameters, Kernel_matrix)
  CBSD_CSS_Dataframe$Total_rate <- CBSD_CSS_Dataframe$Harvest_Rate + CBSD_CSS_Dataframe$Infection
  ## Dataframes to store 
  No_S_C_sequence <- matrix(NA, ncol= 1, nrow =300000)
  No_S_N_sequence <-matrix(NA, ncol= 1, nrow =300000)
  No_I_N_sequence <- matrix(NA, ncol= 1, nrow =300000)
  No_I_C_sequence <- matrix(NA, ncol= 1, nrow =300000)
  Time_sequence <- matrix(NA, ncol= 1, nrow =300000)
  Average_profit <- matrix(NA, ncol= 1, nrow =300000)
  
  ## Initial conditions 
  t <- 0 
  No_S_C <- sum(CBSD_CSS_Dataframe$State=="S_C")     
  No_S_N <- sum(CBSD_CSS_Dataframe$State=="S_N")     
  No_I_C <- sum(CBSD_CSS_Dataframe$State=="I_C") 
  No_I_N <- sum(CBSD_CSS_Dataframe$State=="I_N") 
  Total_time <- t
  Total_rates <- sum(CBSD_CSS_Dataframe$Total_rate)
  i = 0
  No_I <- sum(No_I_C, No_I_N)
  which_I <- as.numeric(c(which(CBSD_CSS_Dataframe$State == "I_C"), which(CBSD_CSS_Dataframe$State == "I_N"))) ### but each individual host will have their own FOI

 
  while(t < tmax ){
    dt <- rexp(1, Total_rates) # decide time until next event
    t = t + dt
    Total_rates <- sum(Total_harvest_rate ,CBSD_CSS_Dataframe$Infection)
    
    Event <- sample(c("Harvest", "Secondary"), 1,replace = F, c(sum(CBSD_CSS_Dataframe$Harvest_Rate), sum(CBSD_CSS_Dataframe$Infection)))
    
    if(Event == "Harvest"){
      Host <- sample(1:parameters$N,1,prob=CBSD_CSS_Dataframe$Harvest_rate)
      Prev_state_host <- CBSD_CSS_Dataframe$Prev_State[Host] <-  CBSD_CSS_Dataframe$State[Host]
      
      if(Prev_state_host == "S_N"){
        No_S_N <- No_S_N - 1
      } else if(Prev_state_host == "S_C"){
        No_S_C <- No_S_C - 1
      }else  if(Prev_state_host == "I_N"){
        No_I_N <- No_I_N - 1
        No_I <- No_I - 1
      } else if(Prev_state_host == "I_C"){
        No_I_C <- No_I_C - 1
        No_I <- No_I - 1 
        
      } 
      CBSD_CSS_Dataframe <- Strategy_profits(CBSD_CSS_Dataframe, parameters) ## get profits for each individual 
      
      State_host <- CBSD_CSS_Dataframe$State[Host]<- next_strategy_compare_grower_to_strategy(CBSD_CSS_Dataframe, parameters, Host, Prev_state_host, no_infected  = No_I, which_I, kernel_matrix = Kernel_Matrix_estimation)
      
      if(State_host == "S_N" ){
        prev_infected <- (No_I)/parameters$N
        infected <- sample(c(T, F),  1,c(parameters$p*(prev_infected), (1-parameters$p*(prev_infected))), replace = F)
        if(infected){
          CBSD_CSS_Dataframe$State[Host] <- State_host <- "I_N"
          CBSD_CSS_Dataframe$Infection[Host] <- 0
          CBSD_CSS_Dataframe$tI[Host] <- t
          No_I <- No_I + 1
          No_I_N <- No_I_N + 1
          if(Prev_state_host == "I_N" | Prev_state_host == "I_C"){
            which_I <- which_I[!(which_I %in% as.numeric(Host))]
          } else{
            inf_pressure <- (Kernel_matrix[,Host])
            CBSD_CSS_Dataframe$Infection <- CBSD_CSS_Dataframe$Infection + parameters$beta* inf_pressure*(CBSD_CSS_Dataframe$State=="S_C" | CBSD_CSS_Dataframe$State=="S_N" ) 
          }
          which_I <- c(as.numeric(which_I), as.numeric(Host))
          
        } else{
          if(Prev_state_host == "I_N" | Prev_state_host == "I_C"){
            inf_pressure <- (Kernel_matrix[,Host])
            CBSD_CSS_Dataframe$Infection <- CBSD_CSS_Dataframe$Infection - parameters$beta* inf_pressure*(CBSD_CSS_Dataframe$State=="S_C" | CBSD_CSS_Dataframe$State=="S_N" ) 
            which_I <- which_I[!(which_I %in% as.numeric(Host))]
          }
          CBSD_CSS_Dataframe$Infection[Host] <-  0 + parameters$beta*sum(Kernel_matrix[Host,which_I])
          No_S_N <- No_S_N + 1
          
        }
      } 
      if(State_host == "S_C") {
        if(Prev_state_host == "I_N" | Prev_state_host == "I_C"){
          inf_pressure <- (Kernel_matrix[,Host])
          CBSD_CSS_Dataframe$Infection <- CBSD_CSS_Dataframe$Infection - parameters$beta*inf_pressure*(CBSD_CSS_Dataframe$State=="S_C" | CBSD_CSS_Dataframe$State=="S_N" ) 
          which_I <- which_I[!(which_I %in% as.numeric(Host))]

        }
        No_S_C <- No_S_C + 1
        CBSD_CSS_Dataframe$Infection[Host] <-  0 + parameters$beta*sum(Kernel_matrix[Host,which_I])
        
      }
      
      
    } else if(Event == "Secondary"){
      Host <- sample(1:parameters$N,1,prob=CBSD_CSS_Dataframe$Infection)
      if(CBSD_CSS_Dataframe$State[Host] == "S_N"){
        No_S_N <- No_S_N - 1
        No_I_N <- No_I_N + 1
        No_I <- No_I + 1
        which_I <- c(as.numeric(which_I), as.numeric(Host))
        
      } else if(CBSD_CSS_Dataframe$State[Host] == "S_C"){
        No_S_C <- No_S_C - 1
        No_I_C <- No_I_C + 1
        No_I <- No_I + 1
        which_I <- c(as.numeric(which_I), as.numeric(Host))
        
      } 
      CBSD_CSS_Dataframe <- epidemic_infection(CBSD_CSS_Dataframe, parameters, Kernel_matrix, Host)
      
    }
    i = i + 1 
    No_S_C_sequence[i,1] <- No_S_C
    No_S_N_sequence[i,1] <- No_S_N
    No_I_N_sequence[i,1] <- No_I_N
    No_I_C_sequence[i,1]<- No_I_C
    Time_sequence[i,1]<- t
    Average_profit[i,1] <- mean(CBSD_CSS_Dataframe$Profit)
    CBSD_CSS_Dataframe$Infection[CBSD_CSS_Dataframe$Infection < 1e-10] <- 0 

  }
  if(sum(CBSD_CSS_Dataframe$State == "I_N",CBSD_CSS_Dataframe$State == "I_C" ) == 0){
    Extinction_events[x] = Extinction_events[x]+ 1
  }
  frame <- data.frame(Time_sequence, No_S_C_sequence, No_S_N_sequence, No_I_C_sequence, No_I_N_sequence, Average_profit) # create dataframe of vectors for each time 
  a = 1
  for(m in t_seq){ # loop through time
    if(length(which_I) == 0){
      Extinction_events[x] = Extinction_events[x]+ 1
    }
    z <- max(which(frame$Time_sequence < m))
    if(is.infinite(z)){
      z <- 1
    }# find the maximum row that is still less than time t 
    S_C_avg[a,x] <- frame$No_S_C_sequence[z] # subset from that row to get the no. S_S etc individuals
    S_N_avg[a,x] <- frame$No_S_N_sequence[z]
    I_C_avg[a,x] <- frame$No_I_C_sequence[z]
    I_N_avg[a,x] <- frame$No_I_N_sequence[z]
    average_profit[a,x] <- frame$Average_profit[z]
    a = a+1
  }
  
  infection <-  I_C_avg[,x] +  I_N_avg[,x]
  control <- I_C_avg[,x] +  S_C_avg[,x]

}


N <- parameters$N

S_C_avg$mean <- rowMeans(S_C_avg, na.rm=TRUE)/N
I_C_avg$mean <- rowMeans(I_C_avg, na.rm=TRUE)/N
S_N_avg$mean <- rowMeans(S_N_avg, na.rm=TRUE)/N
I_N_avg$mean <- rowMeans(I_N_avg, na.rm=TRUE)/N
average_profit$mean <- rowMeans(average_profit, na.rm=TRUE)

S_C_avg$sd <- rowSds(as.matrix(S_C_avg[,1:No_Its]))/N
I_C_avg$sd <- rowSds(as.matrix(I_C_avg[,1:No_Its]))/N
S_N_avg$sd <- rowSds(as.matrix(S_N_avg[,1:No_Its]))/N
I_N_avg$sd <- rowSds(as.matrix(I_N_avg[,1:No_Its]))/N
average_profit$sd <- rowSds(as.matrix( average_profit[,1:No_Its]))

plot(t_seq, S_C_avg$mean , lwd = 4, ty = "l", ylim = c(0,1), col = "coral", ylab = "Proportion of growers", xlab = "Time (seasons)",xaxt = "n", cex.lab = 1.9, cex.axis = 1.6,xaxs = "i", main = "")
axis(side = 1, at = seq(0, max(tmax), max(tmax/5)),labels = seq(0, max(tmax/300), max(tmax/300)/5), cex.axis = 1.6)
polygon(c(t_seq,rev(t_seq)), c(S_C_avg$mean +S_C_avg$sd, rev(S_C_avg$mean - S_C_avg$sd)), col=adjustcolor("coral",alpha.f=0.3), border = NA)
polygon(c(t_seq,rev(t_seq)), c(I_C_avg$mean + I_C_avg$sd , rev(I_C_avg$mean - I_C_avg$sd)), col=adjustcolor("firebrick",alpha.f=0.3), border = NA)
polygon(c(t_seq,rev(t_seq)), c(I_N_avg$mean +I_N_avg$sd, rev(I_N_avg$mean - I_N_avg$sd)), col=adjustcolor("light blue",alpha.f=0.3), border = NA)
polygon(c(t_seq,rev(t_seq)), c(S_N_avg$mean + S_N_avg$sd, rev(S_N_avg$mean - S_N_avg$sd)), col=adjustcolor("dodgerblue3",alpha.f=0.3), border = NA)
lines(t_seq, I_C_avg$mean , lwd = 4, col = "firebrick")
lines(t_seq,S_N_avg$mean , lwd = 4, lty = 1, col = "dodgerblue3")
lines(t_seq, I_N_avg$mean , lwd = 4, lty = 1, col = "light blue")
lines(t_seq, (I_C_avg$mean + S_C_avg$mean), lwd = 4, lty = 2, col = "chartreuse 3")
#lines(t_seq, average_profit$mean , lwd = 2, lty = 2, col = "darkorchid4")
#legend("right", ncol = 2, cex = 1.2,legend = c(expression("S"[CS]*""), expression("S"[NC]*""), expression("I"[CS]*""), expression("I"[NC]*""), "Control"), col = c("coral", " dodgerblue3", "firebrick", "light blue", "chartreuse3"), lwd =2 )


############ Scan over parameters - beta  ########3

parameters <- readParams(beta =110, S_C = 0.1, I_N = 0.01, S_N = .9, eta = 10, phi = 0.25, L = 0.6, p= .8, N = 750, xMax =  3160, yMax =  3160, disp_i_est =  1, beta_s = 1)

host_locs <- generateHosts(parameters)
Initial_states<- Initial_States(parameters)
CBSD_CSS_Dataframe <- CBSD_dataframe(parameters, host_locs, Initial_states)
Distance_matrix <- Distance_Matrix(parameters, CBSD_CSS_Dataframe, host_locs)
Kernel_matrix <- Kernel_Matrix_infection(parameters, Distance_matrix) 

CBSD_CSS_Dataframe <- Infection_pressure(CBSD_CSS_Dataframe, parameters, Kernel_matrix)
No_Its <- 10
tmax <- 10*300
t_seq <- tmax

S_C_avg <-  data.frame(matrix(0, ncol = No_Its, nrow = length(t_seq))) ## Need to be able to store outcomes for each run. 
S_N_avg <- data.frame(matrix(0, ncol = No_Its, nrow = length(t_seq)))  ## Later used to calculate average
I_C_avg <-  data.frame(matrix(0, ncol = No_Its, nrow = length(t_seq)))
I_N_avg <- data.frame(matrix(0, ncol = No_Its, nrow = length(t_seq)))
average_profit <- data.frame(matrix(0, ncol = No_Its, nrow = length(t_seq)))


Extinction_events <- data.frame(matrix(0,nrow = 1, ncol = No_Its ))
P_SC <- parameters$Y - parameters$C
P_IC <- parameters$Y - parameters$C - parameters$L
P_SN <- parameters$Y 
P_IN <- parameters$Y - parameters$L
Total_harvest_rate <- parameters$N*parameters$gamma
prev_infected <- 0
EXPECTED <- F
N <- 750

beta_sequence_stoc <- c(seq(0.1,40.1,20))
control_equilibrium <- c()
infection_equilibrium <- c()
infection_sd <- c()
control_sd <- c()
for(beta in beta_sequence_stoc){
  parameters$beta <- beta
  
  S_C_avg <-  data.frame(matrix(0, ncol = No_Its, nrow = 1)) ## Need to be able to store outcomes for each run. 
  S_N_avg <- data.frame(matrix(0, ncol = No_Its, nrow =1))  ## Later used to calculate average
  I_C_avg <-  data.frame(matrix(0, ncol = No_Its, nrow = 1))
  I_N_avg <- data.frame(matrix(0, ncol = No_Its, nrow = 1))
  Extinction_events <- data.frame(matrix(0,nrow = 1, ncol = No_Its ))
  for(x in 1:No_Its){
    print(x)
    CBSD_CSS_Dataframe <- CBSD_dataframe(parameters, host_locs, Initial_states)
    CBSD_CSS_Dataframe <- Infection_pressure(CBSD_CSS_Dataframe, parameters, Kernel_matrix)
    CBSD_CSS_Dataframe$Total_rate <- CBSD_CSS_Dataframe$Harvest_Rate + CBSD_CSS_Dataframe$Infection
    ## Dataframes to store 
    No_S_C_sequence <- matrix(NA, ncol= 1, nrow =300000)
    No_S_N_sequence <-matrix(NA, ncol= 1, nrow =300000)
    No_I_N_sequence <- matrix(NA, ncol= 1, nrow =300000)
    No_I_C_sequence <- matrix(NA, ncol= 1, nrow =300000)
    Time_sequence <- matrix(NA, ncol= 1, nrow =300000)
    Average_profit <- matrix(NA, ncol= 1, nrow =300000)
    
    ## Initial conditions 
    t <- 0 
    No_S_C <- sum(CBSD_CSS_Dataframe$State=="S_C")     
    No_S_N <- sum(CBSD_CSS_Dataframe$State=="S_N")     
    No_I_C <- sum(CBSD_CSS_Dataframe$State=="I_C") 
    No_I_N <- sum(CBSD_CSS_Dataframe$State=="I_N") 
    Total_time <- t
    Total_rates <- sum(CBSD_CSS_Dataframe$Total_rate)
    i = 0
    No_I <- sum(No_I_C, No_I_N)
    which_I <- as.numeric(c(which(CBSD_CSS_Dataframe$State == "I_C"), which(CBSD_CSS_Dataframe$State == "I_N"))) ### but each individual host will have their own FOI
    
    prev_infected <- 0
    ## Sample run 
    while(t < tmax ){
      dt <- rexp(1, Total_rates) # decide time until next event
      t = t + dt
      Total_rates <- sum(Total_harvest_rate ,CBSD_CSS_Dataframe$Infection)
      
      Event <- sample(c("Harvest", "Secondary"), 1,replace = F, c(sum(CBSD_CSS_Dataframe$Harvest_Rate), sum(CBSD_CSS_Dataframe$Infection)))
      
      if(Event == "Harvest"){
        Host <- sample(1:parameters$N,1,prob=CBSD_CSS_Dataframe$Harvest_rate)
        Prev_state_host <- CBSD_CSS_Dataframe$Prev_State[Host] <-  CBSD_CSS_Dataframe$State[Host]
        
        if(Prev_state_host == "S_N"){
          No_S_N <- No_S_N - 1
        } else if(Prev_state_host == "S_C"){
          No_S_C <- No_S_C - 1
        }else  if(Prev_state_host == "I_N"){
          No_I_N <- No_I_N - 1
          No_I <- No_I - 1
        } else if(Prev_state_host == "I_C"){
          No_I_C <- No_I_C - 1
          No_I <- No_I - 1 
          
        } 
        CBSD_CSS_Dataframe <- Strategy_profits(CBSD_CSS_Dataframe, parameters) ## get profits for each individual 
        
        State_host <- CBSD_CSS_Dataframe$State[Host]<- next_strategy_compare_grower_to_strategy(CBSD_CSS_Dataframe, parameters, Host, Prev_state_host, no_infected  = No_I, which_I, kernel_matrix = Kernel_Matrix_estimation)
        
        if(State_host == "S_N" ){
          prev_infected <- (No_I)/parameters$N
          infected <- sample(c(T, F),  1,c(parameters$p*(prev_infected), (1-parameters$p*(prev_infected))), replace = F)
          if(infected){
            CBSD_CSS_Dataframe$State[Host] <- State_host <- "I_N"
            CBSD_CSS_Dataframe$Infection[Host] <- 0
            CBSD_CSS_Dataframe$tI[Host] <- t
            No_I <- No_I + 1
            No_I_N <- No_I_N + 1
            if(Prev_state_host == "I_N" | Prev_state_host == "I_C"){
              which_I <- which_I[!(which_I %in% as.numeric(Host))]
            } else{
              inf_pressure <- (Kernel_matrix[,Host])
              CBSD_CSS_Dataframe$Infection <- CBSD_CSS_Dataframe$Infection + parameters$beta* inf_pressure*(CBSD_CSS_Dataframe$State=="S_C" | CBSD_CSS_Dataframe$State=="S_N" ) 
            }
            which_I <- c(as.numeric(which_I), as.numeric(Host))
            
          } else{
            if(Prev_state_host == "I_N" | Prev_state_host == "I_C"){
              inf_pressure <- (Kernel_matrix[,Host])
              CBSD_CSS_Dataframe$Infection <- CBSD_CSS_Dataframe$Infection - parameters$beta* inf_pressure*(CBSD_CSS_Dataframe$State=="S_C" | CBSD_CSS_Dataframe$State=="S_N" ) 
              which_I <- which_I[!(which_I %in% as.numeric(Host))]
            }
            CBSD_CSS_Dataframe$Infection[Host] <-  0 + parameters$beta*sum(Kernel_matrix[Host,which_I])
            No_S_N <- No_S_N + 1
            
          }
        } 
        if(State_host == "S_C") {
          if(Prev_state_host == "I_N" | Prev_state_host == "I_C"){
            inf_pressure <- (Kernel_matrix[,Host])
            CBSD_CSS_Dataframe$Infection <- CBSD_CSS_Dataframe$Infection - parameters$beta*inf_pressure*(CBSD_CSS_Dataframe$State=="S_C" | CBSD_CSS_Dataframe$State=="S_N" ) 
            which_I <- which_I[!(which_I %in% as.numeric(Host))]
            
          }
          No_S_C <- No_S_C + 1
          CBSD_CSS_Dataframe$Infection[Host] <-  0 + parameters$beta*sum(Kernel_matrix[Host,which_I])
          
        }
        
        
      } else if(Event == "Secondary"){
        Host <- sample(1:parameters$N,1,prob=CBSD_CSS_Dataframe$Infection)
        if(CBSD_CSS_Dataframe$State[Host] == "S_N"){
          No_S_N <- No_S_N - 1
          No_I_N <- No_I_N + 1
          No_I <- No_I + 1
          which_I <- c(as.numeric(which_I), as.numeric(Host))
          
        } else if(CBSD_CSS_Dataframe$State[Host] == "S_C"){
          No_S_C <- No_S_C - 1
          No_I_C <- No_I_C + 1
          No_I <- No_I + 1
          which_I <- c(as.numeric(which_I), as.numeric(Host))
          
        } 
        CBSD_CSS_Dataframe <- epidemic_infection(CBSD_CSS_Dataframe, parameters, Kernel_matrix, Host)
        
      }
      i = i + 1 
      No_S_C_sequence[i,1] <- No_S_C
      No_S_N_sequence[i,1] <- No_S_N
      No_I_N_sequence[i,1] <- No_I_N
      No_I_C_sequence[i,1]<- No_I_C
      Time_sequence[i,1]<- t
      Average_profit[i,1] <- mean(CBSD_CSS_Dataframe$Profit)
      CBSD_CSS_Dataframe$Infection[CBSD_CSS_Dataframe$Infection < 1e-10] <- 0 
      
    }
    if(sum(CBSD_CSS_Dataframe$State == "I_N",CBSD_CSS_Dataframe$State == "I_C" ) == 0){
      Extinction_events[x] = Extinction_events[x]+ 1
    }
    frame <- data.frame(Time_sequence, No_S_C_sequence, No_S_N_sequence, No_I_C_sequence, No_I_N_sequence, Average_profit) # create dataframe of vectors for each time 
    a = 1
    for(m in t_seq){ # loop through time
      if(length(which_I) == 0){
        Extinction_events[x] = Extinction_events[x]+ 1
      }
      z <- max(which(frame$Time_sequence < m))
      if(is.infinite(z)){
        z <- 1
      }# find the maximum row that is still less than time t 
      S_C_avg[a,x] <- frame$No_S_C_sequence[z] # subset from that row to get the no. S_S etc individuals
      S_N_avg[a,x] <- frame$No_S_N_sequence[z]
      I_C_avg[a,x] <- frame$No_I_C_sequence[z]
      I_N_avg[a,x] <- frame$No_I_N_sequence[z]
      average_profit[a,x] <- frame$Average_profit[z]
      a = a+1
    }
    
    
    infection <-  I_C_avg[,x] +  I_N_avg[,x]
    control <- I_C_avg[,x] +  S_C_avg[,x]

  }
  which(Extinction_events == 0) ### need to condition for non-extinction at low beta values 
  S_C_avg_no_extinct <- S_C_avg[which(Extinction_events == 0)]
  I_C_avg_no_extinct <- I_C_avg[,which(Extinction_events == 0)]
  S_N_avg_no_extinct <- S_N_avg[,which(Extinction_events == 0)]
  I_N_avg_no_extinct<- I_N_avg[,which(Extinction_events == 0)]
  if(length(which(Extinction_events != 0)) > 1){
    S_C_avg$mean <- mean(as.matrix(S_C_avg_no_extinct), na.rm=TRUE)
    I_C_avg$mean <- mean(as.matrix(I_C_avg_no_extinct), na.rm=TRUE)
    S_N_avg$mean <- mean(as.matrix(S_N_avg_no_extinct), na.rm=TRUE)
    I_N_avg$mean <- mean(as.matrix(I_N_avg_no_extinct), na.rm=TRUE)
    
    S_C_avg$sd <- sd(as.matrix(S_C_avg_no_extinct), na.rm = T)
    I_C_avg$sd <- sd(as.matrix(I_C_avg_no_extinct), na.rm = T)
    S_N_avg$sd <- sd(as.matrix(S_N_avg_no_extinct), na.rm = T)
    I_N_avg$sd <- sd(as.matrix(I_N_avg_no_extinct), na.rm = T)
  }else{
    S_C_avg$mean <- rowMeans(S_C_avg, na.rm=TRUE)
    I_C_avg$mean <- rowMeans(I_C_avg, na.rm=TRUE)
    S_N_avg$mean <- rowMeans(S_N_avg, na.rm=TRUE)
    I_N_avg$mean <- rowMeans(I_N_avg, na.rm=TRUE)
    S_C_avg$sd <- rowSds(as.matrix(S_C_avg))
    I_C_avg$sd <- rowSds(as.matrix(I_C_avg))
    I_N_avg$sd <- sd(as.matrix(I_N_avg))
    S_N_avg$sd <- rowSds(as.matrix(S_N_avg))
  }
   control_equilibrium <- c(control_equilibrium, S_C_avg$mean + I_C_avg$mean)
  control_sd <- c(control_sd, S_C_avg$sd + I_C_avg$sd)
  infection_sd <- c(infection_sd, I_C_avg$sd + I_N_avg$sd)
  infection_equilibrium <- c(infection_equilibrium, I_N_avg$mean + I_C_avg$mean)
  
  
}
N <- 750
plot(beta_sequence_stoc[1:length(unlist(control_equilibrium))], unlist(control_equilibrium)/N, ty = "l", lwd = 3, col = "#80cc99", ylim = c(0, 1),  ylab = expression("Proportion of growers"), xlab = expression(paste("Rate of horizontal transmission (", beta, "s)")))
lines(beta_sequence_stoc[1:length(unlist(infection_equilibrium))], unlist(infection_equilibrium)/N, lwd = 3, col ="#9980cc")
polygon(c(beta_sequence_stoc,rev(beta_sequence_stoc)), c(control_equilibrium/N + (control_sd)/N , rev(control_equilibrium/N - (control_sd)/N)), col=adjustcolor("#80cc99",alpha.f=0.6), border = NA)
polygon(c(beta_sequence_stoc,rev(beta_sequence_stoc)), c(infection_equilibrium/N + infection_sd/N , rev(infection_equilibrium/N -infection_sd/N)), col=adjustcolor("#9980cc",alpha.f=0.6), border = NA)

title(main = expression("Spatial, stochastic model"), cex.main = 2.5)
legend("right", 
       legend = c(expression("Control"),expression("Infected")), 
       col = c("#80cc99", "#9980cc"), lwd = 4, cex = 1.5 )
par(xpd = T)