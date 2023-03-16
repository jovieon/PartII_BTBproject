library(deSolve)
library(rootSolve)
library(tidyverse)

ss_SI_dyn <- function(t, var, par) {
  Sc <- var[1]
  Sb <- var[2]
  Ic <- var[3]
  Ib <- var[4]
  
  ss_u = par[[1]]
  ss_v = par[[2]]
  ss_removal = par[[3]]
  beta = par[[4]]
  
  Nc <- Sc + Ic
  Nb <- Sb + Ib
  
  S <- matrix(c(Sc, Sb))
  I <- matrix(c(Ic, Ib))
  Sdiag <- diag(c(Sc, Sb))
  IdivN <- matrix(c((Ic/Nc), (Ib/Nb)))
  N <- matrix(c(Nc, Nb))
  
  dS <- (ss_u %*% N) + ss_removal %*% I - (Sdiag %*% beta %*% IdivN) - (ss_v %*% S)  
  dI <- (Sdiag %*% beta %*% IdivN) - ((ss_removal+ss_v) %*% I)
  
  return(list(c(dS, dI)))
  
}

beta_scenarios = function(scenario) {
  
  if (scenario == 1) {
    
    betaCC <- 0.94*(0.7+vc)
    betaBC <- 0.05*(vb)
    betaCB <- 0.2*(0.7+vc)
    betaBB <- 0.99*(vb)
    ## default scenario as argued by EBP & Wood (2015)
    
  }
  
  if (scenario == 2) {
    
    betaCC <- 0.99*(0.7+vc)
    betaBC <- 0.05*(vb)
    betaCB <- 0.05*(0.7+vc)
    betaBB <- 1.04*(vb)
    ## badger reservoir, with low inter-host transmission
    
  }
  
  if (scenario == 3) {
    
    betaCC <- 1.04*(0.7+vc)
    betaBC <- 0.05*(vb)
    betaCB <- 0.05*(0.7+vc)
    betaBB <- 0.99*(vb)
    ## cattle reservoir, with low inter-host transmission
  }
  
  if (scenario == 4) {
    
    betaCC <- 0.94*(0.7+vc)
    betaBC <- 0.14*(vb)
    betaCB <- 0.07*(0.7+vc)
    betaBB <- 0.99*(vb)
    ## van Tonder scenario
    
  }
  
  beta <- matrix(c(betaCC, betaBC, betaCB, betaBB), 2, 2, byrow = TRUE)
  
  return(beta)
  
}

init_prevalence = function(scenario) {
  
  #initial prevalence according to steady-state analysis
  
  if (scenario==1) {
    
    # default scenario
    I0c = 0.01717
    I0b = 0.1064
    
  }
  
  if (scenario==2) {
    
    #badger reservoir scenario
    I0c = 0.029
    I0b = 0.093
    
  }
  
  if (scenario==3) {
    
    #cattle reservoir scenario
    I0c = 0.0576
    I0b = 0.0976
    
  }
  
  if (scenario==4) {
    
    #van Tonder scenario
    I0c = 0.0324
    I0b = 0.0866
    
  }
  
  return(c(I0c, I0b))
}

N0c <- 1
N0b <- 1

scenario = 4
ss_I0c = init_prevalence(scenario)[1]
ss_I0b = init_prevalence(scenario)[2]
  
  uc <- 0.1
  ## birth rate of cattle = birth rate (constant population size)
  ub <- 0.2
  ## birth rate of badgers = birth rate (constant population size)
  
  vc <- 0.1
  ## death rate of cattle = birth rate (constant population size)
  vb <- 0.2
  ## death rate of badgers = birth rate (constant population size)
  
  tau <- 0.7
  ## removal rate of infected cattle
  ## assumes DIVA test has equal efficacy in vaccinated and unvaccinated cattle
  
  beta <- beta_scenarios(scenario)
  ## 2x2 matrix for values of beta not including effects of vaccination
  
  ss_u <- diag(c(uc, ub))
  ss_v <- diag(c(vc, vb))
  ss_removal <- diag(c(tau, 0))
  
  ss_SI.par <- list(ss_u, ss_v, ss_removal, beta)
  
  ss_SI.init <- c(N0c - ss_I0c,
                  N0b - ss_I0b,
                  ss_I0c, ss_I0b)
  
  SI.t <- seq(0, 100, by = 1/364)
  
  ss_SI.sol <- lsoda(ss_SI.init, SI.t, ss_SI_dyn, ss_SI.par)[,c(1,4,5)]
  
  colnames(ss_SI.sol) = c("Time", "Cattle", "Badgers")
  
  df_ss_SI.sol = as.data.frame(ss_SI.sol)
  df_ss_SI.sol = pivot_longer(df_ss_SI.sol,
                              cols = c("Cattle","Badgers"),
                              names_to = "Species",
                              values_to = "Proportion")
  
  ggplot(df_ss_SI.sol, aes(x=Time, y =Proportion, colour = Species))+
    geom_line()+
    labs(x= "Time (years)",
         y = "Proportion infected")

RS = steady(y=ss_SI.init,func=ss_SI_dyn,parms = ss_SI.par,method='runsteady')
RS$y[c(3,4)]