library(deSolve)

SI_dyn <- function(t, var, par) {
  Sc <- var[1]
  Sb <- var[2]
  Vc <- var[3]
  Vb <- var[4]
  Ic <- var[5]
  Ib <- var[6]
  Wc <- var[7]
  Wb <- var[8]
  CIc <- var[9]
  CWc <- var[10]
  
  u = par[[1]]
  v = par[[2]]
  removal = par[[3]]
  vbeta = par[[4]]
  vbetaCc = par[[5]]
  
  Nc <- Sc + Vc + Ic + Wc
  ##print(Nc)
  Nb <- Sb + Vb + Ib + Wb
  ##print(Nb)
  
  S <- matrix(c(Sc, Sb, Vc, Vb))
  I <- matrix(c(Ic, Ib, Wc, Wb))
  Sdiag <- diag(c(Sc, Sb, Vc, Vb))
  IdivN <- matrix(c((Ic/Nc), (Ib/Nb), (Wc/Nc), (Wb/Nb)))
  N <- matrix(c(Nc, Nb, Nc, Nb))
  
  ScVcdiag <- diag(c(Sc,Vc))

  #print(N)
  #print(u)
  
  dS <- (u %*% N) - (v %*% S) - (Sdiag %*% vbeta %*% IdivN) + removal %*% I
  dI <- (Sdiag %*% vbeta %*% IdivN) - ((removal+v) %*% I)
  dC <- (ScVcdiag %*% vbetaCc %*% IdivN) + matrix(c(Ic,Wc))
  
  return(list(c(dS, dI, dC)))
  
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

scenario = 2

pcattle <- 0.8
## proportion of cattle vaccinated
pbadgers <- 0
## proportion of badgers vaccinated
N0c <- 1
N0b <- 1
I0c <- init_prevalence(scenario)[1]
I0b = init_prevalence(scenario)[2]

epSus <- 0.4
## relative risk of infection (reduction in susceptibility due to vacc = 0.25)
epInf <- 0.75
## relative risk of transmission (reduction in infectiousness due to vacc = 0.36)

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
vbeta <- matrix(0, 4, 4)

vbeta[1,1] <- beta[1,1]
vbeta[1,2] <- beta[1,2]
vbeta[2,1] <- beta[2,1]
vbeta[2,2] <- beta[2,2]

vbeta[3,1] <- beta[1,1]*epSus
vbeta[3,2] <- beta[1,2]*epSus
vbeta[4,1] <- beta[2,1]
vbeta[4,2] <- beta[2,2]

vbeta[1,3] <- beta[1,1]*epInf
vbeta[1,4] <- beta[1,2]
vbeta[2,3] <- beta[2,1]*epInf
vbeta[2,4] <- beta[2,2]

vbeta[3,3] <- beta[1,1]*epSus*epInf
vbeta[3,4] <- beta[1,2]*epSus
vbeta[4,3] <- beta[2,1]*epInf
vbeta[4,4] <- beta[2,2]
## 4x4 matrix for values of beta (including effects of vaccination)

vbetaCc <- matrix(nrow=2,ncol=4)
vbetaCc[1,1] <- beta[1,1]
vbetaCc[1,2] <- beta[1,2]
vbetaCc[1,3] <- epInf*beta[1,1]
vbetaCc[1,4] <- beta[1,2]
vbetaCc[2,1] <- epSus*beta[1,1]
vbetaCc[2,2] <- epSus*beta[1,2]
vbetaCc[2,3] <- epSus*epInf*beta[1,1]
vbetaCc[2,4] <- epSus*beta[1,2]

u <- diag(c(((1-pcattle)*(uc)), (1-pbadgers)*(ub), pcattle*(uc), pbadgers*(ub)))
v <- diag(c(vc, vb, vc, vb))
removal <- diag(c((tau), 0, (tau), 0))

SI.par <- list(u, v, removal, vbeta, vbetaCc)

SI.init <- c(N0c - I0c, N0b - I0b,
             0, 0,
             I0c, I0b,
             0, 0,
             I0c, 0)

SI.t <- seq(0, 400, by = 1)
## since beta and other parameters are measured p/a, timescale should be in years
## e.g. 50 years

SI.sol <- lsoda(SI.init, SI.t, SI_dyn, SI.par)

TIMES <- SI.sol[,1]
Sc <- SI.sol[,2]
Sb <- SI.sol[,3]
Vc <- SI.sol[,4]
Vb <- SI.sol[,5]
Ic <- SI.sol[,6]
Ib <- SI.sol[,7]
Wc <- SI.sol[,8]
Wb <- SI.sol[,9]
Cc <- SI.sol[,10]+SI.sol[,11]

plot(TIMES, Sc, col = "blue",
     pch = ".", 
     ylab = "Proportion of Individuals", 
     xlab = "Time (years)", ylim=c(0,1))

lines(TIMES, Sc, col = "blue")

lines(TIMES, Sb, col = "purple")

lines(TIMES, Vc, col = "green")

lines(TIMES, Vb, col = "brown")

lines(TIMES, Ic, col = "red")

lines(TIMES, Ib, col = "pink")

lines(TIMES, Wc, col = "yellow")

lines(TIMES, Wb, col = "grey")

legend("center",
       legend=c("Sc", "Sb", "Vc","Vb","Ic","Ib","Wc","Wb"),
       col=c("blue","purple","green","brown","red","pink","yellow","grey"),
       pch=c("-","-","-","-","-","-","-","-"),
       box.lty = 1,
       cex=0.85,
       ncol=4)

#CUMULATIVE CASES PLOT
plot(TIMES, Cc, col = "red",
     pch = ".", 
     ylab = "Cumulative cases in cattle (Cc)", 
     xlab = "Time (years)",
     ylim=c(0,Cc[length(Cc)]))

lines(TIMES, Cc, col = "red")

for (i in 1:length(Cc)) {
  
  if (diff(c(Cc))[i] < 0.0001) {
    
    #print(SI.t[i]/(1/364))
    #prints number of days to incidence > 0.01%
    
    print(SI.t[i])
    
    break
    
  }
  
}