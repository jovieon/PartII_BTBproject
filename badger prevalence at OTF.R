library(tidyverse)
library(deSolve)
library(patchwork)

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
    ##as argued by EBP & Wood (2015)
    
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
    ## intermediate inter-host transmission
    
  }
  
  beta <- matrix(c(betaCC, betaBC, betaCB, betaBB), 2, 2, byrow = TRUE)
  
  return(beta)
  
}

init_prevalence = function(scenario) {
  
  if (scenario==1) {
    
    I0c = 0.01717
    I0b = 0.1064
    
  }
  
  if (scenario==2) {
    
    I0c = 0.029
    I0b = 0.093
    
  }
  
  if (scenario==3) {
    
    I0c = 0.0576
    I0b = 0.0976
    
  }
  
  if (scenario==4) {
    
    I0c = 0.0324
    I0b = 0.0866
    
  }
  
  return(c(I0c, I0b))
}

values <- seq(from = 0, to = 1, by = 0.05)
effective_p_scenarios <- c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

scenario=3
epSus=0.7 ## 1-direct vaccine efficacy

N0c = 1
N0b = 1
##Fixed parameters:
tau = 0.7 ## assumes DIVA test of equal sensitivity to SICCT test
pbadgers=0 ##effective vaccination coverage in badgers is set to 0
uc = 0.1 ## cattle birth rate
vc = 0.1 ## cattle death rate
ub = 0.2 ## badger birth rate
vb = 0.2 ## badger death rate

badger_prevalence = data.frame()

beta = beta_scenarios(scenario)
# 2x2 matrix of transmission rates

vbeta <- matrix(0, 4, 4)

vbeta[1,1] <- beta[1,1]
vbeta[1,2] <- beta[1,2]
vbeta[2,1] <- beta[2,1]
vbeta[2,2] <- beta[2,2]

vbeta[3,1] <- beta[1,1]*epSus
vbeta[3,2] <- beta[1,2]*epSus
vbeta[4,1] <- beta[2,1]
vbeta[4,2] <- beta[2,2]

vbeta[1,4] <- beta[1,2]
vbeta[2,4] <- beta[2,2]
vbeta[3,4] <- beta[1,2]*epSus
vbeta[4,4] <- beta[2,2]

vbetaCc <- matrix(nrow=2,ncol=4)
vbetaCc[1,1] <- beta[1,1]
vbetaCc[1,2] <- beta[1,2]

vbetaCc[1,4] <- beta[1,2]
vbetaCc[2,1] <- epSus*beta[1,1]
vbetaCc[2,2] <- epSus*beta[1,2]

vbetaCc[2,4] <- epSus*beta[1,2]

v <- diag(c(vc, vb, vc, vb))
removal <- diag(c((tau), 0, (tau), 0))

I0c = init_prevalence(scenario)[1]
I0b = init_prevalence(scenario)[2]

SI.init <- c(N0c-I0c, N0b-I0b,
             0,0,
             I0c, I0b,
             0,0,
             I0c,0)

SI.t <- seq(from=0, to=1000, by=1)

for (j in 1:21) {
  
  epInf=values[j]
  
  vbeta[1,3] <- beta[1,1]*epInf
  vbeta[2,3] <- beta[2,1]*epInf
  
  vbeta[3,3] <- beta[1,1]*epSus*epInf
  vbeta[4,3] <- beta[2,1]*epInf
  
  vbetaCc[1,3] <- epInf*beta[1,1]
  vbetaCc[2,3] <- epSus*epInf*beta[1,2]
  
  results = epInf
  
  for (k in 1:9) {
    
    pcattle = effective_p_scenarios[k]
    
    u <- diag(c(((1-pcattle)*(uc)), (1-pbadgers)*(ub), pcattle*(uc), pbadgers*(ub)))
    
    SI.par <- list(u, v, removal, vbeta, vbetaCc)
    SI.sol <- lsoda(SI.init, SI.t, SI_dyn, SI.par)
    
    TIMES <- SI.sol[,1]
    Ib <- SI.sol[,7]+SI.sol[,9]
    Cc <- SI.sol[,10]+SI.sol[,11]
    
    #plot(TIMES, Ib, col = "blue", pch = ".")
    
    if ((diff(c(Cc))[length(Cc)-1]) > 0.00085) {
      
      results = c(results, NA)
      #print("NA")
      
    } else {
      
      results = c(results, Ib[min(which(diff(Cc)<0.00085))])
      #print(Ib[min(which(diff(Cc)<0.00085))])
    }
  }
  
  print(results)
  
  badger_prevalence = rbind(badger_prevalence, results)
  
}

colnames(badger_prevalence) = c("epInf", "20","30","40","50","60","70","80","90","100")

badger_prevalence = pivot_longer(badger_prevalence,
                           cols = c("20","30","40","50","60","70","80","90","100"),
                           names_to="p",
                           values_to = "Prevalence")

badger_prevalence = badger_prevalence %>% na.omit()

bp5 = ggplot(badger_prevalence,aes(x=epInf, y=Prevalence, colour=p))+
  geom_line()+
  scale_colour_discrete(name=expression(rho),
                        breaks=c("100","90","80","70","60","50","40","30","20"),
                        labels=c("100%","90%","80%","70%","60%","50%","40%","30%","20%"))+
  labs(y = "Badger prevalence at OTF",
       x = expression(epsilon[I]))+
  ggtitle("Cattle reservoir",
          subtitle = "Direct vaccine efficacy = 30%")+
  ylim(0,0.071)+
  xlim(0,1)

patchwork4 = (bp1|bp2)/(bp3|bp4)/(bp5|bp6)/(bp7|bp8)+
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12))