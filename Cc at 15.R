library(tidyverse)
library(deSolve)
library(scales)
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

tau_values <- seq(from = 0.4, to = 1, by = 0.05)
epI_values <- seq(from=0,to=1,by=0.1)

scenario = 4
epSus=0.4 ## 1-direct vaccine efficacy


##Fixed parameters:
pcattle=0.8 ## effective vaccination coverage in cattle
pbadgers=0 ##effective vaccination coverage in badgers is set to 0
uc = 0.1 ## cattle birth rate
vc = 0.1 ## cattle death rate
ub = 0.2 ## badger birth rate
vb = 0.2 ## badger death rate

N0c = 1
N0b = 1
I0c = init_prevalence(scenario)[1]
I0b = init_prevalence(scenario)[2]

SI.init <- c(N0c - I0c, N0b-I0b, ## initial proportion of susceptible cattle and badgers, respectively
             0, 0, ## initial number of vaccinated cattle and badgers (set to 0)
             I0c, I0b, ## initial prevalence in cattle and badgers, respectively
             0, 0, ## initial number of infected vaccinates (cattle and badgers)
             I0c, 0)

SI.t = seq(from=0,to=15,by=1)

u <- diag(c(((1-pcattle)*(uc)), (1-pbadgers)*(ub), pcattle*(uc), pbadgers*(ub)))
v <- diag(c(vc, vb, vc, vb))

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

Cc_at_OTF = data.frame()

for (j in 1:length(tau_values)) {
  
  tau=tau_values[j]
  
  removal <- diag(c((tau), 0, (tau), 0))
  
  results = tau

  for (k in 1:length(epI_values)) {
    
    epInf = epI_values[k]

    vbeta[1,3] <- beta[1,1]*epInf
    vbeta[2,3] <- beta[2,1]*epInf
    
    vbeta[3,3] <- beta[1,1]*epSus*epInf
    vbeta[4,3] <- beta[2,1]*epInf
    
    vbetaCc[1,3] <- epInf*beta[1,1]
    vbetaCc[2,3] <- epSus*epInf*beta[1,1]
    
    ## 4x4 matrix for values of beta (including effects of vaccination)
    
    SI.par <- list(u, v, removal, vbeta, vbetaCc)
    SI.sol <- lsoda(SI.init, SI.t, SI_dyn, SI.par)
    Cc <- SI.sol[,10]+SI.sol[,11]
    
    results = c(results, Cc[length(Cc)])
    
  }
  Cc_at_OTF = rbind(Cc_at_OTF, results)
  ##print(results)
  ##print(eigenK(epSus,epInf,tau,pcattle))
  
}

colnames(Cc_at_OTF) = c("tau", "0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1")

Cc_at_OTF = pivot_longer(Cc_at_OTF,
                               cols = c("0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1"),
                               names_to="epInf",
                               values_to = "Cc")

Cc_at_OTF = Cc_at_OTF %>% na.omit()

Cc_at_15 = c(0.481,0.813,1.6132,0.9074)

Cc_at_OTF[,3] = (Cc_at_OTF[,3])*(33000/I0c)

cp8 = ggplot(Cc_at_OTF,aes(x=tau, y=Cc, colour=epInf))+
  geom_line()+
  scale_colour_discrete(name="Indirect vaccine efficacy",
                        breaks=seq(from=0,to=1,by=0.1),
                        labels=c("100%","90%","80%","70%","60%","50%","40%","30%","20%","10%","0"))+
  labs(y = "Total cases over 15 years",
       x = expression(tau))+
  ggtitle(expression(paste("van Tonder scenario, ", beta[BC], " > ", beta[CB])),
          subtitle = "Direct vaccine efficacy = 60%")+
  geom_vline(xintercept = 0.6, linetype = "dashed")+
  geom_vline(xintercept = 0.8, linetype = "dashed")+
  geom_hline(yintercept = (Cc_at_15[scenario]*33000)/I0c)+
  scale_y_continuous(labels = scientific, limits = c(0,2e+06))

patchworkCc = (cp1|cp2)/(cp3|cp4)/(cp5|cp6)/(cp7|cp8)+
  plot_annotation(tag_levels = 'a')+ 
  theme(plot.tag = element_text(size = 12))