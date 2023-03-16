library(tidyverse)
library(ggplot2)
library(patchwork)

beta_scenarios = function(scenario) {
  
  if (scenario == 1) {
    
    betaCC <- 0.94*(0.7+vc)
    betaBC <- 0.05*(vb)
    betaCB <- 0.2*(0.7+vc)
    betaBB <- 0.99*(vb)
    ##default scenario, as argued by EBP & Wood (2015)
    
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

eigenK <- function(epSus,epInf,tau,pcattle,beta) {
  
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
  
  R <- matrix(0, 4, 4)
  uc = 0.1
  vc = 0.1
  ub = 0.2
  vb = 0.2
  R[1,] <- vbeta[1,]/(tau + vc)
  R[2,] <- vbeta[2,]/(vb)
  R[3,] <- vbeta[3,]/(tau + vc)
  R[4,] <- vbeta[4,]/(vb)
  ## creating NGM by dividing elements of vbeta matrix
  
  pbadgers = 0
  
  matrixP = diag(c((1-pcattle),(1-pbadgers),pcattle,pbadgers))
  
  K = matrixP %*% R
  
  ##print(K)
  
  return(max(Re(eigen(K)$values)))
  ## reproduction ratio, R = dominant eigenvalue of Next Generation Matrix, K
  
}

values = seq(from=0, to=1, by=0.05)

scenario =3
epSus = 0.4 ## 1 - direct vaccine efficacy

##Fixed parameters:
tau = 0.7 ## assumes DIVA test of equal sensitivity to SICCT test
pbadgers=0 ##effective vaccination coverage in badgers is set to 0
uc = 0.1 ## cattle birth rate
vc = 0.1 ## cattle death rate

ub = 0.2 ## badger birth rate
vb = 0.2 ## badger death rate

vary_p_epI_R = data.frame()

for (i in 1:length(values)) {
  
  epInf = values[i]
  
  results = epInf
  print(results)
  
  for (j in 1:length(values)) {
    
    pcattle = values[j]
    
    results = c(results, eigenK(epSus,epInf, tau, pcattle, beta_scenarios(scenario)))
    
  }
  
  vary_p_epI_R = rbind(vary_p_epI_R, results)
  
}

colnames(vary_p_epI_R) = c("epInf", values)
vary_p_epI_R = pivot_longer(vary_p_epI_R, cols = c("0","0.05","0.1","0.15","0.2","0.25","0.3","0.35","0.4","0.45","0.5","0.55","0.6","0.65","0.7","0.75","0.8","0.85","0.9","0.95","1"),
                            names_to = "p",
                            values_to = "R")

rp6 = ggplot(vary_p_epI_R, aes(epInf, p, fill = R)) +
  geom_raster() + scale_fill_binned(type = "viridis", breaks = seq(0.95,1.4, by = 0.01))+
  labs(x = expression(epsilon[I]),
       y = "Effective vaccination coverage")+
  scale_y_discrete(breaks = seq(from=0,to=1,by=0.1),
                   labels = c("0","10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"))+
  ggtitle("Cattle reservoir",
          subtitle = "Direct vaccine efficacy = 60%")

patchworkR2 = (rp1|rp2)/(rp3|rp4)/(rp5|rp6)/(rp7|rp8) +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 12))

rp2 = rp2+
  ggtitle(expression(paste("Default scenario, ", beta[BC], " < ", beta[CB])),
          subtitle = "Direct vaccine efficacy = 60%")
