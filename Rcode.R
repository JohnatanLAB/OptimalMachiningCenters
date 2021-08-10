library(ggplot2)
library(rstan)

path = "./data.csv"
dataS <- read.csv(path, sep = ",", dec = ".")
str(dataS)
#ITM: Leadwell V-40iT®
#PB: Haas UMC-750

Depth_cut <- ifelse(dataS$depth_of_cut=="1 mm", 1,
                      ifelse(dataS$depth_of_cut=="1.5 mm", 1.5,
                             ifelse(dataS$depth_of_cut=="2 mm", 2,
                                    ifelse(dataS$depth_of_cut=="2.5 mm", 2.5, 3))))
#Transforming the machining centers in dummy variables
MACHINE <- ifelse(dataS$Machine=="PB", 1, 0)

Y <- cbind(log(dataS$ROUGHNESS), log(dataS$POWER))

#Desing matrix X. All regressors are standardized
X <- cbind(1, 
           (Depth_cut - mean(Depth_cut))/sd(Depth_cut),
           (dataS$feed_rate - mean(dataS$feed_rate))/sd(dataS$feed_rate),
           (dataS$spindle_rate -mean(dataS$spindle_rate))/sd(dataS$spindle_rate), 
           ((Depth_cut - mean(Depth_cut))/sd(Depth_cut))^2,
           ((dataS$feed_rate - mean(dataS$feed_rate))/sd(dataS$feed_rate))^2,
           ((dataS$spindle_rate -mean(dataS$spindle_rate))/sd(dataS$spindle_rate))^2, 
           MACHINE, 
           MACHINE*((Depth_cut - mean(Depth_cut))/sd(Depth_cut)), 
           MACHINE*((dataS$feed_rate - mean(dataS$feed_rate))/sd(dataS$feed_rate)),
           MACHINE*((dataS$spindle_rate -mean(dataS$spindle_rate))/sd(dataS$spindle_rate)), 
           ((Depth_cut - mean(Depth_cut))/sd(Depth_cut))*((dataS$feed_rate - mean(dataS$feed_rate))/sd(dataS$feed_rate)), 
           ((Depth_cut - mean(Depth_cut))/sd(Depth_cut))*((dataS$spindle_rate -mean(dataS$spindle_rate))/sd(dataS$spindle_rate)), 
           ((dataS$feed_rate - mean(dataS$feed_rate))/sd(dataS$feed_rate))*((dataS$spindle_rate -mean(dataS$spindle_rate))/sd(dataS$spindle_rate)))


Data.stan <- list(K = nlevels(as.factor(dataS$Machine)), 
                  J = dim(X)[2],
                  N = dim(Y)[1],
                  #N_mis = dim(X[(101:125),])[1],
                  x = X,
                  y = Y
                  #x_mis = X[(101:125),]
)

#Fitting the Bayesian SUR model using STAN
ranIntFit <- stan(file = "~/Downloads/2021/SURmodel.stan", data = Data.stan,
                  iter = 20000, chains = 4, cores = 4, thin = 1, control = list(max_treedepth = 15))

print(ranIntFit)
summary(ranIntFit)


Parameters0 <- extract(ranIntFit, pars = c("beta", "Sigma"))
Parameters <- Parameters0$beta

Beta.power<- Parameters[,2,]
Beta.roughness <- Parameters[,1,]
#HDI for the parameters of the SUR model
apply(Beta.power, 2, hdi)
apply(Beta.roughness, 2, hdi)

#########################
#Optimization process ###
#########################
library("rgenoud")
library("parallel")
f_1 <- function(beta, x){
  y = sum(beta*x)
  return(y)
}


##############################
## ITM: Leadwell V-40iT® #####
##############################

design.ITM <- function(x){
  res_power <- exp(median(apply(Beta.power, 1, f_1, c(1,
                                                          (x[1] - mean(Depth_cut))/sd(Depth_cut),
                                                          (x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate),
                                                          (x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate),
                                                          ((x[1] - mean(Depth_cut))/sd(Depth_cut))^2,
                                                          ((x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate))^2,
                                                          ((x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate))^2,
                                                        
                                                          0, 0, 0, 0,
                                                          ((x[1] - mean(Depth_cut))/sd(Depth_cut))*((x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate)),
                                                          ((x[1] - mean(Depth_cut))/sd(Depth_cut))*((x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate)),
                                                          ((x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate))*((x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate))
  ))))
  res_roughness <- exp(median(apply(Beta.roughness, 1, f_1, c(1,
                                                              (x[1] - mean(Depth_cut))/sd(Depth_cut),
                                                              (x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate),
                                                              (x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate),
                                                              ((x[1] - mean(Depth_cut))/sd(Depth_cut))^2,
                                                              ((x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate))^2,
                                                              ((x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate))^2,
                                                              0, 0, 0, 0,
                                                              ((x[1] - mean(Depth_cut))/sd(Depth_cut))*((x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate)),
                                                              ((x[1] - mean(Depth_cut))/sd(Depth_cut))*((x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate)),
                                                              ((x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate))*((x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate))
  ))))
  
  Y1 <- c(res_roughness, res_power)
  Ymin <- c(0.33, 40) 
  DIST1 <- as.numeric(dist(rbind(Y1, Ymin)))
  return(DIST1)
}

MIN1 <- apply(cbind(Depth_cut, dataS$feed_rate, dataS$spindle_rate), 2, min)
MAX1 <- apply(cbind(Depth_cut, dataS$feed_rate, dataS$spindle_rate), 2, max)
region <- cbind(MIN1, MAX1)

cl <- parallel::makeCluster(4, setup_strategy = "sequential")
clusterExport(cl, c("f_1", "Beta.power", "Beta.roughness", "Depth_cut", "dataS"), )
optim.design.ITM <- genoud(design.ITM, nvars = 3,
                               Domains = region, cluster = cl, boundary.enforcement = 2)

optim.design.ITM


#########################
## PB: Haas UMC-750 #####
#########################

design.PB <- function(x){
  res_power <- exp(median(apply(Beta.power, 1, f_1, c(1,
                                                      (x[1] - mean(Depth_cut))/sd(Depth_cut),
                                                      (x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate),
                                                      (x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate),
                                                      ((x[1] - mean(Depth_cut))/sd(Depth_cut))^2,
                                                      ((x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate))^2,
                                                      ((x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate))^2,
                                                      
                                                      0, 0, 0, 0,
                                                      ((x[1] - mean(Depth_cut))/sd(Depth_cut))*((x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate)),
                                                      ((x[1] - mean(Depth_cut))/sd(Depth_cut))*((x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate)),
                                                      ((x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate))*((x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate))
  ))))
  res_roughness <- exp(median(apply(Beta.roughness, 1, f_1, c(1,
                                                              (x[1] - mean(Depth_cut))/sd(Depth_cut),
                                                              (x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate),
                                                              (x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate),
                                                              ((x[1] - mean(Depth_cut))/sd(Depth_cut))^2,
                                                              ((x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate))^2,
                                                              ((x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate))^2,
                                                              1, 
                                                              (x[1] - mean(Depth_cut))/sd(Depth_cut),
                                                              (x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate),
                                                              (x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate),
                                                              ((x[1] - mean(Depth_cut))/sd(Depth_cut))*((x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate)),
                                                              ((x[1] - mean(Depth_cut))/sd(Depth_cut))*((x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate)),
                                                              ((x[2] - mean(dataS$feed_rate))/sd(dataS$feed_rate))*((x[3] - mean(dataS$spindle_rate))/sd(dataS$spindle_rate))
  ))))
  
  Y1 <- c(res_roughness, res_power)
  Ymin <- c(0.28, 24) 
  DIST1 <- as.numeric(dist(rbind(Y1, Ymin)))
  return(DIST1)
}

MIN1 <- apply(cbind(Depth_cut, dataS$feed_rate, dataS$spindle_rate), 2, min)
MAX1 <- apply(cbind(Depth_cut, dataS$feed_rate, dataS$spindle_rate), 2, max)
region <- cbind(MIN1, MAX1)

cl <- parallel::makeCluster(4, setup_strategy = "sequential")
clusterExport(cl, c("f_1", "Beta.power", "Beta.roughness", "Depth_cut", "dataS"), )
optim.design.PB <- genoud(design.PB, nvars = 3,
                           Domains = region, cluster = cl, boundary.enforcement = 2)

optim.design.PB


