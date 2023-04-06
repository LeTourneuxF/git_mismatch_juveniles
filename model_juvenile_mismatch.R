
if(file.exists("C:\\Users\\Laboratoire\\Documents\\FL")){
  
  setwd('~/FL/juv_mismatch/joint_analysis')
  wd.save<-'~/FL/juv_mismatch/joint_analysis/parallel_saving/'
  
  
}else if(file.exists("C:\\Users\\Frederic Letourneaux\\Documents")){
  
  setwd("~/Doc/juvmismatch")
  wd.save<-'~/Doc/juvmismatch/test_parallel_saving/'
  
  ######ROGER: tu peux ajouter un wd pour ton ordinateur ici, et wd.save est l'endroit où le modèle enregistre les résultats des 3 chaines 
}else if(file.exists("~/WD_de_ROGER")){
  
  #Ajouter le chemin ici aussi, le reste marche tout seul
  setwd("~/WD_de_ROGER/juvmismatch")
  
  
  wd.save<-paste0(getwd(),'/save_results_parallel/')
  
  if(file.exists(paste0(wd.save,'result_chain1'))==F){
    
    dir.create(path = paste0(wd.save,"result_chain1"), recursive = T)
    dir.create(path = paste0(wd.save,"result_chain2"), recursive = T)
    dir.create(path = paste0(wd.save,"result_chain3"), recursive = T)
    
  }
}

#Simulated parameters:
#phiA <- 2        #Intercept survival adults
#phiJ <- 0.35     #Intercept survival young

#Random year effect on survival 
#sd_eps_t_phiA <-0.25
#sd_eps_t_phiJ <-0.35

#Hunt effect on survival
#beta.hunt.ad<- -0.62
#beta.hunt.juv<- -0.3

#s_beta.ms<- -0.045 #beta mismatch on Sj; but this is not the effect when scaling mismatch values!!

#Heterogeneity effect
#beta.p.low<- -1    # Difference between low and high capture probs (real scale)
#p.high <-0.1       # mean capture probabilities for adults (high capture prob group, probability scale)

#sd_eps_t_p<-0.4    # Standard deviation year random effect on capture prob

#s_piaH<-0.4        # proportion of adults in high observability group

#Transitions of juveniles to adults with high/lo/no capturability
#s_pjH<-0.3          # proportion of surviving young that are in the high recapture group as adults
#s_pjE<-0.5          # proportion of surviving young to emigrate
#s_pjL = 1-(pjH+pjE) # proportion of surviving young that are in the low recapture group as adults


#Recoveries:
#int.r <- qlogis(0.3)
#r.eps.t <-rnorm(n=n.occasions-1, mean=0, sd=0.4)
#r <- plogis(int.r + r.eps.t)


library(nimble)
library(nimbleEcology)


#################### USEFUL FUNCTIONS ########################

#Functions to get and set model and MCMC state and variables
getStateVariableNames <- function(samplerDef) {
  resetMethod <- body(samplerDef$reset)
  stateVars <- character()
  if(resetMethod[[1]] != '{') stop('something wrong')
  numLines <- length(resetMethod)
  for(i in 1:numLines) {
    if(i == 1) next
    thisLine <- resetMethod[[i]]
    if(thisLine[[1]] == '<<-') {
      LHS <- thisLine[[2]]
      if(!is.name(LHS)) stop('haven\'t dealt with non-name-LHS case yet')
      stateVars <- c(stateVars, as.character(LHS))
    }
    if('my_calcAdaptationFactor' %in% all.names(thisLine)) {
      stateVars <- c(stateVars, 'my_calcAdaptationFactor')
    }
  }
  setupMethod <- body(samplerDef$setup)
  if('empirSamp' %in% all.names(setupMethod)) stateVars <- c(stateVars, 'empirSamp')
  return(stateVars)
}

getModelState <- function(model) {
  modelVarNames <- model$getVarNames()
  modelVarValuesList <- vector('list', length(modelVarNames))
  names(modelVarValuesList) <- modelVarNames
  for(var in modelVarNames) {
    modelVarValuesList[[var]] <- model[[var]]
  }
  return(modelVarValuesList)
}

getMCMCstate <- function(conf, mcmc) {
  stateVarNamesList <- vector('list', length(conf$samplerConfs))
  mcmcStateValuesList <- vector('list', length(conf$samplerConfs))
  for(i in seq_along(conf$samplerConfs)) {
    samplerDef <- conf$getSamplerDefinition(i)
    theseStateNames <- getStateVariableNames(samplerDef)
    theseStateValuesList <- vector('list', length(theseStateNames))
    names(theseStateValuesList) <- theseStateNames
    for(j in seq_along(theseStateNames)) {
      if(is.nf(mcmc)) {
        if(theseStateNames[j] == 'my_calcAdaptationFactor') {
          theseStateValuesList[[j]] <- list(timesAdapted = mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]$timesAdapted,
                                            gamma1 = mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]$gamma1)
        } else
          theseStateValuesList[[j]] <- mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]
      }
      if(is.Cnf(mcmc)) {
        if(theseStateNames[j] == 'my_calcAdaptationFactor') {
          theseStateValuesList[[j]] <- list(timesAdapted = valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'timesAdapted'),
                                            gamma1 = valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'gamma1'))
        } else
          theseStateValuesList[[j]] <- valueInCompiledNimbleFunction(mcmc$samplerFunctions[[i]], theseStateNames[j])
      }
    }
    mcmcStateValuesList[[i]] <- theseStateValuesList
  }
  return(mcmcStateValuesList)
}

#Function for 'not in'
'%!in%' <- function(x,y)!('%in%'(x,y))

lu<-function(x){length(unique(x))}



n.chains<-3
n.thin<-3       #Thinning des paramètres
n.thin2<-30     #Thinning des prédictions de l'âge car ça prend trop de mémoire
n.iter = 1111
n.burnin = 0


#sim.data<-readRDS('~/Doc/juvmismatch/data_sim_joint_analysis_54k_27occ.RDS')
sim.data<-readRDS("./data_sim_joint_analysis_2_54k_27occ.RDS")

#Information to extract from dataset to feed to model
s_ms<-sim.data$mismatch                          # Individual mismatch
CH<-sim.data$CH                                  # Capture history
CH.TRUE<-sim.data$CH.TRUE                        # 'true' capture history (perfect detection)
samp.mask<-sim.data$samp.mask                    # Individuals for which we masked the age and mismatch value to be estimated
df_sim_juvs<-as.data.frame(sim.data$df_sim_juvs) # Simulated dataset of juveniles with mismatch, age and feather length
samp_mean_covars<-sim.data$samp_mean_covars      # Data used for simulation (peak nitrogen to calculate mismatch, and mean mismatch from which to simulate)
n.occasions <- dim(CH)[2]                        # Number of occasions      
period1<-sim.data$period1                        # Time index for period 1
period2<-sim.data$period2                        # Time index for period 2

#Number of individuals in capture-recapture dataset
nind<-dim(CH)[1]

# Compute time of first capture for each individual
get.first <- function(x) min(which(x!=0))
f <- apply(CH, 1, get.first)

#Get initial latent state (z) for each individual
get.z.y<-function(x) x[min(which(x!=0))]
IS<-apply(CH.TRUE,1,get.z.y) 

#Get ages of individuals
age<-IS
age[age>1]<-2 #individuals in initial state 2 or 3 are adults, coded as age=2, no individual initially marked in state 4 (only juveniles that emigrate)

#Get position of juveniles in dataset
juv <- which(age==1)

#Get position of adults in dataset
ad  <- which(age==2)

#Mismatch data
ind.mismatch<-s_ms

#Remove mismatch info for some individuals
ind.mismatch[samp.mask]<-NA                 

# Determine the position of juveniles with known mismatch and age
ID.ms.kwn <- which(is.na(ind.mismatch)==F)
# Determine the position of juveniles with unknown mismatch and age
ID.ms.unk<-c(which(is.na(ind.mismatch)),ad)


# Recode CH matrix: 0 is not allowed!
# 1 = alive and in study area, 2 = recovered dead, 3 = not seen or recovered
rCH <- CH # Recoded CH
rCH[rCH==0] <- 3


#Compute mean annual mismatch (for predictions of mean annual juvenile survival)
average_yearly_mismatch<-rep(NA, times=max(unique(f)))
for(t in unique(f)){
  ID_yr_t<-which(f[ID.ms.kwn]==t)
  average_yearly_mismatch[t]<-mean(ind.mismatch[ID.ms.kwn][ID_yr_t])
}

#Check simulated data resembles input values
cbind(average_yearly_mismatch,samp_mean_covars$mean_mismatch[-n.occasions])




#### Discretization of mismatch values #####

#Define a mean and sd for scaling mismatch
scale.ms.mean<-12.40132
scale.ms.sd<-6.25861

#Scale mean average mismatch for predictions of mean annual juvenile survival
scaled_ms_average<-(average_yearly_mismatch-scale.ms.mean)/scale.ms.sd

#Create some potential integer mismatch values (we go from -35 to +75 to be certain that no random draw from the model falls out of range)
pot.ms.vals<-seq(from=-35, to=75)

#Scale the potential values with pre-defined mean and sd
scaled.pot.ms<-(pot.ms.vals-scale.ms.mean)/scale.ms.sd

#Create index for each of those values
index.ms.values<-pot.ms.vals-min(pot.ms.vals)+1

#Create vector to classify individual mismatch values (ranging from -35 to 75) into bins of 3-5 days
#This vector is used in the model to reclass mismatch values in 1 of 11 classes
ms.class.vec<-rep(NA, times=length(scaled.pot.ms))

#Because mismatch values < -5 are never observed in the dataset, all mismatch between -35 and -6 classified into the first bin
ms.class.vec[1:30]<-1

#Classify the range of possible mismatch by hand
range.possible.mismatch<-c(-5:33)
df.range.possible.mismatch<- data.frame(real.values=range.possible.mismatch,
                                        reclass=c(2,2,2,2,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,7,8,8,8,8,8,9,9,9,9,9,10,10,10,10,10))

print(df.range.possible.mismatch)

#Associate the right classes to the classifying vector
ms.class.vec[which(pot.ms.vals %in% range.possible.mismatch)]<-df.range.possible.mismatch$reclass

#Beyond 33 days of mismatch (never observed in data set), put everything in the 11th bin
ms.class.vec[(which(pot.ms.vals == max(range.possible.mismatch))+1):length(ms.class.vec)]<-11


#Compute mean of scaled ms value for each mismatch bin, to feed to the model as the mismatch data value for each mismatch class
info.bin.ms<-data.frame(real.ms=pot.ms.vals,scaled_ms=scaled.pot.ms, reclass.val=ms.class.vec)

average.bin.ms.scaled<-rep(NA, times=length(unique(ms.class.vec)))

for(i in unique(ms.class.vec)){
  average.bin.ms.scaled[i]<-mean(info.bin.ms$scaled_ms[which(info.bin.ms$reclass.val==i)])
}

# mean of mcmc samples for mismatch values >33 = 36, so we use scaled value of mismatch = 36 to plug in for all mismatches >33
# Similarly, we use mismatch for -7 for all mismatch values < -5
average.bin.ms.scaled[c(1,11)]<-c(-3.10,3.77)

#Make sure all potential values correspond to an index
length(index.ms.values)==length(scaled.pot.ms)

#Vector containing a sequence of scaled mismatch values that correspond to the observed mismatch values in the data
#Used to make predictions of the mismatch effect
range.mismatch.yr<-seq(from=(min(s_ms)-scale.ms.mean)/scale.ms.sd, to=(max(s_ms)-scale.ms.mean)/scale.ms.sd, length.out=100)

#Index for gamma matrix, where all adult individuals share the same transition matrix slice for a given time
gi.ad<-length(unique(ms.class.vec))+1




#Data, monitors and constants

#Create object to modify with dataframe containing mismatch information for all juveniles
wbtg.dat.miss<-df_sim_juvs

#Hide real age and mismatch value for individuals determined by samp.mask
wbtg.dat.miss$age.at.B[samp.mask]<-NA
wbtg.dat.miss$mismatch.all[samp.mask]<-NA

#Get position in data of juveniles with real age known
b.age.kwn.vec<-which(is.na(wbtg.dat.miss$age.at.B)==F)

#Get position in data of juveniles with real age unknown
b.age.unk.vec<-which(is.na(wbtg.dat.miss$age.at.B))

#Sample of birds for which to keep age and mismatch estimations to compare with real values
age.keep<-sort(sample(b.age.unk.vec, size=round(min(100,length(b.age.unk.vec)/2)), replace=F))


#Data
known_ages <- wbtg.dat.miss$age.at.B[-c(samp.mask)] #Vector with ages of birds with known ages (data)
prim9_lghts <- wbtg.dat.miss$prim9.B                #Vector with length of 9th primary feather
scaled_prim9_lghts <- as.vector(scale(prim9_lghts)) #Scaled length of 9th primary
mismatch.all<-wbtg.dat.miss$mismatch.all            #Real mismatch value for all juveniles
mismatch.known<-mismatch.all[-c(samp.mask)]         #Real mismatch value for only juveniles of known age (not masked)

#Pack data and constants
#Data
my.data.joint <- list(age.at.B = known_ages,
                      prim9 = as.vector(scale(prim9_lghts)),
                      mismatch.known=mismatch.known,
                      peakN=samp_mean_covars$peakN,                                   #Value of peak plant quality (annual) to compute mismatch with hatch date
                      julian_date_B_unk = wbtg.dat.miss$julian.date.B[b.age.unk.vec], #Date of banding of all juveniles from which to substract age to obtain hatch date
                      ms.class.vec = ms.class.vec,                                    #Index to reclassify mismatch values into 11 categories
                      y=rCH)                                                          #CR data

#Constants
my.constants.joint <- list(
  #### FOR AGE MODEL ##### #
  n.wbtg=length(b.age.kwn.vec),                    # Number of birds with known age
  n.not.wbtg=length(b.age.unk.vec),                # Number of birds with unknown age
  b.age.kwn.vec=b.age.kwn.vec,                     # Position (ID) of birds with known age
  b.age.unk.vec=b.age.unk.vec,                     # Position (ID) of birds with unknown age
  yr.kwn = (wbtg.dat.miss$an_B[b.age.kwn.vec])-(min(wbtg.dat.miss$an_B[b.age.kwn.vec])-1), # year of capture of each bird of known age
  yr.unk = (wbtg.dat.miss$an_B[b.age.unk.vec])-(min(wbtg.dat.miss$an_B[b.age.unk.vec])-1), # year of capture of each bird of unknown age
  
  n.yrs = (length(unique(wbtg.dat.miss$an_B))),    # Number of years for estimation of age
  age.keep = age.keep,                             #ID of juveniles whose age we want to monitor
  
  #### FOR SURVIVAL MODEL ##### #
  #Numbers to feed to for-loops of different lengths
  n.phunt=1,                                  # nombre de périodes de chasse (1 seule pour l'instant)
  n.occasions.p1=length(period1),             # nombre d'occasions avant changements de règlements de chasse (i.e. de période 1)
  n.juv=length(juv),                          # Total number of juveniles
  n.age.save = length(age.keep),              # Number of juveniles for which we want to monitor the estimated age (to compare with original data)
  n.ms.class = length(unique(ms.class.vec)),  # Number of classes for mismatch effect estimation
  lpred=length(range.mismatch.yr),            # Number of values to predict mismatch effect in derived parameters
  
  #Time indices
  f = f,                                  # occasion of marking (first capture) for all individuals
  n.occasions = dim(rCH)[2],              # Number of capture occasions
  
  hunt=c(rep(0, times=length(period1)),   # index du beta.p.hunt en vigueur à l'occasion en cours (vaut 1 partout pour modèle test, et 1 ou 2 pour modèle à 3 périodes), 
         rep(1, times=length(period2))),  # longueur = n.occasions; vaut 0 pour période avant règlements de chasse et ensuite 1 (et ensuite 2 pour vrai modèle)
  
  #ID indices
  nind = dim(rCH)[1],                     #Total number of individuals
  gi.ad=gi.ad,                            #Index for transition matrix (all adults share the matrix slice indexed at n.juv+1)
  
  #Other constants:
  min.ms = min(pot.ms.vals),               # Smallest mismatch value to substract from estimated mismatch, so that index of mismatch classes can be used as an index (go from 1 to 11)
  mismatch.scaled = average.bin.ms.scaled, # Scaled mismatch values to put into the corresponding transition matrix slices
  scale.ms.mean=scale.ms.mean,             # Mean mismatch value used for scaling  
  scale.ms.sd=scale.ms.sd,                 # Standard deviation of mismatch used for scaling
  age=age,                                 # Age class for CR model (1: juvenile, 2: adult)
  range.mismatch=range.mismatch.yr,        # Range of mismatch values used for predictions of mismatch effect
  mean.ms.t=scaled_ms_average,             # average scaled mismatch on annual basis to predict average annual juvenile survival
)

#Parameters to track
my.parameters.joint<- c(##### AGE MODEL #### #
                        'int_age','beta_prim9','sd.prim','sd.year.prim', # Age prediction model parameters
                        
                        ##### SURVIVAL MODEL #### #
                        "phi.A", 'sigma.phiA', 'eps.phiA.t',  # Adult survival intercept and random year effect
                        "phi.J", 'sigma.phiJ', 'eps.phiJ.t',  # Juvenile survival intercept and random year effect
                        "int.r", 'sigma.r', 'eps.r.t',        # Recovery prob intercept and random year effect
                        'int.p', 'sigma.p', 'eps.p.t',        # Capture prob intercept and random year effect
                        "beta.p.L",                           # Beta for effect of low capture prob vs. high capture prob
                        "beta.ms",                            # Beta for effect of mismatch on juvenile survival
                        'p.low', 'p.high', 'piaH',            # high and low capture probabilities and heterogeneity for adult capture
                        'lpjH', 'lpjE', 'pjH', 'pjE', 'pjL',  # parameters for transition of juveniles into adult states with heterogeneous capture
                        'Sj.t','Sa',                          # Annual survival predictions for juveniles and adults
                        'Sa.mean.p1', 'Sa.mean.p2',           # Average survival for adults for each period
                        'Sj.mean.p1', 'Sj.mean.p2',           # Average survival for juveniles for each period 
                        'mean_ms',                            # Average estimated mismatch for each time step (to predict annual juvenile survival)
                        'beta.hunt.ad','beta.hunt.juv',       # Hunt effect for 2nd period for adults and juveniles
                        'r',                                  # Predicted recovery values
                        'S.pred.ms')                          # Juvenile survival prediction against range of observed mismatch val

my.parameters.joint.thin2<- c('ms0.keep')                     # Mismatch class for sample of individuals (with different thinning to save memory)

#Initial values
the.inits.joint <- function(){list(
  
  #### AGE MODEL ##### #
  int_age = runif(1, 0, 50),
  beta_prim9 = runif(1, 0, 10),
  sd.prim = runif(1,0,4),
  sd.year.prim = runif(1,0,4),
  age.at.B.est = rep(30, times=length(b.age.unk.vec)),
  
  ##### SURVIVAL MODEL ##### #
  phi.A= runif(1, -3, 3),
  phi.J= runif(1, -3, 3),
  beta.p.L = runif(1, -2, 0),
  beta.ms = runif(1,-2,2),
  beta.hunt.ad = runif(1,-2,2),
  beta.hunt.juv = runif(1,-2,2),
  int.p = runif(1,-5,5),
  int.r = runif(1, -5, 5),
  sigma.phiA = runif(1,0,3),
  sigma.phiJ = runif(1,0,3),
  sigma.p = runif(1,0,3),
  sigma.r = runif(1,0,3),
  piaH = runif(1,0,1),
  lpjH = runif(1,-5,5),
  lpjE = runif(1,-5,5))}

###################################################################### #
######################### MODEL DEFINITION ###########################
###################################################################### #

joint_mismatch_surv_analysis_optz<-nimbleCode({
  
  #################### AGE MODEL ######################
  
  #prior on beta 9th primary to determine age
  beta_prim9~dunif(0,50)  # beta of feather to explain age
  int_age~dunif(0,50)     # intercept for age estimation
  sd.prim~dunif(0,5)      # residual variance for age estimation from feather
  sd.year.prim~dunif(0,5) # variance for annual random effect of year on age estimation from feather
  
  #Random year intercepts on age estimation from feather
  for(t in 1:n.yrs){
    eps.ageB.t[t] ~ dnorm(mean=0, sd=sd.year.prim)
  }
  
  
  #Estimate beta feather and year intercepts for juveniles of known age:
  for(i in 1:n.wbtg){
    
    # linear predictor for relationship between age and 9prim lgth
    mu[i] <- int_age + beta_prim9 * prim9[b.age.kwn.vec[i]] + eps.ageB.t[yr.kwn[i]]
    
    #Likelihood
    age.at.B[i] ~ dnorm(mu[i], sd=sd.prim) 
    
    #Ceci pourrait être calculé à l'extérieur du modèle
    #Mismatch class calculated from known mismatch values for birds of known age
    ms0[b.age.kwn.vec[i]]<- mismatch.known[i] - min.ms + 1 # '-min.ms + 1' to have indices that go from 1 to 100
    
  }
  
  
  #Predict mismatch from feather length for juveniles of unknown age
  for(i in 1:n.not.wbtg){
    
    #Linear predictor with parameters calculated in previous loop
    mu.pred[i]<- int_age + beta_prim9 * prim9[b.age.unk.vec[i]] + eps.ageB.t[yr.unk[i]]# Linear predictor for age of bird i based on its prim9 lgth
    
    #Age at banding prediction
    age.at.B.est[i] ~ dnorm(mu.pred[i], sd=sd.prim) # draw age of bird from normal distribution around linear predictor mean 
    
    #Mismatch class calculated from predicted age for birds of unknown age
    ms0[b.age.unk.vec[i]] <- round((julian_date_B_unk[i]-age.at.B.est[i]) - peakN[yr.unk[i]]) - min.ms + 1
    
    # basically: (date of capture - age at banding) = hatching date
    # then: (date of capture - age at banding) - peakN[yr] = mismatch in number days (rounded to have integer values)
    # then: round((date of capture - age at banding) - peakN[yr]) - min.ms + 1 = mismatch brought back to go between 1 and 100
    
  }
  
  for(i in 1:n.juv){
    ms.class[i]<-ms.class.vec[ms0[i]]  # Classify mismatch in reduced number of classes
  }
  
  #Keep a sample of estimated ages to see how well they are estimated but not the whole dataset because it takes too much memory
  for(i in 1:n.age.save){
    
    ms0.keep[i]<-ms0[age.keep[i]]
    
  }
  
  ############################################################################################# #
  ########################################## SURVIVAL MODEL  #################################### 
  ############################################################################################# #
  
  # ------------------------------------------------
  # Parameters:
  # s: true survival probability
  # r: recovery probability
  # p: recapture probability
  # ------------------------------------------------
  # States (S):
  # 1 alive young
  # 2 alive adult high capture probability
  # 3 alive adult low capture probability
  # 4 alive adult emigrated, capture probability = 0
  # 5 recently dead and recovered
  # 6 recently dead, but not recovered, and dead (absorbing)
  
  # Observations (O):
  # 1 seen alive
  # 2 recovered dead (equal to 1 as recovery probability is coded in the transition matrix)
  # 3 neither seen nor recovered
  # ------------------------------------------------
  
  
  ########################################################## #
  ###################  CONSTRAINTS  ########################
  ########################################################## #
  
  #Constraints for survival:
  
  #First period: no beta.hunt (so it is not sampled uselessly)
  for (t in 1:n.occasions.p1){
    logit(Sa[t]) <- phi.A + eps.phiA.t[t] # Annual adult survival for occasions in period 1
  }
  
  #2nd (and 3rd) periods: go from last occasion of first period +1 until end of dataset
  for (t in (n.occasions.p1+1):(n.occasions-1)){
    logit(Sa[t]) <- phi.A + eps.phiA.t[t] + beta.hunt.ad[hunt[t]] # Annual adult survival for occasions after period 1
  }
  
  
  
  #Index juvenile survival by mismatch class rather than by individual
  for(c in 1:n.ms.class){
    for (t in 1:n.occasions.p1){#période sans chasse
      logit(Sj[c,t])<- phi.J + beta.ms * mismatch.scaled[c] + eps.phiJ.t[t] # Annual juvenile survival for occasions in period 1
    }
    
    for (t in (n.occasions.p1+1):(n.occasions-1)){#période(s) avec chasse
      logit(Sj[c,t])<- phi.J + beta.ms * mismatch.scaled[c] + beta.hunt.juv[hunt[t]] + eps.phiJ.t[t] # Annual juvenile survival for occasions after period 1
    }
  }
  
  
  #Constraints for events
  for (t in 2:n.occasions){
    #Capture probability (here variable through time (random effect))
    logit(p.low[t]) <- int.p + beta.p.L + eps.p.t[t]   # P.cap for low
    logit(p.high[t]) <- int.p + eps.p.t[t]             # P.cap for high
    
    #Old parameterization with beta.p.L on 0-1 scale:
    #p.low[t] <- p.high[t] * beta.p.L 
    
    #Recovery probability 
    logit(r[t]) <- int.r + eps.r.t[t]               
  }
  
  ########################################################## #
  #####################  PRIORS  ###########################
  ########################################################## #

  ## heterogeneity parameters ##
  piaH ~ dunif(0, 1) #Probability of adult being in high recapture group (initial states)
  
  # multinomial logit for young transitions towards alive and adult with high-p, Low-p, Emigration (transition)
  # normal priors on logit of all but one transition probs
  lpjH ~ dnorm(0, sd = 1.6) #High-p
  lpjE ~ dnorm(0, sd = 1.6) #Emigration
  
  # constrain the transitions such that their sum is < 1
  pjH <- exp(lpjH) / (1 + exp(lpjH) + exp(lpjE)) #High-p
  pjE <- exp(lpjE) / (1 + exp(lpjH) + exp(lpjE)) #Emigration
  
  # last transition probability
  pjL <- 1 - pjH - pjE                          #Low-p, constrain so that all 3 add up to 1
  
  
  #Betas
  beta.p.L ~ dunif(-2.5,0) #Prior for difference between high and low recap probability, constrained to be negative so that p.low<p.high
  #Old parameterization on probability scale:
  #beta.p.L ~ dnorm(0,sd=0.8) #Prior for mean effect of heterogeneity on p
  
  #mismatch effect
  beta.ms ~ dnorm(0, sd=1.6)
  
  #hunt effect
  for(ph in 1:n.phunt){
    beta.hunt.ad[ph] ~ dnorm(0, sd=1.6)  # Hunt effect for adults
    beta.hunt.juv[ph] ~ dnorm(0, sd=1.6) # Hunt effect for juveniles
  }
  
  #CR parameters
  phi.A ~ dnorm(0, sd=1.6) # Intercept survival adults
  sigma.phiA ~ dunif(0,5)  # Variance survival adults for random year effect
  
  phi.J ~ dnorm(0, sd=1.6) # Intercept survival juveniles
  sigma.phiJ ~ dunif(0,5)  # Variance survival juveniles for random year effect
  
  int.p ~ dnorm (0, sd=1.6) # Intercept recap
  sigma.p ~ dunif(0,5)     # Variance recap for random year effect
  
  int.r ~ dnorm(0, sd=1.6)  # Intercept recoveries
  sigma.r ~ dunif(0,5)     # Variance recoveries for random year effect
  
  # Random time intercepts
  for(t in 2:n.occasions){
    eps.phiA.t[t-1] ~ dnorm (0, sd=sigma.phiA)    # Adult survival random year intercept
    eps.phiJ.t[t-1] ~ dnorm (0, sd=sigma.phiJ)    # Juvenile survival random year intercept
    
    eps.p.t[t] ~ dnorm (0, sd=sigma.p)         # recapture random year intercept
    eps.r.t[t] ~ dnorm (0, sd=sigma.r)         # recovery random year intercept
  }
  
  
  ########################################################## #
  #####################  MATRICES  #########################
  ########################################################## #
  
  ###################  Initial states  ######################
  
  #1st row: juveniles
  delta[1,1] <- 1            #Probability of being young upon first capture = 1 
  delta[1,2] <- 0            #Probability of being adult high-p upon first capture = 0
  delta[1,3] <- 0            #Probability of being adult low-p upon first capture = 0
  delta[1,4] <- 0            #Probability of being adult emigrated upon first capture = 0
  delta[1,5] <- 0
  delta[1,6] <- 0
  
  #2nd row: adults
  delta[2,1] <- 0           #Probability of being young upon first capture = 0 for ad
  delta[2,2] <- piaH        #Probability of being adult high-p upon first capture
  delta[2,3] <- 1-piaH      #Probability of being adult low-p upon first capture
  delta[2,4] <- 0           #Probability of being adult emigrated upon first capture = 0
  delta[2,5] <- 0
  delta[2,6] <- 0
  
  
  ###################  TRANSITIONS  #########################
  
  ##################  SURVIVAL - Juveniles  #######################
  
  #Transition matrix for juveniles (contains juvenile survival and adult survival probabilities)
  
  for(c in 1:n.ms.class){ #Indexing by classes of mismatch rather than by individual to save memory and computing time
    for(t in 1:n.occasions-1){
      #Alive young
      gamma[1,1,c,t] <- 0                           # Never stays alive young
      gamma[1,2,c,t] <- Sj[c,t] * pjH               # Can survive and become high-p adult
      gamma[1,3,c,t] <- Sj[c,t] * pjL               # Can survive and become low-p adult
      gamma[1,4,c,t] <- Sj[c,t] * pjE               # Can survive and emigrate
      gamma[1,5,c,t] <- (1-Sj[c,t]) * r[t+1]         # Can die (and be recovered at occasion t+1)
      gamma[1,6,c,t] <- (1-Sj[c,t]) * (1-r[t+1])     # Can die (and not be recovered at occasion t+1)
      
      #Alive adult (high-p)
      gamma[2,1,c,t] <- 0                      # Never becomes alive young
      gamma[2,2,c,t] <- Sa[t]                  # Survival (stays high-p)
      gamma[2,3,c,t] <- 0                      # Survival (no tr. to low-p)
      gamma[2,4,c,t] <- 0                      # Survival (no emigration for ad)
      gamma[2,5,c,t] <- (1-Sa[t]) * r[t+1]     # Can die (and be recovered at occasion t+1)
      gamma[2,6,c,t] <- (1-Sa[t]) * (1-r[t+1]) # Can die (and not be recovered at occasion t+1)
      
      #Alive adult (low-p)
      gamma[3,1,c,t] <- 0                      # Never becomes alive young
      gamma[3,2,c,t] <- 0                      # Survival (no tr. to high-p)
      gamma[3,3,c,t] <- Sa[t]                  # Survival (stays low-p)
      gamma[3,4,c,t] <- 0                      # Survival (no emigration for ad)
      gamma[3,5,c,t] <- (1-Sa[t]) * r[t+1]     # Can die (and be recovered at occasion t+1)
      gamma[3,6,c,t] <- (1-Sa[t]) * (1-r[t+1]) # Can die (and not be recovered at occasion t+1)
      
      #Alive adult (emigrated)
      gamma[4,1,c,t] <- 0                      # Never becomes alive young
      gamma[4,2,c,t] <- 0                      # Survival (no tr. to high-p)
      gamma[4,3,c,t] <- 0                      # Survival (no tr. to low-p)
      gamma[4,4,c,t] <- Sa[t]                  # Survival (stays emigrated)
      gamma[4,5,c,t] <- (1-Sa[t]) * r[t+1]     # Can die (and be recovered at occasion t+1)
      gamma[4,6,c,t] <- (1-Sa[t]) * (1-r[t+1]) # Can die (and not be recovered at occasion t+1)
      
      #newly dead (all ages), always becomes 'permanently dead'
      gamma[5,1,c,t] <- 0
      gamma[5,2,c,t] <- 0
      gamma[5,3,c,t] <- 0
      gamma[5,4,c,t] <- 0
      gamma[5,5,c,t] <- 0
      gamma[5,6,c,t] <- 1
      
      #"Old" dead (all ages), stays 'permanently dead'
      gamma[6,1,c,t] <- 0
      gamma[6,2,c,t] <- 0
      gamma[6,3,c,t] <- 0
      gamma[6,4,c,t] <- 0
      gamma[6,5,c,t] <- 0
      gamma[6,6,c,t] <- 1
      
    }# for t in 1:n.occasions-1
  }# for c in 1:n.ms.classes
  
  
  ##################  SURVIVAL - Adults  #######################
  
  #for birds marked as adults, only need to index gamma by time and not individual since they share the same matrix slice for a given time step (no individual covariates)
  #gi.ad takes the value of number of juveniles +1,
  #gi[i] in the likelihood below takes the value of juvenile ID for juveniles and number of juveniles +1 for all adults
  
  for (t in 1:(n.occasions-1)){
    
    #Alive young
    gamma[1,1,gi.ad,t] <- 0                      #Never stays alive young
    gamma[1,2,gi.ad,t] <- 0                      #No juvenile survival for adult birds  
    gamma[1,3,gi.ad,t] <- 0                      #No juvenile survival for adult birds  
    gamma[1,4,gi.ad,t] <- 0                      #No juvenile survival for adult birds  
    gamma[1,5,gi.ad,t] <- 0                      #No juvenile survival for adult birds  
    gamma[1,6,gi.ad,t] <- 0                      #No juvenile survival for adult birds  
    
    #Alive adults (high-p)
    gamma[2,1,gi.ad,t] <- 0                     #Never becomes alive young
    gamma[2,2,gi.ad,t] <- Sa[t]                 #Survival (stays high-p) 
    gamma[2,3,gi.ad,t] <- 0                     #Survival (no tr. to low-p)
    gamma[2,4,gi.ad,t] <- 0                     #Survival (no emigration for ad)
    gamma[2,5,gi.ad,t] <- (1-Sa[t]) *r[t+1]     #Can die and be recovered 
    gamma[2,6,gi.ad,t] <- (1-Sa[t]) *(1-r[t+1]) #Can die and not be recovered 
    
    #Alive adults (low-p)
    gamma[3,1,gi.ad,t] <- 0                     #Never becomes alive young
    gamma[3,2,gi.ad,t] <- 0                     #Survival (no tr. to high-p) 
    gamma[3,3,gi.ad,t] <- Sa[t]                 #Survival (stays low-p) 
    gamma[3,4,gi.ad,t] <- 0                     #Survival (no emigration for ad)
    gamma[3,5,gi.ad,t] <- (1-Sa[t]) *r[t+1]     #Can die and be recovered 
    gamma[3,6,gi.ad,t] <- (1-Sa[t]) *(1-r[t+1]) #Can die and not be recovered 
    
    #Alive adults (low-p)
    gamma[4,1,gi.ad,t] <- 0                     #Never becomes alive young
    gamma[4,2,gi.ad,t] <- 0                     #Survival (no tr. to high-p) 
    gamma[4,3,gi.ad,t] <- 0                     #Survival (no tr. to low-p) 
    gamma[4,4,gi.ad,t] <- Sa[t]                 #Survival (stays emigrated) 
    gamma[4,5,gi.ad,t] <- (1-Sa[t]) *r[t+1]     #Can die and be recovered 
    gamma[4,6,gi.ad,t] <- (1-Sa[t]) *(1-r[t+1]) #Can die and not be recovered 
    
    #newly dead 
    gamma[5,1,gi.ad,t] <- 0
    gamma[5,2,gi.ad,t] <- 0
    gamma[5,3,gi.ad,t] <- 0
    gamma[5,4,gi.ad,t] <- 0
    gamma[5,5,gi.ad,t] <- 0
    gamma[5,6,gi.ad,t] <- 1
    
    #"Old" dead 
    gamma[6,1,gi.ad,t] <- 0
    gamma[6,2,gi.ad,t] <- 0
    gamma[6,3,gi.ad,t] <- 0
    gamma[6,4,gi.ad,t] <- 0
    gamma[6,5,gi.ad,t] <- 0
    gamma[6,6,gi.ad,t] <- 1
    
  }# for t in 1:(n.occasions-1)
  
  ##################  EVENTS  ####################
  
  # To reduce number of nodes sampled, we are indexing the observation matrix per time twice.
  # indexing by m serves to index by marking occasion, so that only m:n.occasions matrix slices are created
  # This allows constraining first capture probability to 1
  for(m in 1:n.occasions){ # m takes the value of the marking occasion
    
    # values for initial capture (at time = m)
    
    #Alive young observed as young first capture = 1
    omega[1,1,m,m] <- 1  
    omega[1,2,m,m] <- 0
    omega[1,3,m,m] <- 0  
    
    #Alive and adult high-p observed on first capture = 1
    omega[2,1,m,m] <- 1
    omega[2,2,m,m] <- 0
    omega[2,3,m,m] <- 0
    
    #Alive and adult low-p observed on first capture = 1
    omega[3,1,m,m] <- 1
    omega[3,2,m,m] <- 0
    omega[3,3,m,m] <- 0
    
    #Alive and adult emigrated impossible on first capture
    omega[4,1,m,m] <- 0
    omega[4,2,m,m] <- 0
    omega[4,3,m,m] <- 0
    
    #Newly dead impossible on first capture
    omega[5,1,m,m] <- 0
    omega[5,2,m,m] <- 0 
    omega[5,3,m,m] <- 0 
    
    #Dead impossible on first capture
    omega[6,1,m,m] <- 0
    omega[6,2,m,m] <- 0 
    omega[6,3,m,m] <- 0 
    
    
    #Then we define event matrices for captures in occasions (m+1):n.occasions:
    for (t in (m+1):n.occasions){
      # Define probabilities of O(t) given S(t)
      #Alive young never observed as young beyond 1st capture
      omega[1,1,m,t] <- 0  
      omega[1,2,m,t] <- 0
      omega[1,3,m,t] <- 1  
      
      #Alive and adult high.p can be observed 
      omega[2,1,m,t] <- p.high[t]
      omega[2,2,m,t] <- 0
      omega[2,3,m,t] <- 1 - p.high[t]
      
      #Alive and adult low.p can be observed
      omega[3,1,m,t] <- p.low[t]
      omega[3,2,m,t] <- 0
      omega[3,3,m,t] <- 1 - p.low[t]
      
      #Alive and adult emigrated, can never be observed; p = 0
      omega[4,1,m,t] <- 0
      omega[4,2,m,t] <- 0
      omega[4,3,m,t] <- 1
      
      #Newly dead and recovered, =1 because recoveries are estimated in transition matrix
      omega[5,1,m,t] <- 0
      omega[5,2,m,t] <- 1
      omega[5,3,m,t] <- 0
      
      omega[6,1,m,t] <- 0
      omega[6,2,m,t] <- 0
      omega[6,3,m,t] <- 1
      
    
    } # for t in (m+1):(n.occasions)
  }#for m in 1:n.occasions
  
  #For all adults, ms.class takes the same value, i.e. gi.ad, or number of mismatch classes + 1
  ms.class[(n.juv+1):nind]<-n.ms.class+1
  
  
  ########################################################## #
  ###################### LIKELIHOOD ########################
  ########################################################## #
  
  
  # Likelihood
  for (i in 1:nind) {
    
    y[i, f[i]:n.occasions] ~ dDHMMo(init = delta[age[i],1:6], 
                                    probTrans = gamma[1:6, 1:6, ms.class[i], f[i]:(n.occasions-1)],
                                    probObs = omega[1:6, 1:3, f[i], f[i]:n.occasions],
                                    len = length(f[i]:n.occasions),
                                    checkRowSums = 0) 
    
  }
  
  ###################  Derived parameters  ######################

  #Predict mean annual juvenile survival from random year effect and average mismatch in year t
  for(t in 1:(n.occasions-1)){
    logit(Sj.t[t])<-  phi.J + beta.ms * mean.ms.t[t] + eps.phiJ.t[t]
  }
  
  #mean survival per period
  Sa.mean.p1 <- mean(Sa[1:n.occasions.p1])
  Sj.mean.p1 <- mean(Sj.t[1:n.occasions.p1])
  
  Sa.mean.p2 <- mean(Sa[(n.occasions.p1+1):(n.occasions-1)])
  Sj.mean.p2 <- mean(Sj.t[(n.occasions.p1+1):(n.occasions-1)])
  
  
  
  #Predictions of mismatch effect on juvenile survival
  for(i in 1:lpred){
    logit(S.pred.ms[i])<- phi.J + beta.ms * range.mismatch[i]
  }
  
})# end model

############################################################### #
######################### RUN MCMC ###########################
############################################################## #


n.chains<-3
n.thin<-3       #Thinning des paramètres
n.thin2<-30     #Thinning des prédictions de l'âge car ça prend trop de mémoire
n.iter = 1111
n.burnin = 0


myConstants = my.constants.joint
myData = my.data.joint
myParameters = my.parameters.joint
myParameters2 = my.parameters.joint.thin2
myCode = joint_mismatch_surv_analysis_optz



#Ce code roule 3 chaînes en parallèle, et enregistre les résultas au fur et à mesure (à toutes les 1000 itérations) dans les fichiers créés en haut du script
run_MCMC_allcode_save <- function(i,nimportequoi,myCode,myConstants,myData,myParameters,myParameters2,n.iter,n.burnin,n.thin,n.thin2,b.age.unk.vec=b.age.unk.vec, wd.save) {
  library(nimble)
  library(nimbleEcology)
  
  #Functions to get and set model and MCMC state and variables (allows re-starting MCMC from where it ended later if needed)
  getStateVariableNames <- function(samplerDef) {
    resetMethod <- body(samplerDef$reset)
    stateVars <- character()
    if(resetMethod[[1]] != '{') stop('something wrong')
    numLines <- length(resetMethod)
    for(i in 1:numLines) {
      if(i == 1) next
      thisLine <- resetMethod[[i]]
      if(thisLine[[1]] == '<<-') {
        LHS <- thisLine[[2]]
        if(!is.name(LHS)) stop('haven\'t dealt with non-name-LHS case yet')
        stateVars <- c(stateVars, as.character(LHS))
      }
      if('my_calcAdaptationFactor' %in% all.names(thisLine)) {
        stateVars <- c(stateVars, 'my_calcAdaptationFactor')
      }
    }
    setupMethod <- body(samplerDef$setup)
    if('empirSamp' %in% all.names(setupMethod)) stateVars <- c(stateVars, 'empirSamp')
    return(stateVars)
  }
  
  getModelState <- function(model) {
    modelVarNames <- model$getVarNames()
    modelVarValuesList <- vector('list', length(modelVarNames))
    names(modelVarValuesList) <- modelVarNames
    for(var in modelVarNames) {
      modelVarValuesList[[var]] <- model[[var]]
    }
    return(modelVarValuesList)
  }
  
  getMCMCstate <- function(conf, mcmc) {
    stateVarNamesList <- vector('list', length(conf$samplerConfs))
    mcmcStateValuesList <- vector('list', length(conf$samplerConfs))
    for(i in seq_along(conf$samplerConfs)) {
      samplerDef <- conf$getSamplerDefinition(i)
      theseStateNames <- getStateVariableNames(samplerDef)
      theseStateValuesList <- vector('list', length(theseStateNames))
      names(theseStateValuesList) <- theseStateNames
      for(j in seq_along(theseStateNames)) {
        if(is.nf(mcmc)) {
          if(theseStateNames[j] == 'my_calcAdaptationFactor') {
            theseStateValuesList[[j]] <- list(timesAdapted = mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]$timesAdapted,
                                              gamma1 = mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]$gamma1)
          } else
            theseStateValuesList[[j]] <- mcmc$samplerFunctions$contentsList[[i]][[theseStateNames[j]]]
        }
        if(is.Cnf(mcmc)) {
          if(theseStateNames[j] == 'my_calcAdaptationFactor') {
            theseStateValuesList[[j]] <- list(timesAdapted = valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'timesAdapted'),
                                              gamma1 = valueInCompiledNimbleFunction(mcmc$Robject$samplerFunctions$contentsList[[i]]$my_calcAdaptationFactor$.CobjectInterface, 'gamma1'))
          } else
            theseStateValuesList[[j]] <- valueInCompiledNimbleFunction(mcmc$samplerFunctions[[i]], theseStateNames[j])
        }
      }
      mcmcStateValuesList[[i]] <- theseStateValuesList
    }
    return(mcmcStateValuesList)
  }
  
  
  the.inits.joint <- function(){list(
    
    #### AGE MODEL ##### #
    int_age = runif(1, 0, 50),
    beta_prim9 = runif(1, 0, 10),
    sd.prim = runif(1,0,4),
    sd.year.prim = runif(1,0,4),
    age.at.B.est = rep(30, times=length(b.age.unk.vec)),
    
    ##### SURVIVAL MODEL ##### #
    phi.A= runif(1, -3, 3),
    phi.J= runif(1, -3, 3),
    beta.p.L = runif(1, -2, 0),
    beta.ms = runif(1,-2,2),
    beta.hunt.ad = runif(1,-2,2),
    beta.hunt.juv = runif(1,-2,2),
    int.p = runif(1,-5,5),
    int.r = runif(1, -5, 5),
    sigma.phiA = runif(1,0,3),
    sigma.phiJ = runif(1,0,3),
    sigma.p = runif(1,0,3),
    sigma.r = runif(1,0,3),
    piaH = runif(1,0,1),
    lpjH = runif(1,-5,5),
    lpjE = runif(1,-5,5))}
  
  
  gsg.mod<-nimbleModel(code = myCode,
                       data = myData,
                       constants = myConstants,
                       inits = the.inits.joint(),
                       calculate=F)
  
  conf <- configureMCMC(gsg.mod,monitors=myParameters,monitors2=myParameters2, useConjugacy = F)
  
  #If we think using block samplers may help with convergence for piaH
  conf$removeSamplers(c('lpjH','lpjE'))
  conf$addSampler(target = c('lpjH','lpjE'),
                  type = "RW_block")
  # conf
  
  
  gsg.MCMC <- buildMCMC(conf)
  
  #Compile model and MCMC 
  C.gsg.mod<-compileNimble(gsg.mod)
  
  Cgsg.MCMC <- compileNimble(gsg.MCMC, project=gsg.mod)
  
  start.chain=Sys.time()
  #First run, keep all samples for the first tests, can always discard them later
  Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T)
  end.chain=Sys.time()
  (time.chain<-end.chain-start.chain)
  
  #Extract and save samples
  chain1_run1<-as.matrix(Cgsg.MCMC$mvSamples)
  saveRDS(chain1_run1,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run1.rds'))
  
  chain1_run1_age<-as.matrix(Cgsg.MCMC$mvSamples2)
  saveRDS(chain1_run1_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run1_age.rds'))
  
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(stateList)
  
  #Run 1.2
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run1<-rbind(chain1_run1,samples)
  saveRDS(chain1_run1,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run1.rds'))
  
  
  chain1_run1_age<-rbind(chain1_run1_age,samples2)
  saveRDS(chain1_run1_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run1_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  
  #Run 1.3
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run1<-rbind(chain1_run1,samples)
  saveRDS(chain1_run1,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run1.rds'))
  
  
  chain1_run1_age<-rbind(chain1_run1_age,samples2)
  saveRDS(chain1_run1_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run1_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 1.4
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run1<-rbind(chain1_run1,samples)
  saveRDS(chain1_run1,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run1.rds'))
  
  
  chain1_run1_age<-rbind(chain1_run1_age,samples2)
  saveRDS(chain1_run1_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run1_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 1.5
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run1<-rbind(chain1_run1,samples)
  saveRDS(chain1_run1,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run1.rds'))
  
  
  chain1_run1_age<-rbind(chain1_run1_age,samples2)
  saveRDS(chain1_run1_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run1_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(chain1_run1, chain1_run1_age, samples, samples2, stateList)
  end.chain_r1=Sys.time()
  (time.chain<-end.chain_r1-start.chain)
  
  
  
  
  
  start.chain=Sys.time()
  #First run, keep all samples for the first tests, can always discard them later
  Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T)
  end.chain=Sys.time()
  (time.chain<-end.chain-start.chain)
  
  #Extract and save samples
  chain1_run2<-as.matrix(Cgsg.MCMC$mvSamples)
  saveRDS(chain1_run2,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run2.rds'))
  
  chain1_run2_age<-as.matrix(Cgsg.MCMC$mvSamples2)
  saveRDS(chain1_run2_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run2_age.rds'))
  
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(stateList)
  
  #Run 2.2
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run2<-rbind(chain1_run2,samples)
  saveRDS(chain1_run2,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run2.rds'))
  
  
  chain1_run2_age<-rbind(chain1_run2_age,samples2)
  saveRDS(chain1_run2_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run2_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  
  #Run 2.3
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run2<-rbind(chain1_run2,samples)
  saveRDS(chain1_run2,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run2.rds'))
  
  
  chain1_run2_age<-rbind(chain1_run2_age,samples2)
  saveRDS(chain1_run2_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run2_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 2.4
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run2<-rbind(chain1_run2,samples)
  saveRDS(chain1_run2,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run2.rds'))
  
  
  chain1_run2_age<-rbind(chain1_run2_age,samples2)
  saveRDS(chain1_run2_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run2_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 2.5
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run2<-rbind(chain1_run2,samples)
  saveRDS(chain1_run2,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run2.rds'))
  
  
  chain1_run2_age<-rbind(chain1_run2_age,samples2)
  saveRDS(chain1_run2_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run2_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(chain1_run2, chain1_run2_age, samples, samples2, stateList)
  end.chain_r1=Sys.time()
  (time.chain<-end.chain_r1-start.chain)
  
  
  
  
  
  
  
  
  
  
  
  
  start.chain=Sys.time()
  #First run, keep all samples for the first tests, can always discard them later
  Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T)
  end.chain=Sys.time()
  (time.chain<-end.chain-start.chain)
  
  #Extract and save samples
  chain1_run3<-as.matrix(Cgsg.MCMC$mvSamples)
  saveRDS(chain1_run3,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run3.rds'))
  
  chain1_run3_age<-as.matrix(Cgsg.MCMC$mvSamples2)
  saveRDS(chain1_run3_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run3_age.rds'))
  
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(stateList)
  
  #Run 3.2
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run3<-rbind(chain1_run3,samples)
  saveRDS(chain1_run3,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run3.rds'))
  
  
  chain1_run3_age<-rbind(chain1_run3_age,samples2)
  saveRDS(chain1_run3_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run3_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  
  #Run 3.3
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run3<-rbind(chain1_run3,samples)
  saveRDS(chain1_run3,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run3.rds'))
  
  
  chain1_run3_age<-rbind(chain1_run3_age,samples2)
  saveRDS(chain1_run3_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run3_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 3.4
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run3<-rbind(chain1_run3,samples)
  saveRDS(chain1_run3,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run3.rds'))
  
  
  chain1_run3_age<-rbind(chain1_run3_age,samples2)
  saveRDS(chain1_run3_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run3_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 3.5
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run3<-rbind(chain1_run3,samples)
  saveRDS(chain1_run3,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run3.rds'))
  
  
  chain1_run3_age<-rbind(chain1_run3_age,samples2)
  saveRDS(chain1_run3_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run3_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(chain1_run3, chain1_run3_age, samples, samples2, stateList)
  end.chain_r1=Sys.time()
  (time.chain<-end.chain_r1-start.chain)
  
  
  
  start.chain=Sys.time()
  #First run, keep all samples for the first tests, can always discard them later
  Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T)
  end.chain=Sys.time()
  (time.chain<-end.chain-start.chain)
  
  #Extract and save samples
  chain1_run4<-as.matrix(Cgsg.MCMC$mvSamples)
  saveRDS(chain1_run4,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run4.rds'))
  
  chain1_run4_age<-as.matrix(Cgsg.MCMC$mvSamples2)
  saveRDS(chain1_run4_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run4_age.rds'))
  
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(stateList)
  
  #Run 4.2
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run4<-rbind(chain1_run4,samples)
  saveRDS(chain1_run4,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run4.rds'))
  
  
  chain1_run4_age<-rbind(chain1_run4_age,samples2)
  saveRDS(chain1_run4_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run4_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  
  #Run 4.3
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run4<-rbind(chain1_run4,samples)
  saveRDS(chain1_run4,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run4.rds'))
  
  
  chain1_run4_age<-rbind(chain1_run4_age,samples2)
  saveRDS(chain1_run4_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run4_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 4.4
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run4<-rbind(chain1_run4,samples)
  saveRDS(chain1_run4,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run4.rds'))
  
  
  chain1_run4_age<-rbind(chain1_run4_age,samples2)
  saveRDS(chain1_run4_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run4_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 4.5
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run4<-rbind(chain1_run4,samples)
  saveRDS(chain1_run4,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run4.rds'))
  
  
  chain1_run4_age<-rbind(chain1_run4_age,samples2)
  saveRDS(chain1_run4_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run4_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(chain1_run4, chain1_run4_age, samples, samples2, stateList)
  end.chain_r1=Sys.time()
  (time.chain<-end.chain_r1-start.chain)
  
  
  
  # library(MCMCvis)
  # 
  # mod.out<-(rbind(chain1_run1, chain1_run2, chain1_run3, chain1_run4))
  # no_est<-as.numeric(which(mod.out[2,]==0))
  # mod.out.nomiss<-mod.out[,-c(no_est)]
  # MCMCtrace(mod.out.nomiss, params= c('Sa.mean.p1','Sa.mean.p2','Sj.mean.p1','Sj.mean.p2','beta.ms','beta.p.L','beta_prim9','int.p','int.r','int_age','phi.A','phi.J','piaH','pjE','pjH','pjL','sd.prim','sd.year.prim','sigma.p','sigma.phiA','sigma.phiJ','sigma.r'),
  #           pdf=F,iter=13000)
  # MCMCtrace(mod.out.nomiss, params= c('beta.hunt.ad','beta.hunt.juv'),pdf=F, iter=13000)
  # 
  
  
  
  
  
  
  
  
  start.chain=Sys.time()
  #First run, keep all samples for the first tests, can always discard them later
  Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T)
  end.chain=Sys.time()
  (time.chain<-end.chain-start.chain)
  
  #Extract and save samples
  chain1_run5<-as.matrix(Cgsg.MCMC$mvSamples)
  saveRDS(chain1_run5,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run5.rds'))
  
  chain1_run5_age<-as.matrix(Cgsg.MCMC$mvSamples2)
  saveRDS(chain1_run5_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run5_age.rds'))
  
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(stateList)
  
  #Run 5.2
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run5<-rbind(chain1_run5,samples)
  saveRDS(chain1_run5,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run5.rds'))
  
  
  chain1_run5_age<-rbind(chain1_run5_age,samples2)
  saveRDS(chain1_run5_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run5_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  
  #Run 5.3
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run5<-rbind(chain1_run5,samples)
  saveRDS(chain1_run5,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run5.rds'))
  
  
  chain1_run5_age<-rbind(chain1_run5_age,samples2)
  saveRDS(chain1_run5_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run5_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 5.4
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run5<-rbind(chain1_run5,samples)
  saveRDS(chain1_run5,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run5.rds'))
  
  
  chain1_run5_age<-rbind(chain1_run5_age,samples2)
  saveRDS(chain1_run5_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run5_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 5.5
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run5<-rbind(chain1_run5,samples)
  saveRDS(chain1_run5,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run5.rds'))
  
  
  chain1_run5_age<-rbind(chain1_run5_age,samples2)
  saveRDS(chain1_run5_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run5_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  
  rm(chain1_run5, chain1_run5_age, samples, samples2, stateList)
  end.chain_r1=Sys.time()
  (time.chain<-end.chain_r1-start.chain)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  start.chain=Sys.time()
  #First run, keep all samples for the first tests, can always discard them later
  Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T)
  end.chain=Sys.time()
  (time.chain<-end.chain-start.chain)
  
  #Extract and save samples
  chain1_run6<-as.matrix(Cgsg.MCMC$mvSamples)
  saveRDS(chain1_run6,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run6.rds'))
  
  chain1_run6_age<-as.matrix(Cgsg.MCMC$mvSamples2)
  saveRDS(chain1_run6_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run6_age.rds'))
  
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(stateList)
  
  #Run 6.2
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run6<-rbind(chain1_run6,samples)
  saveRDS(chain1_run6,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run6.rds'))
  
  
  chain1_run6_age<-rbind(chain1_run6_age,samples2)
  saveRDS(chain1_run6_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run6_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  
  #Run 6.3
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run6<-rbind(chain1_run6,samples)
  saveRDS(chain1_run6,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run6.rds'))
  
  
  chain1_run6_age<-rbind(chain1_run6_age,samples2)
  saveRDS(chain1_run6_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run6_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 6.4
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run6<-rbind(chain1_run6,samples)
  saveRDS(chain1_run6,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run6.rds'))
  
  
  chain1_run6_age<-rbind(chain1_run6_age,samples2)
  saveRDS(chain1_run6_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run6_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 6.5
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run6<-rbind(chain1_run6,samples)
  saveRDS(chain1_run6,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run6.rds'))
  
  
  chain1_run6_age<-rbind(chain1_run6_age,samples2)
  saveRDS(chain1_run6_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run6_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  
  rm(chain1_run6, chain1_run6_age, samples, samples2, stateList)
  end.chain_r1=Sys.time()
  (time.chain<-end.chain_r1-start.chain)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  start.chain=Sys.time()
  #First run, keep all samples for the first tests, can always discard them later
  Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T)
  end.chain=Sys.time()
  (time.chain<-end.chain-start.chain)
  
  #Extract and save samples
  chain1_run7<-as.matrix(Cgsg.MCMC$mvSamples)
  saveRDS(chain1_run7,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run7.rds'))
  
  chain1_run7_age<-as.matrix(Cgsg.MCMC$mvSamples2)
  saveRDS(chain1_run7_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run7_age.rds'))
  
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(stateList)
  
  #Run 7.2
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run7<-rbind(chain1_run7,samples)
  saveRDS(chain1_run7,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run7.rds'))
  
  
  chain1_run7_age<-rbind(chain1_run7_age,samples2)
  saveRDS(chain1_run7_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run7_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  
  #Run 7.3
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run7<-rbind(chain1_run7,samples)
  saveRDS(chain1_run7,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run7.rds'))
  
  
  chain1_run7_age<-rbind(chain1_run7_age,samples2)
  saveRDS(chain1_run7_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run7_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 7.4
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run7<-rbind(chain1_run7,samples)
  saveRDS(chain1_run7,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run7.rds'))
  
  
  chain1_run7_age<-rbind(chain1_run7_age,samples2)
  saveRDS(chain1_run7_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run7_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 7.5
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run7<-rbind(chain1_run7,samples)
  saveRDS(chain1_run7,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run7.rds'))
  
  
  chain1_run7_age<-rbind(chain1_run7_age,samples2)
  saveRDS(chain1_run7_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run7_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  
  rm(chain1_run7, chain1_run7_age, samples, samples2, stateList)
  end.chain_r1=Sys.time()
  (time.chain<-end.chain_r1-start.chain)
  
  
  
  
  
  
  
  
  
  
  
  start.chain=Sys.time()
  #First run, keep all samples for the first tests, can always discard them later
  Cgsg.MCMC$run(niter=n.iter, nburnin=n.burnin, thin=n.thin, thin2=n.thin2, progressBar=T)
  end.chain=Sys.time()
  (time.chain<-end.chain-start.chain)
  
  #Extract and save samples
  chain1_run8<-as.matrix(Cgsg.MCMC$mvSamples)
  saveRDS(chain1_run8,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run8.rds'))
  
  chain1_run8_age<-as.matrix(Cgsg.MCMC$mvSamples2)
  saveRDS(chain1_run8_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run8_age.rds'))
  
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(stateList)
  
  #Run 8.2
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run8<-rbind(chain1_run8,samples)
  saveRDS(chain1_run8,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run8.rds'))
  
  
  chain1_run8_age<-rbind(chain1_run8_age,samples2)
  saveRDS(chain1_run8_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run8_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  
  #Run 8.3
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run8<-rbind(chain1_run8,samples)
  saveRDS(chain1_run8,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run8.rds'))
  
  
  chain1_run8_age<-rbind(chain1_run8_age,samples2)
  saveRDS(chain1_run8_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run8_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 8.4
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run8<-rbind(chain1_run8,samples)
  saveRDS(chain1_run8,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run8.rds'))
  
  
  chain1_run8_age<-rbind(chain1_run8_age,samples2)
  saveRDS(chain1_run8_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run8_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  rm(samples, samples2, stateList)
  
  #Run 8.5
  Cgsg.MCMC$run(niter=n.iter, nburnin=0, thin=n.thin, thin2=n.thin2, progressBar=T, reset=F, resetMV=T)
  
  #Extract and save samples
  samples<-as.matrix(Cgsg.MCMC$mvSamples)
  samples2<-as.matrix(Cgsg.MCMC$mvSamples2)
  
  #Extract and save samples
  chain1_run8<-rbind(chain1_run8,samples)
  saveRDS(chain1_run8,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run8.rds'))
  
  
  chain1_run8_age<-rbind(chain1_run8_age,samples2)
  saveRDS(chain1_run8_age,file=paste0(wd.save,'result_chain',i,'/chain',i,'_run8_age.rds'))
  
  #Extract model and MCMC state and save them
  stateList <- list(modelState = getModelState(C.gsg.mod),
                    mcmcState = getMCMCstate(conf=conf, mcmc=Cgsg.MCMC),
                    rs = .Random.seed)
  
  saveRDS(stateList, file=paste0(wd.save,'result_chain',i,'/model_state.rds'))
  
  
  rm(chain1_run8, chain1_run8_age, samples, samples2, stateList)
  end.chain_r1=Sys.time()
  (time.chain<-end.chain_r1-start.chain)
  
  
  
  
  
  return(Sys.time())
}


this.cluster <- makeCluster(mc <- getOption("cl.cores", n.chains))
registerDoParallel(this.cluster) # register the cluster
whole.run<-Sys.time()
chain_output <- parLapply(cl = this.cluster, 1:n.chains, 
                          fun = run_MCMC_allcode_save,
                          myConstants = my.constants.joint,
                          myData = my.data.joint,
                          myParameters = my.parameters.joint,
                          myParameters2 = my.parameters.joint.thin2,
                          myCode = joint_mismatch_surv_analysis_optz,
                          b.age.unk.vec=b.age.unk.vec,
                          n.iter = n.iter,
                          n.burnin = 0,
                          n.thin = n.thin,
                          n.thin2= n.thin2,
                          wd.save = wd.save)

stopCluster(this.cluster)
fin<-Sys.time()
fin-whole.run




#################################################################################### #
################################ READ AND COMPILE RESULTS ##########################
#################################################################################### #

#Load results from 3 chains
c1r1<-readRDS(paste0(wd.save,"result_chain1/chain1_run1.RDS"))
c1r2<-readRDS(paste0(wd.save,"result_chain1/chain1_run2.RDS"))
c1r3<-readRDS(paste0(wd.save,"result_chain1/chain1_run3.RDS"))
c1r4<-readRDS(paste0(wd.save,"result_chain1/chain1_run4.RDS"))
c1r5<-readRDS(paste0(wd.save,"result_chain1/chain1_run5.RDS"))
c1r6<-readRDS(paste0(wd.save,"result_chain1/chain1_run6.RDS"))
c1r7<-readRDS(paste0(wd.save,"result_chain1/chain1_run7.RDS"))
c1r8<-readRDS(paste0(wd.save,"result_chain1/chain1_run8.RDS"))

c2r1<-readRDS(paste0(wd.save,"result_chain2/chain2_run1.RDS"))
c2r2<-readRDS(paste0(wd.save,"result_chain2/chain2_run2.RDS"))
c2r3<-readRDS(paste0(wd.save,"result_chain2/chain2_run3.RDS"))
c2r4<-readRDS(paste0(wd.save,"result_chain2/chain2_run4.RDS"))
c2r5<-readRDS(paste0(wd.save,"result_chain2/chain2_run5.RDS"))
c2r6<-readRDS(paste0(wd.save,"result_chain2/chain2_run6.RDS"))
c2r7<-readRDS(paste0(wd.save,"result_chain2/chain2_run7.RDS"))
c2r8<-readRDS(paste0(wd.save,"result_chain2/chain2_run8.RDS"))

c3r1<-readRDS(paste0(wd.save,"result_chain3/chain3_run1.RDS"))
c3r2<-readRDS(paste0(wd.save,"result_chain3/chain3_run2.RDS"))
c3r3<-readRDS(paste0(wd.save,"result_chain3/chain3_run3.RDS"))
c3r4<-readRDS(paste0(wd.save,"result_chain3/chain3_run4.RDS"))
c3r5<-readRDS(paste0(wd.save,"result_chain3/chain3_run5.RDS"))
c3r6<-readRDS(paste0(wd.save,"result_chain3/chain3_run6.RDS"))
c3r7<-readRDS(paste0(wd.save,"result_chain3/chain3_run7.RDS"))
c3r8<-readRDS(paste0(wd.save,"result_chain3/chain3_run8.RDS"))

#Bind chain runs together
chain1<-rbind(c1r1,c1r2,c1r3,c1r4,c1r5,c1r6,c1r7,c1r8)#,c1r9,c1r10,c1r11,c1r12,c1r13,c1r14,c1r15)
chain2<-rbind(c2r1,c2r2,c2r3,c2r4,c2r5,c2r6,c2r7,c2r8)#,c2r9,c2r10,c2r11,c2r12,c2r13,c2r14)
chain3<-rbind(c3r1,c3r2,c3r3,c3r4,c3r5,c3r6,c3r7,c3r8)#,c3r9,c3r10,c3r11,c3r12,c3r13,c3r14)

#If chains didn't run at the same speed, resize to smallest chain size
dim1<-min(dim(chain1)[1],dim(chain2)[1],dim(chain3)[1])

# Remove columns with no estimations (i.e. t=1 for capture, recovery, and their random time effects)
no_est<-as.numeric(which(c1r1[2,]==0))


#Put in a list for MCMCvis package
mcmc_out<-list(chain1=chain1[1:dim1,-c(no_est)],chain2=chain2[1:dim1,-c(no_est)],chain3=chain3[1:dim1,-c(no_est)])



#Read the data if necessary
library(MCMCvis)

sim.data<-readRDS('~/Doc/juvmismatch/data_sim_joint_analysis_2_54k_27occ.RDS')
T.recap.H<-sim.data$T.recap.H             #True recapture high
T.recap.L<-sim.data$T.recap.L             #True recapture low
actual.survival<-sim.data$sim.survival    #True survival
T.emigr<-sim.data$T.emigr                 #True emigration probability
T.p.high<-sim.data$T.p.high               #True Probability of juveniles transitioning to high capture group               
T.p.low<-sim.data$T.p.low                 #True Probability of juveniles transitioning to low capture group    
T.piaH<-sim.data$T.piaH                   #True Probability of adults being in the high capture group
T.recov.all<-sim.data$T.recov.all         #True Recovery probabilities
df_sim_juvs<-sim.data$df_sim_juvs         # Simulated mismatch and feather length data
period1<-sim.data$period1                 # Index for period 1
period2<-sim.data$period2[1:13]           # Index for period 2
sim_juvs<-sim.data$df_sim_juvs     
s_ms<-sim.data$mismatch                   #Mismatch values

###Input parameters for the simulation
library(lme4)
lm_gr<-lmer(data=df_sim_juvs, age.at.B~scale(prim9.B)+(1|an_B))
lm_gr_sum<-summary(lm_gr)

s_int.age<-lm_gr_sum$coefficients[1,1]      # Age model intercept feather 
s_beta.prim9<-lm_gr_sum$coefficients[2,1]   # Age model beta feather
s_sd.prim<-lm_gr_sum$sigma                  # Age model residual variance
s_sd.yr.prim<-sqrt(var(ranef(lm_gr)$an_B))  # Age model random year effect variance
phiA<-2                                     # true Phi adults, 
sigma.phiA<-0.25                            # etc... the rest is intuitive
phiJ<-0.35
sigma.phiJ<-0.35
int.r<-qlogis(0.3)
sigma.r<-0.4
int.p<-qlogis(0.1)
sigma.p<-0.4


s_Sa.mean.p1<-mean(actual.survival[period1,3])
s_Sa.mean.p2<-mean(actual.survival[period2,3])

s_Sj.mean.p1<-mean(actual.survival[period1,2])
s_Sj.mean.p2<-mean(actual.survival[period2,2])

beta.p.low<--1
hunt.effect.ad<- -0.62
hunt.effect.juv<- -0.3
s_beta.ms<- -0.045


s_piaH<-0.4 # proportion of adults in high observability group

#Transitions of juveniles to adults with high/lo/no capturability
s_pjH<-0.3 # proportion of surviving young that are in the high recapture group as adults
s_pjE<-0.5 # proportion of surviving young to emigrate
s_pjL <-1-s_pjH-s_pjE #Proportion of surviving young to transition to 'low'

#Pack these values to put into traceplot function
packed.gvalues<-c(s_int.age,s_beta.prim9,s_sd.prim,s_sd.yr.prim,
                  phiA,sigma.phiA,
                  phiJ,sigma.phiJ,
                  int.r,sigma.r,
                  int.p,sigma.p,
                  beta.p.low,s_beta.ms,hunt.effect.ad,hunt.effect.juv,
                  s_piaH,s_pjH,s_pjE,s_pjL,
                  s_Sa.mean.p1,s_Sa.mean.p2,s_Sj.mean.p1,s_Sj.mean.p2)

#each sample burned removes 3 iterations
burn<-1000

MCMCtrace(mcmc_out,params = c('int_age','beta_prim9','sd.prim','sd.year.prim',
                                                 "phi.A", 'sigma.phiA',
                                                 "phi.J", 'sigma.phiJ', 
                                                 "int.r", 'sigma.r', 
                                                 'int.p', 'sigma.p', 
                                                 "beta.p.L","beta.ms",'beta.hunt.ad','beta.hunt.juv',
                                                 'piaH', 'pjH', 'pjE', 'pjL',
                                                 'Sa.mean.p1', 'Sa.mean.p2', 'Sj.mean.p1', 'Sj.mean.p2'),
          pdf=F,iter = (dim1-burn),
          gvals = packed.gvalues)


#Range of iterations to keep for visualizing results
range<-c(burn:dim1)

mcmc_results_burned<-list(chain1=chain1[burn:dim1,-c(no_est)],chain2=chain2[burn:dim1,-c(no_est)],chain3=chain3[burn:dim1,-c(no_est)])


#Graphs

# Recapture high P
#summarize results from MCMC
(summary.recap.H<-MCMCsummary(mcmc_results_burned, params=c('p.high')))

#Create data frame with results and true input values
comp.recap.H<-data.frame(time=1:nrow(summary.recap.H),sim=T.recap.H,mean.est=summary.recap.H$mean,lci.est=summary.recap.H$"2.5%",uci.est=summary.recap.H$"97.5%")

#Plot both (red is input)
ggplot(data=comp.recap.H,
       aes(y=mean.est,
           x=time))+
  geom_point()+
  geom_errorbar(aes(ymax=uci.est,ymin=lci.est), width=0.5)+
  geom_point(aes(y=sim), col='red')+
  ggtitle(label='recaptures - high')


# Recapture low P
(summary.recap.L<-MCMCsummary(mcmc_results_burned, params=c('p.low')))

comp.recap.L<-data.frame(time=1:nrow(summary.recap.L),sim=T.recap.L,mean.est=summary.recap.L$mean,lci.est=summary.recap.L$"2.5%",uci.est=summary.recap.L$"97.5%")

ggplot(data=comp.recap.L,
       aes(y=mean.est,
           x=time))+
  geom_point()+
  geom_errorbar(aes(ymax=uci.est,ymin=lci.est), width=0.5)+
  geom_point(aes(y=sim), col='red')+
  ggtitle(label='recaptures - low')


#Recoveries
(summary.recov<-MCMCsummary(mcmc_results_burned, params=c('r')))

comp.recov<-data.frame(time=1:nrow(summary.recov),sim=T.recov.all,mean.est=summary.recov$mean,lci.est=summary.recov$"2.5%",uci.est=summary.recov$"97.5%")

ggplot(data=comp.recov,
       aes(y=mean.est,
           x=time))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymax=uci.est,ymin=lci.est), width=0.5)+
  geom_point(aes(y=sim), col='red')+
  geom_line(aes(y=sim), col='red')+
  ggtitle(label='recoveries')

#Survie adulte ####
actual.survival<-as.data.frame(actual.survival)

(summary.survie<-MCMCsummary(mcmc_results_burned, params=c('Sa','Sj.t')))


comp.survival<-data.frame(parameter=rep(c('adultes','jeunes'), each=nrow(actual.survival)),
                          time=rep(1:nrow(actual.survival),times=2),
                          sim=c(actual.survival$'Survie adultes',actual.survival$'Survie jeunes'),
                          mean.est=summary.survie$mean,
                          lci.est=summary.survie$"2.5%",
                          uci.est=summary.survie$"97.5%")

#Graph adult survival
ggplot(data=comp.survival%>%filter(parameter=='adultes'),
       aes(y=mean.est,
           x=time))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymax=uci.est,ymin=lci.est), width=0.5)+
  geom_point(aes(y=sim), col='red')+
  geom_line(aes(y=sim), col='red')

#Graph juvenile survival
ggplot(data=comp.survival%>%filter(parameter=='jeunes'),
       aes(y=mean.est,
           x=time))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymax=uci.est,ymin=lci.est), width=0.5)+
  geom_point(aes(y=sim), col='red')+
  geom_line(aes(y=sim), col='red')


# Heterogeneity parameters
(summary.hetero<-MCMCsummary(mcmc_results_burned, params=c('pjE','pjH','pjL','piaH')))
n.occasions<-nrow(actual.survival)+1
comp.hetero<-data.frame(parameter=rep(c('pjE','pjH','pjL','piaH'), each=n.occasions-1),
                        sim=c(T.emigr,T.p.high,T.p.low,T.piaH),
                        mean.sim=rep(c(mean(T.emigr),mean(T.p.high),mean(T.p.low),mean(T.piaH)), each=n.occasions-1),
                        sd.sim=rep(c(sd(T.emigr),sd(T.p.high),sd(T.p.low),sd(T.piaH)), each=n.occasions-1),
                        mean.est=rep(summary.hetero$mean, each=n.occasions-1),
                        lci.est=rep(summary.hetero$"2.5%",each=n.occasions-1),
                        uci.est=rep(summary.hetero$"97.5%", each=n.occasions-1))


ggplot(data=comp.hetero,
       aes(y=mean.est,
           x=parameter))+
  geom_point()+
  geom_errorbar(aes(ymax=uci.est,ymin=lci.est), width=0.1)+
  geom_jitter(aes(y=sim), col='red', width=0.1, pch=4)+
  geom_point(aes(x=parameter, y=mean.sim), col='red',size=2)+
  geom_errorbar(aes(ymax=mean.sim+1.96*sd.sim, ymin=mean.sim-1.96*sd.sim), width=0.1, col='red')+
  scale_y_continuous(lim=c(0,1))




#Predictions mismatch effect
(summary.pred.ms<-MCMCsummary(mcmc_results_burned, params=c('S.pred.ms')))

#Compute mean mismatch by year
sim_juvs%>%dplyr::rename(time=an_B)%>%group_by(time)%>%dplyr::summarize(average_mismatch=mean(mismatch.all))->mean_mismatch

#Scaling parameters
scale.ms.mean<-12.40132
scale.ms.sd<-6.25861

#Join mismatch data and juvenile survival estimations
check_ms_effect<-mean_mismatch%>%
  left_join(comp.survival%>%filter(parameter=='jeunes'), by='time')

#Predictions directly from MCMC chains
preds_mismatch_effect<-rbind(mcmc_results_burned$chain1[,1:100],mcmc_test_age_model_no_miss$chain2[,1:100],mcmc_test_age_model_no_miss$chain2[,1:100])

#Keep only 95% of samples
vals<-apply(preds_mismatch_effect, MARGIN = 2, FUN = function(x){sort(x)[c(round(0.025*length(x)),round(0.5*length(x)),round(0.975*length(x)))]} )

#Compute confidence intervals for mismatch predictions
CI_mismatch<-pivot_longer(as.data.frame(vals),names_to='xvalue', cols=c(1:100))
CI_mismatch<-data.frame(LCI=CI_mismatch[1:100,2],mean=CI_mismatch[101:200,2],UCI=CI_mismatch[201:300,2])

#Compute range of mismatch for predictions used in model
range.mismatch<-seq(from=(min(s_ms)-scale.ms.mean)/scale.ms.sd, to=(max(s_ms)-scale.ms.mean)/scale.ms.sd, length.out=100)

#Pack all information in a dataframe
effect_mismatch<-cbind(CI_mismatch, range.mismatch)
names(effect_mismatch)<-c('LCI','mean','UCI','range_mismatch.scaled')

#Compute the simulated effect based on the parameter used and the range of mismatch values fed to model
effect_mismatch<-effect_mismatch%>%
  mutate(mismatch_original=range_mismatch.scaled*scale.ms.sd+scale.ms.mean)%>%
  mutate(sim_effect=plogis(mismatch_original*-0.045))

#Plot mismatch effect estimated vs the one simulated
ggplot(effect_mismatch,
       aes(x=mismatch_original,
           y=(mean)))+
  geom_line()+
  geom_line(aes(y=sim_effect), col='red', lty=2)+
  #geom_line(aes(y=sim_effect2), col='red', lty=3)+
  geom_ribbon(aes(ymax=UCI, ymin=LCI), alpha=0.3)+
  geom_point(data=check_ms_effect, aes(x=average_mismatch,y=mean.est), col='green')













