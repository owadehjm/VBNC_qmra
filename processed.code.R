
#==================================================================================
### QMRA OF ECO157:H7 IN LETTUCE ### PROCESSED LETTUCE ###
#==================================================================================

# The QMRA tracks the produce from production to consumption
# The produce is whole and nut cut

# The model depicts the risk of infection with EcO157:H7 among the general public after consumption of lettuce
# Load libraries that will be used in the code
library(truncnorm)#Generating truncated normal distribution data
library(EnvStats) #Generate triangle distributions
library(TruncExpFam) #Generate truncated exponential values
library(mc2d) # Pert distribution and truncation
library(e1071) # Discrete distributions
library(dplyr)# For data manipulation
library(ggplot2) # For visualization of the data

#===================================================================================
### EXPOSURE ASSESSMENT
## A. INFIELD PRODUCTION
#===================================================================================
# NAMING CONVENTION
# t- Time
# T- temperature
# C- Concentration
# p- prevalence
# a - annual

# Create a dataframe to save outputs
P.ill<-runif(1e4)

processed <- data.frame(P.ill_0_0 = P.ill)

a_processed <- data.frame(P.ill_0_0 = P.ill)

# Function for processed lettuce
processed_scenario<-function(col1,col2){
  # Set seed for reproducibility (optional)
  set.seed(1234)
  
  # Set the number of iterations to be conducted using Monte Carlo techniques
  iter=1e4
  mypath<- getwd()
  
  path<-paste(mypath,"/Plots",sep="")
  #-------------------------------------------------------------------------------
  ### Specify the scenarios
  #-------------------------------------------------------------------------------
  # Scenario 2- Harvest temp
  #Specify the harvest temperature as "baseline" or 0, 17 or 9
  harvest.T= col1 # The value can be baseline or 0, 9 or 17
  
  # Scenario 3- Cold storage time
  # Specify the cold storage time as either "None", 0, 24, 48, 72, 96, 120 hours
  cold.storage.time= col2
  
  ## 1. Irrigation water quality (Rw)
  # Quality of the irrigation water used on the romaine lettuce in the field
  # Uniform distribution of values is assumed
  # Source Citation: LGMA 2010
  # E coli level for foliar application should not be higher than 235 and the limit of detection for culturable counts is 1 cfu/100ml
  # Ottoson 2011 defines water quality using normal distribution
  # The unit value for the (Cw) is cfu/100 ml
  Cw = runif(n=iter,min=1,max=235)
  
  # EcO157:H7 in the irrigated water (VRw)
  # This is determined by obtaining the ratio
  # Guidelines are based on VTEC ratio in the irrigation water (Ottoson 2011)
  # A normal distribution with truncation with exponent truncated to zero so that the ratio cannot be 1 
  # Citations: Ottoson et al. 2011
  VRw = 10^(rtruncnorm(n=iter, b=0, mean = -1.9, sd = 0.6))
  
  # EcO157:H7 concentration in water in cfu/ml (Wc)
  # This is calculated as CFU/ml
  Wc = Cw/100*VRw
  
  # Water holding capacity (W)
  # A normal distribution truncated at minimum of 0. 
  # Source citation: Hamilton et al. 2006
  W= mc2d::rtrunc(distr = rnorm,n=iter,linf = 0,mean = 0.108, sd = 0.019)
  #W=rtruncnorm(n=iter, a=0, mean = 0.108, sd = 0.019)
  
  # 2. Holding time (t.hold) in days
  # Inactivation during holding time (time.hold)
  # A triangle distribution  is assumed. FDA 2012
  t.hold= rtri(n=iter, min = 2, max = 8, mode = 4)
  
  ## B. HARVESTING
  ## 1. Harvest temperatures and time
  # This data was fitted from the time-temperature profile from industry
  # Temperature of lettuce during harvest (T.har)
  #-------------------------------------------------------------------------------
  # Change the growth temperature
  #-------------------------------------------------------------------------------
  ## 1. Harvest temperatures and time
  # This data was fitted from the time-temperature profile from industry
  # Total harvest time of lettuce (t.har)
  t.Har= mc2d::rtrunc(distr = rexp,n=iter,lsup = 8,rate = 1.113)
  # Temperature of lettuce during harvest (T.har)
  if (harvest.T == 9) {
    T.Har <- 9
  } else if (harvest.T == 17) {
    T.Har <- 17
  } else {
    T.Har <- rnorm(n = iter, mean = 13.458, sd = 4.443)
  } 
  
  # Soil Concentration (Cs)
  # A normal distribution of values with truncated values
  # Source citation: Lonehan et al. 2006
  # mean and standard deviation computed from the average and SD of values in the study by Lonehan 2006
  Cs= 10^(rtruncnorm(n=iter, a=0,b=3.67, mean = 0.928, sd = 1.11))
  Cs= 10^(mc2d::rtrunc(distr=rnorm,n=iter, linf=0,lsup=3.67, mean = 0.928, sd = 1.11))
  # VTEC in the soil
  # Risk truncate normal distribution- Ottoson et al. 2011
  # Truncated to ensure the value is a fraction less than 1
  VR.s= 10^(rtruncnorm(n=iter, b=0, mean = -1.9, sd = 0.6))
  
  # Concentration of STEC O157:H7 in the soils (Conc.s)
  Conc.s= Cs*VR.s
  
  # Attached soils on harvesting tools (Ts)
  # No distinct distribution was identified from the data sourced for Yang et al. 2012
  # Discrete distribution of the attached soil was done
  # Citation: Yang et al. 2012
  soil.size= c(31.21,31.26,30.58,22.68,21.23,22.14,10.12,9.87,10.66,0.81,0.99,1.02,0.06,0.05,0.05)
  
  p.soil.size=rep(1,length(soil.size))
  
  Ts= rdiscrete(n=iter,probs = p.soil.size,values = soil.size)
  
  # Number of E coli cells in the blades (Nb)
  Nb= Conc.s*Ts
  
  # Transfer rate from harvesting tools to lettuce (TrT)
  # The transfer rate from blade to lettuce was calculated as a proportion (%) retained in the blades
  # Distribution of both worst case and real life scenarios were combined from the data by Yang et al. 2012
  # The fraction cannot be higher than 100%
  # Nondetected values were replaced with
  # Log normal distribution of the rates (%) was selected
  # Citation: Yang et al. 2012
  TrT= 10^ (mc2d::rtrunc(distr = rnorm,n=iter,lsup = 0,mean = -2.459, sd = 2.014))
  
  ## C. TRANPORTATION TO COOLING
  # Temperature at the end of transportation of produce to processing
  # Normal distribution of temperature profile data
  # Temperature was curtailed at 20 degrees as quality declined (Tian 2014)
  T.bTran= rtruncnorm(n=iter,mean = 15.379,sd = 4.353,b=20)
  
  # Transportation temperature in degrees celsius (T.trans2)
  T.Tran1= 1/2*(T.Har+T.bTran)
  
  # Time taken for the transportation of the produce (t.tran1)
  # The lower limit for transportation time is 0
  # Maximum is 48 hours. Guidelines for Head Lettuce Production in Arizona
  t.Tran1= rtrunclnorm(n=iter,meanlog = 1.128,sdlog = 0.459,a=0,b=48)
  
  ## D. COOLING OF THE PRODUCE
  ## Time of cooling the produce
  # The minimum time is 0
  # The truncated Cauchy lies between low=0 and high =400/60
  # Use the trunction of mc2d package
  t.cool = mc2d::rtrunc(distr = rcauchy,n=iter,linf=0,lsup = (400/60),location = 0.530, scale = 0.084)
  
  # Temperature at the start of cooling (T.scool)
  # Same as the harvest temperature
  T.scool = rtruncnorm(n=iter,mean=15.403,sd=4.318,a=0)
  
  # Temperature at the end of cooling (T.ecool)
  # Log normal distribution is assumed
  T.ecool = rtrunclnorm(n=iter,meanlog = 0.380, sdlog = 0.484, a=0)  
  
  # Temperature of cooling is the average of the two (T.cool)
  T.cool = 1/2 *(T.scool+T.ecool)
  
  ## E. POSTCOOLING TRANSPORTATION AND STORAGE
  # Temperature profile data of commercial transit was used. Zeng 2011
  # Recommended transportation temperatures of 4 degrees was used here
  
  ## Post cooling storage time
  t.str1=runif(n=iter,min = 0, max = 48)
  
  ## Pre-processing storage time
  t.str2=runif(n=iter,min = 0, max = 48)
  
  ## Transportation time after cooling
  # This is assumed to added to storage time to be a maximum of 5 days based on data from partner
  t.iter= runif(n=iter,min = 0, max = 48)
  
  # Calculate the total time
  t.total= t.str1 + t.str2 + t.iter
  
  # Conditionally calculate t.tran2
  t.tran2 <- ifelse(t.total > 120, 120 - (t.str1 + t.str2), t.iter)
  
  # Transportation temperature
  T.tran2 = rlnorm(n=iter,meanlog=0.678,sdlog=0.256)
  
  #-------------------------------------------------------------------------------
  ## Specify the cold storage time
  #-------------------------------------------------------------------------------
  # Microbial die-off during post-cooling storage and transportation
  # # Microbial die-off during post-cooling storage and transportation
  if(cold.storage.time=="None"){
    t.pcool= t.tran2+t.str1+t.str2
  }else if(cold.storage.time==0){
    t.pcool=0.00001
  }else {
    t.pcool=cold.storage.time
  }
  
  ## F. PROCESSING
  # Initial prevalence (Prev0)
  # Discrete distributions from the Versar (2021)
  prev= c(0.02, 0.05, 0.075, 0.1)
  
  p.prev=rep(1,length(prev))
  
  prev0= rdiscrete(n=iter,values=prev, probs = p.prev)
  
  #Log reduction due to washing with water (Dw)
  LR.w= rpert(n=iter,min = 0.6, mode = 1, max = 1.4 )
  
  # Transfer (%) from contaminated lettuce to Flume (TR1)
  # Citation Perez 2011
  TR1= rtriang(n=iter,min=0,mode=0.01,max=0.02)
  
  # Transfer (%) from contaminated lettuce to shredder (TR2)
  # Citation Perez 2011
  TR2= rtriang(n=iter,min=0,mode=0.02,max=0.02)
  
  # Transfer (%) from contaminated lettuce to shaker (TR3)
  # Citation Perez 2011
  TR3= rtriang(n=iter,min=0,mode=0.01,max=0.02)
  
  # Transfer (%) from contaminated lettuce to centrifuge (TR4)
  # Citation Perez 2011
  TR4= rtriang(n=iter,min=0.01,mode=0.04,max=0.08)
  
  # Transfer (%) from contaminated lettuce to conveyor (TR5)
  # Citation Perez 2011
  TR5= rtriang(n=iter,min=0,mode=0.1,max=0.24)
  
  # Total transfer (%) from contaminated lettuce to facilities
  total.TR= TR1+TR2+TR3+TR4+TR5
  
  # Overall transfer coefficient (%) from facility to uncontaminated lettuce ((Ofu))
  # Citation Perez 2011
  ofu= rtri(n=iter,min=9.9,mode=15.33,max=18.83)
  
  # Percentage (%) of pathogen cells remain on facilities
  per.fac= 100- ofu
  
  ## CHLORINE WASHING
  # Die-off due to chlorine washing of the prioduce
  # An assumption that washing in chlorine works against cross-contamination of produce
  # Time taken during chlorine washing
  # cl.t= c(0.5, 0.75, 1.0)
  
  # p.cl=rep(1,length(cl.t))
  
  # t.Cl= rdiscrete(n=iter,values=cl.t, probs = p.cl)
  t.Cl= runif(n=iter, min=1, max=5)
  
  
  ## G. RETAIL STORAGE
  # Retail storage time in hours (Time.R[tR] truncated normal distribution)
  # Assumption was made that minimal growth and decay occurred during the post processing transportation
  # Citation: Ecosure 2008
  
  t.Rs= rtri(n=iter,min=0.5, mode=4, max=7)*24
  
  # Retail temperature in degrees celsius (temp.R[TR], truncated normal distribution)
  # Citation: Ecosure 2008
  # Normal distribution with values truncated at 0 was used as refrigeration at retail seldom goes below zero and beyond 20.56
  T.Rs= rtruncnorm(n=iter, a=0, b=20.56, mean = 4.4441, sd = 2.9642)
  
  ## H. RETAIL DISPLAY
  # Retail Display temperature
  # Retail display time (t.Rd in hours- Uniform distribution)
  # Random values generated using cumulative distribution function. The values are from Zeng 2011
  t.Rd= runif(n=iter, min=0,max = 72)
  
  # Retail display temperature (degrees celsius- )
  Rd.dat<-read.csv("Retail.temp.csv")
  
  # Generate random temperatures based on the observed distribution
  T.Rd <- sample(Rd.dat$Temperature, size = iter, replace = TRUE)
  
  ## I. TRANSPORTATION- RETAIL TO HOME
  # Transportation time (t.Tran2). A distribution
  # Citation: Ecosure 2008
  t.Tran4= rtrunclnorm(n=iter,meanlog = 1.421,sdlog = 0.46478,a=0.1833,b=3.8667)-0.24609
  
  # Temperature before putting in home storage (TbH)- A distribution
  # Citation: Ecosure 2008
  T.bH= rtruncnorm(n=iter,mean=8.386,sd=3.8314,a=0,b=20)
  
  # Transportation temperature in degrees celsius (temp.trans2)
  T.Tran4=1/2*(T.Rd+T.bH)
  
  ## Home storage
  # Time to first (home storage)- A distribution (time.F)
  # Convert this into hours from days by multiplying by 24
  # Citation Puillot 2010
  t.F= rweibull(n=iter,shape = 1.13,scale=2.84)*24
  
  # Time to last (Home storage)- a distribution (time.L)
  # Citation Puillot 2010
  # A truncation of consumption within 7 days was set- Arienzo 2020
  t.L= rtrunc(rweibull,n=iter, shape = 1.7,scale=7.96,lsup = 7)*24
  
  # Time selected for home storage (time.Home)
  t.H = ifelse(1/2*(t.F+t.L)>168,168,1/2*(t.F+t.L))
  
  # Home storage temperature- A distribution (temp.Home)
  # Citation: Ecosure 2008
  T.H=rtruncnorm(n=iter, a=-5, b=17.22, mean = 3.4517, sd = 2.4442)
  
  # Home storage temperature- A distribution (temp.Home)
  # Citation: Ecosure 2008
  T.H=rtruncnorm(n=iter, a=-5, b=17.22, mean = 3.4517, sd = 2.4442)
  
  ## EXPOSURE ASSESSMENT
  # Quantifies the amount of pathogen consumed
  # 1. EXPOSURE DURING HARVESTING
  # Growth parameters
  b= 0.033
  Tmin=  1.335-5.766*b
  C = 2.303
  
  # Inactivation
  # Juneja and marks model
  # Derives the values as log reduction values showing Nt/N0
  pred_jm2 = function(k1,k2,t,LRV){	#Always tails. Neg. k1 is acceptable. Decreasing k2 leads to more rapid plateauing. Juneja VK 200, eq. 1?
    #k2=abs(k2)
    out = 1/(1+exp(k1+k2*log(t)))
    if (LRV==TRUE) out=log10(out)
    return(out)
  }
  
  # Concentration after irrigation (Ic)
  
  Ic= Wc*W
  
  # Log concentration after irrigation (log.Ci)
  # Concentration in log cfu/g
  log.Ic=log10(Ic)
  
  # Die-off during holding of produce post irrigation
  # Log reduction during holding time (LR.hold)
  # Computed from the meta-analysis from primary production data
  # The log-logistic juneja and marks model
  # Read the parameters for irrigation
  hold.par= read.csv("par_irrigation.csv")
  
  list.hold= lapply(1:nrow(hold.par),function(i)#The lapply function shall be looped across the 35 pairs of k1 and k2 parameters
  {LR= pred_jm2(k1=hold.par[i,]$k1,k2=hold.par[i,]$k2,t=t.hold,LRV=TRUE)
  
  return(LR)#All the LR values are returned
  })
  
  # Convert the list all the LRs to a dataframe
  df.hold <- as.data.frame(do.call(cbind, list.hold))
  
  # Calculate the mean of LR for each
  LR.hold <- rowMeans(df.hold)
  
  # Log concentration after holding time (logC.hold) in log cfu/g
  logC.hold= log.Ic+LR.hold
  
  # Concentration after holding (C.hold) in cfu/g
  C.hold= 10^logC.hold
  
  # Transfer from harvesting tools to lettuce (Tsl)
  # Only a fraction of the cells in the harvesting tools (Nb) are transferred to 500 g of lettuce head
  # The blade does the cutting of three lettuce heads (Yang et al. 2012)
  # The counts are in cfu/g
  # It is assumed three lettuce heads weighing 500g were cut
  Tsl=Nb*TrT/1500
  
  # Concentration of E coli O157 in lettuce (C.har) after harvest (cfu/g)
  C.Har=C.hold + Tsl
  
  # log Concentration of E coli O157 in lettuce (Ch) after harvest (log cfu/g)
  logC.Har= log10(C.Har)
  
  ## 3. Microbial growth and die-off during harvesting (logD.har)
  # The log-logistic juneja and marks model was used to generate decay (LR.har)
  # A cut off of 8 degrees celsius was utilized to indicate growth and decay (Dinu 2011)
  # Growth parameters
  R.har= (b*(T.Har-Tmin))^2/C
  
  # Die-off during harvesting
  # Make the time to be in days as the computation was done with the time in Days
  # Read the parameters for irrigation
  # Die-off during harvesting
  # Make the time to be in days as the computation was done with the time in Days
  # Read the parameters for irrigation
  Har.par= read.csv("par_harvest.csv")
  
  list.Har= lapply(1:nrow(Har.par),function(i)#The lapply function shall be looped across the 35 pairs of k1 and k2 parameters
  {LR= pred_jm2(k1=Har.par[i,]$k1,k2= Har.par[i,]$k2,t=t.Har/24,LRV=TRUE)
  
  return(LR)#All the LR values are returned
  })
  
  # Convert the list all the LRs to a dataframe
  df.Har <- as.data.frame(do.call(cbind, list.Har))
  
  # Calculate the mean of LR for each
  LR.Har <- rowMeans(df.Har)
  
  # Determine the change in population during harvesting of the produce
  
  logD.Har= ifelse(T.Har>8,R.har*t.Har,LR.Har)
  
  # Final concentration of E coli in harvested produce taken to processing
  # This account for growth and decay during harvesting of the produce
  logC.HarP = logC.Har + logD.Har
  
  # 2. EXPOSURE DURING TRANSPORTATION-COOLING
  # Microbial growth and decay during transportation
  # The log-logistic juneja and marks model was used to generate decay (LR.har)
  # A cut off of 8 degrees celsius was utilized to indicate growth and decay (Dinu 2011)
  # Growth
  R.Tran1= (b*(T.Tran1-Tmin))^2/C
  
  # Die-off
  # A function of jm2 has been pasted at the beginning of the exposure assessment section
  tran.par= read.csv("par_sto.csv")
  
  list.Tran1= lapply(1:nrow(tran.par),function(i)#The lapply function shall be looped across the 35 pairs of k1 and k2 parameters
  {LR= pred_jm2(k1=tran.par[i,]$k1,k2= tran.par[i,]$k2,t=t.cool/24,LRV=TRUE)
  
  return(LR)#All the LR values are returned
  })
  
  # Convert the list all the LRs to a dataframe
  df.Tran1 <- as.data.frame(do.call(cbind, list.Tran1))
  
  # Calculate the mean of LR for each
  LR.Tran1 <- rowMeans(df.Tran1)
  
  LogD.tran1= ifelse(T.Tran1>5,R.Tran1*t.Tran1,LR.Tran1)
  
  # Log Concentration of E coli in transported produce
  # Accounts for growth or decay of e coli on the harvested produce under transport
  logC.Tran1= logC.HarP + LogD.tran1
  
  # 3. DIE-OFF DURING COOLING OF PRODUCE
  # An assumption is made that no growth occurs during cooling
  # Die-off during cooling
  # Parameters or cooling
  # The average parameters from the meta-analysis have been used
  # Use the lapply function fior the temperature ranges of 0 to 15 degrees celsius
  # The data has 35 values for each k1 and k2 parameters
  # A function of jm2 has been pasted at the beginning of the exposure assessment section
  cool.par= read.csv("par_cool.csv")
  
  list.cool= lapply(1:nrow(cool.par),function(i)#The lapply function shall be looped across the 35 pairs of k1 and k2 parameters
  {LR= pred_jm2(k1=cool.par[i,]$k1,k2= cool.par[i,]$k2,t=t.cool/24,LRV=TRUE)
  
  return(LR)#All the LR values are returned
  })
  
  # Convert the list all the LRs to a dataframe
  df.cool <- as.data.frame(do.call(cbind, list.cool))
  
  # Calculate the mean of LR for each
  LR.cool <- rowMeans(df.cool)
  
  #Growth during cooling of the produce
  # Concentration of e coli in the produce after cooling
  # An assumption was made that during cooling no growth took place
  D.cool= LR.cool
  
  #Log concentration of E coli on produce after cooling (logC.cool)
  logC.cool= logC.Tran1 + D.cool
  
  # 4. DIE-OFF DURING POST-COOLING STORAGE AND TRANSPORTATION
  # Die-off during post cooling
  tran.par= read.csv("par_sto.csv")
  
  list.Tran2= lapply(1:nrow(tran.par),function(i)#The lapply function shall be looped across the 35 pairs of k1 and k2 parameters
  {LR= pred_jm2(k1=tran.par[i,]$k1,k2= tran.par[i,]$k2,t=t.cool/24,LRV=TRUE)
  
  return(LR)#All the LR values are returned
  })
  
  # Convert the list all the LRs to a dataframe
  df.Tran2 <- as.data.frame(do.call(cbind, list.Tran2))
  
  # Calculate the mean of LR for each
  LR.Tran2 <- rowMeans(df.Tran2)
  
  #Change in population during transporation to after cooling
  D.Tran2= LR.Tran2
  
  #Concentration of  STEC in the produce at after decay post cooling
  logC.tran2= logC.cool + D.Tran2
  #-------------------------------------------------------------------------------
  ## VBNC switch rates
  #-------------------------------------------------------------------------------
  
  # Changes in the concentration due to the VBNC fraction
  # Add the VBNC fraction to the culturable fraction
  # Concentration of postcooled produce after accounting for VBNC
  # Add the fraction of VBNC to the counts
  if(harvest.T==0 & cold.storage.time==0){
    f.vbnc = (10^logC.tran2)*0.000
  }else if(harvest.T==9 & cold.storage.time==0){
    f.vbnc = (10^logC.tran2)*0.000355
  }else if(harvest.T==9 & cold.storage.time==24){
    f.vbnc = (10^logC.tran2)*0.0431
  }else if(harvest.T==9 & cold.storage.time==48){
    f.vbnc = (10^logC.tran2)*0.00169
  }else if(harvest.T==9 & cold.storage.time==72){
    f.vbnc = (10^logC.tran2)*0.0534
  }else if(harvest.T==9 & cold.storage.time==96){
    f.vbnc = (10^logC.tran2)*0.000442
  }else if(harvest.T==9 & cold.storage.time==120){
    f.vbnc = (10^logC.tran2)*0.00166
  }else if(harvest.T==17 & cold.storage.time==0){
    f.vbnc = (10^logC.tran2)*0.0282
  }else if(harvest.T==17 & cold.storage.time==24){
    f.vbnc = (10^logC.tran2)*0.0202
  }else if(harvest.T==17 & cold.storage.time==48){
    f.vbnc = (10^logC.tran2)*0.0447
  }else if(harvest.T==17 & cold.storage.time==72){
    f.vbnc = (10^logC.tran2)*0.00531
  }else if(harvest.T==17 & cold.storage.time==96){
    f.vbnc = (10^logC.tran2)*0.0196
  }else if(harvest.T==17 & cold.storage.time==120){
    f.vbnc = (10^logC.tran2)*0.0142
  }
  
  # Concentration of postcooled produce after accounting for VBNC
  # Add the fraction of VBNC to the counts
  
  LogC.cool2= log10((10^logC.tran2) + f.vbnc)
  
  # 5. DIE-OFF DURING PROCESSING
  # I. WASHING OF THE PRODUCE
  # Log Concentration of STEC after washing
  logC.w= LogC.cool2 - LR.w
  
  # Concentration after washing (Cw)
  C.w= 10^logC.w
  
  # Concentration of E coli j unit batch after washing (N.int)
  # Make the Prev0 prevalence a proportion for it is a percentage
  N.int= C.w * (prev0 * 0.01)
  
  # II. CROSS-CONTAMINATION DURING PROCESSING
  # CFU left on originally contaminated portion in (CFU/batch)-res.frac
  # Transform the total.TR into a proportion since it is a percentage
  res.frac= N.int*(1-(total.TR*0.01))
  
  # CFU getting into newly contaminated produce (new.cont)
  new.cont= N.int* (total.TR * 0.01) * (ofu * 0.01)  
  
  # CFU in the final batch after cross contamination (N.final)
  N.final=res.frac + new.cont
  
  # Spread of contamination due to processing (S)
  # Citation FDA 2012
  S= rpert(n=iter,min=1,mode=1.2,max=2)
  
  # Prevalence after cross contamination (Prev.f %)-it is a percentage (%)
  Prev.f=prev0*S
  
  # Concentration of E coli on lettuce due to cross contamination
  # Compute the % to proportion
  C.pw= N.final/(Prev.f * 0.01)
  
  # Log concentration of E coli in washed produce + coss contamination
  logC.pw= log10(C.pw)
  
  # III DIE-OFF DUE TO CHLORINE WASHING
  # Die-off during chlorine washing
  # Source citation. Madamba 2022b.
  # Critical point for deactivation due to chlorine was established as >5ppm
  # The log-logistic juneja and marks model used to show inactivation due to chlorine
  # The log-logistic juneja and marks model used to show inactivation due to chlorine
  Cl.par<-read.csv("par_Cl.csv")
  
  # An assumption that temperature was maintained at refrigerated conditions
  list.Cl= lapply(1:nrow(Cl.par),function(i)#The lapply function shall be looped across the 35 pairs of k1 and k2 parameters
  {LR= pred_jm2(k1=Cl.par[i,]$k1,k2= Cl.par[i,]$k2,t=t.Cl/24,LRV=TRUE)
  
  return(LR)#All the LR values are returned
  })
  
  # Convert the list all the LRs to a dataframe
  df.Cl <- as.data.frame(do.call(cbind, list.Cl))
  
  # Calculate the mean of LR for each
  LR.Cl <- rowMeans(df.Cl)
  
  # Change in STEC in leafy greens after washing
  logD.Cl= LR.Cl
  
  # Level of E coli in lettuce after washing (Clw)
  logC.p= logC.pw + logD.Cl
  
  # 6. DIE-OFF AND GROWTH DURING STORAGE AND TRANSPORTATION
  # I. RETAIL STORAGE
  # Storage and transportation experience growth and decay
  # Growth
  R.Rs= (b*(T.Rs-Tmin))^2/C
  
  # Die-off
  # We import all the values of the k1 and k2 parameters that range between
  sto.par= read.csv("par_sto.csv")
  
  # Use the lapply function fior the temperature ranges of 0 to 15 degrees celsius
  # The data has 35 values for each k1 and k2 parameters
  # A function of jm2 has been pasted at the beginning of the exposure assessment section
  list.Rs= lapply(1:nrow(sto.par),function(i)#The lapply function shall be looped across the 35 pairs of k1 and k2 parameters
  {LR= pred_jm2(k1=sto.par[i,]$k1,k2=sto.par[i,]$k2,t=t.Rs/24,LRV=TRUE)
  
  return(LR)#All the LR values are returned
  })
  
  # Convert the list all the LRs to a dataframe
  df.Rs <- as.data.frame(do.call(cbind, list.Rs))
  
  # Calculate the mean of LR for each
  LR.Rs <- rowMeans(df.Rs)
  
  # Change of microbial growth during retail storage
  LogD.Rs= ifelse(T.Rs>5,R.Rs * t.Rs,LR.Rs)
  
  # Log concentration after retail storage
  LogC.Rs= logC.p + LogD.Rs
  
  # II. RETAIL DISPLAY
  # Growth and die-off during retail display
  # Juneja and Marks log logistic decay model was used to get the die-off
  # Growth
  R.Rd= (b*(T.Rd-Tmin))^2/C
  
  # Die-off during retail display
  # A function of jm2 has been pasted at the beginning of the exposure assessment section
  list.Rd= lapply(1:nrow(sto.par),function(i)#The lapply function shall be looped across the 35 pairs of k1 and k2 parameters
  {LR= pred_jm2(k1=sto.par[i,]$k1,k2=sto.par[i,]$k2,t=t.Rd/24,LRV=TRUE)
  
  return(LR)#All the LR values are returned
  })
  
  # Convert the list all the LRs to a dataframe
  df.Rd <- as.data.frame(do.call(cbind, list.Rd))
  
  # Calculate the mean of LR for each
  LR.Rd <- rowMeans(df.Rd)
  
  # Indicate inactivation or growth
  LogD.Rd= ifelse(T.Rd>5,R.Rd*t.Rd,LR.Rd)
  
  # Log concentration after retail display
  LogC.Rd= LogC.Rs + LogD.Rd
  
  # III. TRANSPORTATION TO HOME
  # Growth and die-off during transportation to home
  # Juneja and Marks log logistic decay model was used to get the die-off
  # Growth
  R.Tran4= (b*(T.Tran4-Tmin))^2/C
  
  # Die-off
  tran.par= read.csv("par_sto.csv")
  
  list.Tran4= lapply(1:nrow(tran.par),function(i)#The lapply function shall be looped across the 35 pairs of k1 and k2 parameters
  {LR= pred_jm2(k1=tran.par[i,]$k1,k2= tran.par[i,]$k2,t=t.cool/24,LRV=TRUE)
  
  return(LR)#All the LR values are returned
  })
  
  # Convert the list all the LRs to a dataframe
  df.Tran4 <- as.data.frame(do.call(cbind, list.Tran4))
  
  # Calculate the mean of LR for each
  LR.Tran4 <- rowMeans(df.Tran4)
  
  #Change in population during transporation to after cooling
  LogD.Tran4= ifelse(T.Tran4>5,R.Tran4*t.Tran4,LR.Tran4)
  
  # Log concentration after transportation to Home
  LogC.Tran4= LogC.Rd + LogD.Tran4
  
  # IV. HOME STORAGE
  # Microbial die-off or growth in distribution and on transit and Home storage
  # Growth and die-off parameters
  # Growth
  R.H= (b*(T.H-Tmin))^2/C
  
  # Use the lapply function fior the temperature ranges of 0 to 15 degrees celsius
  # The data has 35 values for each k1 and k2 parameters
  # A function of jm2 has been pasted at the beginning of the exposure assessment section
  list.H= lapply(1:nrow(sto.par),function(i)#The lapply function shall be looped across the 35 pairs of k1 and k2 parameters
  {LR= pred_jm2(k1=sto.par[i,]$k1,k2=sto.par[i,]$k2,t=t.H/24,LRV=TRUE)
  
  return(LR)#All the LR values are returned
  })
  
  # Convert the list all the LRs to a dataframe
  df.H <- as.data.frame(do.call(cbind, list.H))
  
  # Calculate the mean of LR for each
  LR.H <- rowMeans(df.H)
  
  # Change of microbial growth during retail storage
  LogD.H= ifelse(T.H>5,R.H*t.H,LR.H)
  
  # Log concentration after home storage (LogC.Hc)
  LogC.Hc= LogC.Tran4 + LogD.H
  
  # Limit the level of contamination to <log 7
  LogC.H= ifelse(LogC.Hc<7,LogC.Hc,7)
  
  # Concentration of STEC in food after home storage
  C.H=10^LogC.H
  
  # 7. SERVING
  # Serving sizes (ser)
  # The size of the serving in grams
  # The serving sizes are based on the recommended reference amounts
  # Source citation: USFDA 2002. Reference amounts customarily consumed per eating occasion
  ser= 70
  #Source Citation:
  # Dose of pathogen taken in a serving (CFU/ serving)
  # Obtained by multiplying the concentration after home storage with the serving sizes
  # The dose is in cfu/serving
  D= C.H*ser
  
  ## K. DOSE RESPONSE MODEL
  # beta poisson dose response model was the most fitted
  # The estimated risk of illness (Danyluk and Schaffner 2011)
  B= 229.2928 #Dose response parameter.
  
  #Obtain this value from literature. This is the beta value
  
  a= 0.267#Dose response parameter alpha
  
  ## L. RISK CHARACTERIZATION
  # Risk of illness per serving without accounting for prevalence
  # The risk is calculated risk of illness/serving/person/day
  P.ill=(1-(1+D/B)^(-a))
  
  # Risk of illness accounting for prevalence
  P=(1-(1+D/B)^(-a))
  
  processed[[paste0("P.ill_", harvest.T, "_", cold.storage.time)]] <- P.ill
  ## Consumption
  # Compute the annual risk of illness
  # Total US population (N.pop)
  # Total population of persons in the United States
  # The estimate was as at January 2024
  # Source citation: US Census Beuareau
  N.pop= 335893089 ## Find out how many eat lettuce
  
  # Annual lettuce consumption in the United States in grams/year
  # Calculated from annual leetice availability per capita (A)
  # The data is from 2020:
  # Source citation: USDA 2023: Food Availability (Per Capita) Data System
  A= 5195.123
  
  # Number of consumed servings per person per year
  # This is determined by dividing the per capita consumption by the serving size
  NP= A/ser 
  
  ## ANNUAL RISK
  # The risk outcomes are independent of each other
  ## ANNUAL RISK
  # The risk outcomes are independent of each other
  P.A.ill=1-((1-P.ill)^NP)
  
  a_processed[[paste0("P.A.ill_",harvest.T, "_", cold.storage.time)]] <- P.A.ill
  
  # Sensitivity analysis
  # Create a data frame of the distributions in the model
  sensitivity.dat<-data.frame(as.numeric(Cw),as.numeric(VRw),as.numeric(W),as.numeric(t.hold),as.numeric(Cs),
                              as.numeric(VR.s),as.numeric(TrT),as.numeric(T.Tran1),as.numeric(t.Tran1),as.numeric(t.cool),
                              as.numeric(T.cool),as.numeric(t.pcool),as.numeric(prev0),as.numeric(LR.w),as.numeric(TR1),as.numeric(TR2),
                              as.numeric(TR3),as.numeric(TR4),as.numeric(TR5),as.numeric(ofu),as.numeric(t.Cl),
                              as.numeric(t.Rs),as.numeric(T.Rs),as.numeric(t.Rd),as.numeric(T.Rd),as.numeric(t.Tran4),
                              as.numeric(T.Tran4),as.numeric(t.H),as.numeric(T.H),as.numeric(T.Har),as.numeric(t.Har),
                              as.numeric(S),as.numeric(P.A.ill))
  
  names(sensitivity.dat)<-c("Ec_Water","STEC_ratio_water","Water_holding_capacity","Post_irrigation_Holding_time",
                            "Ec_soil","STEC_ratio_soil","Soil_trasfer.rt","Tranport_temp_cool",
                            "Tranport_time_cool","Cooling_time","Cooling_temp","Postcooling_time",
                            "Prevalence","LR_washing","Flume_transfer.rt","Shredder_transfer.rt","Shaker_transfer.rt",
                            "Centrifuge_transfer.rt","conceyor.transfer","Facility_transfer.rt","Chlorine_washing_time",
                            "Retail_storage_time","Retail_storage_temp","Retail_display_time",
                            "Retail_display_temp","Transport_home_time","Transport_home_temp",
                            "Home_storage_time","Home_storage_temp","Harvesting_temp",
                            "Harvesting_holding_time","Spreading","Prob.Risk")
  
  # Conduct random forest
  rf <- randomForest::randomForest(Prob.Risk~., data=sensitivity.dat, proximity=TRUE,importance=T)#Create the random forest
  
  #Plot variable imporance
  rf.data<-as.data.frame(rf$importance)
  names(rf.data)<-c("%IncMSE","IncNodePurity")
  rf.data$condition<-row.names(rf.data)
  rf.data$`%IncMSE`<-rf.data$`%IncMSE`
  
  # Save the random forest plot
  rf.plot<-ggplot(rf.data[rf.data$`%IncMSE`>4.00e-5,], aes(x = reorder(condition, `%IncMSE`), y = `%IncMSE`)) +
    geom_bar(stat = "identity", fill = "grey") +  # Adjust fill color
    geom_text(label="A",y=0.20,x=1,size=5,family = "serif",col = "black")+
    labs(x = "Conditions", y = "% MSE",title = "Shredded and packaged lettuce") +        # Label axes
    coord_flip() +                                # Flip the coordinates for horizontal bars
    theme_light()+                               # Use a minimal theme (you can customize as needed)
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=20,family = "serif",colour = "black"),
          axis.text =  element_text(size=20,family = "serif",colour = "black"),
          text = element_text(size=20,family = "serif",colour = "black"))
  
  # Save the sensitivity analysis plot
  # Create directory if it doesn't exist
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  
  # Filename with extension
  filename <- paste0(path, "/processed_", harvest.T, "_", cold.storage.time, ".png")
  ggsave(filename = filename,plot = rf.plot, width = 8, 
         height = 6, units = "in", dpi = 300)

}
