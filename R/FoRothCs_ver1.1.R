##########################################################################################
## ForRothCs model v1.0
## Simulation for Cs (& C) cycling in forest ecosystem (espesially for artificial forest)
## Coded by K. Nishina (National Institute for Environmental Studies)
## Contact adresss; nishina.kazuya@nies.go.jp
## Detail infromation in Nishina & Hayashi (2015) in Frontiers in Environmental Sciences
##########################################################################################

# This model is based on Roth-C model (Jenkinson et al., 1996) and DDPS (Hayashi et al., 2002&2006)
# Roth-C model is carbon cycling model for Soil, and DDPS is self-thinning growth model for artificial forest.
# ForRothCs was coupled these two models and Cs dynamics is incorporated using transfer parameters obtained by previous literatures.


RothCs <- function( 
  #Simulation setup
  Nyear = 30, #Input simulation duration
  Depth = 20, #Input soil depth (cm)
  
  # Site location
  LAT = 37.6, # Latitude degree
  
  #Input monthly average temperature, precipitation, evaporation, soil surface condition: bare (0) or covered (1)   
  ### For Soil Carbon and Cs dynamics ###
  Temp = c(2.7,  3.6,  7,	12.6,	17.3,	20.6,	24.1,	25.5,	22,	16.1,	10.1, 4.9), 
  Prep = c(35.5,  44.9,	85,	101.1,	121.8,	131.1, 140.4,	141.8, 176,	155.8, 68.2, 39.2),
  bare = c(1, 1, 1, 1, 1, 1.0, 1.0, 1.0, 1, 1, 1, 1),
  Clay = 20.4, # %clay contents
  BD = 0.6, # Bulk density (Mg/m3)
  Till = 1, # %tillage factor (Nishina et al., in prep)
  inputC = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), #Input carbon in each month from external management
  
  # Initial for Soil Carbon fraction
  init_c = rbind(0.05, 0.44, 0.05, 1.4), #Input initial SOC value (DPM, RPM, BIO, HUM)
  
  ### For Forest Dynamics
  init_cs = rbind(10, 10, 0, 10, 100, 100), #Input initial 137Cs (RPM, DPM, BIO, HUM, Mineral soil, and Tree Leave)
  Cs_dep = rbind(1500, 600, 0.5, 100, 0, 1000, 300, 0), #Input 137Cs deposition (RPM, DPM, BIO, HUM, Mineral soil, Tree Leave, Branch, Stem)
  
  Season = c(0.03,0.01,0.01,0,0,0,0,0.15,0.1,0.1,0.01,0.01), #Litterfall fraction to Leaf biomass in each month 
  
  Transfac = 0.00005, #7.2e-2, ## Aggregated Transfactor of Cs for Tree (Bq/kg(Tree)/Bq/m2(Soil))
  Litterfac = 0.0042, ## Trans factor from soil to litter # 0.0042 from Fukuyama & Takenaka 2004
  Trlocate = 0.3, ## Pullback Cs in Leaves before litter falling
  Relocfac = 0.2, ## Translocation Cs from stem to Leaves
  
  init_BA = 50, # m2/ha or C/ha 
  init_age = 40, # Tree age in 2011 (2011 for Fukushima)  
  init_rho = 1600, # Number of tree / ha
  
  R_s_f = 0.22,  #Slow to fast ratio (slow/fast) in leacheable fraction by TF Loffredo et al. 2014
  
  ## Management options
  
  thin_year = 0, # Thinning year after deosition
  M_thin = 10, # Thinning month
  thin_rate = 0,  # Thinning rate 0-1
  
  
  litter_year = 0,  # Litter Removal year after 2011
  M_rm = 10,
  lit_removal_rate = 0,  # Litter Removal rate 0-1
  
  int_start =FALSE,
  start_year = 0,  # Litter Removal year after 2011
  start_month = 9,
  
  ## plot 
  autoplot = TRUE
){  
  
  ####################################
  ##Simulation for tree growth model##
  ####################################
  
  qden <- numeric(Nyear)
  qden[1] <- qden_t <- init_rho
  baha <- baha2 <- numeric(Nyear)
  baha[1] <- baha2[1] <- init_BA
  qden[1] <- init_rho
  thin <- rep(1, Nyear)
  thin[thin_year] <- thin_rate
  thin[thin_year] <- (1 - thin_rate)
  
  ##### Tree growth model loop #####
  
  
  for(i in 1:Nyear){
    k <- i+init_age
    logq_i <- log(2582.2) - 1/2.92*log((baha[i]/85.15)^(2.92*10.82) + (2582.2/init_rho)^2.92)   
    qden_t <- qden[i]
    
    while(qden_t > exp(logq_i)){ 
      #qden_t <- print(qden_t)* thin[i] - 1  
      qden_t <- qden_t - 1
    }  
    
    fp <- exp((1/(k+1) - 1/k)*(-8.76))*((k+1)/k)^0.908*(qden_t/qden[i])^0.709
    baha[i+1] <- baha[i] * fp * thin[i]
    #qden[i+1] <- qden_t
    qden_t2 <- qden_t * thin[i]
    qden[i+1] <-  qden_t2  
    
    #print(c(i, fp, baha[i], qden[i], exp(logq_i), qden_t), digit=4)
  }
  
  ##### Conversion from basal area to "height (no use)", "Leaf_mass", "Stem mass", ""　by allometry equation
  
  a_height <- -0.7331 + 0.01676 * c(init_age:(init_age+Nyear)) + 0.000218 * init_rho + 0.6708 * mean(baha) # Hayashi 2006 
  height <- 1/(1/(a_height*10) + 1/45)  # Hmax 40m from 資料 in 福島スギ２ 
  
  #Leaf_mass <- exp(-3.6153)*(sqrt(baha/qden/3.14)*2*100)^2.0929 * qden/1000  #葉量変換 from 根岸1988　東大演習林報   t/ha
  #Leaf_mass <- 0.0019*(sqrt(baha/qden/3.14)*2*100)^2.6924 * qden/1000 #葉量変換 from 渡邉2002　岐阜県スギ   t/ha
  #Leaf_mass <- 0.0840*((sqrt(baha/qden/3.14)*2*100)^2*height)^0.5791 * qden/1000/10 #葉量変換 from  田内・字都木　全国スギ   t/ha -> kg/m2
  Leaf_mass <- 0.004327*((sqrt(baha/qden/3.14)*2*100))^2.61 * qden/10000 #葉量変換 from  梶本2014    kg/tree -> kg/m2
  
  #Branch_mass <- exp(-2.2678)*(sqrt(baha/qden/3.14)*2*100)^2.4822 * qden/1000  #枝重量変換 from 根岸1988　東大演習林報 
  #Branch_mass <- 0.0069*(sqrt(baha/qden/3.14)*2*100)^2.4838 * qden/1000  #枝重量変換 from 渡邉2002　岐阜県スギ t/ha
  #Branch_mass <- 0.0032*((sqrt(baha/qden/3.14)*2*100)^2*height)^0.900 * qden/1000/10  #枝重量変換 from  田内・字都木　全国スギ   t/ha -> kg/m2
  Branch_mass <- 0.000436*((sqrt(baha/qden/3.14)*2*100))^3.17 * qden/10000  #枝重量変換 from  梶本2014  kg/tree -> kg/m2
  
  #Stem_mass <-  exp(-1.4278)*(sqrt(baha/qden/3.14)*2*100)^2.6188 * qden/1000  #ステム変換 from 根岸1988　東大演習林報 
  #Stem_mass <- 0.01633*((sqrt(baha/qden/3.14)*2*100)^2*height)^0.8577 * qden/1000  #ステム変換 from 渡邉2002　岐阜県スギ t/ha
  Stem_mass <- 0.0308*((sqrt(baha/qden/3.14)*2*100)^2*height)^0.9106 * qden/1000/10  #ステム変換 from 田内・字都木　全国スギ   t/ha -> kg/m2
  
  #qdiff_loop <- rep(1- c(qden[2:(length(qden))], qden[(length(qden))])/qden[1:(length(qden))], each=12)/12  # Montly mortality rate (year to month by 12)  
  
  #Leaf_m_loop_t <- rep(Leaf_mass[1:(length(Stem_mass)-1)], each=12)
  #Leaf_m_loop_t[((1:length(Leaf_m_loop_t))-1)%%12 != 0] <- NA
  #Branch_m_loop_t <- rep(Branch_mass[1:(length(Stem_mass)-1)], each=12)
  #Branch_m_loop_t[((1:length(Branch_m_loop_t))-1)%%12 != 0] <- NA
  #Stem_m_loop_t <- rep(Stem_mass[1:(length(Stem_mass)-1)], each=12)
  #Stem_m_loop_t[((1:length(Stem_m_loop_t))-1)%%12 != 0] <- NA
  
  #Leaf_m_loop <- c(Leaf_mass[1], rep(Leaf_mass[1:(length(Stem_mass)-1)], each=12))    # Check later
  #Branch_m_loop <- c(Branch_mass[1], rep(Branch_mass[1:(length(Stem_mass)-1)], each=12))
  #Stem_m_loop <- c(Stem_mass[1], rep(Stem_mass[1:(length(Stem_mass)-1)], each=12))
  qdiff_loop <- rep(1- c(qden[2:(length(qden))], qden[(length(qden))])/qden[1:(length(qden))], each=12)/12  # Montly mortality rate (year to month by 12)
  
  w_growth <- (Temp/sum(Temp))  # Seasonality for tree growth
  w_growth_leaf <- c(0,0,Temp[3:6]/sum(Temp[3:6]),0,0,0,0,0,0)  # Seasonality for leaf growth (according to Miyaura et al. 1995 in Bulletin of the Nagoya University Forests )
  wm_litterfall <- colSums(sapply(Leaf_mass[-1], FUN = "*", Season, simplify = TRUE))
  litterfall_season <- c(0, as.numeric(sapply(Leaf_mass[-1], FUN = "*", Season, simplify = TRUE)))
  
  Leaf_m_loop <- c(Leaf_mass[1], Leaf_mass[1] + cumsum(as.vector(sapply(diff(Leaf_mass) + wm_litterfall, FUN = "*", w_growth_leaf, simplify = TRUE)))) - cumsum(litterfall_season)
  Branch_m_loop <- c(Branch_mass[1], Branch_mass[1] + cumsum(as.vector(sapply(diff(Branch_mass), FUN = "*", w_growth, simplify = TRUE))))
  Stem_m_loop <- c(Stem_mass[1], Stem_mass[1] + cumsum(as.vector(sapply(diff(Stem_mass), FUN = "*", w_growth, simplify = TRUE))))
  
  
  for(i in 2:(12*Nyear)){
    if(((i-2)%/%12+1) ==  thin_year){
      qdiff_loop[i-1] <- 0  
    }
    if(((i-2)%/%12+1) ==  thin_year && (i-2)%%12+1 >= M_thin){#(i%%12+1) ==  M_thin
      
      Leaf_m_loop[i] <- Leaf_mass[(i-2)%/%12+2]   # Check later
      Branch_m_loop[i] <- Branch_mass[(i-2)%/%12+2]
      Stem_m_loop[i] <- Stem_mass[(i-2)%/%12+2]
      
      print(c(i, (i-2)%/%12+1, (i-2)%%12+1, (i-2)%/%12+2))
    }
  }
  
  
  for(i in 1:(12*Nyear)){
    if(((i-1)%/%12+1) ==  thin_year && (i-2)%%12+1 == M_thin){       #(i%%12+1) ==  M_thin
      Volume_thin_l <- Leaf_mass[(i-1)%/%12+2]/Leaf_mass[(i-1)%/%12+1]   # Check later
      Volume_thin_b <- Branch_mass[(i-1)%/%12+2]/Branch_mass[(i-1)%/%12+1]
      Volume_thin_s <- Stem_mass[(i-1)%/%12+2]/Stem_mass[(i-1)%/%12+1]
      
      if(Volume_thin_l > 1){
        Volume_thin_l <- Volume_thin_b <- Volume_thin_s <- 1
      }
      
      print(c(i, ((i-1)%/%12+1),Volume_thin_l, Volume_thin_b, Volume_thin_s))    
    }
  }
  
  
  StBr_m_loop <- Branch_m_loop + Stem_m_loop 
  Total_m_loop <- Leaf_m_loop + Branch_m_loop + Stem_m_loop
  
  
  ############################################################################
  #### For Potential Evaporation estimation (Thornthwaite equation (1948)) ####
  ############################################################################
  
  init <- init_c*10 # Conversion from kg-C/m2 to t/ha (<- model internal unit)
  
  ND <- 0:364 # Julius day
  ND_ve <- 80 # vernal equinox
  monthday <- c(1,31,28,31,30,31,30,31,31,30,31,30,31)
  
  ### for caluculation of day length ###
  
  SLD <- 23.4*sin((ND-ND_ve)/364*pi*2) # SLD Solar declination: degree ## Pending!!
  
  SE_md <- sin(LAT/180*pi)*sin(SLD/180*pi) + cos(LAT/180*pi)*cos(SLD/180*pi) # solar elevation; degree
  SE_md <- asin(SE_md)
  HA <- acos(-tan(LAT/180*pi)*tan(SLD/180*pi)) # hour angle from sunrise to midday; degree
  HA <- HA/2/pi*360
  DL <- 24*HA/180 # day length; hour
  
  L_month <- numeric(12)  ## monthly mean day length
  for(i in 1:12) L_month[i] <- mean(DL[monthday[i]:monthday[i+1]])
  
  ### Thornthwaite equation for PE ###
  Temp2 <- Temp
  Temp2[Temp2<0] <- 0
  Iheat <- sum((Temp2/5)^1.514)
  alphaI <- (6.75e-7)*Iheat^3 - (7.71e-5)*Iheat^2 + (1.792e-2)*Iheat + 0.49239
  Evapo <- 16*(L_month/12)*(monthday[2:13]/30)*((10*Temp2/Iheat)^alphaI)
  
  
  ######################################
  ##Simulation for Roth-C and Cs model##
  ######################################
  ## For Roth-C model and parameters, please see Jenkinson et al (1996)
  
  # DPM : RPM ratio
  f_dpm = 0.2      # From Jenkinson et al (1996)
  f_rpm = 1 - f_dpm
  
  ##### Definition of limitting functions #####
  
  #calculation of IOM by faloon et al.
  TC <- sum(init)
  IOM <- 0.049*TC^1.139
  
  #Input C rep
  Resid <- rep(inputC, Nyear)
  
  #calculation of partitioning between CO2 evolved and (BIO+HUM) formed
  
  CBH <- 1.67*(1.85+1.6*exp(-0.0786*Clay))
  C <- 1/(1+CBH)
  
  #Related to moisture deficit
  
  maxTSMD <- -(20.0 + 1.3*Clay - 0.01*Clay^2)*Depth/23
  accTSMD <- Prep - 0.75*Evapo
  accTSMD <-  replace(accTSMD, which(accTSMD > 0), 0)
  accTSMD <-  replace(accTSMD, which(accTSMD < maxTSMD), maxTSMD)
  
  #calculation of modifying factors
  
  ifelse(maxTSMD > 0.444*accTSMD,
         modB <- 1,
         modB <- 0.2 + (1-0.2)*(maxTSMD - accTSMD)/(maxTSMD - 0.444*accTSMD)
  )
  
  modA <- 47.91/(1 + exp(106.06/(Temp+18.27)))
  modC <- bare
  
  #Modified function from Roth C model
  
  modfac <- modA * modB * modC
  modfac <- rep(modfac, Nyear)
  
  
  #Plant dynamics for Evergreen and Prep
  
  fac_seo <- rep(Season, Nyear)
  Prep_long <- rep(Prep, Nyear)
  
  #Prepareing Matrix for simulation results
  
  Yn1 <- matrix(NA, nrow=4, ncol=(1+12*Nyear))
  Yn2 <- matrix(NA, nrow=8, ncol=(1+12*Nyear))
  Yn1[,1] <- rbind(init)
  Yn2[,1] <- rbind(init_cs, 0.1, 0.1)
  
  
  ###########################
  ##### Simulation loop #####
  ###########################
  
  for(i in 1:(12*Nyear)){
    
    # monthly decomposition limitting functions for each soil fraction
    
    fa <- (exp(-10/12*modfac[i]+log(Till))) 
    fb <- (exp(-0.3/12*modfac[i]+log(Till)))
    fc <- (exp(-0.66/12*modfac[i])+log(Till))
    fd <- (exp(-0.02/12*modfac[i])+log(Till))
    
    # For litterfall fraction  
    fs <- fac_seo[i]
    
    ##Transition matrix##
    
    X1 <-  c(fa, 0, 0, 0)
    X2 <-  c(0,  fb,  0,	0)
    X3 <-  c(C*(1-fa)*0.46,	C*(1-fb)*0.46,	fc + C*(1-fc)*0.46,	C*(1-fd)*0.46)
    X4 <-  c(C*(1-fa)*0.54, 	C*(1-fb)*0.54, 	C*(1-fc)*0.54, 	fd + C*(1-fd)*0.54)
    
    
    # For translocation from stem (or branch) to leaves
    
    Reloc_leaf  <- Relocfac*(Yn2[7,i]/Branch_m_loop[i])*Leaf_m_loop[i]/Yn2[7,i]
    
    Reloc_branch  <- (Yn2[8,i]/Stem_m_loop[i] - Yn2[7,i]/Branch_m_loop[i])*Branch_m_loop[i]/Yn2[8,i]
    #if(Reloc_branch < 0) Reloc_branch  <- Branch_m_loop[i]/Stem_m_loop[i] #0.1*(1/(Yn2[7,i]/Branch_m_loop[i])-1/(Yn2[8,i]/Stem_m_loop[i]))*Branch_m_loop[i]/Stem_m_loop[i]
    #if(Reloc_branch < 0) Reloc_branch  <- 10*((Yn2[8,i]/Stem_m_loop[i])/(Yn2[7,i]/Branch_m_loop[i]))*Branch_m_loop[i]/Yn2[8,i]
    if(Reloc_branch < 0) Reloc_branch  <- 0.5*((Yn2[8,i]/Stem_m_loop[i]))*Branch_m_loop[i]/Yn2[8,i]
    #Reloc_branch  <- 5*(Yn2[8,i]/Stem_m_loop[i] - Yn2[7,i]/Branch_m_loop[i])*Branch_m_loop[i]/Yn2[8,i]
    #if(Reloc_branch < 0) Reloc_branch  <- 0.1*(Yn2[8,i]/Stem_m_loop[i])*Branch_m_loop[i]/Stem_m_loop[i]
    #if(Reloc_branch < 0) Reloc_branch  <- Relocfac #0.1*(1/(Yn2[7,i]/Branch_m_loop[i])-1/(Yn2[8,i]/Stem_m_loop[i]))*Branch_m_loop[i]/Stem_m_loop[i]
    #if(Reloc_branch < 0) Reloc_branch  <- 10*(Yn2[8,i]/Stem_m_loop[i]/(Yn2[7,i]/Branch_m_loop[i]))*Branch_m_loop[i]/Yn2[8,i]
    Reloc_stem  <- 1 - qdiff_loop[i]  - Reloc_branch   
    
    #print(c(Reloc_branch, 0.1*(Yn2[8,i]/Stem_m_loop[i])*Branch_m_loop[i]/Stem_m_loop[i]))
    #print(c(Reloc_leaf, Reloc_branch, Reloc_stem))
    
    ### Throughfall model proposed by Loffredo et al. 2014 in Science of tghe Total Environment ###
    # Simple Monthly Throughfall function 
    
    k1 <- 5.0E-04*31 # kinetic for slow fraction
    k2 <- 1.2E-02*31 # kinetic for fast fraction
    # R_f_s <- 0.22 # e.g. Slow 52.6 kBq/m2; Fast 237.8 kBq/m2
    b1 <- 0.0172
    TFmonth <- b1* Prep_long[i]/(1+1/R_s_f) *(k1/R_s_f*exp(-k1*i) + k2*exp(-k2*i)) # Deriv; A1*(1-exp(-k1*i)) + A2*(1-exp(-k2*i))
    
    
    ## Transition matrix for SOC dynamics (Xn1) and Cs dynamics (Xn2)  ##
    
    Cs1 <- c(fa, 0, 0, 0, Litterfac*f_dpm, fs*(1-Trlocate)*f_dpm + qdiff_loop[i]*f_dpm + TFmonth*f_dpm, fs*f_dpm + qdiff_loop[i]*f_dpm, qdiff_loop[i]*f_dpm)
    Cs2 <- c(0, fb, 0, 0, Litterfac*f_rpm, fs*(1-Trlocate)*f_rpm + qdiff_loop[i]*f_rpm + TFmonth*f_rpm, fs*f_rpm + qdiff_loop[i]*f_rpm, qdiff_loop[i]*f_rpm)
    Cs3 <- c(C*(1-fa)*0.46,  C*(1-fb)*0.46,  fc + C*(1-fc)*0.46,  C*(1-fd)*0.46-Transfac*(Leaf_m_loop[i]+Branch_m_loop[i]+Stem_m_loop[i]), 0, 0, 0, 0)  
    Cs4 <- c(C*(1-fa)*0.54,   C*(1-fb)*0.54,   C*(1-fc)*0.54,   fd + C*(1-fd)*0.54, 0, 0, 0, 0)
    CsM <- c((1-C)*(1-fa),  (1-C)*(1-fb),  (1-C)*(1-fc), (1-C)*(1-fd), 1-Transfac*(Leaf_m_loop[i]+Branch_m_loop[i]+Stem_m_loop[i]) - Litterfac, 0, 0, 0)
    L_Cs <- c(0, 0, 0, Transfac*(Leaf_m_loop[i]), Transfac*(Leaf_m_loop[i]), 1 - fs - qdiff_loop[i] - TFmonth, Reloc_leaf, 0)
    B_Cs <- c(0, 0, 0, Transfac*(Branch_m_loop[i]), Transfac*(Branch_m_loop[i]), fs*Trlocate*Branch_m_loop[i]/StBr_m_loop[i], 1 - fs - qdiff_loop[i] - Reloc_leaf, Reloc_branch)
    S_Cs <- c(0, 0, 0, Transfac*(Stem_m_loop[i]), Transfac*(Stem_m_loop[i]), fs*Trlocate*Stem_m_loop[i]/StBr_m_loop[i], 0, Reloc_stem)
    
    
    Xn1 <- as.matrix(rbind(X1, X2, X3, X4))
    Xn2 <- as.matrix(rbind(Cs1, Cs2, Cs3, Cs4, CsM, L_Cs, B_Cs, S_Cs))
    
    #print(colSums(Xn2)) # check conservation of matter
    
    ## Monthly Carbon input from external dose and tree litter fall ## 
    
    Cinput <- rbind((Resid[i] + 0.4*Leaf_m_loop[i]*(fs + qdiff_loop[i]) + 0.4*Branch_m_loop[i]*(fs + qdiff_loop[i])) * f_dpm, 
                    (Resid[i] + 0.4*Leaf_m_loop[i]*(fs + qdiff_loop[i]) + 0.4*Branch_m_loop[i]*(fs + qdiff_loop[i])) * f_rpm + 0.4*Stem_m_loop[i]*qdiff_loop[i], 
                    0,
                    0)
    
    ## Thinning and Litter removal ##
    
    if(((i-2)%/%12+1) ==  thin_year && (i-2)%%12+1 == M_thin){       #(i%%12+1) ==  M_thin
      Yn2[6,i] <- Yn2[6,i] * Volume_thin_l
      Yn2[7,i] <- Yn2[7,i] * Volume_thin_b
      Yn2[8,i] <- Yn2[8,i] * Volume_thin_s
    }
    
    if(((i-2)%/%12+1) ==  thin_year && (i-2)%%12+1 == M_thin){      
      Yn1[c(1,2),i] <- Yn1[c(1,2),i] * (1 - lit_removal_rate)
      Yn2[c(1,2),i] <- Yn2[c(1,2),i] * (1 - lit_removal_rate)
    }    
    
    
    
    
    #### Simulation transition matrix ####
    Yn1[,i+1] <- Xn1 %*% Yn1[,i]  + Cinput
    Yn2[,i+1] <- Xn2 %*% Yn2[,i]
    
    
    #### Deposition in 3/11 ####    
    if(int_start == FALSE){    
      if(i == 2){
        Yn2[,i+1] <- Yn2[,i+1] + Cs_dep
      }
    }
    
    #### Start from arbitrary date ####    
    
    if(int_start == TRUE){    
      if(((i-1)%/%12+1) ==  start_year && (i-2)%%12+2 == start_month){
        Yn2[,i+1] <- Cs_dep + Yn2[,i]
      }
    }
    
    ## 137Cs decay function
    Yn2[,i+1] <- Yn2[,i+1]*exp(-log(2)/12/30.167)
    
    ## 134Cs decay function
    #if(Cs134 == 1){
    #Yn2[,i+1] <- Yn2[,i+1]*exp(-log(2)/12/2.06)
    #}    
    
  }
  
  #For output data format 
  Yn1 <- as.data.frame(t(Yn1))*1000/10000 # convert from t/ha to kg/m2
  Yn2 <- as.data.frame(t(Yn2))
  #Yn <- rbind(Yn1, Yn2)  
  #Yres <- as.data.frame(t(Yn))
  #names(Yres) <- c("DPM", "RPM", "BIO", "HUM", "DPM_Cs", "RPM_Cs", "BIO_Cs", "HUM_Cs", "Min_Cs", "L_Cs", "B_Cs", "S_Cs")
  names(Yn1) <- c("DPM", "RPM", "BIO", "HUM")
  names(Yn2) <- c("DPM_Cs", "RPM_Cs", "BIO_Cs", "HUM_Cs", "Min_Cs", "L_Cs", "B_Cs", "S_Cs")
  
  #Yres$Total <- (Yn[4, ] + Yn[3, ] + Yn [2, ] + Yn [1, ]) + IOM
  #Yres$Total_Cs <- (Yn[5, ] + Yn[6, ] + Yn [7, ] + Yn [8, ] + Yn [9, ] + Yn [10, ] + Yn [11, ] + Yn [12, ])
  
  Yn1$Soil_TC <- (Yn1[, 1] + Yn1[, 2] + Yn1[, 3] + Yn1[, 4]) + IOM/10
  Yn2$Total_Cs <- (Yn2[, 1] + Yn2[, 2] + Yn2[, 3] + Yn2[, 4] + Yn2[, 5] + Yn2[, 6] + Yn2[, 7] + Yn2[, 8])
  
  
  Tres <- as.data.frame(cbind(baha, qden, Leaf_mass, Stem_mass, Branch_mass, height), 
                        row.names=c("Basal_area_m2_ha", "Tree_Density_ha", "Leaf_mass_kg_m2", "Stem_mass_kg_m2", "Branch_mass_kg_m2", "Height_m"))
  
  Yn3 <- data.frame(Leaf_act = Yn2$L_Cs/Leaf_m_loop, Branch_act = Yn2$B_Cs/Branch_m_loop, Stem_act = Yn2$S_Cs/Stem_m_loop, 
                    Litter_act = (Yn2$DPM + Yn2$RPM)/((Yn1$DPM + Yn1$RPM)/0.3), Soil_act = (Yn2$HUM_Cs + Yn2$Min_Cs)/(BD*1000*Depth/100))
  
  
  #ResAll <- list(Tree = Tres, C_kg_m2 = Yn1, Cs_bq_m2 = Yn2)
  ResAll <- list(Tree = Tres, C_kg_m2 = Yn1, Cs_bq_m2 = Yn2, Cs_bq_kg = Yn3)
  
  
  
  #######################################################
  #######   Auto plot functions for 137Cs ###############
  #######################################################
  
  if(autoplot == 1){
    
    start_xmin <- 0
    if(int_start == TRUE) start_xmin <- start_month/12 + (start_year-1)
    
    x_tmp <- c(seq(start_xmin, (12*Nyear+1)/12, 1/12), rev(seq(start_xmin, (12*Nyear+1)/12, 1/12)))
    x_start_tmp <- start_xmin/(1/12)
    x_num_tmp <- length(x_tmp)/2
    end_tmp <- (12*Nyear+1)/12
    end_num_tmp <- length(1:(12*Nyear+1)/12)
    
    
    #par(mfrow=c(2,2), mar=c(3,4,0.5,0), oma=c(0.5,0.5,0.5,0.5), mgp=c(1.75,0.75,0))
    #par(mar=c(3,4,0.5,0), oma=c(0.5,0.5,0.5,0.5), mgp=c(1.75,0.75,0))
    # Tree dynamics
    #plot(0:Nyear, Tres$Leaf_mass, type="n", xlim=c(start_xmin, Nyear+0.1), ylim=c(0, max(Tres$Stem_mass+Tres$Branch_mass+Tres$Leaf_mass)), 
    #     ylab= expression(paste("Above-ground biomass", " [kg ", m^{-2}, "]")),
    #     xaxs = "i", xlab="Year")
    #polygon(c(0:Nyear, Nyear:0), c(Tres$Stem_mass+Tres$Branch_mass+Tres$Leaf_mass, rep(0, Nyear+1)), col="darkolivegreen")
    #polygon(c(0:Nyear, Nyear:0), c(Tres$Stem_mass+Tres$Branch_mass, rep(0, Nyear+1)), col="darkolivegreen3")
    #polygon(c(0:Nyear, Nyear:0), c(Tres$Stem_mass, rep(0, Nyear+1)), col="darkolivegreen1")
    
    #legend("bottomright", pch=15, y.intersp=0.8, bg="white", 
    #       legend=c("Leaves", "Branch", "Stem"),
    #       col=c("darkolivegreen", "darkolivegreen3", "darkolivegreen1"))
    
    #mtext(side=3, "(a)", adj=0, line=0.1)
    
    #SoilC
    #plot(1:(12*Nyear+1)/12, Yn1$Soil_TC, type="n", xlim=c(start_xmin, Nyear+0.1), ylim=c(0, max(Yn1$Soil_TC)+0.3), 
    #     ylab=expression(paste("Litter & Soil C", " [kg-C ", m^{-2}, "]")), 
    #     xaxs = "i", xlab="Year")
    #polygon(c(1:(12*Nyear+1)/12, ((12*Nyear+1):1)/12), c(Yn1$DPM + Yn1$RPM + Yn1$HUM + Yn1$BIO, rep(0, length((12*Nyear+1):1))), col="orange")
    #polygon(c(1:(12*Nyear+1)/12, ((12*Nyear+1):1)/12), c(Yn1$RPM + Yn1$HUM + Yn1$BIO, rep(0, length((12*Nyear+1):1))), col="orange")
    #polygon(c(1:(12*Nyear+1)/12, ((12*Nyear+1):1)/12), c(Yn1$HUM + Yn1$BIO, rep(0, length((12*Nyear+1):1))), col="tan3")
    #polygon(c(1:(12*Nyear+1)/12, ((12*Nyear+1):1)/12), c(Yn1$BIO, rep(0, length((12*Nyear+1):1))), col="grey")  
    
    #legend("topright", pch=15, y.intersp=0.8, bty="n",
    #       legend=c("Litter (DPM+RPM)", "Humus", "Microbe"),
    #       col=c("orange", "tan3", "grey"))
    #mtext(side=3, "(b)", adj=0, line=0.1)
    
    # All Cs
    #plot(1:(12*Nyear+1)/12, Yn2$Total_Cs, type="n", xlim=c(start_xmin, Nyear+0.1), ylim=c(0, sum(Cs_dep)+300), 
    #     ylab=expression(paste({}^{137}, "Cs inventory", " [Bq ", m^{-2}, "]")), #"137Cs inventory [Bq/m2]",
    #     xaxs = "i", xlab="Year")
    #polygon(x_tmp, c(Yn2$Total_Cs[x_start_tmp:end_num_tmp], rep(0, x_num_tmp)), col="darkolivegreen")
    ##polygon(x_tmp, c((Yn2$HUM_Cs + Yn2$Min_Cs + Yn2$RPM_Cs + Yn2$DPM_Cs + Yn2$L_Cs + Yn2$B_Cs + Yn2$S_Cs + Yn2$BIO_Cs)[x_start_tmp:end_num_tmp], rep(0, x_num_tmp)), col="darkolivegreen2")
    #polygon(x_tmp, c((Yn2$HUM_Cs + Yn2$Min_Cs + Yn2$RPM_Cs + Yn2$DPM_Cs + Yn2$BIO_Cs + Yn2$S_Cs + Yn2$B_Cs)[x_start_tmp:end_num_tmp], rep(0, x_num_tmp)), col="darkolivegreen3")
    #polygon(x_tmp, c((Yn2$HUM_Cs + Yn2$Min_Cs + Yn2$RPM_Cs + Yn2$DPM_Cs + Yn2$BIO_Cs + Yn2$S_Cs)[x_start_tmp:end_num_tmp], rep(0, x_num_tmp)), col="darkolivegreen1")
    #polygon(x_tmp, c((Yn2$HUM_Cs + Yn2$Min_Cs + Yn2$RPM_Cs + Yn2$DPM_Cs + Yn2$BIO_Cs)[x_start_tmp:end_num_tmp], rep(0, x_num_tmp)), col="orange")
    #polygon(x_tmp, c((Yn2$HUM_Cs + Yn2$Min_Cs)[x_start_tmp:end_num_tmp], rep(0, x_num_tmp)), col="tan3")
    #polygon(x_tmp, c((Yn2$BIO_Cs + Yn2$Min_Cs)[x_start_tmp:end_num_tmp], rep(0, x_num_tmp)), col="grey")
    #polygon(x_tmp, c(Yn2$Min_Cs[x_start_tmp:end_num_tmp], rep(0, x_num_tmp)), col="burlywood2")
    
    #mtext(side=3, "(c)", adj=0, line=0.1)  
    #abline(h=c(sum(Cs_dep))/2, col=2)
    #legend("topright", pch=15, 
    #       legend=c("Stem", "Branch", "Leaves", "Litter", "Humus", "Mineral", "Microbe"),
    #       col=c("darkolivegreen1", "darkolivegreen3", "darkolivegreen", "orange", "tan3", "burlywood2", "grey"))
    
    #legend("topright", pch=15, y.intersp=0.8, bty="n",
    #       legend=c("Leaves", "Branch", "Stem", "Litter", "Humus", "Microbe", "Mineral"),
    #       col=c("darkolivegreen", "darkolivegreen3", "darkolivegreen1", "orange", "tan3", "grey", "burlywood2"))
    
    #legend("topright", pch=15, ncol=2, cex=1.1, y.intersp=0.8, 
    #       legend=c("Stem", "Branch", "Leaves", "Litter", "Humus", "Mineral", "Microbe"),
    #       col=c("darkolivegreen1", "darkolivegreen3", "darkolivegreen", "orange", "tan3", "burlywood2", "grey"))
    
    # Tree Cs inventory
    #plot(1:(12*Nyear+1)/12, Yn2$Total_Cs, type="n", xlim=c(start_xmin, Nyear+0.1), ylim=c(0, max(Yn2$L_Cs)), ylab="137Cs inventory [Bq/m2]", xaxs = "i", xlab="Year")
    #polygon(c(1:(12*Nyear+1)/12, ((12*Nyear+1):1)/12), c(Yn2$B_Cs + Yn2$S_Cs + Yn2$L_Cs, rep(0, length((12*Nyear+1):1))), col="darkolivegreen")
    #polygon(c(1:(12*Nyear+1)/12, ((12*Nyear+1):1)/12), c(Yn2$B_Cs + Yn2$S_Cs, rep(0, length((12*Nyear+1):1))), col="darkolivegreen3")
    #polygon(c(1:(12*Nyear+1)/12, ((12*Nyear+1):1)/12), c(Yn2$S_Cs, rep(0, length((12*Nyear+1):1))), col="darkolivegreen1")
    
    
    # Tree Cs conc
    plot(1:(12*Nyear+1)/12, Yn2$Total_Cs, type="n", xlim=c(start_xmin, Nyear+0.1), ylim=c(0.5, max(Yn3$Litter_act, Yn3$Leaf_act)), 
         ylab=expression(paste({}^{137}, "Cs activity", " [Bq ", kg^{-1}, "]")), 
         xaxs = "i", xlab="Year", log="y")
    #plot(1:(12*Nyear+1)/12, Yn2$Total_Cs, type="n", xlim=c(start_xmin, Nyear+0.1), ylim=c(0.1, 10000), ylab="137Cs activity [Bq/kg]", xaxs = "i", xlab="Year", log="y")
    points(0:1000, 10000*(exp(-log(2)/30.167))^(0:1000), type="l", col="grey80")
    points(0:1000, 1000*(exp(-log(2)/30.167))^(0:1000), type="l", col="grey80")
    points(0:1000, 100*(exp(-log(2)/30.167))^(0:1000), type="l", col="grey80")
    points(0:1000, 10*(exp(-log(2)/30.167))^(0:1000), type="l", col="grey80")
    points(0:1000, 1*(exp(-log(2)/30.167))^(0:1000), type="l", col="grey80")
    points(0:1000, 0.1*(exp(-log(2)/30.167))^(0:1000), type="l", col="grey80")
    
    points((1:(12*Nyear+1)/12)[(x_start_tmp+1):end_num_tmp], Yn3$Leaf_act[(x_start_tmp+1):end_num_tmp], col="darkolivegreen", type="l", lwd=2)
    points((1:(12*Nyear+1)/12)[(x_start_tmp+1):end_num_tmp], Yn3$Branch_act[(x_start_tmp+1):end_num_tmp], col="darkolivegreen3", type="l", lwd=2)
    points((1:(12*Nyear+1)/12)[(x_start_tmp+1):end_num_tmp], Yn3$Stem_act[(x_start_tmp+1):end_num_tmp], col="darkolivegreen1", type="l", lwd=2)
    points((1:(12*Nyear+1)/12)[(x_start_tmp+1):end_num_tmp], Yn3$Litter_act[(x_start_tmp+1):end_num_tmp], col="orange", type="l", lwd=2, lty=1)
    points((1:(12*Nyear+1)/12)[(x_start_tmp+1):end_num_tmp], Yn3$Soil_act[(x_start_tmp+1):end_num_tmp], col="brown", type="l", lwd=2, lty=5)
    
    legend("bottomright", pch=-1, lty=c(1,1,1,1,5), lwd=2, y.intersp=0.8, ncol=3, bty="n",
           legend=c("Leaves", "Branch", "Stem", "Litter", "Soil"),
           col=c("darkolivegreen", "darkolivegreen3", "darkolivegreen1", "orange", "brown"))
    
    #mtext(side=3, "(d)", adj=0, line=0.1)
    
  }
  return(ResAll)
}



##### Usage demo #####

#RothCs()  ## Default setting 

#  test2 <- RothCs(Nyear=50, init_age=45, init_rho=1600, Temp = c(c(0,0.2,3.3,9.5,14.6,18.1,21.6,23.4,19.1,13.1,7.2,2.4)),
#                  init_BA=75.4, Cs_dep= rbind(5000*1.6/(1.6+2.71), 5000*2.71/(1.6+2.71), 0.5, 100, 0, 6000, 0, 0), Litterfac=0.004,
#                  thin_year=1, thin_rate=0, M_thin=10)


#test2 <- RothCs(Nyear=50, init_age=45, init_rho=1600, Temp = c(c(0,0.2,3.3,9.5,14.6,18.1,21.6,23.4,19.1,13.1,7.2,2.4)),
#                init_BA=75.4, Cs_dep= rbind(5000*1.6/(1.6+2.71), 5000*2.71/(1.6+2.71), 0.5, 100, 0, 6000, 0, 0), Litterfac=0.004,
#                thin_year=10, thin_rate=0.7, M_thin=8)


#test2 <- RothCs(Nyear=30, init_age=45, init_rho=836, Temp = c(c(0,0.2,3.3,9.5,14.6,18.1,21.6,23.4,19.1,13.1,7.2,2.4)),
#                init_BA=iniBA, Cs_dep= rbind(3000*1.6/(1.6+2.71), 3000*2.71/(1.6+2.71), 0.5, 100, 0, 5000, 0, 0), Litterfac=0.004,
#                thin_year=2, thin_rate=0.5, M_thin=3)
