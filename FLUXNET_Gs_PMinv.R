# Calculate hourly surface conductance
#####
rm(list = ls())
library(lubridate)
###############
# define functions
###############
# slope of the Ta-VPD curve (KPa K-1)
delta_fun <- function(Ta){
  # input - air temperature, (C)
  # output - slope of T-es curve (KPa*K-1)
  h2osat_fun <- function(x){
    h2o_base = 0.61121
    h2o_mult = 17.502 
    h2o_add = 240.97
    esT = h2o_base*exp((h2o_mult*x)/(h2o_add+x))
    return(esT)
  }
  psi_fun <- function(x,y){
    psi_mult = 240.97*17.502
    psi_add = 240.97
    delta = (psi_mult*x)/((y+psi_add)^2)
    return(delta)
  }
  #
  esT = h2osat_fun(Ta) # h2o saturation
  delta = psi_fun(esT, Ta)
  return(delta)
}

# Calculate aerodynamic conductance (m s-1)
Ga_calfromflx_fun <- function(u, ustar){
  ra = u/(ustar^2) + 6.2*ustar^(-2/3)
  Ga = 1/ra;
  return(Ga)
}

# Calculate surface conductance from FLUXNET data (m s-1)
Gs_calfromflx_fun <- function(VPD, A, ro, Cp, Ga, delta, gamma, LEobs){
  # calculate surface conductance by using Penman-Monteith equation
  # Estimate the Gcx by inverting P-M equation
  # VPD - vapour pressure deficit kPa
  # A - available energy W*m-2
  # ro - air density
  # Cp - 
  # Ga - aerodynamic conductance 
  # delta - 
  # gamma - 
  # LEobs - LE observations W*m-2
  #
  ratio = delta/gamma
  Gs = LEobs*Ga/(ratio*A - (ratio+1)*LEobs + Ga*ro*Cp*VPD/gamma)
  return(Gs)
}


################
# end of definitions of functions
################

#######
# define constant
# filling values
mfill = -9999
# psychrometric constant
gamma = 0.066 # kPa*K-1
# air density
ro = 1.15 # kg*m-3
# specific heat capacity of air
Cp = 1012 # J*kg-1*K-1
# gas constant 
Rgas <- 8.3144598 # J mol-1 K-1

###
# file paths and names
infilepath1 = 'F:/FLUXNET2015_Nov/'
infilepath4 = 'F:/FLUXNET2015_PATCH1/' # patch for flux data quality controls
infilename2 <- 'G:/FLUXNET2015_siteinfo.csv'

outfilepath1 <- 'G:/FLUXNET2015_Gs/'

###
# site information
# read site information in a file
siteinfo <- read.table(infilename2, header=T, sep=',', stringsAsFactors=F)
sitepft <- siteinfo$pft
siteid <-siteinfo$site
nsite = length(siteid)

###

time.start <- Sys.time()
###
# save peak-gpp month data
df.out <- data.frame()
for(i in 2:nsite){
  
  ######
  # surface conductance under dry periods
  ######
  # read hourly flux data
  flxfiledir1 = list.files(infilepath1, pattern = paste0('FLX_',siteid[i],'.*'))
  # no such site file, exit and go to next step
  if (isTRUE(all.equal(flxfiledir1, character(0)))) next 
  flxfiledir2hr = list.files(paste0(infilepath1, flxfiledir1), pattern = '.*_FULLSET_H{1}')
  flxT_hr = read.csv(paste0(infilepath1,flxfiledir1,'/',flxfiledir2hr))
  ## get timestamp of hourly data
  ## First get this because 'flx_year' will be used
  time_hr = flxT_hr$TIMESTAMP_START
  flx_year_hr = floor(time_hr/100000000)
  flx_month_hr = floor((time_hr-flx_year_hr*100000000)/1000000)
  flx_day_hr <- floor((time_hr-flx_year_hr*100000000-flx_month_hr*1000000)/10000)
  flx_year = unique(flx_year_hr)  
  nyear = length(flx_year) # number of years
  # get hourly steps
  time_step = (time_hr[2]-time_hr[1])/100
  n_steps_day = ifelse(time_step==1, 24, 48)
  # get site PFT
  
  ######
  # read daily flux data
  flxfiledir2d <- list.files(paste0(infilepath1, flxfiledir1), pattern = '.*_FULLSET_DD')
  flxT_d <- read.csv(paste0(infilepath1,flxfiledir1,'/',flxfiledir2d))
  # get timestamp of daily data
  time_d <- flxT_d$TIMESTAMP
  flx_year_d <- floor(time_d/10000)
  flx_month_d <- floor((time_d - flx_year_d*10000)/100)
  
  ######
  # read monthly flux data
  flxfiledir2mo <- list.files(paste0(infilepath1, flxfiledir1), pattern = '.*_FULLSET_MM')
  flxT_mo <- read.csv(paste0(infilepath1,flxfiledir1,'/',flxfiledir2mo))
  # timestamp of monthly data
  time_mo <- flxT_mo$TIMESTAMP
  flx_year_mo <- floor(time_mo/100)
  flx_month_mo <- (time_mo-flx_year_mo*100)
  
  #######
  # Deal with QC problem in heat flux data
  # Heat flux QC of hourly data
  hflxdir1 = list.files(infilepath4, pattern = paste0('FLX_',siteid[i],'_FLUXNET2015_PATCH1_H{1}'))
  hflxqc_hr = read.csv(paste0(infilepath4, hflxdir1))
  timeqc_hr <- hflxqc_hr$TIMESTAMP_START
  Hqc_hr = hflxqc_hr$H_F_MDS_QC
  LEqc_hr = hflxqc_hr$LE_F_MDS_QC
  # years of hourly QC
  qc_year_hr <- floor(timeqc_hr/100000000)
  
  # Heat flux QC data in daily scale
  hflxdir2 <- list.files(infilepath4,pattern=paste0('FLX_',siteid[i],'_FLUXNET2015_PATCH1_DD'))
  hflxqc_d <- read.csv(paste0(infilepath4, hflxdir2))
  timeqc_d <- hflxqc_d$TIMESTAMP
  Hqc_d <- hflxqc_d$H_F_MDS_QC
  LEqc_d <- hflxqc_d$LE_F_MDS_QC
  # years of daily QC
  qc_year_d <- floor(time_d/10000)
  
  ######
  # 'clean' data 
  # e.g. site 'NL-Loo' data lacks the last record
  # A ROBUST way to clean data
  # assume: problem with the last year data
  # hourly data: 24 or 48 data in the last day
  # daily data: 31 days in December
  # monthly data: December is the last month
  # annual data: QC data and original data have the same size
  ndif_hr <- n_steps_day-
    length(which(flx_year_hr==flx_year[nyear]&flx_month_hr==12&flx_day_hr==31))#
  if (ndif_hr>0) {
    # less than 24 or 48 
    # add rows
    newrows <- array(dim = c(ndif_hr,dim(flxT_hr)[2]))
    newrows <- as.data.frame(newrows)
    names(newrows) <- names(flxT_hr) # have the same names
    flxT_hr <- rbind(flxT_hr, newrows)
  }
  #
  ndif_d <- 31 - 
    length(which(flx_year_d==flx_year[nyear]&flx_month_d==12))
  if (ndif_d>0) {
    newrows <- as.data.frame(array(dim = c(ndif_d, dim(flxT_d)[2])))
    names(newrows) <- names(flxT_d)
    flxT_d <- rbind(flxT_d, newrows)
  }
  #
  ndif_mo <- 12 - 
    length(which(flx_year_mo==flx_year[nyear]))
  if (ndif_mo>0) {
    newrows <- as.data.frame(array(dim = c(ndif_mo, dim(flxT_mo)[2])))
    names(newrows) <- names(flxT_mo)
    flxT_mo <- rbind(flxT_mo, newrows)
  }
  #
  # same as length of FLUXNET hourly data
  Hqc_hr <- Hqc_hr[1:length(flxT_hr[,1])]
  LEqc_hr <- LEqc_hr[1:length(flxT_hr[,1])]
  # same as length of FLUXNET daily data
  Hqc_d <- Hqc_d[1:length(flxT_d[,1])]
  LEqc_d <- LEqc_d[1:length(flxT_d[,1])]
  
  ######
  # read yearly flux data

  
  #####
  # hourly data: obtain measured and good-quality data
  # carbon fluxes
  # NEE/NEP
  NEE_hr = flxT_hr$NEE_VUT_REF
  NEEqc_hr = flxT_hr$NEE_VUT_REF_QC
  NEE_hr[NEE_hr==mfill|NEEqc_hr>1] = NA
  NEP_hr = -NEE_hr
  # GPP
  GPP_hr = flxT_hr$GPP_NT_VUT_REF
  GPPqc_hr = flxT_hr$NEE_VUT_REF_QC
  GPP_hr[GPP_hr==mfill|GPPqc_hr>1] = NA
  # RE
  RE_hr = flxT_hr$RECO_NT_VUT_REF
  REqc_hr = flxT_hr$NEE_VUT_REF_QC
  RE_hr[RE_hr==mfill|REqc_hr>1] = NA
  
  # water fluxes
  # LE
  LE_hr = flxT_hr$LE_F_MDS
  LE_hr[LE_hr==mfill|LEqc_hr>1] = NA
  
  # sensible heat fluxes (H)
  H_hr = flxT_hr$H_F_MDS
  H_hr[H_hr==mfill|Hqc_hr>1] = NA
  
  # radiation fluxes
  SW_hr = flxT_hr$SW_IN_F_MDS
  SWqc_hr = flxT_hr$SW_IN_F_MDS_QC
  SW_hr[SW_hr==mfill|SWqc_hr>1] = NA
  
  # soil heat flux
  G_hr = flxT_hr$G_F_MDS
  Gqc_hr = flxT_hr$G_F_MDS_QC
  G_hr[G_hr==mfill|Gqc_hr>1] = NA
  
  # net radiation
  Rn_hr = flxT_hr$NETRAD
  Rn_hr[Rn_hr==mfill] = NA
  if (isTRUE(all.equal(Rn_hr,logical(0)))) {
    print(paste0(siteid[i],' no net radiation'))
    next
  }
  
  # air temperature, C
  Ta_hr = flxT_hr$TA_F_MDS
  Taqc_hr = flxT_hr$TA_F_MDS_QC
  Ta_hr[Ta_hr==mfill|Taqc_hr>1] = NA
  
  # VPD, hPa
  VPD_hr = flxT_hr$VPD_F_MDS
  VPDqc_hr = flxT_hr$VPD_F_MDS_QC
  VPD_hr[VPD_hr==mfill|VPDqc_hr>1] = NA
  VPD_hr = VPD_hr/10 # KPa
  
  # precipitation
  P_hr = flxT_hr$P_F
  Pqc_hr = flxT_hr$P_F_QC
  
  # soil moisture (NOTE: this could be NULL)
  SWC_hr = flxT_hr$SWC_F_MDS_1
  SWCqc_hr = flxT_hr$SWC_F_MDS_1_QC
  SWC_hr[SWC_hr==mfill|SWCqc_hr>1] = NA
  ###
  # Sites have no SWC records
  if (isTRUE(all.equal(SWC_hr,logical(0)))) {
    # create a vector
    SWC_hr <- rep(NA, length(NEP_hr))
  }
  
  # wind speed
  u_hr <- flxT_hr$WS_F
  uqc_hr <- flxT_hr$WS_F_QC
  u_hr[u_hr==mfill|uqc_hr>0] = NA
  
  # friction velocity
  ustar_hr = flxT_hr$USTAR
  ustar_hr[ustar_hr==mfill] = NA
  
  # CO2 concentration
  CO2_hr <- flxT_hr$CO2_F_MDS
  co2qc_hr <- flxT_hr$CO2_F_MDS_QC
  CO2_hr[CO2_hr==mfill|co2qc_hr>1] <- NA
  
  # nighttime NEE
  
  
  ##########
  # daily precipitation, mm
  P_d <- flxT_d$P_F
  P_d[P_d==mfill] <- NA
  
  # daily soil water content, %
  SWC_d <- flxT_d$SWC_F_MDS_1
  SWCqc_d <- flxT_d$SWC_F_MDS_1_QC
  SWC_d[SWC_d==mfill|SWCqc_d<0.7] <- NA
  
  # daily carbon flux
  # GPP
  GPP_d <- flxT_d$GPP_NT_VUT_REF
  GPPqc_d <- flxT_d$NEE_VUT_REF_QC
  GPP_d[GPP_d==mfill] <- NA
  
  # daily heat fluxes
  # latent heat flux
  LE_d <- flxT_d$LE_F_MDS
  LE_d[LE_d==mfill|LEqc_d<0.7] <- NA
  
  
  ##########
  # monthly CO2 concentration from FLUXNET
  CO2_mo <- flxT_mo$CO2_F_MDS
  CO2_mo[CO2_mo==mfill] <- NA
  
  # monthly GPP
  GPP_mo <- flxT_mo$GPP_NT_VUT_REF
  GPP_mo[GPP_mo==mfill] <- NA
  
  # monthly ET/LE
  LE_mo <- flxT_mo$LE_F_MDS
  LEqc_mo <- flxT_mo$LE_F_MDS_QC
  LE_mo[LE_mo==mfill|LEqc_mo<0.7] <- NA
  
  # monthly SWC
  SWC_mo <- flxT_mo$SWC_F_MDS_1
  SWCqc_mo <- flxT_mo$SWC_F_MDS_1_QC
  SWC_mo[SWC_mo==mfill] <- NA
  
  #####
  # Calculate surface conductance
  delta <- delta_fun(Ta_hr)
  Ga_hr <- Ga_calfromflx_fun(u_hr,ustar_hr)
  Gs_hr <- Gs_calfromflx_fun(VPD_hr,Rn_hr,ro,Cp,Ga_hr,delta,gamma,LE_hr)
  # post-processing, data selection
  # Gs_hr_neg <- Gs_hr # negative Gs
  # Gs_hr_neg[Gs_hr_neg>=0] = NA
  # Gs_hr_pos <- Gs_hr # positive Gs
  # Gs_hr_pos[Gs_hr_pos<0] = NA
  Gs_hr_filter <- Gs_hr # data selection
  Gs_hr_filter[SW_hr<100|is.na(SW_hr)] <- NA # daytime data
  Gs_hr_filter[GPP_hr<5|is.na(GPP_hr)] <- NA # positvie carbon uptake
  Gs_hr_filter[Ta_hr<5|is.na(Ta_hr)] <- NA # extreme low air temperature
  Gs_hr_filter[Gs_hr_filter<=0.000001] <- NA # unrealistic values
  Gs_hr_filter[(Rn_hr-LE_hr)<5] <- NA # Rn generally larger than LE
  
  #########
  # Unit transformation m s-1 to mol m-2 s-1
  # GPP - umolCO2 m-2 s-1
  # VPD - KPa
  # Ca - umol mol-1
  # Gs - m s-1, turn into mol m-2 s-1
  # P - Pa
  # Rgas - J mol-1 K-1
  # Tk - K
  Tk_hr <- Ta_hr+273.15 # from C to K
  PA_hr <- flxT_hr$PA_F*1000 # from kPa to Pa
  PA_hr[PA_hr==mfill] <- NA
  Ca_hr <- flxT_hr$CO2_F_MDS
  Ca_hr[Ca_hr==mfill] <- NA
  Gs_mol <- Gs_hr_filter*PA_hr/(Rgas*Tk_hr) 
  
  ########
  # remove impact of Precipitation
  P_flag <- P_d
  P_flag[] <- 0 # initiate an array for flagging
  P1d_index <- which(P_d>0.1&!is.na(P_d))
  P_flag[P1d_index] <- 1
  P2d_index <- P1d_index+1 # remove the day after the rainy day
  P2d_index[P2d_index>length(P_d)] <- length(P_d) # in case beyond the boundary
  P_flag[P2d_index] <- 1
  # from Daily to Hourly
  P_flag_hr <- rep(P_flag, each=n_steps_day)
  Gs_mol[P_flag_hr==1] <- NA

  
  ############
  ## save results
  ## save estimated surface conductances for each site
  ## create a data frame
  time_start <- flxT_hr$TIMESTAMP_START
  # time_end <- flxT_hr$TIMESTAMP_END
  Gs_df <- data.frame(time_start,Ta_hr,u_hr,ustar_hr,VPD_hr,Rn_hr,ro,Cp,Ga_hr,delta,gamma,LE_hr, # variables for Gs estimation
                      Gs_mol,NEE_hr,GPP_hr,RE_hr,SWC_hr) # Gs and other variables
  outfilename1 <- paste0(outfilepath1, siteid[i], '_Gs_hr.csv')
  write.table(Gs_df, outfilename1, row.names = F, sep = ',')

}
time.end <- Sys.time()