#'@title Compute the water leaving reflectance from ASD data
#'
#'
#'This function compute the water reflectance (rho_w = pi * Rrs) for a water surface.
#'The sky glint is removed using various methods as detailed in \code{\link{ASD.go}.
#'
#'
#' @param raw.asd is a long list returned by \code{\link{merge.ASD.radiances.for.rhow}} containing all
#'the data necessary to compute the diffuse water reflectance. That include the raw radiance data,
#'the sun-viewing geometry, the location, the time and the wind speed.
#' @param rho.panel is the reflectivity of the spectralon panel. The Default is 0.985
#' @param quantile.prob is a value (between 0.25 to 1) for the maximum quantile probability
#'for which the surface radiance (Lt) values will be discarded.
#'The quantile probability is the value at which the probability
#'of the random variable being less than or equal to this value.
#'For example, a value of 0.5 means that every Lt value above the
#'50\% quantile probability will be discarded, while a value of 0.8 means
#'that only the values above the 80\% quantile probability will be discarded.
#'The rational is to eliminate outliers resulting from sun glitters, foam, sea ice, etc.
#'The default value is 0.5 to eliminate Lt value above the 50\% quantile probability.
#' @param VERBOSE is a bolean to output the processing steps in the terminal (Default is TRUE).
#' @param COPS is logical parameter to force the water reflectance to pass through the COPS
#' reflectance measurements made a priori. It must be turn on only if COPS data have been
#' processed and validated
#'
#'@details See User's Guide (in french) for details
#'
#'
#'@author Simon Bélanger

compute.ASD.rhow <- function(raw.asd,
                            rho.panel = 0.985,
                            quantile.prob = 0.5,
                            COPS = FALSE,
                            VERBOSE=TRUE) {

  asd.path <- getwd()
  Ltot <-raw.asd$Ltot
  Lpanel <-raw.asd$Lpanel
  Lsky <- raw.asd$Lsky
  anc <- raw.asd$anc
  waves <- Ltot$waves
  nb.waves <- length(waves)

  # Trim data to remove high Ltot spectra which may be affected by sun glint
  ix350 = which.min(abs(waves - 350))

  # use this to avoid negativa values of (ix350-5),May 22,2020, Yanqun
  if (ix350<6) ix350 = 6

  ix380 = which.min(abs(waves - 380))
  ix490 = which.min(abs(waves - 490))
  ix720 = which.min(abs(waves - 720))
  ix750 = which.min(abs(waves - 750))
  ix780 = which.min(abs(waves - 780))
  ix810 = which.min(abs(waves - 810))
  ix840 = which.min(abs(waves - 840))
  ix870 = which.min(abs(waves - 870))
  ix890 = which.min(abs(waves - 890))
  ix900 = which.min(abs(waves - 900))


  if (VERBOSE) print("Averaging Radiance spcetra...")
  ix.Lt.good = which(Ltot$Lum[ix490,] < quantile(Ltot$Lum[ix490,], probs = quantile.prob) &
                    Ltot$Lum[ix490,] > quantile(Ltot$Lum[ix490,], probs = 0.10))

  # Interquantile mean for Lsky and Lpanel
  ix.Lsky.good = which(Lsky$Lum[ix490,] < quantile(Lsky$Lum[ix490,], probs = 0.90) &
                         Lsky$Lum[ix490,] > quantile(Lsky$Lum[ix490,], probs = 0.10))

  ix.Lpanel.good = which(Lpanel$Lum[ix490,] < quantile(Lpanel$Lum[ix490,], probs = 0.90) &
                           Lpanel$Lum[ix490,] > quantile(Lpanel$Lum[ix490,], probs = 0.10))

  Lt    = Ltot$Lum[,ix.Lt.good]
  Li  = Lsky$Lum[,ix.Lsky.good]
  Lpanel= Lpanel$Lum[,ix.Lpanel.good]

  Lt.time.mean = mean.POSIXct(Ltot$DateTime[ix.Lt.good])


  ######### Average data
  Ed.mean = pi*apply(Lpanel, 1, mean)/rho.panel # See eq 5 i Mobley (1999)
  Li.mean = apply(Li, 1, mean)
  if (ncol(as.matrix(Lt))==1) {
    Lt.mean = Lt
  } else {
    Lt.mean = apply(Lt, 1, mean)
  }
  #Lt.mean = apply(Lt, 1, mean)

  Ed.sd = pi*apply(Lpanel, 1, sd)/rho.panel
  Li.sd = apply(Li, 1, sd)
  if (ncol(as.matrix(Lt)) == 1) {
    print("Only one observation, can't calculate standard deviation")
  } else {
    Lt.sd = apply(Lt, 1, sd)
  }

  ######## Compute Sky reflectance and rho and apply a smoothing function
  if (VERBOSE) print("Compute Sky and Surface reflectance and apply a smoothing function")
  x = range(waves)

  sky = Li.mean/Ed.mean
  mod= loess(sky~waves, data=data.frame(waves=Lsky$waves, sky=sky), span=20/(x[2]-x[1])) # 20 nm window
  sky.smooth = predict(mod, Lsky$waves)

  sea =Lt.mean/Ed.mean
  mod= loess(sea~waves, data=data.frame(waves=Ltot$waves, sea=sea), span=10/(x[2]-x[1])) # 10 nm window
  sea.smooth = predict(mod, Ltot$waves)

  ######### Get rho fresnel from MOBLEY LUT or use a constant rho of 0.0256 if cloudy
  if (sky.smooth[ix750] >= 0.05){
    #Then  CLOUDY SKY (Ruddick et al L&O 2006, eq. 23, 24)
    if (VERBOSE) print("Cloudy sky, Use rho = 0.0256")
    rho = 0.0256
    CLEARSKY = FALSE

  }  else {
    if (VERBOSE) print("Interpolate Mobley LUT for Fresnel reflectance")
      CLEARSKY = TRUE
      rho = get.rho550(anc$ThetaV, anc$Dphi, anc$Windspeed,anc$ThetaS)
  }


  ################### remove the reflected sky
  if (VERBOSE) print("Apply NIR corrections")
  
  #### Compute Rrs by removing sky reflectance to total reflectance.
  #### 10 methods are implemented
  nb.methods <- 10
  methods=c("Mobley+NONE", "Mobley+NULL", "Mobley+SIMILARITY1", "Mobley+SIMILARITY2",
            "NIR", "UV", "UV+NIR",  "COPS", "Kutser", "Jiang")
  Rrs <- matrix(NA, nrow=nb.methods, ncol=nb.waves)
  ##### store the QWIP parameters computed using QWIP.Rrs.R
  AVW <- rep(NA,nb.methods)
  NDI <- rep(NA,nb.methods)
  QWIP <- rep(NA,nb.methods)
  QWIP.score <- rep(NA,nb.methods)
  FU <- rep(NA,nb.methods)
  
  ##### Method 1. Mobley LUT standard method with no white correction ("NONE")
  i=1
  Rrs[i,] <- sea.smooth - (rho*sky.smooth) ### Mobley LUT standard method
  tmp=QWIP.Rrs(waves, Rrs[i,])
  AVW[i] <- tmp$AVW
  NDI[i] <- tmp$NDI
  QWIP[i] <-tmp$QWIP
  QWIP.score[i] <- tmp$QWIP.score
  FU[i] <- Rrs2FU(waves, Rrs[i,])$FU
  
  #####
  ###### Method 2.  A standard NULL correction
  i=2
  offset = mean(Rrs[1,ix890:ix900],na.rm = T)
  Rrs[i,] <- Rrs[1,] - offset # Apply NIR correction
  tmp=QWIP.Rrs(waves, Rrs[i,])
  AVW[i] <- tmp$AVW
  NDI[i] <- tmp$NDI
  QWIP[i] <-tmp$QWIP
  QWIP.score[i] <- tmp$QWIP.score
  FU[i] <- Rrs2FU(waves, Rrs[i,])$FU
  
  ###### Methods 3 & 4. Estimation of the NIR Rrs offset correction based on
  #      Ruddick et al L&O 2006, SPIE 2005
  i=3
  offset.simil1 = 2.35*Rrs[1,ix780] - Rrs[1,ix720]/(2.35-1)
  Rrs[i,] <- Rrs[1,] - offset.simil1 # Apply NIR correction
  tmp=QWIP.Rrs(waves, Rrs[i,])
  AVW[i] <- tmp$AVW
  NDI[i] <- tmp$NDI
  QWIP[i] <-tmp$QWIP
  QWIP.score[i] <- tmp$QWIP.score
  FU[i] <- Rrs2FU(waves, Rrs[i,])$FU
  Rrs.SIMILARITY1 = Rrs - offset # Apply NIR correction
  
  
  i=4
  offset.simil2 = 1.91*Rrs[1,ix870] - Rrs[1,ix780]/(1.91-1)
  Rrs[i,] <- Rrs[1,] - offset.simil2 # Apply NIR correction
  tmp=QWIP.Rrs(waves, Rrs[i,])
  AVW[i] <- tmp$AVW
  NDI[i] <- tmp$NDI
  QWIP[i] <-tmp$QWIP
  QWIP.score[i] <- tmp$QWIP.score
  FU[i] <- Rrs2FU(waves, Rrs[i,])$FU
  Rrs.SIMILARITY1 = Rrs - offset # Apply NIR correction
  
  #####
  ###### Method 5. Estimation of the rho.fresnel assuming BLACK Pixel assumption in the NIR
  i=5
  rho.sky.NIR =  (mean(sea.smooth[ix890:ix900], na.rm = T) /
                    mean(sky.smooth[ix890:ix900], na.rm = T))
  Rrs[i,] <-  sea.smooth - (rho.sky.NIR*sky.smooth)
  tmp=QWIP.Rrs(waves, Rrs[i,])
  AVW[i] <- tmp$AVW
  NDI[i] <- tmp$NDI
  QWIP[i] <-tmp$QWIP
  QWIP.score[i] <- tmp$QWIP.score
  FU[i] <- Rrs2FU(waves, Rrs[i,])$FU
  
  #####
  ###### Method 6. Estimation of the rho.fresnel assuming BLACK Pixel assumption in the UV
  i=6
  rho.sky.UV =  (mean(sea.smooth[ix350:ix380], na.rm = T) /
                   mean(sky.smooth[ix350:ix380], na.rm = T))
  Rrs[i,] <-  sea.smooth - (rho.sky.UV*sky.smooth)
  tmp=QWIP.Rrs(waves, Rrs[i,])
  AVW[i] <- tmp$AVW
  NDI[i] <- tmp$NDI
  QWIP[i] <-tmp$QWIP
  QWIP.score[i] <- tmp$QWIP.score
  FU[i] <- Rrs2FU(waves, Rrs[i,])$FU
  
  #####
  ###### Method 7. Estimation of the rho.fresnel assuming BLACK Pixel assumption in both UV and NIR (spectrally dependent)
  i=7
  rho.sky.UV.NIR = spline(c(mean(waves[ix350:ix380], na.rm = T),
                            mean(waves[ix890:ix900], na.rm = T)),
                          c(rho.sky.UV, rho.sky.NIR),
                          xout = waves)$y
  Rrs[i,] <-  sea.smooth - (rho.sky.UV.NIR*sky.smooth)
  tmp=QWIP.Rrs(waves, Rrs[i,])
  AVW[i] <- tmp$AVW
  NDI[i] <- tmp$NDI
  QWIP[i] <-tmp$QWIP
  QWIP.score[i] <- tmp$QWIP.score
  FU[i] <- Rrs2FU(waves, Rrs[i,])$FU
  
  #####
  ###### Method 8. (OPTIONAL). Only if COPS is available
  # this method forces the ASD-derived Rrs to pass through
  # the cops Rrs at two wavelenghts (second shortest and longest respectively)
  Rrs.COPS = NA
  rho.sky.COPS = NA
  if (COPS) {
    
    # Finding the COPS files avaiblable
    # Check for COPS folder in the parent directory
    ld = list.dirs("..", recursive = F)
    ix.d = grep("COPS", ld)
    if (length(ix.d) >= 1) {
      if (length(ix.d) > 1) {
        print("More than one COPS folder found")
        cops.path = ld[ix.d[1]]
      } else {
        cops.path = ld[ix.d]
      }
      
      setwd(cops.path)
      
      remove.file <- "remove.cops.dat"
      select.file <- "select.cops.dat"
      
      select.file.exists <- FALSE
      
      if (file.exists(remove.file)) {
        remove.tab <- read.table(remove.file, header = FALSE, colClasses = "character", sep = ";")
        kept.cast <- remove.tab[[2]] == "1"
      }
      if (file.exists(select.file)) {
        select.file.exists <- TRUE
        remove.tab <- read.table(select.file, header = FALSE, colClasses = "character", sep = ";")
        kept.cast <- remove.tab[[2]] == "1"
        Rrs_method <- remove.tab[kept.cast, 3]
      }
      listfile  <- remove.tab[kept.cast, 1]
      
      
      setwd("./BIN/")
      
      nf = length(listfile)
      print(listfile)
      
      if (nf > 1) {
        
        mRrs = matrix(ncol=19, nrow = nf)
        
        for (j in 1:nf) {
          
          load(paste(listfile[j], ".RData", sep=""))
          waves.COPS = cops$LuZ.waves
          
          # extract Rrs
          if (select.file.exists) {
            #mRrs[j,] = eval(parse(text=paste0("cops$",Rrs_method[j],"[xw]")))
            mRrs[j,] = eval(parse(text=paste0("cops$",Rrs_method[j])))
            
          } else {
            mRrs[j,] = cops$Rrs.0p.linear
          }
          
          
        }
        
        cops.Rrs.m = apply(mRrs, 2, mean, na.rm=T)
      } else {
        load(paste(listfile, ".RData", sep=""))
        waves.COPS = cops$LuZ.waves
        if (select.file.exists) {
          #cops.Rrs.m = eval(parse(text=paste0("cops$",Rrs_method,"[xw]")))
          cops.Rrs.m = eval(parse(text=paste0("cops$",Rrs_method)))
          
        } else {
          cops.Rrs.m = cops$Rrs.0p.linear
        }
      }
      
      # Find the shortest and longest valid wavelength for the Rrs
      ix.good.Rrs = which(cops.Rrs.m > 0)
      ix.waves.min.cops = min(ix.good.Rrs)
      ix.waves.max.cops = max(ix.good.Rrs)
      waves.max = waves.COPS[ix.waves.max.cops]
      waves.min = waves.COPS[ix.waves.min.cops]
      
      
      ix.waves.min.ASD = which.min(abs(waves - waves.min))
      if (ix.waves.min.ASD<6) ix.waves.min.ASD=6
      
      ix.waves.max.ASD = which.min(abs(waves - waves.max))
      
      # Estimate rho.shy at the two wavelenghts selected
      rho.sky.min  = ((mean(sea.smooth[(ix.waves.min.ASD-5):(ix.waves.min.ASD+5)],na.rm = T) - cops.Rrs.m[ix.waves.min.cops])
                      / mean(sky.smooth[(ix.waves.min.ASD-5):(ix.waves.min.ASD+5)], na.rm = T))
      rho.sky.max  = ((mean(sea.smooth[(ix.waves.max.ASD-5):(ix.waves.max.ASD+5)],na.rm = T) - cops.Rrs.m[ix.waves.max.cops])
                      / mean(sky.smooth[(ix.waves.max.ASD-5):(ix.waves.max.ASD+5)], na.rm = T))
      
      rho.sky.COPS = spline(c(waves.min, waves.max), c(rho.sky.min, rho.sky.max),
                            xout = waves)$y
      
      #Rrs.COPS = sea.smooth - (rho.sky.COPS*sky.smooth)
      
      i=8
      Rrs[i,] <- sea.smooth - (rho.sky.COPS*sky.smooth)
      tmp=QWIP.Rrs(waves, Rrs[i,])
      AVW[i] <- tmp$AVW
      NDI[i] <- tmp$NDI
      QWIP[i] <-tmp$QWIP
      QWIP.score[i] <- tmp$QWIP.score
      FU[i] <- Rrs2FU(waves, Rrs[i,])$FU
      
      setwd(asd.path)
      
    } else {
      print("No COPS folder found!!!")
      print("Stop processing")
      return(0)
    }
    
  }
  
  #####
  ##### Method 9: Implementation of Kutser et al. 2013 for removal of sky glint
  if (VERBOSE) print("Begining the Kutser correction for glint")
  UVdata <- sea.smooth[ix350:ix380]
  NIRdata <- sea.smooth[ix890:ix900]
  
  UVwave <- waves[ix350:ix380]
  NIRwave <- waves[ix890:ix900]
  if (VERBOSE) print("Wavelength and data at UV and NIR binned")
  
  UV.NIR.data <- c(UVdata,NIRdata)
  UV.NIR.wave <- c(UVwave,NIRwave)
  Kutserdata <- data.frame(UV.NIR.wave, UV.NIR.data)
  names(Kutserdata) <- c("waves", "urhow")
  if (VERBOSE) print("Starting the NLS")
  
  try(glint.fit <- nls(urhow ~ b*waves^z,start = list(b = 1, z = -1),data=Kutserdata,
                       control = nls.control(maxiter = 300)))
  
  
  if (exists("glint.fit")) {
    p <- coef(glint.fit)
    Kutserestimate <- p[1]*(waves)^p[2]
    i=9
    Rrs[i,]  <- sea.smooth - Kutserestimate
    tmp=QWIP.Rrs(waves, Rrs[i,])
    AVW[i] <- tmp$AVW
    NDI[i] <- tmp$NDI
    QWIP[i] <-tmp$QWIP
    QWIP.score[i] <- tmp$QWIP.score
    FU[i] <- Rrs2FU(waves, Rrs[i,])$FU
    if (VERBOSE) print("Kutser correction finished")
  } else {
    print("******** NLS fit failed for Kutser method")
    #Rrs.Kutser <- NA
  }
  
  #####
  # Method 10: from Jiang et al 2020 (https://doi.org/10.1016/j.isprsjprs.2020.05.003)
  #####
  Rrs.M99 = sea.smooth - (0.028*sky.smooth)
  md_750_780<-median(as.numeric(Rrs.M99[ix750:ix780]))
  Rrs780<-median(Rrs.M99[(ix780-5):(ix780+5)])
  Rrs810<-median(Rrs.M99[(ix810-5):(ix810+5)])
  Rrs840<-median(Rrs.M99[(ix840-5):(ix840+5)])
  RHW<-Rrs810-Rrs780-(Rrs840-Rrs780)*(810.0-780.0)/(840.0-780.0)
  #
  est_md_750_780<-18267.884*RHW^3-129.158*RHW^2+3.072*RHW
  if (RHW>0){
    delta<-md_750_780-est_md_750_780
  }else{
    delta<-md_750_780
  }
  
  i=10
  Rrs[i,] <-  Rrs.M99-delta 
  tmp=QWIP.Rrs(waves, Rrs[i,])
  AVW[i] <- tmp$AVW
  NDI[i] <- tmp$NDI
  QWIP[i] <-tmp$QWIP
  QWIP.score[i] <- tmp$QWIP.score
  FU[i] <- Rrs2FU(waves, Rrs[i,])$FU
  
  #####
  list.rho = list(
    waves = waves,
    methods = methods,
    Rrs = Rrs,
    rho.sky = rho,
    rho.sky.NIR = rho.sky.NIR,
    rho.sky.UV = rho.sky.UV,
    rho.sky.UV.NIR = rho.sky.UV.NIR,
    offset = offset,
    offset.simi1 = offset.simil1,
    offset.simi2 = offset.simil2,
    AVW = AVW,
    NDI = NDI,
    QWIP = QWIP,
    QWIP.score = QWIP.score,
    FU = FU,
    rhow = Rrs*pi,
    rho.sky = rho,
    rho.sky.NIR = rho.sky.NIR,
    rho.sky.UV = rho.sky.UV,
    rho.sky.UV.NIR = rho.sky.UV.NIR,
    rho.sky.COPS = rho.sky.COPS,
    Lpanel=Lpanel,
    Lt=Lt,
    Li=Li,
    Ed.mean = Ed.mean,
    Ed.sd = Ed.sd,
    Lt.mean = Lt.mean,
    Lt.sd = Lt.sd,
    Li.mean = Li.mean,
    Li.sd = Li.sd,
    ix.Lt.good=ix.Lt.good,
    ix.Lsky.good=ix.Lsky.good,
    ix.Lpanel.good=ix.Lpanel.good,
    DateTime = Lt.time.mean,
    anc=anc,
    CLEARSKY = CLEARSKY)

  if (VERBOSE) str(list.rho)
  return(list.rho)

}

