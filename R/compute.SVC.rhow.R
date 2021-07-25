#'@title Compute the water leaving reflectance from SVC data
#'
#'@description
#'This function compute the water reflectance (rho_w = pi * Rrs) for a water surface.
#'The sky glint is removed using various methods as detailed in \code{\link{SVC.go}.
#'
#'
#' @param raw.SVC is a  list returned by \code{\link{merge.SVC.radiances.for.rhow}} 
#' containing all the data necessary 
#' to compute the diffuse above-water reflectance (pi*Rrs). 
#' That include the raw radiance data,
#'the sun-viewing geometry, 
#'the location, 
#'the time and the wind speed.
#' @param rho.panel is the reflectivity of the spectralon panel. The Default is 0.985
#' @param VERBOSE is a bolean to output the processing steps in the terminal (Default is TRUE).
#' @param COPS is logical parameter to force the water reflectance to pass through the COPS
#' reflectance measurements made a priori. It must be turn on only if COPS data have been
#' processed and validated
#'
#'@details See User's Guide (in french) for details
#'
#'
#'@author Simon Bélanger

compute.SVC.rhow <- function(raw.SVC,
                             rho.panel = 0.985,
                             COPS = FALSE,
                             VERBOSE=TRUE) {
  
  SVC.path <- getwd()
  waves <- raw.SVC$Ltot$waves
  
  anc <- raw.SVC$anc
  
  # Ltot <-raw.SVC$Ltot$Target$L
  # Lsky <- raw.SVC$Lsky$Target$L
  # Lpannel <-
  #anc <- raw.SVC$Ltot$Target
  
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
  
  
  # Lt    = Ltot$Lum[,ix.Lt.good]
  # Li  = Lsky$Lum[,ix.Lsky.good]
  # Lpanel= Lpanel$Lum[,ix.Lpanel.good]
  # 
  # Lt.time.mean = mean.POSIXct(Ltot$DateTime[ix.Lt.good])
  
  
  ######### Average data
  # Ed.mean = pi*apply(Lpanel, 1, mean)/rho.panel # See eq 5 i Mobley (1999)
  # Li.mean = apply(Li, 1, mean)
  # if (ncol(as.matrix(Lt))==1) {
  #   Lt.mean = Lt
  # } else {
  #   Lt.mean = apply(Lt, 1, mean)
  # }
  # #Lt.mean = apply(Lt, 1, mean)
  # 
  # Ed.sd = pi*apply(Lpanel, 1, sd)/rho.panel
  # Li.sd = apply(Li, 1, sd)
  # if (ncol(as.matrix(Lt)) == 1) {
  #   print("Only one observation, can't calculate standard deviation")
  # } else {
  #   Lt.sd = apply(Lt, 1, sd)
  # }
  Ed = pi*raw.SVC$Ltot$Reference$L/rho.panel 
  
  ######## Compute Sky reflectance and rho and apply a smoothing function
  if (VERBOSE) print("Compute Sky and Surface reflectance and apply a smoothing function")
  x = range(waves)
  
  #sky = raw.SVC$Lsky$reflectance
  sky = raw.SVC$Lsky$Target$L / Ed
  mod= loess(sky~waves, data=data.frame(waves=waves, sky=sky), span=20/(x[2]-x[1])) # 20 nm window
  sky.smooth = predict(mod, waves)
  
  #sea = raw.SVC$Ltot$reflectance
  sea = raw.SVC$Ltot$Target$L / Ed
  mod= loess(sea~waves, data=data.frame(waves=waves, sea=sea), span=10/(x[2]-x[1])) # 10 nm window
  sea.smooth = predict(mod, waves)
  
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
  Rrs = sea.smooth - (rho*sky.smooth)
  
  # Apply a correction
  
  #####
  ###### Method 1.  A standard NULL correction
  
  offset = Rrs[ix900]
  Rrs.NULL = Rrs - offset # Apply NIR correction
  
  #####
  ###### Methods 2 and 3. Estimation of the NIR Rrs offset correction based on
  #      Ruddick et al L&O 2006, SPIE 2005
  offset = 2.35*Rrs[ix780] - Rrs[ix720]/(2.35-1)
  Rrs.SIMILARITY1 = Rrs - offset # Apply NIR correction
  offset = 1.91*Rrs[ix870] - Rrs[ix780]/(1.91-1)
  Rrs.SIMILARITY2 = Rrs - offset # Apply NIR correction
  
  #####
  ###### Method 4. Estimation of the rho.fresnel assuming BLACK Pixel assumption in the NIR
  rho.sky.NIR =  (mean(sea.smooth[ix890:ix900], na.rm = T) /
                    mean(sky.smooth[ix890:ix900], na.rm = T))
  Rrs.BP = sea.smooth - (rho.sky.NIR*sky.smooth)
  
  #####
  ###### Method 5. Estimation of the rho.fresnel assuming BLACK Pixel assumption in the UV
  rho.sky.UV =  (mean(sea.smooth[ix350:ix380], na.rm = T) /
                   mean(sky.smooth[ix350:ix380], na.rm = T))
  Rrs.UV = sea.smooth - (rho.sky.UV*sky.smooth)
  
  #####
  ###### Method 6. Estimation of the rho.fresnel assuming BLACK Pixel assumption in both UV and NIR (spectrally dependent)
  rho.sky.UV.NIR = spline(c(mean(waves[ix350:ix380], na.rm = T),
                            mean(waves[ix890:ix900], na.rm = T)),
                          c(rho.sky.UV, rho.sky.NIR),
                          xout = waves)$y
  Rrs.UV.NIR = sea.smooth - (rho.sky.UV.NIR*sky.smooth)
  
  #####
  ###### Method 7. (OPTIONAL). Only if COPS is available
  # this method forces the SVC-derived Rrs to pass through
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
      
      
      ix.waves.min.SVC = which.min(abs(waves - waves.min))
      if (ix.waves.min.SVC<6) ix.waves.min.SVC=6
      
      ix.waves.max.SVC = which.min(abs(waves - waves.max))
      
      # Estimate rho.shy at the two wavelenghts selected
      rho.sky.min  = ((mean(sea.smooth[(ix.waves.min.SVC-5):(ix.waves.min.SVC+5)],na.rm = T) - cops.Rrs.m[ix.waves.min.cops])
                      / mean(sky.smooth[(ix.waves.min.SVC-5):(ix.waves.min.SVC+5)], na.rm = T))
      rho.sky.max  = ((mean(sea.smooth[(ix.waves.max.SVC-5):(ix.waves.max.SVC+5)],na.rm = T) - cops.Rrs.m[ix.waves.max.cops])
                      / mean(sky.smooth[(ix.waves.max.SVC-5):(ix.waves.max.SVC+5)], na.rm = T))
      
      rho.sky.COPS = spline(c(waves.min, waves.max), c(rho.sky.min, rho.sky.max),
                            xout = waves)$y
      
      Rrs.COPS = sea.smooth - (rho.sky.COPS*sky.smooth)
      
      setwd(SVC.path)
      
    } else {
      print("No COPS folder found!!!")
      print("Stop processing")
      return(0)
    }
    
  }
  
  #####
  ##### Method 8: Implementation of Kutser et al. 2013 for removal of sky glint
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
  
  try(glint.fit <- nls(urhow ~ b*waves^z,start = list(b = 1, z = 0),data=Kutserdata))
  
  if (exists("glint.fit")) {
    p <- coef(glint.fit)
    Kutserestimate <- p[1]*(waves)^p[2]
    Rrs.Kutser <- sea.smooth - Kutserestimate
    if (VERBOSE) print("Kutser correction finished")
    # TO BE REMOVED EVENTUALLY
    #plot(Kutserdata$waves, Kutserdata$urhow)
    #summary(glint.fit)
    #plot(sea.smooth)
    #lines(Kutserestimate)
  } else {
    print("******** NLS fit failed for Kutser method")
    Rrs.Kutser <- NA
  }
  
  
  
  #####
  # Method 9: from Jiang et al 2020 (https://doi.org/10.1016/j.isprsjprs.2020.05.003)
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
  Rrs.Jiang <- Rrs.M99-delta
  
  #####
  list.rho = list(
    waves = waves,
    rhow = Rrs*pi,
    rhow.NULL = Rrs.NULL*pi,
    rhow.SIMILARITY1 = Rrs.SIMILARITY1*pi,
    rhow.SIMILARITY2 = Rrs.SIMILARITY2*pi,
    rhow.NIR = Rrs.BP * pi,
    rhow.UV = Rrs.UV * pi,
    rhow.UV.NIR = Rrs.UV.NIR * pi,
    rhow.COPS = Rrs.COPS * pi,
    rhow.Kutser = Rrs.Kutser * pi,
    rhow.Jiang = Rrs.Jiang * pi,
    rho.sky = rho,
    rho.sky.NIR = rho.sky.NIR,
    rho.sky.UV = rho.sky.UV,
    rho.sky.UV.NIR = rho.sky.UV.NIR,
    rho.sky.COPS = rho.sky.COPS,
    #Lpanel=raw.SVC$Ltot$Reference$L,
    #Lt=raw.SVC$Ltot$Target$L,
    #Li=raw.SVC$Lsky$Target$L,
    Ed.mean = Ed,
    #Ed.sd = Ed.sd,
    Lt.mean = raw.SVC$Ltot$Target$L,
    #Lt.sd = Lt.sd,
    Li.mean = raw.SVC$Lsky$Target$L,
    #Li.sd = Li.sd,
    #ix.Lt.good=ix.Lt.good,
    #ix.Lsky.good=ix.Lsky.good,
    #ix.Lpanel.good=ix.Lpanel.good,
    DateTime = raw.SVC$Ltot$Target$DateTime,
    anc=anc,
    CLEARSKY = CLEARSKY)
  if (VERBOSE) str(list.rho)
  return(list.rho)
  
}

