#' Merge SVC radiances files from the surface, the sky and the spectralon
#' panel taken to determine the diffuse above-water reflectance
#'
#' @param Ltot is a list returned by the function \code{\link{read.SVC}}
#' containing the radiances of the surface (Ltot)
#' @param Lsky is a list returned by the function \code{\link{read.SVC}}
#' containing the radiances of the sky (Lsky)
#' @param StationID is the station ID
#' @param Dphi is the sensor azimuth angle relative to the sun
#' @param Windspeed is the wind speed in m/s.
#'
#'
#' @return Returns a list containing all the radiance data of the three
#' type of measurements required for the Rrs calculation of the water
#' diffuse reflectance. In addition, ancillary data list needed to compute the Rrs
#' is also provided (lat,lon, DateTime, ThetaS, ThetaV and DeltaPhi)
#'
#' @author Simon Bélanger


merge.SVC.radiances.for.rhow <- function(Ltot,
                                              Lsky,
                                              StationID="StationX",
                                              Dphi,
                                              Windspeed) {

  # calculate the sun zenith angle from time and position.
  # Derive data from above information
  # sun-zenith angle
  
  lat <- Ltot$Target$Latitude
  lon <- Ltot$Target$Longitude
  if (is.na(Ltot$Target$DateTime)) {
    print("************No GPS available for TARGET***********************")
    print("Edit SIG file fields with format:")
    print("longitude= DDDMM.mmmmW, DDDMM.mmmmW")
    print("latitude = DDMM.mmmmN, DDMM.mmmmN")
    print("gpstime  = HHMMSS.sss, HHMMSS.sss")
    print("**************************************************************")
    stop()
  } 
  
  DateTime <- Ltot$Target$DateTime
  ThetaV <-   Ltot$Target$Tilt
  
  day <- as.numeric(format(DateTime, format = "%d"))
  month <- as.numeric(format(DateTime, format = "%m"))
  year <- as.numeric(format(DateTime, format = "%Y"))
  hour <- as.numeric(format(DateTime, format = "%H"))
  minute <- as.numeric(format(DateTime, format = "%M"))
  second <- as.numeric(format(DateTime, format = "%S"))
  ah <- hour + minute / 60 + second / 3600

  sunpos <- possol(month, day, ah, lon, lat)
  ThetaS <- sunpos[1]
  PhiS <-sunpos[2]
  anc <- list(StationID=as.character(StationID),
              lat=lat,
              lon=lon,
              DateTime=DateTime,
              ThetaV=ThetaV,
              Dphi=Dphi,
              Windspeed=Windspeed,
              ThetaS= ThetaS,
              PhiS=PhiS)
  SVC = list(Ltot=Ltot,
             Lsky=Lsky,
             anc=anc)
  return(SVC)

}
