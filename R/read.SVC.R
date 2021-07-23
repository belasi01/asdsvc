#' Read SVC files in sig format
#'
#' @param filen is the file name (*.sig)
#' 
#' 
#' @author Simon Bélanger
#'
#'@export
#'@name read.SVC
read.SVC <- function(filen){

  print(paste("Reading:", filen))

  # Count the number of header lines
  id = file(filen, "r")
  line = strsplit(readLines(con=id, n =1), " ") # Reads the first header line
  nrec = 1
  while (line != "data="){
    x = readLines(con=id, n =1)
    line = unlist(strsplit(x, " "))
    nrec <- nrec+1
    #print(line)
    ##### Extract information from header
    if (line[1] == "integration=") {
      y=unlist(strsplit(x, "="))
      z=as.numeric(unlist(strsplit(y[2], ",")))
      IntTime.Ref = z[1]
      IntTime.Target = z[4]
    }
    if (line[1] == "temp=") {
      y=unlist(strsplit(x, "="))
      z=as.numeric(unlist(strsplit(y[2], ",")))
      Temp.Ref = z[1]
      Temp.Target = z[4]
    }
    if (line[1] == "battery=") {
      y=unlist(strsplit(x, "="))
      z=as.numeric(unlist(strsplit(y[2], ",")))
      Battery.Ref = z[1]
      Battery.Target = z[2]
    }
    if (line[1] == "time=") {
      y=unlist(strsplit(x, "="))
      z=as_date(mdy_hms(str_trim(unlist(strsplit(y[2], ",")))))
      Date=z
    }
    if (line[1] == "longitude=") {
      y=unlist(strsplit(x, "="))
      z=str_trim(unlist(strsplit(y[2], ",")))
      DDD=as.numeric(str_sub(z,1,3))
      mm=as.numeric(str_sub(z,4,10))
      C = str_sub(z,11,11)
      fact=c(1,1)
      fact[C == "W"] = -1*fact[C=="W"]
      Long = (DDD + (mm/60)) * fact
      Longitude.Ref = Long[1]
      Longitude.Target = Long[2]
    }
    if (line[1] == "latitude=") {
      y=unlist(strsplit(x, "="))
      z=str_trim(unlist(strsplit(y[2], ",")))
      DD=as.numeric(str_sub(z,1,2))
      mm=as.numeric(str_sub(z,3,9))
      C = str_sub(z,10,10)
      fact=c(1,1)
      fact[C == "S"] = -1*fact[C=="S"]
      Lat = (DD + (mm/60)) * fact
      Latitude.Ref = Lat[1]
      Latitude.Target = Lat[2]
    }
    if (line[1] == "gpstime=") {
      y=unlist(strsplit(x, "="))
      z=str_trim(unlist(strsplit(y[2], ",")))
      HH=as.numeric(str_sub(z,1,2))
      MM=as.numeric(str_sub(z,3,4))
      SS=as.numeric(str_sub(z,5,8))
      GPSt = paste(HH,MM,SS, sep=":")
    }
    if (length(line) > 1) {
      if (line[2] == "method=") {
        y=unlist(strsplit(x, "="))
        z=(unlist(strsplit(y[2], ",")))
        ScanMethod = str_trim(z[1])
      }
      if (line[2] == "coadds=") {
        y=unlist(strsplit(x, "="))
        z=as.numeric(unlist(strsplit(y[2], ",")))
        Nspectra.Ref = z[1]
        Nspectra.Target = z[4]
      }
      if (line[2] == "time=") {
        y=unlist(strsplit(x, "="))
        z=as.numeric(unlist(strsplit(y[2], ",")))
        ScanTime.Ref = z[1]
        ScanTime.Target = z[2]
      }
      if (line[1] == "inclinometer"  & line[2] == "x") {
        y=unlist(strsplit(x, "="))
        z=as.numeric(unlist(strsplit(y[2], ",")))
        pitch = z
      }
      if (line[1] == "inclinometer"  & line[2] == "y") {
        y=unlist(strsplit(x, "="))
        z=as.numeric(unlist(strsplit(y[2], ",")))
        roll = z
      }
    }


  }
  
  #### Compute Tilt from Roll and Pitch for each sensor
  d2r <- pi / 180
  tilt <- atan(sqrt(tan(roll*d2r)^2+tan(pitch*d2r)^2))/d2r
  Tilt.Ref <- tilt[1]
  Tilt.Target <- tilt[2]
  
  ##### Generate UTC Date Time object
  DateTime = ymd_hms(paste(Date,GPSt))
  DateTime.Ref =  DateTime[1]
  DateTime.Target =  DateTime[2]
  reftarget.diff = difftime(DateTime.Ref, DateTime.Target)
  print(reftarget.diff)
  print("between reference and target measurements...")

  #### Close file and read the data
  print(paste("Number of header lines:", nrec))
  close(id)
  df<-fread(file = filen, skip = nrec)
  
  if (is.na(DateTime.Target)) {
    print("************No GPS available for TARGET***********************")
  }
  
  
  return(list(ScanMethod = ScanMethod,
              Reference  = list(IntTime=IntTime.Ref,
                       Temperature = Temp.Ref,
                       Battery = Battery.Ref,
                       DateTime =  DateTime.Ref,
                       Longitude = Longitude.Ref,
                       Latitude  = Latitude.Ref,
                       Tilt = Tilt.Ref,
                       Nspectra = Nspectra.Ref,
                       ScanTime = ScanTime.Ref,
                       L = df$V2), 
              Target =    list(IntTime=IntTime.Target,
                          Temperature = Temp.Target,
                          Battery = Battery.Target,
                          DateTime =  DateTime.Target,
                          Longitude = Longitude.Target,
                          Latitude  = Latitude.Target,
                          Tilt = Tilt.Target,
                          Nspectra = Nspectra.Target,
                          ScanTime = ScanTime.Target,
                          L = df$V3),
              waves = df$V1,
              reflectance = df$V4))

  


}
