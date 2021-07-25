#'
#'@title Process a SVC data folder
#'
#'@description
#'This function is called by the higher level function \code{\link{SVC.go}}.
#'It can be also called in the command line to process a given data directory.
#'
#'
#' @param dirdat is the current directory path contaning the files to process.
#' @param  PNG is a logical parameter indicating whether or not diagnostic plots are saved in PNG format.
#'Two types of plot are produced by the function \code{\link{plot.SVC.rhow}} and saved in a sub-folder
#'/PNG/ that will be created in working directory. Default is PNG=FALSE.
#' @param COPS is logical parameter to force the water reflectance to pass through the COPS
#' reflectance measurements made a priori. It must be turn on only if COPS data have been
#' processed and validated
#'
#'@return The function returns a list containing the rhow data.
#'The same list (rhow) is saved into a rhow.RData file in dirdat.
#'
#'@details The function first looks into the cast.info.dat file found
#'in the dirdat and check whether all the mandatory fields are present.
#'If not, the processing may be stop or a default value is taken.
#'Next the SVC files are merged using \code{\link{merge.SVC.radiances.for.rhow}}.
#'Finally, the function calls \code{\link{compute.SVC.rhow}} to
#'compute rhow spectra using various methods (see \code{\link{SVC.go} for available methods).
#'
#'@seealso
#'\code{\link{compute.SVC.rhow}} and \code{\link{SVC.go}}
#'
#'@author Simon Bélanger
process.SVC<- function(dirdat, PNG=FALSE,
                       COPS=FALSE) {

  # Cast info file
  default.cast.info.file <- paste( Sys.getenv("R_SVC_DATA_DIR"), "cast.info.SVC.dat", sep = "/")

  cast.info.file <- paste(dirdat, "cast.info.SVC.dat", sep = "/")
  if(!file.exists(cast.info.file)) {
    file.copy(from = default.cast.info.file, to = cast.info.file)
    cat("EDIT file", cast.info.file, "and CUSTOMIZE IT\n")
    cat("This file contains the information on:\n")
    cat("   SVC.sky, SVC.surface,	ID,	\n")
    cat(" 	Dphi,	Windspeed, Wind.units,	rho.panel, rhow.Method\n")
    cat("   Read the User's Guide for more details.\n")
    stop("Abort processing...")
  }

  cast.info <- read.table(cast.info.file, header=T, comment.char = "#")

  if (is.null(cast.info$ID)) {
    print("ID (the station ID) not found in cast.info.SVC.dat")
    print("Abort processing.")
    stop()
  }

  ################ Loop on each lines found in cast.info.dat
  experiments <- nrow(cast.info)
  print(paste("There are ", experiments, " cast(s) in the folder ", dirdat))
  for(experiment in 1:experiments) {
    if (is.null(cast.info$Dphi[experiment])) {
      print("Dphi (delta azimuth betwen sun and sensor viewing aximuth) not found in cast.info.dat")
      print("Assuming 90 degrees")
      cast.info$Dphi[experiment] = 90
    }
    if (is.null(cast.info$Windspeed[experiment])) {
      print("Windspeed not found in cast.info.dat")
      print("Assuming wind = 4 m/s")
      cast.info$Windspeed[experiment] = 4
    } else {
      if (is.null(cast.info$Wind.units[experiment])) {
        print("Wind.units (i.e., Kts or m.s or km.h) not found in cast.info.dat")
        print("Assuming wind units is knots!!!!")
        cast.info$Wind.units[experiment] = "Kts"
      } else
      {
        #        for (i in 1:nreplicates) {
        if (cast.info$Wind.units[experiment] == "Kts") {
          print("Convert wind speed from Knots to m/s" )
          cast.info$Windspeed[experiment] = cast.info$Windspeed[experiment]/1.9426
        }
        if (cast.info$Wind.units[experiment] == "Km.h" | cast.info$Wind.units[experiment] == "Km/h") {
          print("Convert wind speed from Km/h to m/s" )
          cast.info$Windspeed[experiment] = cast.info$Windspeed[experiment]/3.6
        }
        #        }
      }
    }

    if (is.null(cast.info$rho.panel[experiment])) {
      print("rho.panel not found in cast.info.dat")
      print("Assumes 0.985")
      cast.info$rho.panel[experiment] = 0.985#rep(0.5,nreplicates)
    }

    if (is.null(cast.info$rhow.Method)) {
      print("rhow.Method not found in cast.info.dat")
      print("Assumes no NIR correction (=0)")
      cast.info$rhow.Method[experiment] = 0#rep(0.5,nreplicates)
    }

    #####
    ##### Launch the processing
    # Lecture des données SVC du ciel et de la surface
    Lsky = read.SVC(cast.info$SVC.sky[experiment])
    Ltot = read.SVC(cast.info$SVC.surface[experiment])
      # Les données SVC et les informations auxiliaires sont passées à la fonction suivante
    SVCtot = merge.SVC.radiances.for.rhow(Ltot,
                                          Lsky,
                                          StationID=cast.info$ID[experiment],
                                          Replicate=cast.info$Replicate[experiment],
                                          Dphi=cast.info$Dphi[experiment],
                                          Windspeed=cast.info$Windspeed[experiment])


    ####
    rhow = compute.SVC.rhow(SVCtot, 
                            rho.panel = cast.info$rho.panel[experiment],
                            COPS)

    plot.SVC.rhow(rhow, PNG, RADIANCES = TRUE)
    plot.SVC.rhow(rhow, PNG)

    if (!dir.exists("RData")) dir.create("RData")
    save(file = paste(dirdat,"/RData/",cast.info$ID[experiment],"_",cast.info$Replicate[experiment],".SVC.rhow.RData", sep=""), rhow)
  }
}

