#'
#'@title SVC data processing launcher
#'
#'@description This function launches the SVC data processing,
#'which constists in calculating the hyperspectral marine reflectance (rho_w=pi*Rrs)
#'from radiometric measurements performed using a SVC system. It processes
#'data files found in each directories specify in the file named directories.for.SVC.dat.
#'
#'For each directory, the logging information must be stored in an
#'ASCII file named cast.info.dat. Note that the SVC files were previously
#'processed from raw to L2 using ViewSpecPro and exported to ASCII format (*.SVC.txt).
#'
#' @param  PNG is a logical parameter indicating whether or not diagnostic plots are saved in PNG format.
#'Two types of plot are produced by the function \code{\link{plot.SVC.rhow}} and saved in a sub-folder
#'/PNG/ that will be created in working directory. Default is PNG=FALSE.
#' @param COPS is logical parameter to force the water reflectance to pass through the COPS
#' reflectance measurements made a priori. It must be turn on only if COPS data have been
#' processed and validated
#'
#'
#'@details First when SVC.go is executed, it reads a file named directories.for.SVC.dat
#'in the working directory from which \code{\link{SVC.go}} was launched.
#'
#'Second, in each folder found in directories.for.SVC.dat, the programm will look for a
#'file named cast.info.SVC.dat. This file contains the logging information need to process each station
#'The cast.info.SVC.dat file contains the information on viewing geometry, windspeed,
#'the "white correction method" to eliminate sun glint, foam,
#'ocean spray etc. In details, the following fields should be found in cast.info.dat:
#'
#' * SVC.sky           is the file name for the sky reflectane (e.g., "YYMMDD_HHMM_Rxxx_Tyyy.sig")
#' * SVC.surface        is the file name for the surface reflectane
#' * ID           is acharacter string corresponding to the Station or Transect ID
#' *  Dphi         is the diffirence in azimuth between sun and sensor
#' *  Windspeed    is the wind speed
#' *  Wind.units   is the units of the wind speed ("Kts", "Km.h" or "m.s")
#' *  rhow.Method is an integer (0 to 7 or 999) indicating the best method for the specular sky reflectance removal
#' (see User Guide). 0 is for Mobley rho_sky with no NIR correction.
#' 1 is for Mobley rho_sky with NIR correction based on black pixel assumption.
#' 2 is for Mobley rho_sky with NIR correction based on similarity spectrum using 720 and 780 nm.
#' 3 is for Mobley rho_sky with NIR correction based on similarity spectrum using 780 and 870 nm.
#' 4 is for rho_sky estimated using the black pixel assumption in NIR (900 nm).
#' 5 is for rho_sky estimated using the black pixel assumption in UV (350 nm).
#' 6 is for rho_sky estimated using the black pixel assumption in both UV and NIR (spectrally dependent).
#' 7 is for rho_sky estimated using the COPS reflectance (optional processing).
#' 8 is for Kutser et al RSE 2013 method.
#' 9 if for Jiang et al. 2020 correction method. 
#' 999 is for BAD data.
#'
#'
#'Finally, the function \code{\link{process.SVC}} will be called to process each data folder.
#' @md
#'
#'@seealso See \code{\link{process.SVC}}, \code{\link{compute.SVC.rhow}} and \code{\link{plot.SVC.rhow}} for more details
#'about the processing parameters of the cast.info.dat file listed above.
#'
#'@author Simon Bélanger
SVC.go <- function(PNG=FALSE, 
                   COPS=FALSE) {

  parent.dir <- getwd()
  if(!file.exists("directories.for.SVC.dat")) {
    cat("CREATE a file named directories.for.SVC.dat in current directory (where R is launched)\n")
    cat("  and put in it the names of the directories where data files can be found (one by line)\n")
    stop()
  } else {
    dirdats <- scan(file = "directories.for.SVC.dat", "", sep = "\n", comment.char = "#")
    starting.dir <- getwd()
    for(dirdat in dirdats) {
      if(!file.exists(dirdat)) {
        cat(dirdat, "does not exist")
        stop()
      } else setwd(dirdat)
      #mymessage(paste("PROCESSING DIRECTORY", dirdat), head = "@", tail = "@")
      print(paste("PROCESSING DIRECTORY", dirdat, "@@@@@@@@@@"))
      rhow = process.SVC(dirdat, PNG, COPS)

    }
    setwd(starting.dir)
  }
}

