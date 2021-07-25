#'  Produit des figures pour la réflectance marine
#'
#'  @param SVC est une liste produite par la fonction \code{\link{compute.SVC.rhow}}
#'  @param PNG est une variable booléenne (TRUE ou FALSE) qui permet de produire un fichier png.
#'  Par défaut PNG=FALSE
#'  @param RADIANCES est une variable booléenne (TRUE ou FALSE) qui permet de produire une figure
#'  avec les mesures de luminances et la réflectance de la surface et du ciel. Par défaut RADIANCES=FALSE
#'
#'  @author Simon Bélanger

plot.SVC.rhow <- function (SVC, PNG=FALSE, RADIANCES=FALSE) {


  ix.wl = which(SVC$waves > 350 & SVC$waves <900)

  if (PNG & !dir.exists("PNG")) dir.create("PNG")

  if (RADIANCES) {
    if (PNG) {
      png(paste("PNG/",SVC$anc$StationID,"_",SVC$anc$Replicate,"_Radiances.png",sep=""), units = "in",
          width = 5, height = 4, res = 300)

      # first plot of the raw radiances 

      Df = as.data.frame(cbind(wavelength=SVC$waves[ix.wl],
                               Li=SVC$Li.mean[ix.wl], 
                               Lt=SVC$Lt.mean[ix.wl], 
                               Ed=SVC$Ed.mean[ix.wl], 
                               rhosky=pi*SVC$rho.sky*(SVC$Li.mean[ix.wl]/SVC$Ed.mean[ix.wl]),
                               rhosurf=pi*SVC$Lt.mean[ix.wl]/SVC$Ed.mean[ix.wl]))

      p1=ggplot(data=Df, aes(x=wavelength, y=Ed))  + geom_line()  + scale_x_continuous(limits = c(350, 950))#+ geom_ribbon(aes(ymin=Ed-Ed.sd, ymax=Ed+Ed.sd, x=wavelength), alpha = 0.5)
      p2=ggplot(data=Df, aes(x=wavelength, y=Li))  + geom_line()  + scale_x_continuous(limits = c(350, 950))#+ geom_ribbon(aes(ymin=Li-Li.sd, ymax=Li+Li.sd, x=wavelength), alpha = 0.5)
      p3=ggplot(data=Df, aes(x=wavelength, y=Lt))  + geom_line()  + scale_x_continuous(limits = c(350, 950))#+ geom_ribbon(aes(ymin=Lt-Lt.sd, ymax=Lt+Lt.sd, x=wavelength), alpha = 0.5)
      p4=ggplot(data=Df, aes(x=wavelength, y=rhosky))  + geom_line()
      p4 = p4  + geom_line(aes(x=wavelength, y=rhosurf), linetype=2) + labs(x=expression(lambda),
                                                                            y=expression(paste(rho[sky],rho[surf])))


      pushViewport(viewport(layout = grid.layout(4, 1)))
      print(p1, vp = viewport(layout.pos.row = 1))
      print(p2, vp = viewport(layout.pos.row = 2))
      print(p3, vp = viewport(layout.pos.row = 3))
      print(p4, vp = viewport(layout.pos.row = 4))

      dev.off()
    } else {
      Df = as.data.frame(cbind(wavelength=SVC$waves[ix.wl],
                               Li=SVC$Li.mean[ix.wl], 
                               Lt=SVC$Lt.mean[ix.wl], 
                               Ed=SVC$Ed.mean[ix.wl], 
                               rhosky=pi*SVC$rho.sky*(SVC$Li.mean[ix.wl]/SVC$Ed.mean[ix.wl]),
                               rhosurf=pi*SVC$Lt.mean[ix.wl]/SVC$Ed.mean[ix.wl]))

      p1=ggplot(data=Df, aes(x=wavelength, y=Ed))  + geom_line()  + scale_x_continuous(limits = c(350, 950))#+ geom_ribbon(aes(ymin=Ed-Ed.sd, ymax=Ed+Ed.sd, x=wavelength), alpha = 0.5)
      p2=ggplot(data=Df, aes(x=wavelength, y=Li))  + geom_line()  + scale_x_continuous(limits = c(350, 950))#+ geom_ribbon(aes(ymin=Li-Li.sd, ymax=Li+Li.sd, x=wavelength), alpha = 0.5)
      p3=ggplot(data=Df, aes(x=wavelength, y=Lt))  + geom_line()  + scale_x_continuous(limits = c(350, 950))#+ geom_ribbon(aes(ymin=Lt-Lt.sd, ymax=Lt+Lt.sd, x=wavelength), alpha = 0.5)
      p4=ggplot(data=Df, aes(x=wavelength, y=rhosky))  + geom_line()
      p4 = p4  + geom_line(aes(x=wavelength, y=rhosurf), linetype=2) + labs(x=expression(lambda),
                                                                            y=expression(paste(rho[sky],rho[surf])))


      pushViewport(viewport(layout = grid.layout(4, 1)))
      print(p1, vp = viewport(layout.pos.row = 1))
      print(p2, vp = viewport(layout.pos.row = 2))
      print(p3, vp = viewport(layout.pos.row = 3))
      print(p4, vp = viewport(layout.pos.row = 4))
    }


  } else {

    if (!is.na(SVC$rhow.COPS[1])) {
      Df = as.data.frame(cbind(wavelength=SVC$waves[ix.wl],
                               None=SVC$rhow[ix.wl],
                               Null_900=SVC$rhow.NULL[ix.wl],
                               Similarity_720_780=SVC$rhow.SIMILARITY1[ix.wl],
                               Similarity_780_870=SVC$rhow.SIMILARITY2[ix.wl],
                               NIR = SVC$rhow.NIR[ix.wl],
                               UV  = SVC$rhow.UV[ix.wl],
                               UV.NIR = SVC$rhow.UV.NIR[ix.wl],
                               COPS= SVC$rhow.COPS[ix.wl],
                               Kutser13 = SVC$rhow.Kutser[ix.wl],
                               Jiang20 = SVC$rhow.Jiang[ix.wl]))
    } else {
      Df = as.data.frame(cbind(wavelength=SVC$waves[ix.wl],
                               None=SVC$rhow[ix.wl],
                               Null_900=SVC$rhow.NULL[ix.wl],
                               Similarity_720_780=SVC$rhow.SIMILARITY1[ix.wl],
                               Similarity_780_870=SVC$rhow.SIMILARITY2[ix.wl],
                               NIR = SVC$rhow.NIR[ix.wl],
                               UV = SVC$rhow.UV[ix.wl],
                               UV.NIR = SVC$rhow.UV.NIR[ix.wl],
                               Kutser13 = SVC$rhow.Kutser[ix.wl],
                               Jiang20 = SVC$rhow.Jiang[ix.wl]))
    }
    
    
    Dfm = melt(Df, id.vars = c("wavelength"))
    names(Dfm) = c("wavelength", "rho_w", "value" )

    if (PNG) {
      png(paste("PNG/",SVC$anc$StationID,"_",SVC$anc$Replicate,"_rhow.png",sep=""), units = "in",
          width = 5, height = 4, res = 300)

      p1 <- ggplot(data=Dfm, aes(x=wavelength, y=value, colour=rho_w)) + geom_line()
      p1 <- p1 + scale_x_continuous(limits = c(350, 950))
      p1 <- p1 + labs(x=expression(lambda), y=expression(paste(rho[w])), colour="Correction method")
      p1 <- p1 + ggtitle(paste(SVC$anc$StationID,SVC$anc$Replicate, SVC$DateTime, "Lat:", SVC$anc$lat, "Lon:", SVC$anc$lon),
                         subtitle = bquote(rho[Fresnel]^Mobley2015 == .(SVC$rho.sky) ~
                                             "   "~ rho[Fresnel]^NIR == .(SVC$rho.sky.NIR)~
                                             "   "~ rho[Fresnel]^UV == .(SVC$rho.sky.UV)))
      print(p1)
      dev.off()
    } else {
      p1 <- ggplot(data=Dfm, aes(x=wavelength, y=value, colour=rho_w)) + geom_line()
      p1 <- p1 + scale_x_continuous(limits = c(350, 950))
      p1 <- p1 + labs(x=expression(lambda), y=expression(paste(rho[w])), colour="Correction method")
      p1 <- p1 + ggtitle(paste(SVC$anc$StationID, SVC$anc$Replicate, SVC$DateTime, "Lat:", SVC$anc$lat, "Lon:", SVC$anc$lon),
                         subtitle = bquote(rho[Fresnel]^Mobley2015 == .(SVC$rho.sky) ~
                                             "   "~ rho[Fresnel]^NIR == .(SVC$rho.sky.NIR)~
                                             "   "~ rho[Fresnel]^UV == .(SVC$rho.sky.UV)))
      print(p1)
    }
  }
}
