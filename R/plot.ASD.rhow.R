#'  Produit des figures pour la réflectance marine
#'
#'  @param ASD est une liste produite par la fonction \code{\link{compute.ASD.rhow}}
#'  @param PNG est une variable booléenne (TRUE ou FALSE) qui permet de produire un fichier png.
#'  Par défaut PNG=FALSE
#'  @param RADIANCES est une variable booléenne (TRUE ou FALSE) qui permet de produire une figure
#'  avec les mesures de luminances et la réflectance de la surface et du ciel. Par défaut RADIANCES=FALSE
#'
#'  @author Simon Bélanger

plot.ASD.rhow <- function (ASD, PNG=FALSE, RADIANCES=FALSE) {


  ix.wl = which(ASD$waves > 350 & ASD$waves <900)

  if (PNG & !dir.exists("PNG")) dir.create("PNG")

  if (RADIANCES) {
    if (PNG) {
      png(paste("PNG/",ASD$anc$StationID,"_Radiances.png",sep=""), units = "in",
          width = 5, height = 4, res = 300)

      # first plot of the raw radiances with error bars

      Df = as.data.frame(cbind(wavelength=ASD$waves[ix.wl],
                               Li=ASD$Li.mean[ix.wl], Li.sd=ASD$Li.sd[ix.wl],
                               Lt=ASD$Lt.mean[ix.wl], Lt.sd=ASD$Lt.sd[ix.wl],
                               Ed=ASD$Ed.mean[ix.wl], Ed.sd=ASD$Ed.sd[ix.wl],
                               rhosky=pi*ASD$rho.sky*(ASD$Li.mean[ix.wl]/ASD$Ed.mean[ix.wl]),
                               rhosurf=pi*ASD$Lt.mean[ix.wl]/ASD$Ed.mean[ix.wl]))

      p1=ggplot(data=Df, aes(x=wavelength, y=Ed))  + geom_line()  + scale_x_continuous(limits = c(350, 950))  + geom_ribbon(aes(ymin=Ed-Ed.sd, ymax=Ed+Ed.sd, x=wavelength), alpha = 0.5)
      p2=ggplot(data=Df, aes(x=wavelength, y=Li))  + geom_line()  + scale_x_continuous(limits = c(350, 950))+ geom_ribbon(aes(ymin=Li-Li.sd, ymax=Li+Li.sd, x=wavelength), alpha = 0.5)
      p3=ggplot(data=Df, aes(x=wavelength, y=Lt))  + geom_line()  + scale_x_continuous(limits = c(350, 950))+ geom_ribbon(aes(ymin=Lt-Lt.sd, ymax=Lt+Lt.sd, x=wavelength), alpha = 0.5)
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
      Df = as.data.frame(cbind(wavelength=ASD$waves[ix.wl],
                               Li=ASD$Li.mean[ix.wl], Li.sd=ASD$Li.sd[ix.wl],
                               Lt=ASD$Lt.mean[ix.wl], Lt.sd=ASD$Lt.sd[ix.wl],
                               Ed=ASD$Ed.mean[ix.wl], Ed.sd=ASD$Ed.sd[ix.wl],
                               rhosky=pi*ASD$rho.sky*(ASD$Li.mean[ix.wl]/ASD$Ed.mean[ix.wl]),
                               rhosurf=pi*ASD$Lt.mean[ix.wl]/ASD$Ed.mean[ix.wl]))

      p1=ggplot(data=Df, aes(x=wavelength, y=Ed))  + geom_line()  + scale_x_continuous(limits = c(350, 950))  + geom_ribbon(aes(ymin=Ed-Ed.sd, ymax=Ed+Ed.sd, x=wavelength), alpha = 0.5)
      p2=ggplot(data=Df, aes(x=wavelength, y=Li))  + geom_line()  + scale_x_continuous(limits = c(350, 950))+ geom_ribbon(aes(ymin=Li-Li.sd, ymax=Li+Li.sd, x=wavelength), alpha = 0.5)
      p3=ggplot(data=Df, aes(x=wavelength, y=Lt))  + geom_line()  + scale_x_continuous(limits = c(350, 950))+ geom_ribbon(aes(ymin=Lt-Lt.sd, ymax=Lt+Lt.sd, x=wavelength), alpha = 0.5)
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
    
    
    rrs.df <- as.data.frame(t(ASD$Rrs[,ix.wl]))
    names(rrs.df) <- ASD$methods
    Df = as.data.frame(cbind(wavelength=ASD$waves[ix.wl],rrs.df))
    
    # remove COPS column if all NA
    if (all(is.na(Df$COPS))) {
      Df <- subset(Df, select = -COPS)
      methods <-   ASD$methods[-8]
    } else {
      methods <-   ASD$methods
    }
    
    Dfm = melt(Df, id.vars = c("wavelength"))
    names(Dfm) = c("wavelength", "Methods", "value" )
    
    
    # Define meaningful colors for the points and match them to the levels of the Methods variable
    method_colors <- viridis(ncol(Df)-1)
    names(method_colors) <- methods
    
    plot.rrs <- ggplot(data=Dfm, aes(x=wavelength, y=value, colour=Methods)) +
      geom_line() + scale_x_continuous(limits = c(350, 900)) +
      labs(x=expression(lambda), y=expression(paste(R[rs])), colour="Methods") +
      scale_color_manual(name = "Methods",
                         values = method_colors,
                         labels = methods) +
      ggtitle(paste( "Lat:", signif(ASD$anc$lat,5), "Lon:", signif(ASD$anc$lon,6)),
              subtitle = bquote(rho[Fresnel]^Mobley2015 == .(signif(ASD$rho.sky,3)) ~
                                  "   "~ rho[Fresnel]^NIR == .(signif(ASD$rho.sky.NIR, 3))~
                                  "   "~ rho[Fresnel]^UV == .(signif(ASD$rho.sky.UV, 3))))
    
    
    
    ##### generate the QWIP plot
    # QWIP coefficients
    p1 <- -8.399885e-9
    p2 <- 1.715532e-5
    p3 <- -1.301670e-2
    p4 <- 4.357838e0
    p5 <- -5.449532e2
    
    predicted.AVW <- 440:600
    predicted.NDI <- p1*(predicted.AVW^4) +
      p2*(predicted.AVW^3) +
      p3*(predicted.AVW^2) +
      p4*predicted.AVW   + p5
    
    # My line data frame
    df <- data.frame(AVW = predicted.AVW,
                     NDI = predicted.NDI,
                     NDI.minus.0.1 = predicted.NDI-0.1,
                     NDI.plus.0.1 = predicted.NDI+0.1,
                     NDI.minus.0.2 = predicted.NDI-0.2,
                     NDI.plus.0.2 = predicted.NDI+0.2)
    
    # Reshaping
    dfm <- melt(df,id.vars = "AVW")
    names(dfm) <- c("AVW", "Predicted", "NDI")
    
    
    # my point data frame
    df.qwip <- data.frame(AVW = ASD$AVW,
                          NDI = ASD$NDI,
                          Methods = ASD$methods,
                          FU = ASD$FU)
    
    # remove COPS line if all NA
    if (all(is.na(Df$COPS))) df.qwip <- df.qwip[-8,]
    
    
    # Plotting
    plot.QWIP <- ggplot() +
      geom_line(data = dfm, aes(x = AVW, y = NDI, color = Predicted, linetype = Predicted)) +
      geom_point(data = df.qwip, aes(x = AVW, y = NDI, fill = Methods), shape = 21, size = 3, color = "black") +
      geom_text(data = df.qwip, aes(x = AVW, y = NDI, label = FU), vjust = -0.5) + # Add labels
      scale_color_manual(name = "Lines",
                         labels = c("Predicted", "-0.1", "+0.1", "-0.2", "+0.2"),
                         values = c("black", "orange", "orange", "red", "red")) +
      scale_fill_manual(name = "Methods",
                        values = method_colors,
                        labels = methods) +
      scale_linetype_manual(name = "Lines",
                            labels = c("Predicted", "-0.1", "+0.1", "-0.2", "+0.2"),
                            values = c("solid", "dashed", "dashed", "dotted", "dotted"))
    
    fullplot <- plot.rrs / plot.QWIP +
      plot_annotation(title = paste(ASD$anc$StationID, ASD$DateTime),
                      theme = theme(plot.title = element_text(hjust = 0.5))) +
      plot_layout(guides = "collect")
    suppressMessages(plot(fullplot))
    
    
    if (PNG) ggsave(paste("PNG/Rrs_",ASD$anc$StationID,"_",ASD$anc$Replicate,".png",sep=""), units = "in",
                    width = 8, height = 7)
    

    # if (!is.na(ASD$rhow.COPS[1])) {
    #   Df = as.data.frame(cbind(wavelength=ASD$waves[ix.wl],
    #                            None=ASD$rhow[ix.wl],
    #                            Null_900=ASD$rhow.NULL[ix.wl],
    #                            Similarity_720_780=ASD$rhow.SIMILARITY1[ix.wl],
    #                            Similarity_780_870=ASD$rhow.SIMILARITY2[ix.wl],
    #                            NIR = ASD$rhow.NIR[ix.wl],
    #                            UV  = ASD$rhow.UV[ix.wl],
    #                            UV.NIR = ASD$rhow.UV.NIR[ix.wl],
    #                            COPS= ASD$rhow.COPS[ix.wl],
    #                            Kutser13 = ASD$rhow.Kutser[ix.wl],
    #                            Jiang20 = ASD$rhow.Jiang[ix.wl]))
    # } else {
    #   Df = as.data.frame(cbind(wavelength=ASD$waves[ix.wl],
    #                            None=ASD$rhow[ix.wl],
    #                            Null_900=ASD$rhow.NULL[ix.wl],
    #                            Similarity_720_780=ASD$rhow.SIMILARITY1[ix.wl],
    #                            Similarity_780_870=ASD$rhow.SIMILARITY2[ix.wl],
    #                            NIR = ASD$rhow.NIR[ix.wl],
    #                            UV = ASD$rhow.UV[ix.wl],
    #                            UV.NIR = ASD$rhow.UV.NIR[ix.wl],
    #                            Kutser13 = ASD$rhow.Kutser[ix.wl],
    #                            Jiang20 = ASD$rhow.Jiang[ix.wl]))
    # }
    # Dfm = melt(Df, id.vars = c("wavelength"))
    # names(Dfm) = c("wavelength", "rho_w", "value" )
    # 
    # if (PNG) {
    #   png(paste("PNG/",ASD$anc$StationID,"_rhow.png",sep=""), units = "in",
    #       width = 5, height = 4, res = 300)
    # 
    #   p1 <- ggplot(data=Dfm, aes(x=wavelength, y=value, colour=rho_w)) + geom_line()
    #   p1 <- p1 + scale_x_continuous(limits = c(350, 950))
    #   p1 <- p1 + labs(x=expression(lambda), y=expression(paste(rho[w])), colour="Correction method")
    #   p1 <- p1 + ggtitle(paste(ASD$anc$StationID, ASD$DateTime, "Lat:", signif(ASD$anc$lat,5), "Lon:", signif(ASD$anc$lon,5)),
    #                      subtitle = bquote(rho[Fresnel]^Mobley == .(signif(ASD$rho.sky,3)) ~
    #                                          "   "~ rho[Fresnel]^NIR == .(signif(ASD$rho.sky.NIR,3))~
    #                                          "   "~ rho[Fresnel]^UV == .(signif(ASD$rho.sky.UV,3))))
    #   print(p1)
    #   dev.off()
    # } else {
    #   p1 <- ggplot(data=Dfm, aes(x=wavelength, y=value, colour=rho_w)) + geom_line()
    #   p1 <- p1 + scale_x_continuous(limits = c(350, 950))
    #   p1 <- p1 + labs(x=expression(lambda), y=expression(paste(rho[w])), colour="Correction method")
    #   p1 <- p1 + ggtitle(paste(ASD$anc$StationID, ASD$DateTime, "Lat:", signif(ASD$anc$lat,5), "Lon:", signif(ASD$anc$lon,5)),
    #                      subtitle = bquote(rho[Fresnel]^Mobley ==  .(signif(ASD$rho.sky,3)) ~
    #                                          "   "~ rho[Fresnel]^NIR == .(signif(ASD$rho.sky.NIR,3))~
    #                                          "   "~ rho[Fresnel]^UV == .(signif(ASD$rho.sky.UV,3))))
    #   print(p1)
    #}
  }
}
