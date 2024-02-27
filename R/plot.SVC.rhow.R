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

    rrs.df <- as.data.frame(t(SVC$Rrs[,ix.wl]))
    names(rrs.df) <- SVC$methods
    Df = as.data.frame(cbind(wavelength=SVC$waves[ix.wl],rrs.df))
    
    # remove COPS column if all NA
    if (all(is.na(Df$COPS))) {
      Df <- subset(Df, select = -COPS)
      methods <-   SVC$methods[-8]
    } else methods <-   SVC$methods

    
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
        ggtitle(paste( "Lat:", signif(SVC$anc$lat,5), "Lon:", signif(SVC$anc$lon,6)),
                subtitle = bquote(rho[Fresnel]^Mobley2015 == .(signif(SVC$rho.sky,3)) ~
                                    "   "~ rho[Fresnel]^NIR == .(signif(SVC$rho.sky.NIR, 3))~
                                    "   "~ rho[Fresnel]^UV == .(signif(SVC$rho.sky.UV, 3))))
      
      
      
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
      df.qwip <- data.frame(AVW = SVC$AVW,
                            NDI = SVC$NDI,
                            Methods = SVC$methods,
                            FU = SVC$FU)
      
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
        plot_annotation(title = paste(SVC$anc$StationID, SVC$DateTime),
                        theme = theme(plot.title = element_text(hjust = 0.5))) +
        plot_layout(guides = "collect")
      suppressMessages(plot(fullplot))
      
      
      if (PNG) ggsave(paste("PNG/Rrs_",SVC$anc$StationID,"_",SVC$anc$Replicate,".png",sep=""), units = "in",
                      width = 8, height = 7)

  }
}
