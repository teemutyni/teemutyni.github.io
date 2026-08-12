library(animation)
ani.options(interval=.1) #.2


k0<-5
pts<-seq(-5,5,by=.01)
usc<-u(pts,k0,M1)-u_j(pts,k0,0)

saveGIF(for(time in seq(0,2*pi/(k0^2),by=(.025*2*pi/(k0^2)))){
  #pdf(paste0(c("wave_k",k_0,".pdf"), collapse = ""),paper="USr")
  plot(pts,Re(usc*exp(-1i*k0^2*time)),type="l",col="blue",ylim=c(-1.5e-3,1.5e-3),ylab=expression('u'[sc]),cex=1.2)
  lines(pts,Im(usc*exp(-1i*k0^2*time)),type="l",col="red")
  polygon(c(0,0,.5,.5),c(-10,10,10,-10),col = adjustcolor("gray", alpha=0.3), border = NA)
  
  par(new=T)
  plot(function(x) V(x),-5,5,
       type="l",
       #xlab="x",
       col="black",
       ylim=c(-1,1),
       n=600,
       axes=F, xlab=NA, ylab=NA
  )
  axis(side = 4)
  mtext(side = 4, line = 3, 'Coeffs.')
  title(main=paste("Scattered wave, k=", k0))
  legend( x="topleft", 
          legend=c(expression('Re(u'[sc]*')'), expression('Im(u'[sc]*')'),"coeffs."),
          col=c("blue","red","black"),
          lty=c(1,1) 
  )
  #dev.off()
})