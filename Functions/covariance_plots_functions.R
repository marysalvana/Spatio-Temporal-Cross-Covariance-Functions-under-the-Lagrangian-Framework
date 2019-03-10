covariance_plots_functions <- function(filename,simcov,simcov2,simcov3){
  
  pdf(file=filename, width=20, height=13)
  split.screen( rbind(c(0.05,0.95,0.652,0.95),c(0.05,0.95,0.35,0.70),c(0.05,0.95,0.048,0.45), c(.95,0.99,0,0.98)))
  split.screen( figs = c( 1, 5 ), screen = 1 )
  split.screen( figs = c( 1, 5 ), screen = 2 )
  split.screen( figs = c( 1, 5 ), screen = 3 )
  for(bb in 1:5){
    screen(bb+4)
    if(bb==1){
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(simcov[[bb]],col=tim.colors(),zlim=c(0,1),axes=FALSE,main = "")
      axis(2, at=seq(0,1,length.out = 5), labels = seq(0,1,length.out=5))
      mtext(expression(s[y]), side = 2, line = 1.75, adj = 0.5, cex = 1, font=2)
      mtext(expression(C[11]), side = 2, line = 3, adj = 0.5, cex = 1.5, font=2,col="#0086FF")
      mtext(paste("u=",bb-1,sep=""), side = 3, line = 1, adj = 0.5, cex = 1.5, font=2,col="#4EC1DE")
    }else{
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(simcov[[bb]],col=tim.colors(),zlim=c(0,1),axes=FALSE,main = "")
      mtext(paste("u=",bb-1,sep=""), side = 3, line = 1, adj = 0.5, cex = 1.5, font=2,col="#4EC1DE")
    }
  }
  for(bb in 1:5){
    screen(bb+9)
    if(bb==1){
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(simcov2[[bb]],col=tim.colors(),zlim=c(0,1),axes=FALSE)
      axis(2, at=seq(0,1,length.out = 5), labels = seq(0,1,length.out=5))
      mtext(expression(s[y]), side = 2, line = 1.75, adj = 0.5, cex = 1, font=2)
      mtext(expression(C[22]), side = 2, line = 3, adj = 0.5, cex = 1.5, font=2,col="#0086FF")
    }else{
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(simcov2[[bb]],col=tim.colors(),zlim=c(0,1),axes=FALSE)
    }
  }
  for(bb in 1:5){
    screen(bb+14)
    if(bb==1){
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(simcov3[[bb]],col=tim.colors(),zlim=c(0,1),axes=FALSE)
      axis(1, at=seq(0,1,length.out = 5), labels = seq(0,1,length.out=5))
      axis(2, at=seq(0,1,length.out = 5), labels = seq(0,1,length.out=5))
      mtext(expression(s[x]), side = 1, line = 2, adj = 0.5, cex = 1, font=2)
      mtext(expression(s[y]), side = 2, line = 1.75, adj = 0.5, cex = 1, font=2)
      mtext(expression(C[12]), side = 2, line = 3, adj = 0.5, cex = 1.5, font=3,col="#0086FF")
    }else{
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(simcov3[[bb]],col=tim.colors(),zlim=c(0,1),axes=FALSE)
      axis(1, at=seq(0,1,length.out = 5), labels = seq(-0.5,0.5,length.out=5))
      mtext(expression(s[x]), side = 1, line = 2, adj = 0.5, cex = 1, font=2)
    }
  }
  screen(4)
  x1 <- c(0.3,0.4,0.4,0.3)
  y1 <- c(0.3,0.3,0.75,0.75)
  legend.gradient2(cbind(x1,y1),cols = colorRampPalette(colors) (colsteps), title = "", 
                   limits = seq(0,1,by=0.25), cex=1)
  
  close.screen( all=TRUE)
  dev.off()
}

simulation_plots_functions <- function(filename,var1sim,var2sim){
  
  pdf(file=filename, width=20, height=9)
  
  zr1=range(c(var1_sim,var2_sim))
  split.screen( rbind(c(0.05,0.95,0.4,1),c(0.05,0.95,0,0.6), c(.95,0.99,0,0.98)))
  split.screen( figs = c( 1, 5 ), screen = 1 )
  split.screen( figs = c( 1, 5 ), screen = 2 )
  for(bb in 1:5){
    screen(bb+3)
    if(bb==1){
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(matrix(var1sim[,bb],nrow=11,ncol=11,byrow=T),col=tim.colors(),zlim=zr1,main = "",axes=FALSE,col.main="#4EC1DE")
      mtext(expression(s[y]), side = 2, line = 1.75, adj = 0.5, cex = 1, font=2)
      axis(2, at=seq(0,1,by=0.2), labels = seq(0,1,by=0.2))
      mtext(expression(Z[1]), side = 2, line = 3, adj = 0.5, cex = 1.5, font=3,col="#0086FF")
      mtext(paste("t=",bb,sep=""), side = 3, line = 1, adj = 0.5, cex = 1.5, font=2,col="#4EC1DE")
    }else{
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(matrix(var1sim[,bb],nrow=11,ncol=11,byrow=T),col=tim.colors(),zlim=zr1,main = "",axes=FALSE)
      mtext(paste("t=",bb,sep=""), side = 3, line = 1, adj = 0.5, cex = 1.5, font=2,col="#4EC1DE")
    }
    
  }
  for(bb in 1:5){
    screen(bb+8)
    if(bb==1){
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(matrix(var2sim[,bb],nrow=11,ncol=11,byrow=T),col=tim.colors(),zlim=zr1,
            ylab=expression(bold(Z[2])),font.main = 2,font.lab = 2,axes=FALSE)
      mtext(expression(s[y]), side = 2, line = 1.75, adj = 0.5, cex = 1, font=2)
      mtext(expression(s[x]), side = 1, line = 2, adj = 0.5, cex = 1, font=2)
      axis(1, at=seq(0,1,by=0.2), labels = seq(0,1,by=0.2))
      axis(2, at=seq(0,1,by=0.2), labels = seq(0,1,by=0.2))
      mtext(expression(Z[2]), side = 2, line = 3, adj = 0.5, cex = 1.5, font=3,col="#0086FF")
    }else{
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(matrix(var2sim[,bb],nrow=11,ncol=11,byrow=T),col=tim.colors(),zlim=zr1,
            font.main = 1,axes=FALSE)
      mtext(expression(s[x]), side = 1, line = 2, adj = 0.5, cex = 1, font=2)
      axis(1, at=seq(0,1,by=0.2), labels = seq(0,1,by=0.2))
    }
  }
  screen(3)
  x1 <- c(0.5,0.6,0.6,0.5)
  y1 <- c(0.3,0.3,0.75,0.75)
  legend.gradient2(cbind(x1,y1),cols = colorRampPalette(colors) (colsteps), title = "", 
                   limits = seq(round(zr1[1],0),round(zr1[2],0),length.out = 5), cex=1)
  #title("Example Realizations from the model in Example 2", line = -1, outer = TRUE)
  close.screen( all=TRUE)
  dev.off()
}

covariance_plots_functions_section5 <- function(filename, simcov, rotated){
  pdf(paste(filename,sep=''), width=8.7, height=8)
  
  split.screen( rbind(c(0.1,0.93,0.1,1), c(.93,0.99,0.1,1)))
  split.screen( figs = c( 3, 1 ), screen = 1 )
  split.screen( figs = c( 1, dt ), screen = 3 )
  split.screen( figs = c( 1, dt ), screen = 4 )
  split.screen( figs = c( 1, dt ), screen = 5 )
  
  #par(mfrow=c(3,3), tcl=-0.5,omi=c(0.2,0.2,0,0))
  # Top left panel
  for(bb in 1:dt){
    screen(bb+5)
    if(bb==1){
      par(pty="s") 
      par(mai=c(0.1,0.15,0.2,0.1))
      
      if(rotated){
        plot((simcov[[1]][,1:2]%*%R)/1000000, col=colorRampPalette(colors) (colsteps) [ findInterval(simcov[[bb]][,4], seq(zr1[1],zr1[2], length.out=colsteps)) ]
             , pch = 20, xaxt="n", cex.axis=0.8,
             main=paste('u =',bb-1,sep=' '),xlab='', ylab='',col.main= "#4EC1DE")
      }else{
        plot((simcov[[1]][,1:2])/1000000, col=colorRampPalette(colors) (colsteps) [ findInterval(simcov[[bb]][,4], seq(zr1[1],zr1[2], length.out=colsteps)) ]
             , pch = 20, xaxt="n", cex.axis=0.8,
             main=paste('u =',bb-1,sep=' '),xlab='', ylab='',col.main= "#4EC1DE")
      }
      mtext(expression(hat(C)[11]), side = 2, line = 3, adj = 0.5, font=3,col="#0086FF",cex=1.3)
      
    }else{
      par(pty="s") 
      par(mai=c(0.1,0.15,0.2,0.1))
      if(rotated){
        plot((simcov[[1]][,1:2]%*%R)/1000000, col=colorRampPalette(colors) (colsteps) [ findInterval(simcov[[bb]][,4], seq(zr1[1],zr1[2], length.out=colsteps)) ]
             , pch = 20, yaxt="n", xaxt="n", cex.axis=0.8,
             main=paste('u =',bb-1,sep=' '),xlab='', ylab='',col.main= "#4EC1DE")
      }else{
        plot((simcov[[1]][,1:2])/1000000, col=colorRampPalette(colors) (colsteps) [ findInterval(simcov[[bb]][,4], seq(zr1[1],zr1[2], length.out=colsteps)) ]
             , pch = 20, yaxt="n", xaxt="n", cex.axis=0.8,
             main=paste('u =',bb-1,sep=' '),xlab='', ylab='',col.main= "#4EC1DE")
      }
    }
    abline(h=0,v=0,lwd=1,lty=2,col=3)
  }
  
  for(bb in 1:dt){
    screen(bb+5+dt)
    if(bb==1){
      par(pty="s") 
      par(mai=c(0.1,0.15,0.2,0.1))
      if(rotated){
        plot((simcov[[1]][,1:2]%*%R)/1000000, col=colorRampPalette(colors) (colsteps) [ findInterval(simcov[[bb]][,5], seq(zr1[1],zr1[2], length.out=colsteps)) ]
             ,cex.axis=0.8, pch = 20,xlab='',ylab='',col.main= "blue", xaxt="n")
      }else{
        plot((simcov[[1]][,1:2])/1000000, col=colorRampPalette(colors) (colsteps) [ findInterval(simcov[[bb]][,5], seq(zr1[1],zr1[2], length.out=colsteps)) ]
             ,cex.axis=0.8, pch = 20,xlab='',ylab='',col.main= "blue", xaxt="n")
      }
      mtext(expression(hat(C)[22]), side = 2, line = 3, adj = 0.5, font=3,col="#0086FF",cex=1.3)
      
    }else{
      par(pty="s") 
      par(mai=c(0.1,0.15,0.2,0.1))
      if(rotated){
        plot((simcov[[1]][,1:2]%*%R)/1000000, col=colorRampPalette(colors) (colsteps) [ findInterval(simcov[[bb]][,5], seq(zr1[1],zr1[2], length.out=colsteps)) ]
             ,cex.axis=0.8, pch = 20, xlab='',ylab='',col.main= "blue", yaxt="n",xaxt="n")
      }else{
        plot((simcov[[1]][,1:2])/1000000, col=colorRampPalette(colors) (colsteps) [ findInterval(simcov[[bb]][,5], seq(zr1[1],zr1[2], length.out=colsteps)) ]
             ,cex.axis=0.8, pch = 20, xlab='',ylab='',col.main= "blue", yaxt="n",xaxt="n")
      }
    }
    abline(h=0,v=0,lwd=1,lty=2,col=3)
    
  }
  
  for(bb in 1:dt){
    screen(bb+5+dt*2)
    if(bb==1){
      par(pty="s") 
      par(mai=c(0.1,0.15,0.2,0.1))
      if(rotated){
        plot((simcov[[1]][,1:2]%*%R)/1000000, col=colorRampPalette(colors) (colsteps) [ findInterval(simcov[[bb]][,6], seq(zr1[1],zr1[2], length.out=colsteps)) ]
             ,cex.axis=0.8, pch = 20, xlab='',ylab='',col.main= "blue")
      }else{
        plot((simcov[[1]][,1:2])/1000000, col=colorRampPalette(colors) (colsteps) [ findInterval(simcov[[bb]][,6], seq(zr1[1],zr1[2], length.out=colsteps)) ]
             ,cex.axis=0.8, pch = 20, xlab='',ylab='',col.main= "blue")
      }
      mtext(expression(hat(C)[12]), side = 2, line = 3, adj = 0.5, font=3,col="#0086FF",cex=1.3)
    }else{
      par(pty="s") 
      par(mai=c(0.1,0.15,0.2,0.1))
      if(rotated){
        plot((simcov[[1]][,1:2]%*%R)/1000000, col=colorRampPalette(colors) (colsteps) [ findInterval(simcov[[bb]][,6], seq(zr1[1],zr1[2], length.out=colsteps)) ]
             , cex.axis=0.8,pch = 20, xlab='',ylab='',col.main= "blue", yaxt="n")
      }else{
        plot((simcov[[1]][,1:2])/1000000, col=colorRampPalette(colors) (colsteps) [ findInterval(simcov[[bb]][,6], seq(zr1[1],zr1[2], length.out=colsteps)) ]
             , cex.axis=0.8,pch = 20, xlab='',ylab='',col.main= "blue", yaxt="n")
      }
    }
    abline(h=0,v=0,lwd=1,lty=2,col=3)
    
  }
  
  screen(2)
  x1 <- c(0,0.15,0.15,0)
  y1 <- c(0.25,0.25,0.7,0.7)
  legend.gradient2(cbind(x1,y1),cols = colorRampPalette(colors) (colsteps), title = "", 
                   limits = seq(0,1,length.out = 5), cex=1)
  mtext("E-W distance (x 1000 km)", side=1, line = 7, adj = 3, font=2)
  mtext("N-S distance (x 1000 km)", side = 2, line = 41.3, adj = 0.5, font=2)
  
  close.screen( all=TRUE)
  dev.off()
}

sigma_plots_functions <- function(filename,simcov,simcov2,simcov3,var,
                                  simcov.2,simcov2.2,simcov3.2){
  
  lim1 <- range(exp(simcov[[1]]),exp(simcov[[2]]),exp(simcov[[3]]),exp(simcov[[4]]),exp(simcov[[5]]),
                exp(simcov.2[[1]]),exp(simcov.2[[2]]),exp(simcov.2[[3]]),exp(simcov.2[[4]]),exp(simcov.2[[5]]))
  lim2 <- range(exp(simcov2[[1]]),exp(simcov2[[2]]),exp(simcov2[[3]]),exp(simcov2[[4]]),exp(simcov2[[5]]),
                exp(simcov2.2[[1]]),exp(simcov2.2[[2]]),exp(simcov2.2[[3]]),exp(simcov2.2[[4]]),exp(simcov2.2[[5]]))
  lim3 <- range(pi/2 * exp(simcov3[[1]])/(1 + exp(simcov3[[1]])),
                pi/2 * exp(simcov3[[2]])/(1 + exp(simcov3[[2]])),
                pi/2 * exp(simcov3[[3]])/(1 + exp(simcov3[[3]])),
                pi/2 * exp(simcov3[[4]])/(1 + exp(simcov3[[4]])),
                pi/2 * exp(simcov3[[5]])/(1 + exp(simcov3[[5]])),
                pi/2 * exp(simcov3.2[[1]])/(1 + exp(simcov3.2[[1]])),
                pi/2 * exp(simcov3.2[[2]])/(1 + exp(simcov3.2[[2]])),
                pi/2 * exp(simcov3.2[[3]])/(1 + exp(simcov3.2[[3]])),
                pi/2 * exp(simcov3.2[[4]])/(1 + exp(simcov3.2[[4]])),
                pi/2 * exp(simcov3.2[[5]])/(1 + exp(simcov3.2[[5]])))
  
  pdf(file=filename, width=20, height=13)
  split.screen( rbind(c(0.05,0.95,0.652,0.95),c(0.05,0.95,0.35,0.70),c(0.05,0.95,0.048,0.45), 
                      c(.95,0.99,0.652,0.95),c(.95,0.99,0.35,0.70),c(.95,0.99,0.048,0.45)))
  split.screen( figs = c( 1, 5 ), screen = 1 )
  split.screen( figs = c( 1, 5 ), screen = 2 )
  split.screen( figs = c( 1, 5 ), screen = 3 )
  for(bb in 1:5){
    screen(bb+6)
    if(bb==1){
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(exp(simcov[[bb]]),col=tim.colors(),zlim=lim1,axes=FALSE,main = "")
      axis(2, at=seq(0,1,length.out = 5), labels = seq(0,1,length.out=5))
      mtext(expression(s[y]), side = 2, line = 1.75, adj = 0.5, cex = 1, font=2)
      if(var==1){
        mtext(expression(lambda[1.1]), side = 2, line = 3, adj = 0.5, cex = 2, font=2,col="#0086FF")
      }else{
        mtext(expression(lambda[2.1]), side = 2, line = 3, adj = 0.5, cex = 2, font=2,col="#0086FF")
      }
      mtext(paste("t=",bb,sep=""), side = 3, line = 1, adj = 0.5, cex = 1.5, font=2,col="#4EC1DE")
    }else{
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(exp(simcov[[bb]]),col=tim.colors(),zlim=lim1,axes=FALSE,main = "")
      mtext(paste("t=",bb,sep=""), side = 3, line = 1, adj = 0.5, cex = 1.5, font=2,col="#4EC1DE")
    }
  }
  for(bb in 1:5){
    screen(bb+11)
    if(bb==1){
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(exp(simcov2[[bb]]),col=tim.colors(),zlim=lim2,axes=FALSE)
      axis(2, at=seq(0,1,length.out = 5), labels = seq(0,1,length.out=5))
      mtext(expression(s[y]), side = 2, line = 1.75, adj = 0.5, cex = 1, font=2)
      if(var==1){
        mtext(expression(lambda[1.2]), side = 2, line = 3, adj = 0.5, cex = 2, font=2,col="#0086FF")
      }else{
        mtext(expression(lambda[2.2]), side = 2, line = 3, adj = 0.5, cex = 2, font=2,col="#0086FF")
      }
    }else{
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(exp(simcov2[[bb]]),col=tim.colors(),zlim=lim2,axes=FALSE)
    }
  }
  for(bb in 1:5){
    screen(bb+16)
    if(bb==1){
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(pi/2 * exp(simcov3[[bb]])/(1 + exp(simcov3[[bb]])),col=tim.colors(),zlim=lim3,axes=FALSE)
      axis(1, at=seq(0,1,length.out = 5), labels = seq(0,1,length.out=5))
      axis(2, at=seq(0,1,length.out = 5), labels = seq(0,1,length.out=5))
      mtext(expression(s[x]), side = 1, line = 2, adj = 0.5, cex = 1, font=2)
      mtext(expression(s[y]), side = 2, line = 1.75, adj = 0.5, cex = 1, font=2)
      if(var==1){
        mtext(expression(phi[1]), side = 2, line = 3, adj = 0.5, cex = 2, font=3,col="#0086FF")
      }else{
        mtext(expression(phi[2]), side = 2, line = 3, adj = 0.5, cex = 2, font=3,col="#0086FF")
      }
    }else{
      par(pty="s") 
      par(mai=c(0.2,0.2,0.2,0.2))
      image(pi/2 * exp(simcov3[[bb]])/(1 + exp(simcov3[[bb]])),col=tim.colors(),zlim=lim3,axes=FALSE)
      axis(1, at=seq(0,1,length.out = 5), labels = seq(-0.5,0.5,length.out=5))
      mtext(expression(s[x]), side = 1, line = 2, adj = 0.5, cex = 1, font=2)
    }
  }
  screen(4)
  x1 <- c(0.3,0.4,0.4,0.3)
  y1 <- c(0.3,0.3,0.75,0.75)
  legend.gradient2(cbind(x1,y1),cols = colorRampPalette(colors) (colsteps), title = "", 
                   limits = seq(lim1[1],lim1[2],length.out = 5), cex=1)
  
  screen(5)
  legend.gradient3(cbind(x1,y1),cols = colorRampPalette(colors) (colsteps), title = "", 
                   limits = seq(lim2[1],lim2[2],length.out = 3), cex=1)
  
  screen(6)
  legend.gradient2(cbind(x1,y1),cols = colorRampPalette(colors) (colsteps), title = "", 
                   limits = seq(lim3[1],lim3[2],length.out = 5), cex=1)
  
  close.screen( all=TRUE)
  dev.off()
}
