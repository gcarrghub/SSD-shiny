genericPLL.AddOne <- function(
  pllData,
  speciesLabels,
  distName,
  effectLevel=0.05,
  gridSize=50,
  titleString=NULL,
  italicFont=2,
  xlimVals=NULL,
  xlabString="",
  showUpper=FALSE,
  roundTo=0,
  doPlots=FALSE,
  par1.LB=-Inf,
  par2.LB=-Inf,
  confLevel=0.95,
  logTransform=TRUE,
  startVals=NULL,
  quietTF=TRUE){
  
  if(logTransform)pllData <- log(pllData)
  #reorder
  newOrder <- order(pllData)
  pllData <- pllData[newOrder]
  speciesLabels <- speciesLabels[newOrder]
  
  dFUN <- get(x=paste("d",distName,sep=""),pos=1)
  qFUN <- get(x=paste("q",distName,sep=""),pos=1)
  pFUN <- get(x=paste("p",distName,sep=""),pos=1)
  
  #here, instead of survreg, etc, simply perform my own maximum likelihood, borrowing the rriskFitdist
  #function from package rriskDistributions (normal, logistic, weibull, gamma, etc are available)
  #the ONE restriction at the moment is that it must be a two-parameter distribution
  generalLL <- function(pars,x){
    sum(dFUN(x,pars[1],pars[2],log=TRUE))
  }
  if( is.null(startVals))MLE.fit <- rriskFitdist.GJC(data=pllData,distr=distName)
  if(!is.null(startVals))MLE.fit <- rriskFitdist.GJC(data=pllData,distr=distName,start=startVals)
  MLE.pars <- MLE.fit$estimate
  if(!is.null(startVals))names(MLE.pars)<-names(formals(dFUN))[2:3]
  MLE.pars.sd <- MLE.fit$sd
  MLE.LL <- MLE.fit$loglik
  #CI.critval <- qchisq(confLevel,1)
  CI.critval <- qf(confLevel,1,MLE.fit$n)
  
  #this is going to be the biggest issue by this approach.  Should probably walk across the parameter space by
  #sd units one at a time (separately for each direction), until the convex hull of the contour is totally enclosed
  #ie, all edges of the parameter space rectangle are OUTSIDE the confidence region
  
  #the figure is also informative; this region is independent of ECx level chosen
  enclosedTF <- FALSE
  par1.down <- par1.up <- par2.down <- par2.up <- 1
  regionSize <- gridSize
  while(!enclosedTF){
    testGrid <- expand.grid(
      par1Vals<-seq(max(par1.LB,MLE.pars[1]-par1.down*MLE.pars.sd[1]),MLE.pars[1]+par1.up*MLE.pars.sd[1],length=regionSize),
      par2Vals<-seq(max(par2.LB,MLE.pars[2]-par2.down*MLE.pars.sd[2]),MLE.pars[2]+par2.up*MLE.pars.sd[2],length=regionSize))
    testGrid.LL <- matrix(2*(apply(testGrid,1,FUN=function(pars){generalLL(pars,pllData)})-MLE.LL),ncol=regionSize)
    testGrid.HC05 <- matrix(apply(testGrid,1,FUN=function(pars){qFUN(0.05,pars[1],pars[2])}),ncol=regionSize)
    edges <- c(par2low<-testGrid.LL[,1],par2high<-testGrid.LL[,regionSize],par1low<-testGrid.LL[1,],par1high<-testGrid.LL[regionSize,])
    if(!any(abs(edges)<CI.critval))enclosedTF <- TRUE
    if(any(abs(edges)<CI.critval)){
      if(any(abs(par2low )<CI.critval))par2.down <- par2.down+0.5
      if(any(abs(par2high)<CI.critval))par2.up   <- par2.up  +0.5
      if(any(abs(par1low )<CI.critval))par1.down <- par1.down+0.5
      if(any(abs(par1high)<CI.critval))par1.up   <- par1.up  +0.5
    }
    if(!quietTF)print(c(par1.down=par1.down,par1.up=par1.up,par2.down=par2.down,par2.up=par2.up))
  }
  testGrid.LL.vec <- as.vector(testGrid.LL)
  #cLine is the bounding line (in x,y pairs) for the confidence region on the distribution
  #these points can then be querried to ask "what is the smallest ECx in this confidence space,
  #and what is the largest, to translate the region to limits on the desired ECx.
  cLine <- contourLines(par1Vals,par2Vals,abs(testGrid.LL),levels=CI.critval)[[1]]
  
  #columns are levels of par2, rows are levels of par1
  if(FALSE)matrix(testGrid.LL,nrow=length(par1Vals),ncol=length(par1Vals))
  
  if(doPlots | FALSE){
    par(omi=c(0,0,0,0))
    
    plot(x=testGrid[,1],y=testGrid[,2],xlab=names(MLE.pars)[1],ylab=names(MLE.pars)[2])
    points(x=testGrid[abs(testGrid.LL.vec)<CI.critval,1],y=testGrid[abs(testGrid.LL.vec)<CI.critval,2],col="blue",pch=16)
    lines(x=cLine$x,y=cLine$y,col="cyan",lwd=3)
    mtext(side=1,outer=TRUE,text=paste("Distribution:",distName),line=-1,adj=0)
    
    #library(rgl)
    open3d()
    baseRadius <- max(diff(range(testGrid[,1])),diff(range(testGrid[,2])),diff(range(as.vector(testGrid.HC05))))/200
    maxPoint <- which.max(as.vector(testGrid.HC05)[abs(testGrid.LL.vec)<CI.critval])
    minPoint <- which.min(as.vector(testGrid.HC05)[abs(testGrid.LL.vec)<CI.critval])
    spheres3d(x=testGrid[abs(testGrid.LL.vec)<CI.critval,1],y=testGrid[abs(testGrid.LL.vec)<CI.critval,2],z=as.vector(testGrid.HC05)[abs(testGrid.LL.vec)<CI.critval],col="blue",radius=baseRadius)
    spheres3d(x=testGrid[abs(testGrid.LL.vec)>CI.critval,1],y=testGrid[abs(testGrid.LL.vec)>CI.critval,2],z=as.vector(testGrid.HC05)[abs(testGrid.LL.vec)>CI.critval],alpha=0.2,radius=baseRadius)    
    spheres3d(x=testGrid[abs(testGrid.LL.vec)<CI.critval,1][maxPoint],y=testGrid[abs(testGrid.LL.vec)<CI.critval,2][maxPoint],z=as.vector(testGrid.HC05)[abs(testGrid.LL.vec)<CI.critval][maxPoint],col="red",radius=baseRadius*1.5)
    spheres3d(x=testGrid[abs(testGrid.LL.vec)<CI.critval,1][minPoint],y=testGrid[abs(testGrid.LL.vec)<CI.critval,2][minPoint],z=as.vector(testGrid.HC05)[abs(testGrid.LL.vec)<CI.critval][minPoint],col="red",radius=baseRadius*1.5)
    
    axes3d(labels=TRUE)
    aspect3d(1,1,1)
    
    open3d()
    surface3d(par1Vals,par2Vals,testGrid.LL,color="red",alpha=0.3)
    baseRadius <- max(diff(range(testGrid[,1])),diff(range(testGrid[,2])),diff(range(testGrid.LL.vec)))/200
    spheres3d(x=testGrid[abs(testGrid.LL.vec)<CI.critval,1],y=testGrid[abs(testGrid.LL.vec)<CI.critval,2],z=testGrid.LL.vec[abs(testGrid.LL.vec)<CI.critval],col="blue",radius=2*baseRadius)
    spheres3d(x=testGrid[abs(testGrid.LL.vec)>CI.critval,1],y=testGrid[abs(testGrid.LL.vec)>CI.critval,2],z=testGrid.LL.vec[abs(testGrid.LL.vec)>CI.critval],radius=baseRadius)  
    spheres3d(x=MLE.pars[1],y=MLE.pars[2],z=0,col="red",radius=baseRadius*3)
    lines3d(cLine$x,cLine$y,rep(-CI.critval,length(cLine$x)))
    axes3d(labels=TRUE)
    aspect3d(1,1,1)
  }
  
  
  #the interval, for a range of quantiles so confidence bands can be plotted
  #10^do.call(qFUN,as.list(c(p=.05,MLE.pars)))
  #10^qFUN(.05,MLE.pars[1],MLE.pars[2])
  pVals <- seq(.001,.999,by=.001)
  #pVals <- seq(.01,.99,by=.01)
  if( logTransform)ciVals <- t(sapply(pVals,FUN=function(ppp)exp(range(apply(cbind(cLine$x,cLine$y),1,FUN=function(pars)do.call(qFUN,as.list(c(p=ppp,pars))))))))
  if(!logTransform)ciVals <- t(sapply(pVals,FUN=function(ppp)   range(apply(cbind(cLine$x,cLine$y),1,FUN=function(pars)do.call(qFUN,as.list(c(p=ppp,pars)))))))
  #ciVals[pVals==0.05,]
  if(!quietTF)print(head(ciVals))
  
  if(doPlots){
    require(ADGofTest)
    print(ad.test(pllData,distr=pFUN,MLE.pars[1],MLE.pars[2]))
  }
  if(!quietTF)print(MLE.pars)
  if( logTransform)modelVals <- exp(do.call(qFUN,c(list(p=pVals),as.list(MLE.pars))))
  if(!logTransform)modelVals <-    do.call(qFUN,c(list(p=pVals),as.list(MLE.pars)))
  if(is.null(xlimVals)){
    if( logTransform)xlims <- range(log(modelVals))
    if(!logTransform)xlims <- range(   (modelVals))
    xlims[1] <- floor(na.omit(c(min(xlims[1],pllData))))
    xlims[2] <-  ceiling(na.omit(c(max(xlims[2],pllData))))
    if( logTransform)xlims <- exp(xlims)
  }
  if(!is.null(xlimVals))xlims <- xlimVals
  
  if( logTransform)HC05All <- exp(do.call(qFUN,c(list(p=effectLevel),as.list(MLE.pars))))
  if(!logTransform)HC05All <-    do.call(qFUN,c(list(p=effectLevel),as.list(MLE.pars)))
  LCL.HC05All <- ciVals[pVals==effectLevel,1]
  UCL.HC05All <- ciVals[pVals==effectLevel,2]
  #find params that go with ci endpoints
  ciParms05 <- cbind(cLine$x,cLine$y,apply(cbind(cLine$x,cLine$y),1,FUN=function(pars)do.call(qFUN,as.list(c(p=effectLevel,pars)))))
  ciParms05.LCL <- ciParms05[which.min(ciParms05[,3]),1:2]
  ciParms05.UCL <- ciParms05[which.max(ciParms05[,3]),1:2]
  if(!quietTF)print(rbind(ciParms05.LCL,ciParms05.UCL))
  
  
  
  if(doPlots){
    par(omi=c(0,0,0.5,2.5))
    
    if( logTransform)plot(x=modelVals,y=pVals,ylim=c(0,1),log='x',
                          type='n',xlab=xlabString,ylab="Probability",axes=F,main=titleString,xlim=xlims)
    if(!logTransform)plot(x=modelVals,y=pVals,ylim=c(0,1),
                          type='n',xlab=xlabString,ylab="Probability",axes=F,main=titleString,xlim=xlims)
    mtext(side=1,outer=TRUE,text=paste("Distribution:",distName),line=-1,adj=0)
    axis(side=2)
    if(logTransform){
      majorTicks <- 10^seq(ceiling(par("usr")[1]),floor(par("usr")[2]))
      axis(side=1,at=majorTicks,label=majorTicks)
      axis(side=1,at=10^as.vector(sapply(ceiling(par("usr")[1]):(ceiling(par("usr")[2])-1),FUN=function(x)x+log10(1:9))),
           tck=-.01,labels=F)
    }
    if(!logTransform)axis(side=1)
    box()
    #      mtext(side=1,expression(paste("Concentration (",mu,"g/L)")),line=3)
    if( logTransform)points(x=exp(pllData),y=order(pllData)/(max(order(pllData))+1),col='gray',pch=16)
    if(!logTransform)points(x=   pllData,y=order(pllData)/(max(order(pllData))+1),col='gray',pch=16)
    
    lines(x=ciVals[,1],y=pVals,col="blue",lwd=3)
    lines(x=ciVals[,2],y=pVals,col="blue",lwd=3)
    lines(x=modelVals,y=pVals,lwd=3)
    
    if(!quietTF & logTransform){
      lines(x=exp(qFUN(pVals,ciParms05.LCL[1],ciParms05.LCL[2])),y=pVals,col="cyan")
      lines(x=exp(qFUN(pVals,ciParms05.UCL[1],ciParms05.UCL[2])),y=pVals,col="cyan")
    }
    
    #predictedValues <- exp((log(pllData) - interceptParam)/scaleParam)
    mtext(side=4,at=order(pllData)/(max(order(pllData))+1),text=speciesLabels,las=2,font=italicFont,cex=.7,line=1)
    #points(x=(pllData),y=1 - 1/(1 + predictedValues),col=8,pch=16)
    if(!showUpper)abline(h=effectLevel,v=c(LCL.HC05All,HC05All),col='gray')
    if(showUpper)abline(h=effectLevel,v=c(LCL.HC05All,HC05All,UCL.HC05All),col='gray')
    CItexts <- format(sapply(c(LCL.HC05All,HC05All,UCL.HC05All),signif,roundTo))
    mtext(side=3,at=HC05All*.85,text=bquote(HC[5]==.(CItexts[2])),adj=1,las=2,line=-.1)
    mtext(side=3,at=LCL.HC05All*.85,text=bquote(LCL==.(CItexts[1])),adj=1,las=2,line=-.1)
    if(showUpper)mtext(side=3,at=UCL.HC05All*.85,text=bquote(UCL==.(CItexts[3])),adj=1,las=2,line=-.1)
    par(omi=rep(0,4))
    
    #xVals comes from the plot range, which would always be log10
    xVals <- seq(par("usr")[1]-1,par("usr")[2],length=1000)
    ylimMax <- max(c(
      dFUN(log(10^xVals),ciParms05.LCL[1],ciParms05.LCL[2]),
      dFUN(log(10^xVals),ciParms05.UCL[1],ciParms05.UCL[2]),
      dFUN(log(10^xVals),MLE.pars[1],MLE.pars[2])))
    
    if( logTransform){
      plot(x=modelVals,y=rep(0,length(modelVals)),ylim=c(0,ylimMax),log='x',
           type='n',xlab=xlabString,ylab="Distribution Density",axes=F,main=titleString,xlim=xlims)
      majorTicks <- 10^seq(ceiling(par("usr")[1]),floor(par("usr")[2]))
      axis(side=1,at=majorTicks,label=majorTicks)
      axis(side=1,at=10^as.vector(sapply(ceiling(par("usr")[1]):(ceiling(par("usr")[2])-1),FUN=function(x)x+log10(1:9))),
           tck=-.01,labels=F)
      axis(side=2)
      lines(y=dFUN(log(10^xVals),MLE.pars[1],MLE.pars[2]),x=10^(xVals),lwd=2)
      lines(y=dFUN(log(10^xVals),ciParms05.LCL[1],ciParms05.LCL[2]),x=10^(xVals),col="cyan",lwd=2)
      lines(y=dFUN(log(10^xVals),ciParms05.UCL[1],ciParms05.UCL[2]),x=10^(xVals),col="cyan",lwd=2,lty=2)
      rug(side=1,x=exp(pllData),lwd=2,col=rgb(0,0,0,.2))
      abline(v=HC05All,col='gray')
      abline(v=LCL.HC05All,col='cyan')
      if(showUpper)abline(v=UCL.HC05All,col='cyan',lty=2)
      points(x=HC05All,y=dFUN(log(HC05All),MLE.pars[1],MLE.pars[2]))
      points(x=LCL.HC05All,y=dFUN(log(LCL.HC05All),ciParms05.LCL[1],ciParms05.LCL[2]),col="cyan")
      points(x=UCL.HC05All,y=dFUN(log(UCL.HC05All),ciParms05.UCL[1],ciParms05.UCL[2]),col="cyan")
      CItexts <- format(sapply(c(LCL.HC05All,HC05All,UCL.HC05All),signif,roundTo))
      mtext(side=3,at=HC05All*.85,text=bquote(HC[5]==.(CItexts[2])),adj=1,las=2,line=-.1)
      mtext(side=3,at=LCL.HC05All*.85,text=bquote(LCL==.(CItexts[1])),adj=1,las=2,line=-.1)
      if(showUpper)mtext(side=3,at=UCL.HC05All*.85,text=bquote(UCL==.(CItexts[3])),adj=1,las=2,line=-.1)
      box()
    }
    
    
  }
  list(fit=c(X=effectLevel,HCx=HC05All,LowerCL=LCL.HC05All,UpperCL=UCL.HC05All,MLE.pars),CI.line=cbind(pVals,ciVals))
}
#x <- genericPLL.AddOne(foo$Norm.C12.C13,speciesLabels=foo$Species,distName="logis")
#head(x[["CI.line"]])
#plot(x=as.vector(x[["CI.line"]][,2:3]),y=rep(x[["CI.line"]][,1],2),log='x')


#genericPLL.AddOne(foo$Norm.C12.C13,speciesLabels=foo$Species,distName="logis",doPlots=TRUE,roundTo=3)
