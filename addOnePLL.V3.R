#this program ASSUMES that the data are not log-transformed, but the analysis will be on log scale
#in larger data sets, the added dose gets so far away that errors occur in distribution estimation...

addOnePLL <- function(effectValues,speciesLabels,distName,dataTag,xlabString="AS Concentration (mg/L)",effectLevelA1=0.05,debugTF=FALSE){
  #climEnd is for built-in distributions (typically either "norm" or "logis")
  #distName is for the survreg call argument dist (typically either "lognormal" or "loglogistic", paired with climEnd)
  standardPlot <- function(caseData,speciesLabels,distName,climEnd,titleString=NULL,italicFont=2,
                           xlimVals="NULL",showUpper=FALSE,roundTo=NULL,unitStr=xlabString){
    #print(caseData)
    distributionFitPLL <- genericPLL.AddOne(caseData,speciesLabels,distName,logTransform=TRUE,effectLevel=effectLevelA1)
    #require(survival)
    #distributionFitPLL<-survreg(formula = Surv(caseData) ~ 1, na.action = na.exclude, dist = distName)
    #if(debugTF)print(distributionFitPLL)
    interceptParamAll <- interceptParam <- distributionFitPLL[[1]][5]
    scaleParamAll <- scaleParam <- distributionFitPLL[[1]][6]
    par(omi=c(0,0,0.5,2.5))
    probs <- distributionFitPLL[[2]][,1]
    qFUN <- get(paste("q",distName,sep=""))
    pFUN <- get(paste("p",distName,sep=""))
    
    
    #if(require(ADGofTest))print(ad.test(log(caseData),distr=pFUN,interceptParamAll,scaleParamAll))
    
    quantiles <- exp(qFUN(probs,interceptParamAll,scaleParamAll))
    #climTableCase <- get(paste("myCLTable",climEnd,length(caseData),sep="."),pos=1)
    #this to show all confidence limits
    percentileStats.p <- distributionFitPLL[[2]][,1]
    percentileStatsLCLall <- distributionFitPLL[[2]][,2]
    percentileStatsUCLall <- distributionFitPLL[[2]][,3]
    #this for specific HC05
    HCxAll<-distributionFitPLL[[1]][2]
    LCL.HCxAll <- distributionFitPLL[[1]][3]
    UCL.HCxAll <- distributionFitPLL[[1]][4]
    #if(debugTF)print(c(HCxAll,LCL.HCxAll,UCL.HCxAll))
    #print(range(quantiles))
    if(is.character(xlimVals)){
      xlims <- range(quantiles[probs>=.01 & probs<=0.99])
      xlims[1] <- floor(min(na.omit(log10(c(LCL.HCxAll,xlims[1],caseData)))))
      xlims[2] <-  ceiling(max(na.omit(log10(c(xlims[2],caseData)))))
      xlims <- 10^xlims
    }else
      xlims <- xlimVals
    #no matter what, the lower xlim must encompass the LCL:
    if(LCL.HCxAll<xlims[1])xlims[1]<-10^floor(log10(LCL.HCxAll))
    plot(x=quantiles,y=probs,ylim=c(0,1),log='x',
         type='n',xlab=unitStr,ylab="Probability",axes=F,main=titleString,xlim=xlims)
    axis(side=2)
    majorTicks <- 10^seq(ceiling(par("usr")[1]),floor(par("usr")[2]))
    axis(side=1,at=majorTicks,label=majorTicks)
    box()
    #      mtext(side=1,expression(paste("Concentration (",mu,"g/L)")),line=3)
    axis(side=1,at=10^as.vector(sapply(ceiling(par("usr")[1]):(ceiling(par("usr")[2])-1),FUN=function(x)x+log10(1:9))),
         tck=-.01,labels=F)
    points(x=(caseData),y=order(caseData)/(max(order(caseData))+1),col='gray',pch=16)
    
    mySpline<-spline(x=log(percentileStatsLCLall),y=percentileStats.p)
    mySpline$x<-exp(mySpline$x)
    lines(mySpline,col="blue",lwd=3)
    mySpline<-spline(x=log(percentileStatsUCLall),y=percentileStats.p)
    mySpline$x<-exp(mySpline$x)
    lines(mySpline,col="blue",lwd=3)
    lines(x=(quantiles),y=probs,lwd=3)
    
    if(is.null(roundTo)){
      roundTo <- -(floor(log10(LCL.HCxAll))-2)
    }
    #predictedValues <- exp((log(caseData) - interceptParam)/scaleParam)
    mtext(side=4,at=order(caseData)/(max(order(caseData))+1),text=speciesLabels,las=2,font=italicFont,cex=.7,line=1)
    #points(x=(caseData),y=1 - 1/(1 + predictedValues),col=8,pch=16)
    if(!showUpper)abline(h=effectLevelA1,v=c(LCL.HCxAll,HCxAll),col='gray')
    if(showUpper)abline(h=effectLevelA1,v=c(LCL.HCxAll,HCxAll,UCL.HCxAll),col='gray')
    CItexts <- format(sapply(c(LCL.HCxAll,HCxAll,UCL.HCxAll),signif,roundTo))
    mtext(side=3,at=HCxAll*.85,text=bquote(HC[.(round(100*effectLevelA1))]==.(CItexts[2])),adj=1,las=2,line=-.1)
    mtext(side=3,at=LCL.HCxAll*.85,text=bquote(LCL==.(CItexts[1])),adj=1,las=2,line=-.1)
    if(showUpper)mtext(side=3,at=UCL.HCxAll*.85,text=bquote(UCL==.(CItexts[3])),adj=1,las=2,line=-.1)
    
    c(HCx=HCxAll,LowerCL=LCL.HCxAll,UpperCL=UCL.HCxAll)
  }
  
  
  #this program will assume that tox values are given in dose units (NOT log transformed), and that they
  #will be transformed by ln(), or log()
  
  #require(survival)
  caseData <- effectValues
  caseTaxa <- speciesLabels
  newOrder <- order(caseData)
  caseData <- caseData[newOrder]
  caseTaxa <- caseTaxa[newOrder]
  
  dataTable <- data.frame(caseData,caseTaxa)
  row.names(dataTable) <- 1:nrow(dataTable)
  
  qFUN <- get(paste("q",distName,sep=""))
  pFUN <- get(paste("p",distName,sep=""))
  
  
  distributionFit <- rriskFitdist.GJC(data=log(caseData),distr=distName)
  #distributionFit<-survreg(formula = Surv(caseData) ~ 1, na.action= na.exclude, dist = distName)
  interceptParamAll <- interceptParam <- distributionFit$estimate[1]
  scaleParamAll <- scaleParam <- distributionFit$estimate[2]
  HCxAll<-exp(qFUN(effectLevelA1,interceptParamAll,scaleParamAll))
  if(debugTF)print(distributionFit$estimate)
  if(debugTF)print(c(HCxAll=HCxAll))
  
  #solve for the value of a new observation that will shift the HC05 to a new target value
  #input candidate value is on log-scale.  When the point is added, we are looking for the HCx
  #value to shift lower.
  #These were tested because I had a bug in the code, where the second argument to qFUN in each
  #of these referred to distributionFitA1$estimate[2] rather than distributionFit$estimate[2],
  #where distributionFit is for the ORIGINAL data, whereas it is supposed to be changing in here
  funOpt.old<-function(addedPoint,CaseData,HCxTargetValue){
    #print(c(HCxgiven=HCxgiven))
    distributionFitA1<-rriskFitdist.GJC(log(c(exp(addedPoint),CaseData)), dist = distName)
    #log(predict(distributionFit,type="quantile",p=c(0.05),newdata=data.frame(1)))-HCxgiven
    #The HCx percentile, when adding the data value, should equal the target, eq HCx/2
    HCx.new <- qFUN(effectLevelA1,distributionFitA1$estimate[1],distributionFitA1$estimate[2])
    #print(c(HCx.new=HCx.new))
    HCx.new-HCxTargetValue
  }

  #different approach.  When the new point is added, it should result in a fit where the target HC level
  #is returned when we put in the target HCx value, eg, HCx/2
  funOpt<-function(addedPoint,CaseData,HCxTargetValue){
    #print(c(HCxgiven=HCxgiven))
    distributionFitA1<-rriskFitdist.GJC(log(c(exp(addedPoint),CaseData)), dist = distName)
    #log(predict(distributionFit,type="quantile",p=c(0.05),newdata=data.frame(1)))-HCxgiven
    #The HCx percentile, when adding the data value, should equal the target, eq HCx/2
    HCx.new <- pFUN(HCxTargetValue,distributionFitA1$estimate[1],distributionFitA1$estimate[2])
    #print(c(HCx.new=HCx.new))
    effectLevelA1-HCx.new
  }
    
  #this assumes that the added point will be less than the HC5 (upper=), but that is not always the case. An alternative is to take the
  #GREATER of the HC5, and the smallest data value.
  #upperValue <- log(max(HCxAll,min(caseData)))
  print(caseData)
  upperValue <- log(max(HCxAll,max(caseData)))
  #checks on interval -- should be opposite signs
  print(c(effectLevelA1=effectLevelA1,HCxAll=HCxAll))
  for(logVal2 in seq(5,200,by=5)){
    if(funOpt(log(HCxAll)-logVal2,caseData,log(HCxAll/2))<0)break
  }
  for(logVal3 in seq(5,200,by=5)){
    if(funOpt(log(HCxAll)-logVal3,caseData,log(HCxAll/3))<0)break
  }
  for(logVal5 in seq(5,200,by=5)){
    if(funOpt(log(HCxAll)-logVal5,caseData,log(HCxAll/5))<0)break
  }
  for(logVal10 in seq(5,200,by=5)){
       if(funOpt(log(HCxAll)-logVal10,caseData,log(HCxAll/10))<0)break
  }
  print(c(logVal2=logVal2,logVal3=logVal3,logVal5=logVal5,logVal10=logVal10))
  print(
    intervalChecks <- rbind(
      c(lower=try(funOpt(log(HCxAll)-logVal2,caseData,log(HCxAll/2))),upper=try(funOpt(upperValue,caseData,log(HCxAll/2)))),
      c(lower=try(funOpt(log(HCxAll)-logVal3,caseData,log(HCxAll/3))),upper=try(funOpt(upperValue,caseData,log(HCxAll/3)))),
      c(lower=try(funOpt(log(HCxAll)-logVal5,caseData,log(HCxAll/5))),upper=try(funOpt(upperValue,caseData,log(HCxAll/5)))),
      c(lower=try(funOpt(log(HCxAll)-logVal10,caseData,log(HCxAll/10))),upper=try(funOpt(upperValue,caseData,log(HCxAll/10))))
    )
  )
  if(!all(rowSums(sign(intervalChecks))==0))stop("Dataset may require add-one-in values too far from observed data.")
  addPoints <- c(
    HC5.2=try(uniroot(funOpt,lower=log(HCxAll)-logVal2,upper=upperValue,CaseData=caseData,HCxTargetValue=log(HCxAll/2),extendInt="no")$root),
    HC5.3=try(uniroot(funOpt,lower=log(HCxAll)-logVal3,upper=upperValue,CaseData=caseData,HCxTargetValue=log(HCxAll/3),extendInt="no")$root),
    HC5.5=try(uniroot(funOpt,lower=log(HCxAll)-logVal5,upper=upperValue,CaseData=caseData,HCxTargetValue=log(HCxAll/5),extendInt="no")$root),
    HC5.10=try(uniroot(funOpt,lower=log(HCxAll)-logVal10,upper=upperValue,CaseData=caseData,HCxTargetValue=log(HCxAll/10),extendInt="no")$root))
  
  print(addPoints)
  print(origPars <- rriskFitdist.GJC(log(c(caseData)), dist = distName)$estimate)
  print(pFUN(log(HCxAll),origPars[1],origPars[2]))
  print(newPars <- rriskFitdist.GJC(log(c(exp(addPoints[1]),caseData)), dist = distName)$estimate)
  print(pFUN(log(HCxAll/2),newPars[1],newPars[2]))
  percentageSTR <- round(100*effectLevelA1)
  percentageSTR <- formatC(percentageSTR,width=2,format="d",flag="0")
  addProbs <- pFUN(addPoints,interceptParamAll,scaleParamAll)
  names(addProbs) <- paste(paste("HC",percentageSTR,sep=""),c(2,3,5,10),sep=".")
  oneIn <- 1/addProbs
  #print(addPoints)
  #print(exp(addPoints))
  
  ### calculate a goodness of fit -- not used
  #ks.p <- stats::ks.test(x = (log(caseData)+seq(-1,1,along=caseData)/10000), y = "plogis", location = interceptParam, scale = scaleParam)$p.value
  ###define SSD curve
  probs <- c(1e-10,seq(.0001,.9999,by=.0001))
  #quantiles <- (predict(distributionFit,type="quantile",p=probs,newdata=data.frame(1)))
  quantiles <- exp(qFUN(probs,interceptParamAll,scaleParamAll))
  ###PLOT
  xlims <- range(log10(quantiles))
  xlims[1] <- floor(xlims[1])
  xlims[2] <-  ceiling(xlims[2])
  xlims <- 10^xlims    
  plot(x=quantiles,y=probs,ylim=c(0,1),log='x',type='n',
       xlab=xlabString,ylab="Probability",axes=F,xlim=xlims)
  #print(par("usr"))
  axis(side=2)
  majorTicks <- 10^seq(ceiling(par("usr")[1]),floor(par("usr")[2]))
  axis(side=1,at=majorTicks,label=majorTicks)
  box()
  #     mtext(side=1,expression(paste("Concentration (",mu,"g/L)")),line=3)
  axis(side=1,
       at=10^as.vector(sapply(ceiling(par("usr")[1]):(floor(par("usr")[2])-1),
                              FUN=function(x)x+log10(1:9))),
       tck=-.01,labels=F)
  
  points(x=(caseData),y=order(caseData)/(max(order(caseData))+1),col='red',pch=16)
  
  #climTableCase <- get(paste("myCLTable",climEnd,length(caseData),sep="."),pos=1)
  #percentileStatsLCLall <- exp(interceptParam - climTableCase[["lower_c"]] * scaleParam)
  #percentileStatsUCLall <- exp(interceptParam - climTableCase[["upper_c"]] * scaleParam)
  #if(climEnd=="logis")HC05<-qlogis(0.05,interceptParam,scaleParam)#exp(scaleParam*log(0.05/(1-0.05)) + interceptParam)
  #if(climEnd=="norm")HC05<-qnorm(0.05,interceptParam,scaleParam)
  #mySpline<-spline(x=log(percentileStatsLCLall),y=unique(climTable$p))
  #mySpline$x<-exp(mySpline$x)
  #lines(mySpline,col=2,lwd=3)
  #mySpline<-spline(x=log(percentileStatsUCLall),y=unique(climTable$p))
  #mySpline$x<-exp(mySpline$x)
  #lines(mySpline,col=2,lwd=3)
  lines(x=quantiles,y=probs,lwd=3)
  abline(h=effectLevelA1,v=c(HCxAll,HCxAll/2,HCxAll/3,HCxAll/5,HCxAll/10),col='gray')
  #HC05All<-HC05
  
  #for (newCase in exp(seq(log(.001),log(HC05All),length=10))){
  for (newCase in exp(addPoints)){
    newDistributionFit<-rriskFitdist.GJC(log(newCaseData<-c(newCase,caseData)), dist = distName)
    #distributionFit<-survreg(formula = Surv(newCaseData<-c(newCase,caseData)) ~ 1, na.action = na.exclude, dist = distName)
    newQuantiles <- exp(qFUN(probs,newDistributionFit$estimate[1],newDistributionFit$estimate[2]))
    #newQuantiles <- (predict(distributionFit,type="quantile",p=probs,newdata=data.frame(1)))
    interceptParam <- newDistributionFit$estimate[1]
    scaleParam <- newDistributionFit$estimate[2]
    ### calculate a goodness of fit -- not used
    #ks.p <- stats::ks.test(x = log(newCaseData)+seq(-1,1,along=newCaseData)/1000000, y = "plogis", location = interceptParam, scale = scaleParam)$p.value
    #print(c(ks.p=ks.p))
    ###define SSD curve
    #	allEC20 <- exp(seq(from=log(.01),to=log(100),length=100))
    #	logEC20All <- log(allEC20)
    #	expTerm <- exp((logEC20All - interceptParam)/scaleParam)
    #	estimatedCdfAll <- 1 - 1/(1 + expTerm)
    #	tableCDF <- scaleParam*log(unique(climTable$p)/(1-unique(climTable$p))) + interceptParam 
    #	1 - 1/(1 + exp((tableCDF - interceptParam)/scaleParam))
    #	unique(climTable$p)
    #	ord<-order(newCaseData)
    #	phat<-ord/(max(ord)+1)
    
    
    #climTableCase <- get(paste("myCLTable",climEnd,length(newCaseData),sep="."),pos=1)
    #percentileStatsLCLall <- exp(interceptParam - climTableCase[["lower_c"]] * scaleParam)
    #percentileStatsUCLall <- exp(interceptParam - climTableCase[["upper_c"]] * scaleParam)
    #if(climEnd=="logis")HC05<-exp(scaleParam*log(0.05/(1-0.05)) + interceptParam)
    #if(climEnd=="norm")HC05<-qnorm(0.05,interceptParam,scaleParam)
    #mySpline<-spline(x=log(percentileStatsLCLall),y=unique(climTable$p))
    #mySpline$x<-exp(mySpline$x)
    #	lines(mySpline,col=2,lwd=.5)
    #mySpline<-spline(x=log(percentileStatsUCLall),y=unique(climTable$p))
    #mySpline$x<-exp(mySpline$x)
    #	lines(mySpline,col=2,lwd=.5)
    lines(x=newQuantiles,y=probs,lwd=.5)
    points(x=(newCaseData)[1],y=(order(newCaseData)/(max(order(newCaseData))+1))[1],col=8,pch=1)
    probNewCase <- pFUN(log(newCase),interceptParamAll,scaleParamAll)
    #print(format(1/probNewCase,digits=2,scientific=F,big.mark=","))
    
    oneInVal <- 1/probNewCase
    oneInVal <- signif(oneInVal,3)
    oneInVal <- sapply(oneInVal,FUN=function(x){
      ifelse(log10(x)>10,yes = format(x,scientific = TRUE),no = format(x,scientific = FALSE, big.mark=","))
    })
    
    text(
      x=newCaseData[1],
      y=.1,
      label=paste("1 / ",oneInVal,sep=""),
      adj=0,srt=90
    )
  }
  mtext(side=3,at=HCxAll,text=bquote(HC[.(round(100*effectLevelA1))]),adj=0,las=2)
  mtext(side=3,at=HCxAll/2,text=bquote(HC[.(round(100*effectLevelA1))]/2),adj=0,las=2)
  mtext(side=3,at=HCxAll/3,text=bquote(HC[.(round(100*effectLevelA1))]/3),adj=0,las=2)
  mtext(side=3,at=HCxAll/5,text=bquote(HC[.(round(100*effectLevelA1))]/5),adj=0,las=2)
  mtext(side=3,at=HCxAll/10,text=bquote(HC[.(round(100*effectLevelA1))]/10),adj=0,las=2)
  baseSTR <- paste("HC",round(100*effectLevelA1),sep="")
  
  oneInVal <- sapply(oneIn,FUN=function(x){
    ifelse(log10(x)>10,yes = format(signif(x,3),scientific = TRUE),no = format(signif(x,3),scientific = FALSE, big.mark=","))
  })
  results <- data.frame(target=paste(baseSTR,c(2,3,5,10),sep="/"),targetValue=HCxAll/c(2,3,5,10),addOneValue=signif(exp(addPoints),4),addOneProb=signif(addProbs,4),oneIn=oneInVal,
                        oneInFormated=paste(
                          "1 / ",
                          oneInVal,
                          sep=""))
  
  title(paste("Log",switch(distName,logis="Logistic",norm="Normal"),sep="-"),adj=0)
  mtext(side=1,outer=TRUE,line=-1,text=dataTag,adj=.05)
  
  #verification plots
  #print(caseData)
  standardPlot(caseData,as.character(caseTaxa),distName,xlimVals=" ")
  switch(distName,logis=title("Logistic (Original Data)",adj=0),norm=title("Normal (Original Data)",adj=0))
  sapply(1:nrow(results),FUN=function(i){
    x<-results[i,"addOneValue"]
    #print(results)
    #print(x)
    standardPlot(
      c(x,caseData),
      c(paste("Added value (",format(c(x,caseData),scientific=FALSE,digits=5)[1],")",sep=""),as.character(caseTaxa)),
      distName,showUpper=TRUE)
    points(x=x,y=1/(length(caseData)+2))
    
    #rug(side=1,x=x,col="red")
    switch(distName,
           logis=title(paste("Logistic (",as.character(results[i,"target"]),")",sep=""),adj=0),
           norm=title(paste("Normal (",as.character(results[i,"target"]),")",sep=""),adj=0))
  })
  data.frame(responseVar=rep(dataTag,4),distribution=rep(distName,4),results)
}

if(FALSE){
  addOnePLL(effectValues=foo[,2],speciesLabels=foo$Species,distName="logis",dataTag="C12-13 All Species",effectLevelA1=0.05)
  addOnePLL(effectValues=foo[,2],speciesLabels=foo[,1],distName="logis",dataTag="C12-13 All Species",effectLevelA1=0.5,debugTF=TRUE)
  addOnePLL(effectValues=foo$Norm.C12.C13,speciesLabels=foo$Species,distName="norm",dataTag="C12-13 All Species")
}


