#########
#########  Check parallel settings for use by others
#########

#the leave-one-out analysis can be very time consuming for large datasets
#even if the user says yes, we might want to warn for large datasets?
#multiple HCx levels will mulitply the time waiting as well.
doLeaveOneOut <- TRUE
#add-one-in is not as bad.
doAddOneIn <- TRUE

HC.primary <- 0.4
#########  HC5 AND HC50?
library(tcltk)
if(FALSE){
  #if interactive...
  wb  <- loadWorkbook(file.choose())
  foo <- readWorksheet(wb,sheetName <- select.list(getSheets(wb),graphics=TRUE))
  resVar <- names(foo)[2]
  
  tagString <- substring(sheetName,first=1,last=nchar(sheetName)-4)
}

pageBreak <- function(messageString){
  plot(x=0,y=0,type="n",axes=FALSE,xlab="",ylab="")
  text(x=0,y=0,labels=messageString,adj=0.5,cex=3)  
}

#source files that do the specific tasks
source("genericPLL.R")
source("addOne.genericPLL.R")
source("addOnePLL.V3.R")
source("rriskFitdist.GJC.R")

fileFilters <- matrix(c(
  "Excel", ".xls*", "All files", "*"),
  ncol=2, byrow = TRUE)

library(XLConnect)
#dataWB <- "AS data sheets.xlsx"
dataWB <- tcltk::tk_choose.files(caption="Input File",multi=FALSE,default=paste(getwd(),"/*.xls*",sep=""),filters=fileFilters)
pathPart <- substring(dataWB,first = 1,last = tail(gregexpr(pattern = "\\/",text = dataWB)[[1]],1))
baseName <- unlist(strsplit(tail(unlist(strsplit(dataWB,"\\/")),1),"\\."))[1]
#GJC modification -- default file is built on name of input file
outputWB <- paste(pathPart,baseName,paste(" HC",100*HC.primary," Output.xlsx",sep=""))
#now let them change it if they want
outputWB <- filenameDialog("Name of output file",question = "Filename",entryInit = outputWB,returnValOnCancel = NULL)

masterWB <- outputWB
#problem here is that every sheet will have the same label for graphics, so all the sheets have to be somewhat related
concentrationLabel <- "AS Concentration (mg/L)"

if(file.exists(masterWB))unlink(masterWB)

wb.out <- loadWorkbook(masterWB,create=TRUE)
wb.data <- loadWorkbook(dataWB)

#the program will go through each sheet, and analyze each set of data, one per sheet
for(sheetName in getSheets(wb.data)){
  print(sheetName)
  #GJC modification -- check for empty, or nonsense sheets (2 or fewer rows)
  #now, don't try them (errors for empty sheets now won't occur)
  if(getLastRow(wb.data,sheetName)<=2)next
  foo <- readWorksheet(wb.data,sheetName)
  #data are assumed to have two columns, first is species label, second is data
  speciesVar <- names(foo)[1]
  resVar <- names(foo)[2]
  foo[,resVar] <- as.numeric(as.character(foo[,resVar]))
  badRows <- which(rowSums(is.na(foo))>0)
  if(length(badRows)>0){
    cat("\nRow(s)",badRows,"will be dropped from the file.  Please inspect the data.\n")
    foo <- na.omit(foo)
  }
  
  
  tagString <- sheetName
  #tagString <- substring(sheetName,first=1,last=nchar(sheetName)-4)
  
  #use XLConnect::writeWorksheet to put data into a workbook.
  #do this once for each file, like the plots, and let Scott go from there
  createSheet(wb.out,name=tagString)
  
  pdf(file=paste(baseName," - ",tagString,".pdf",sep=""),height=8.5,width=11)
  #full data fits -- plots directly from this
  #if both 50 and 95th are estimated, be sure to add lines to space out results
  #from subsequent tables
  pageBreak("Data Fits")
  #full data fits -- plots directly from this
  pageBreak("A:  Logistic Data Fits")
  fit.out1 <- rbind(
    #genericPLL(foo[,resVar],speciesLabels=foo[,speciesVar],distName="logis",doPlots=TRUE,xlabString=concentrationLabel),
    genericPLL(foo[,resVar],speciesLabels=foo[,speciesVar],distName="logis",doPlots=TRUE,xlabString=concentrationLabel,effectLevel = HC.primary)
    )
  pageBreak("B:  Normal Data Fits")
  fit.out2 <- rbind(
    #genericPLL(foo[,resVar],speciesLabels=foo[,speciesVar],distName="norm",doPlots=TRUE,xlabString=concentrationLabel),
    genericPLL(foo[,resVar],speciesLabels=foo[,speciesVar],distName="norm",doPlots=TRUE,xlabString=concentrationLabel,effectLevel = HC.primary)
  )
  fit.out <- rbind(fit.out1,fit.out2)
  fit.out2 <- data.frame(
    Distribution=rep(c("Logistic","Normal"),times=c(nrow(fit.out1),nrow(fit.out2))),
    fit.out)
  
  head.out <- data.frame(spreadsheet=tagString)
  writeWorksheet(wb.out,data.frame(spreadsheet=c(tagString,"SSD Analysis of Full Data")),sheet=tagString,startRow=1,startCol=1)
  writeWorksheet(wb.out,foo,sheet=tagString,startRow=1,startCol=10)
  writeWorksheet(wb.out,fit.out2,sheet=tagString,startRow=4)
  rowCount <- nrow(fit.out2)+9
  
  #dev.off()
  #saveWorkbook(wb.out)
  
  #leave-one-out analyses
  clusterExport(snowCluster,varlist=c("genericPLL.AddOne","foo", "rriskFitdist.GJC","resVar","speciesVar", "HC.primary"))
  pageBreak("A:  Logistic Leave-One-Out")
  LOO.results.logis <- clusterApplyLB(snowCluster,x=0:nrow(foo),fun=function(i){
       if(i==0)result <- genericPLL.AddOne(foo[,resVar],speciesLabels=foo[,speciesVar],distName="logis",doPlots=FALSE,effectLevel = HC.primary)
       if(i>=1)result <- genericPLL.AddOne(foo[,resVar][-i],speciesLabels=foo[-i,speciesVar],distName="logis",doPlots=FALSE,effectLevel = HC.primary)
       result
  })
  LOO.results.norm <- clusterApplyLB(snowCluster,x=0:nrow(foo),fun=function(i){
       if(i==0)result <- genericPLL.AddOne(foo[,resVar],speciesLabels=foo[,speciesVar],distName="norm",doPlots=FALSE,effectLevel = HC.primary)
       if(i>=1)result <- genericPLL.AddOne(foo[,resVar][-i],speciesLabels=foo[-i,speciesVar],distName="norm",doPlots=FALSE,effectLevel = HC.primary)
       result
  })
  stopCluster(snowCluster)
  
    
    #results table
    do.call(rbind,lapply(LOO.results.logis,FUN=function(resList){
      resList[[1]]
    }))
    do.call(rbind,lapply(LOO.results.norm,FUN=function(resList){
      resList[[1]]
    }))
    
    #logistic LOO
    fit.out <- cbind(structure(data.frame(c("NONE",as.character(foo[,speciesVar])),c(NA,foo[,resVar])),names=c("Species.Out",names(foo)[2])),
                     do.call(rbind,lapply(LOO.results.logis,FUN=function(resList){
                       resList[[1]]
                     })))
    
    writeWorksheet(wb.out,data.frame(Logistic=c(tagString,"Leave One Out Analysis")),sheet=tagString,startRow=rowCount)
    rowCount <- rowCount+3
    writeWorksheet(wb.out,fit.out,sheet=tagString,startRow=rowCount,startCol=1)
    rowCount <- rowCount+nrow(fit.out)+3
    
    #normal LOO
    fit.out <- cbind(structure(data.frame(c("NONE",as.character(foo[,speciesVar])),c(NA,foo[,resVar])),names=c("Species.Out",names(foo)[2])),
                     do.call(rbind,lapply(LOO.results.norm,FUN=function(resList){
                       resList[[1]]
                     })))
    writeWorksheet(wb.out,data.frame(Normal=c(tagString,"Leave One Out Analysis")),sheet=tagString,startRow=rowCount)
    rowCount <- rowCount+3
    writeWorksheet(wb.out,fit.out,sheet=tagString,startRow=rowCount,startCol=1)
    rowCount <- rowCount+nrow(fit.out)+3
    
    xLimLower <- 10^floor(log10(min(c(
      exp(qlogis(.05,LOO.results.logis[[1]][[1]][5],LOO.results.logis[[1]][[1]][6]))/5,
      sort(unlist(lapply(LOO.results.logis,FUN=function(resList){
        linesMat <- resList[[2]]
        #print(head(linesMat[100:200,]))
        linesMat[round(linesMat[,1],3)==0.05,2]
      })))))))
    plot(x=foo[,resVar],y=(1:nrow(foo))/(nrow(foo)+1),log='x',xlim=c(xLimLower,10^(max(ceiling(log10(foo[,resVar]))))),ylim=c(0,1),axes=FALSE,
         xlab=concentrationLabel,ylab="Probability")
    majorTicks <- 10^seq(ceiling(par("usr")[1]),floor(par("usr")[2]))
    axis(side=1,at=majorTicks,label=majorTicks)
    axis(side=1,at=10^as.vector(sapply(ceiling(par("usr")[1]):(ceiling(par("usr")[2])-1),FUN=function(x)x+log10(1:9))),
         tck=-.01,labels=F)
    axis(side=2)
    box(lwd=2)
    lines(x=exp(qlogis(seq(.001,.999,length=200),LOO.results.logis[[1]][[1]][5],LOO.results.logis[[1]][[1]][6])),
          y=seq(.001,.999,length=200),col="black",lwd=3)
    lines(y=LOO.results.logis[[1]][[2]][,1],x=LOO.results.logis[[1]][[2]][,2],col="blue",lwd=3)
    lines(y=LOO.results.logis[[1]][[2]][,1],x=LOO.results.logis[[1]][[2]][,3],col="blue",lwd=3)
    
    invisible(lapply(LOO.results.logis[-1],FUN=function(resList){
      lines(x=exp(qlogis(seq(.001,.999,length=200),resList[[1]][5],resList[[1]][6])),
            y=seq(.001,.999,length=200),col=do.call(rgb,c(as.list(col2rgb("gray")[,1]),maxColorValue=255,alpha=.5*255)))
      lines(y=resList[[2]][,1],x=resList[[2]][,2],col=do.call(rgb,c(as.list(col2rgb("cyan")[,1]),maxColorValue=255,alpha=.5*255)))
      lines(y=resList[[2]][,1],x=resList[[2]][,3],col=do.call(rgb,c(as.list(col2rgb("cyan")[,1]),maxColorValue=255,alpha=.5*255)))
      invisible()
    }))
    abline(h=HC.primary,lty=3)
    mtext(side=1,outer=TRUE,text="Distribution: logis",line=-1,adj=0)
    hcX <- exp(qlogis(HC.primary,LOO.results.logis[[1]][[1]][5],LOO.results.logis[[1]][[1]][6]))
    mtext(side=3,at=hcX,text=bquote(HC[.(round(100*HC.primary))]),adj=0,las=2)
    mtext(side=3,at=hcX/3,text=bquote(HC[.(round(100*HC.primary))]/3),adj=0,las=2)
    mtext(side=3,at=hcX/5,text=bquote(HC[.(round(100*HC.primary))]/5),adj=0,las=2)
    abline(v=hcX/c(1,2,3,5),col="gray")
    
    
    xLimLower <- 10^floor(log10(min(c(
      exp(qlogis(.05,LOO.results.logis[[1]][[1]][5],LOO.results.logis[[1]][[1]][6]))/5,
      sort(unlist(lapply(LOO.results.logis,FUN=function(resList){
        linesMat <- resList[[2]]
        #print(head(linesMat[100:200,]))
        linesMat[round(linesMat[,1],3)==0.05,2]
      })))))))
    yLimUpper <- max(unlist(lapply(LOO.results.logis,FUN=function(resList){
      distPars <- resList[[1]][5:6]
      dlogis(distPars[1],distPars[1],distPars[2])
    })))
    plot(x=xLimLower,y=yLimUpper,log='x',xlim=c(xLimLower,10^(max(ceiling(log10(foo[,resVar]))))),ylim=c(0,yLimUpper),axes=FALSE,
         xlab=concentrationLabel,ylab="Density",type="n")
    majorTicks <- 10^seq(ceiling(par("usr")[1]),floor(par("usr")[2]))
    axis(side=1,at=majorTicks,label=majorTicks)
    axis(side=1,at=10^as.vector(sapply(ceiling(par("usr")[1]):(ceiling(par("usr")[2])-1),FUN=function(x)x+log10(1:9))),
         tck=-.01,labels=F)
    axis(side=2)
    box(lwd=2)
    qVals <- qlogis(seq(.001,.999,length=200),LOO.results.logis[[1]][[1]][5],LOO.results.logis[[1]][[1]][6])
    lines(y=dlogis(qVals,LOO.results.logis[[1]][[1]][5],LOO.results.logis[[1]][[1]][6]),
          x=exp(qVals),col="black",lwd=3)
    
    invisible(lapply(LOO.results.logis[-1],FUN=function(resList){
      modelPars <- resList[[1]][5:6]
      #print(modelPars)
      qVals <- qlogis(seq(.001,.999,length=200),modelPars[1],modelPars[2])
      #print(qVals)
      lines(y=dlogis(qVals,modelPars[1],modelPars[2]),
            x=exp(qVals),col=do.call(rgb,c(as.list(col2rgb("cyan")[,1]),maxColorValue=255,alpha=.5*255)))
      invisible()
    }))
    rug(side=1,x=foo[,resVar],lwd=1,col=rgb(0,0,0,.25))
    mtext(side=1,outer=TRUE,text="Distribution: logis",line=-1,adj=0)
    #hc5 <- exp(qlogis(.05,LOO.results.logis[[1]][[1]][5],LOO.results.logis[[1]][[1]][6]))
    #mtext(side=3,at=hc5,text=expression(HC[5]),adj=0,las=2)
    #mtext(side=3,at=hc5/3,text=expression(HC[5]/3),adj=0,las=2)
    #mtext(side=3,at=hc5/5,text=expression(HC[5]/5),adj=0,las=2)
    #abline(v=hc5/c(1,2,3,5),col="gray")
    hcX <- exp(qlogis(HC.primary,LOO.results.logis[[1]][[1]][5],LOO.results.logis[[1]][[1]][6]))
    mtext(side=3,at=hcX,text=bquote(HC[.(round(100*HC.primary))]),adj=0,las=2)
    mtext(side=3,at=hcX/3,text=bquote(HC[.(round(100*HC.primary))]/3),adj=0,las=2)
    mtext(side=3,at=hcX/5,text=bquote(HC[.(round(100*HC.primary))]/5),adj=0,las=2)
    abline(v=hcX/c(1,2,3,5),col="gray")
    
    
    pageBreak("B:  Normal Leave-One-Out")
    xLimLower <- 10^floor(log10(min(c(
      exp(qnorm(.05,LOO.results.norm[[1]][[1]][5],LOO.results.norm[[1]][[1]][6]))/5,
      sort(unlist(lapply(LOO.results.norm,FUN=function(resList){
        linesMat <- resList[[2]]
        #print(head(linesMat[100:200,]))
        linesMat[round(linesMat[,1],3)==0.05,2]
      })))))))
    plot(x=foo[,resVar],y=(1:nrow(foo))/(nrow(foo)+1),log='x',xlim=c(xLimLower,10^(max(ceiling(log10(foo[,resVar]))))),ylim=c(0,1),axes=FALSE,
         xlab=concentrationLabel,ylab="Probability")
    majorTicks <- 10^seq(ceiling(par("usr")[1]),floor(par("usr")[2]))
    axis(side=1,at=majorTicks,label=majorTicks)
    axis(side=1,at=10^as.vector(sapply(ceiling(par("usr")[1]):(ceiling(par("usr")[2])-1),FUN=function(x)x+log10(1:9))),
         tck=-.01,labels=F)
    axis(side=2)
    box(lwd=2)
    lines(x=exp(qnorm(seq(.001,.999,length=200),LOO.results.norm[[1]][[1]][5],LOO.results.norm[[1]][[1]][6])),
          y=seq(.001,.999,length=200),col="black",lwd=3)
    lines(y=LOO.results.norm[[1]][[2]][,1],x=LOO.results.norm[[1]][[2]][,2],col="blue",lwd=3)
    lines(y=LOO.results.norm[[1]][[2]][,1],x=LOO.results.norm[[1]][[2]][,3],col="blue",lwd=3)
    
    invisible(lapply(LOO.results.norm[-1],FUN=function(resList){
      lines(x=exp(qnorm(seq(.001,.999,length=200),resList[[1]][5],resList[[1]][6])),
            y=seq(.001,.999,length=200),col=do.call(rgb,c(as.list(col2rgb("gray")[,1]),maxColorValue=255,alpha=.5*255)))
      lines(y=resList[[2]][,1],x=resList[[2]][,2],col=do.call(rgb,c(as.list(col2rgb("cyan")[,1]),maxColorValue=255,alpha=.5*255)))
      lines(y=resList[[2]][,1],x=resList[[2]][,3],col=do.call(rgb,c(as.list(col2rgb("cyan")[,1]),maxColorValue=255,alpha=.5*255)))
      invisible()
    }))
    abline(h=HC.primary,lty=3)
    mtext(side=1,outer=TRUE,text="Distribution: norm",line=-1,adj=0)
    #hc5 <- exp(qnorm(.05,LOO.results.norm[[1]][[1]][5],LOO.results.norm[[1]][[1]][6]))
    #mtext(side=3,at=hc5,text=expression(HC[5]),adj=0,las=2)
    #mtext(side=3,at=hc5/3,text=expression(HC[5]/3),adj=0,las=2)
    #mtext(side=3,at=hc5/5,text=expression(HC[5]/5),adj=0,las=2)
    #abline(v=hc5/c(1,2,3,5),col="gray")
    hcX <- exp(qnorm(HC.primary,LOO.results.norm[[1]][[1]][5],LOO.results.norm[[1]][[1]][6]))
    mtext(side=3,at=hcX,text=bquote(HC[.(round(100*HC.primary))]),adj=0,las=2)
    mtext(side=3,at=hcX/3,text=bquote(HC[.(round(100*HC.primary))]/3),adj=0,las=2)
    mtext(side=3,at=hcX/5,text=bquote(HC[.(round(100*HC.primary))]/5),adj=0,las=2)
    abline(v=hcX/c(1,2,3,5),col="gray")
    
    
    xLimLower <- 10^floor(log10(min(c(
      exp(qnorm(.05,LOO.results.norm[[1]][[1]][5],LOO.results.norm[[1]][[1]][6]))/5,
      sort(unlist(lapply(LOO.results.norm,FUN=function(resList){
        linesMat <- resList[[2]]
        #print(head(linesMat[100:200,]))
        linesMat[round(linesMat[,1],3)==0.05,2]
      })))))))
    yLimUpper <- max(unlist(lapply(LOO.results.norm,FUN=function(resList){
      distPars <- resList[[1]][5:6]
      dnorm(distPars[1],distPars[1],distPars[2])
    })))
    plot(x=xLimLower,y=yLimUpper,log='x',xlim=c(xLimLower,10^(max(ceiling(log10(foo[,resVar]))))),ylim=c(0,yLimUpper),axes=FALSE,
         xlab=concentrationLabel,ylab="Density",type="n")
    majorTicks <- 10^seq(ceiling(par("usr")[1]),floor(par("usr")[2]))
    axis(side=1,at=majorTicks,label=majorTicks)
    axis(side=1,at=10^as.vector(sapply(ceiling(par("usr")[1]):(ceiling(par("usr")[2])-1),FUN=function(x)x+log10(1:9))),
         tck=-.01,labels=F)
    axis(side=2)
    box(lwd=2)
    qVals <- qnorm(seq(.001,.999,length=200),LOO.results.norm[[1]][[1]][5],LOO.results.norm[[1]][[1]][6])
    lines(y=dnorm(qVals,LOO.results.norm[[1]][[1]][5],LOO.results.norm[[1]][[1]][6]),
          x=exp(qVals),col="black",lwd=3)
    
    invisible(lapply(LOO.results.norm[-1],FUN=function(resList){
      modelPars <- resList[[1]][5:6]
      #print(modelPars)
      qVals <- qnorm(seq(.001,.999,length=200),modelPars[1],modelPars[2])
      #print(qVals)
      lines(y=dnorm(qVals,modelPars[1],modelPars[2]),
            x=exp(qVals),col=do.call(rgb,c(as.list(col2rgb("cyan")[,1]),maxColorValue=255,alpha=.5*255)))
      invisible()
    }))
    rug(side=1,x=foo[,resVar],lwd=1,col=rgb(0,0,0,.25))
    mtext(side=1,outer=TRUE,text="Distribution: norm",line=-1,adj=0)
    #hc5 <- exp(qnorm(.05,LOO.results.norm[[1]][[1]][5],LOO.results.norm[[1]][[1]][6]))
    #mtext(side=3,at=hc5,text=expression(HC[5]),adj=0,las=2)
    #mtext(side=3,at=hc5/3,text=expression(HC[5]/3),adj=0,las=2)
    #mtext(side=3,at=hc5/5,text=expression(HC[5]/5),adj=0,las=2)
    #abline(v=hc5/c(1,2,3,5),col="gray")
    hcX <- exp(qnorm(HC.primary,LOO.results.norm[[1]][[1]][5],LOO.results.norm[[1]][[1]][6]))
    mtext(side=3,at=hcX,text=bquote(HC[.(round(100*HC.primary))]),adj=0,las=2)
    mtext(side=3,at=hcX/3,text=bquote(HC[.(round(100*HC.primary))]/3),adj=0,las=2)
    mtext(side=3,at=hcX/5,text=bquote(HC[.(round(100*HC.primary))]/5),adj=0,las=2)
    abline(v=hcX/c(1,2,3,5),col="gray")
  }
    
    
    
    
    
    
    
  if(doAddOneIn){
    #add-one-in is much more difficult to program
    #add-one-in analyses -- plots directly from this
    pageBreak("Add-One-In Analysis")
    pageBreak("A:  Logistic Add-One-In")
    fit.out1 <- try(addOnePLL(effectValues = foo[,resVar],speciesLabels=foo[,speciesVar],distName="logis",xlabString=concentrationLabel,dataTag=tagString,effectLevel = HC.primary))
    if(class(fit.out1)==Error){
    pageBreak("B:  Normal Add-One-In")
    fit.out2 <- addOnePLL(effectValues = foo[,resVar],speciesLabels=foo[,speciesVar],distName="norm",xlabString=concentrationLabel,dataTag=tagString,effectLevel = HC.primary)
    fit.out <- rbind(fit.out1,fit.out2)
    print(fit.out)
    }
    writeWorksheet(wb.out,data.frame(AddOne=c(tagString,"Add One In Analysis")),sheet=tagString,startRow=rowCount)
    rowCount <- rowCount+2
    writeWorksheet(wb.out,fit.out[,-1],sheet=tagString,startRow=rowCount,startCol=1)
    setColumnWidth(wb.out, sheet = tagString, column = 1:11, width = -1)
    
    #rgl.quit()#clean up the rgl stuff
  }
  

  dev.off()#close the pdf output
}

saveWorkbook(wb.out)
