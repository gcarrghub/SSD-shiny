library(openxlsx)
library(reshape2)
#library(shinyIncubator)
#library(tcltk)
library(dplyr)

#tablerows <- 0

options(shiny.maxRequestSize=11*1024^2)

gm_mean = function(x, na.rm=TRUE){
     exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

shinyServer(function(input, output, session) {
     output$activeTab <- reactive({
          return(input$tab)
     })
     outputOptions(output, 'activeTab', suspendWhenHidden=FALSE)
     
     getsheets <- function(){
          
          inFile <- input$file1
          
          if (is.null(inFile))
               return(NULL)
          if(substr(inFile$name,nchar(inFile$name)-3,nchar(inFile$name))=="xlsx" | 
                  substr(inFile$name,nchar(inFile$name)-2,nchar(inFile$name))=="xls"){
               dataWB <- inFile$datapath
               wb.data <- openxlsx::loadWorkbook(file = dataWB)
               sheetNames <- openxlsx::getSheetNames(file = dataWB)
          }
          
          if(substr(inFile$name,nchar(inFile$name)-2,nchar(inFile$name))=="txt"){
               sheetNames <- inFile$name
               
          }
          
          return (sheetNames)
     }
     
     getvars <- function(type){
          
          inFile <- input$file1
          
          if (is.null(inFile))
               return(NULL)
          
          if(substr(inFile$name,nchar(inFile$name)-3,nchar(inFile$name))=="xlsx" | 
                  substr(inFile$name,nchar(inFile$name)-2,nchar(inFile$name))=="xls"){
               dataWB <- inFile$datapath
               #wb.data <- loadWorkbook(dataWB)
               #sheetNames <- getSheetNames(wb.data)
               foo <- read.xlsx(dataWB,sheet = input$sheets)
          }
          
          if(substr(inFile$name,nchar(inFile$name)-2,nchar(inFile$name))=="txt"){
               dataTXT <- inFile$datapath
               foo <- read.delim(file=dataTXT, stringsAsFactors=FALSE)
          }
          
          vartypes <- sapply(foo, class)
          
          vartypetotal <- dcast(as.data.frame(table(vartypes)), Freq ~ vartypes)
          
          vars <- names(vartypes[vartypes==type])
          
          return (vars)
     }
     
     loadtable <- reactive({
          
          withProgress({
               setProgress(message = "Loading Data")
               inFile <- input$file1
               
               if (is.null(inFile) | is.null(input$sheets))
                    return(NULL)
               
               if(substr(inFile$name,nchar(inFile$name)-3,nchar(inFile$name))=="xlsx" | 
                  substr(inFile$name,nchar(inFile$name)-2,nchar(inFile$name))=="xls"){
                    dataWB <- inFile$datapath
                    #wb.data <- loadWorkbook(dataWB)
                    df <- openxlsx::read.xlsx(dataWB,sheet=input$sheets)
               }
               
               if(substr(inFile$name,nchar(inFile$name)-2,nchar(inFile$name))=="txt"){
                    dataTXT <- inFile$datapath
                    df <- read.delim(file=dataTXT, stringsAsFactors=FALSE)
               }
               
               tablerows <- nrow(df)
               if(input$calcgeommeans){
                    req(input$groupvar)
                    if((length(input$groupvar) == 1 & input$groupvar[1] != input$charvar) | length(input$groupvar) > 1){
                         groupvars <- c(input$charvar, input$groupvar[which(input$groupvar != input$charvar)])
                         
                         names(df)[which(names(df) == input$numvar)] <- "numvar"
                         
                         groupvars.sym <- lapply(groupvars, as.symbol)
                         df <- df %>% 
                              group_by_(.dots = groupvars.sym) %>% 
                              summarise(geommean = gm_mean(numvar))
                         
                         names(df)[which(names(df) == "geommean")] <- input$numvar
                         
                         
                    } else {
                         names(df)[which(names(df) == input$charvar)] <- "charvar"
                         names(df)[which(names(df) == input$numvar)] <- "numvar"
                         
                         df <- df %>% group_by(charvar) %>% summarise(geommean = gm_mean(numvar))
                         
                         names(df)[which(names(df) == "geommean")] <- input$numvar
                         names(df)[which(names(df) == "charvar")] <- input$charvar
                    }
                    
               }
               #print("Load data completed")
               #print(df)
               return(df)
          })

          
          
     })
     
     analysis <- function() {
     	if(input$ssdttc=="TTC"){
     		doLeaveOneOut <- FALSE
     	} else {
     		doLeaveOneOut <- input$leaveone
     	}
     	#add-one-in is not as bad.
     	if(input$ssdttc=="TTC"){
     		doAddOneIn <- FALSE
     	} else {
     		doAddOneIn <- input$addone
     	}
     	
     	HC.primary <- input$percent

          pageBreak <- function(messageString){
               plot(x=0,y=0,type="n",axes=FALSE,xlab="",ylab="")
               text(x=0,y=0,labels=messageString,adj=0.5,cex=3)  
          }
          

          source("genericPLL.R")
          source("addOne.genericPLL.R")
          source("addOnePLL.V3.R")
          source("rriskFitdist.GJC.R")

           
          setProgress(value = 1, detail="10% Complete")
          
          foo <- loadtable()
          
          
          #tagString <- substring(inFile$name,first=1,last=nchar(inFile$name)-5)
          pdfString <- "tempout"
          
          if(!is.na(as.numeric(substr(isolate(input$sheets),nchar(isolate(input$sheets)),nchar(isolate(input$sheets)))))){
               sheetname <- paste0(isolate(input$sheets)," Output")
          } else {
               sheetname <- isolate(input$sheets)
          }
          
          tagString <- sheetname
          
          resVar <- isolate(input$numvar)
          masterWB <- "temp/tempout.xlsx"
          concentrationLabel <- paste0("Concentration (", input$units, ")")
          
          if(file.exists(masterWB))unlink(masterWB)
          
          wb.out <- createWorkbook()
          
          
          colnames(foo)[which(names(foo) == paste0(isolate(input$charvar)))] <- "Species"
          speciesVar <- "Species"
          
          #str(foo)
          
          foo <- foo[which(!is.na(foo[[input$numvar]])),]

					foo <- foo[order(foo[,resVar]),]
          
          # badRows <- which(rowSums(is.na(foo))>0)
          # 
          # if(length(badRows)>0){
          #      #cat("\nRow(s)",badRows,"will be dropped from the file.  Please inspect the data.\n")
          #      foo <- na.omit(foo)
          # }
          
          
          
          #str(foo)
          #return(as.data.frame(badRows))
          addWorksheet(wb.out,sheetName = sheetname)
          
          setProgress(value = 2, detail="20% Complete")
          printSpeciesLabels <- input$listinmargins
          ital <- ifelse(input$ital, 3, 2)
                    
          pdf(file=paste("temp/",pdfString,".pdf",sep=""),height=8.5,width=11)
          #full data fits -- plots directly from this
          #if both 50 and 95th are estimated, be sure to add lines to space out results
          #from subsequent tables
          pageBreak("Data Fits")
          #full data fits -- plots directly from this
          pageBreak("A:  Logistic Data Fits")
          fit.out1 <- rbind(
               #genericPLL(foo[,resVar],speciesLabels=foo[,speciesVar],distName="logis",doPlots=TRUE,xlabString=concentrationLabel),
               genericPLL(foo[,resVar],speciesLabels=foo[,speciesVar],
                          distName="logis",doPlots=TRUE,xlabString=concentrationLabel,
                          effectLevel = HC.primary, showUpper=input$ucl, 
                          printSpeciesLabels=printSpeciesLabels,
                          italicFont=ital,
                          cex.in = input$textSize)
          )
          pageBreak("B:  Normal Data Fits")
          
          setProgress(value = 3, detail="30% Complete")
          
          fit.out2 <- rbind(
               #genericPLL(foo[,resVar],speciesLabels=foo[,speciesVar],distName="norm",doPlots=TRUE,xlabString=concentrationLabel),
               genericPLL(foo[,resVar],speciesLabels=foo[,speciesVar],
                          distName="norm",doPlots=TRUE,xlabString=concentrationLabel,
                          effectLevel = HC.primary, showUpper=input$ucl, printSpeciesLabels=printSpeciesLabels,
                          italicFont=ital,
                          cex.in = input$textSize)
          )
          
          fit.out <- rbind(fit.out1,fit.out2)
          fit.out2 <- data.frame(
               Distribution=rep(c("Logistic","Normal"),times=c(nrow(fit.out1),nrow(fit.out2))),
               fit.out)
          
          ###################################
          pageBreak("C: Non-Parametric Data Fits")
           normParms <- c(mean(log10(foo[,resVar])),sd(log10(foo[,resVar])))
           xRange <- 10^c(do.call(qnorm,as.list(c(0.001,normParms))),do.call(qnorm,as.list(c(0.999,normParms))))
           plot(x=10^seq(log10(xRange[1]),log10(xRange[2]),length=1000),y=pnorm(seq(log10(xRange[1]),log10(xRange[2]),length=1000),mean=normParms[1],sd=normParms[2]),
                type="l",lwd=3,col="magenta",log="x", ylab="Probability", xlab = concentrationLabel)
           
           HCx <- 10^qnorm(HC.primary,mean=normParms[1],sd=normParms[2])
                
           rug(side=1,col="magenta",lwd=2,x=HCx)               
           
           qvals <- quantile(log10(foo[,resVar]),prob=seq(0,1,by=0.001),type=8)
           #invert quantile function
           pDat <- sapply(log10(foo[,resVar]),FUN = function(x){
                uniroot(f = function(px)quantile(log10(foo[,resVar]),prob=px,type=8)-x,interval = c(0,1))$root
           })
           points(x=foo[,resVar],y=pDat,col="blue",pch=16)
           lines(x=10^qvals,y=seq(0,1,by=0.001),col="cyan",lwd=2)
           rug(side=1,col="cyan",lwd=2,x=10^quantile(log10(foo[,resVar]),prob=0.05,type=8))
           
          fit.out.parametric <- fit.out2[1,]
          fit.out.parametric[1,] <- NA
          fit.out.parametric$Distribution <- "Non-Parametric"
          fit.out.parametric$X <- HC.primary
          fit.out.parametric[1,3] <- HCx
          
          fit.out2 <- rbind(fit.out2, fit.out.parametric)
          ###################################
          
          
          head.out <- data.frame(spreadsheet=tagString)
       
          writeData(wb.out,sheet = tagString,x = data.frame(
            spreadsheet=c(tagString, 
                          paste0(input$ssdttc, 
                                 " Analysis of Full Data (", 
                                 nrow(loadtable()), 
                                 " points)"))))
          footoshow <- foo
          if(input$ssdttc=="TTC"){
          	names(footoshow)[1] <- "Chemical"
          }
          writeData(wb.out,x = footoshow,sheet=tagString,startRow=1,startCol=10)
          #writeWorksheet(wb.out,footoshow,sheet=tagString,startRow=1,startCol=10)
          writeData(wb.out,sheet=tagString,startRow=4,x = fit.out2)
          #writeWorksheet(wb.out,fit.out2,sheet=tagString,startRow=4)
          rowCount <- nrow(fit.out2)+9
          
          setProgress(value = 4, detail="40% Complete")
          
          #dev.off()
          #saveWorkbook(wb.out)
          
          #leave one out is easy, but computationally expensive
          if(doLeaveOneOut){
               pageBreak("Leave-One-Out Analysis")
               library(parallel)
               snowCluster <- makePSOCKcluster(8)
               
               #leave-one-out analyses
               clusterExport(snowCluster,varlist=c("genericPLL.AddOne", "rriskFitdist.GJC"))
               pageBreak("A:  Logistic Leave-One-Out")
               
               setProgress(value = 5, detail="50% Complete")
               
               LOO.results.logis <- clusterApplyLB(snowCluster,x=0:nrow(foo),fun=function(i){
                    if(i==0)result <- genericPLL.AddOne(foo[,resVar],speciesLabels=foo[,speciesVar],distName="logis",doPlots=FALSE,effectLevel = HC.primary)
                    if(i>=1)result <- genericPLL.AddOne(foo[,resVar][-i],speciesLabels=foo[-i,speciesVar],distName="logis",doPlots=FALSE,effectLevel = HC.primary)
                    result
               })
               
               setProgress(value = 6, detail="60% Complete")
               
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
               
               writeData(wb.out,
                         x = data.frame(Logistic=c(tagString,"Leave One Out Analysis")),
                         sheet=tagString,startRow=rowCount)
               #writeWorksheet(wb.out,data.frame(Logistic=c(tagString,"Leave One Out Analysis")),sheet=tagString,startRow=rowCount)
               rowCount <- rowCount+3
               writeData(wb.out,x = fit.out,
                         sheet=tagString,
                         startRow=rowCount,
                         startCol=1)
               #writeWorksheet(wb.out,fit.out,sheet=tagString,startRow=rowCount,startCol=1)
               rowCount <- rowCount+nrow(fit.out)+3
               
               setProgress(value = 7, detail="70% Complete")
               
               #normal LOO
               fit.out <- cbind(structure(data.frame(c("NONE",as.character(foo[,speciesVar])),c(NA,foo[,resVar])),names=c("Species.Out",names(foo)[2])),
                                do.call(rbind,lapply(LOO.results.norm,FUN=function(resList){
                                     resList[[1]]
                                })))
               writeData(wb.out,
                        x=data.frame(Normal=c(tagString,"Leave One Out Analysis")),
                         sheet=tagString,
                         startRow=rowCount)
               #writeWorksheet(wb.out,data.frame(Normal=c(tagString,"Leave One Out Analysis")),sheet=tagString,startRow=rowCount)
               rowCount <- rowCount+3
               writeData(wb.out,
                         x=fit.out,
                         sheet=tagString,
                         startRow=rowCount,
                         startCol=1)
               #writeWorksheet(wb.out,fit.out,sheet=tagString,startRow=rowCount,startCol=1)
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
               
               setProgress(value = 8, detail="80% Complete")
               
               
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
          
          
          setProgress(value = 9, detail="90% Complete")
          
          
          
          
          if(doAddOneIn){
               #add-one-in is much more difficult to program
               #add-one-in analyses -- plots directly from this
               pageBreak("Add-One-In Analysis")
               pageBreak("A:  Logistic Add-One-In")
               fit.out1 <- try(addOnePLL(effectValues = foo[,resVar],speciesLabels=foo[,speciesVar],distName="logis",xlabString=concentrationLabel,dataTag=tagString,effectLevel = HC.primary))
               pageBreak("B:  Normal Add-One-In")
               fit.out2 <- try(addOnePLL(effectValues = foo[,resVar],speciesLabels=foo[,speciesVar],distName="norm",xlabString=concentrationLabel,dataTag=tagString,effectLevel = HC.primary))
               fit.out <- rbind(fit.out1,fit.out2)
               print(fit.out)
               writeData(wb.out,x = data.frame(AddOne=c(tagString,"Add One In Analysis")),
                         sheet=tagString,startRow=rowCount)
               #writeWorksheet(wb.out,data.frame(AddOne=c(tagString,"Add One In Analysis")),sheet=tagString,startRow=rowCount)
               rowCount <- rowCount+2
               writeData(wb.out,x = fit.out[,-1],
                         sheet=tagString,startRow=rowCount,startCol=1)
               #writeWorksheet(wb.out,fit.out[,-1],sheet=tagString,startRow=rowCount,startCol=1)
               
               #rgl.quit()#clean up the rgl stuff
          }
          
          
          setProgress(value = 10, detail="100% Complete")
          
          
          dev.off()#close the pdf output
          #setColumnWidth(wb.out, sheet = tagString, column = 1:11, width = -1)
          saveWorkbook(wb.out,file=masterWB,overwrite = TRUE)
          
          
          
          return()
          
     }
     
     
     output$sheets<-renderUI({
          req(input$file1)
          selectInput("sheets", "2.) Select the Sheet to Analyze", choices=getsheets())
     })
     
     output$charvars<-renderUI({
          req(input$sheets)
     	if(input$ssdttc=='SSD'){
     		label1 <- "3.) Select the Species Name Variable"
     	} else if(input$ssdttc=='TTC'){
     		label1 <- "3.) Select the Chemical Name Variable"
     	}
          tagList(selectInput("charvar", label = label1, choices=getvars('character')))
     	
     })
     
     output$geommeansUI <- renderUI({
          req(input$sheets)
          if(input$calcgeommeans){
               return(selectInput("groupvar", "5.) Select the Variable to Group By", choices=getvars('character'), multiple = T))
          }
     })
     
     output$numvars<-renderUI({
          req(input$sheets)
          selectInput("numvar", "4.) Select the Variable to Analyze", choices=getvars('numeric'))
     })
     
     output$slider<-renderUI({
          req(input$sheets)
          sliderInput("percent", "Choose Percentage of Interest:", min=.05, max=0.5, step=0.05, value=0.05)
     })
     
     output$uclbox<-renderUI({
     	checkboxInput("ucl", "Include UCL on Data Fit Plots", value=FALSE)
     })
     
     output$listinmarginsbox<-renderUI({
     	if(input$ssdttc=='SSD'){
     		label <- "Include list of Species in Margins of Plot"
     		defaultvalue <- TRUE
     		checkboxInput("listinmargins", label = label, value=defaultvalue)
     	} else if(input$ssdttc=='TTC'){
     		label <- "Include list of Chemicals in Margins of Plot"
     		defaultvalue <- FALSE
     	}
          tagList(
               checkboxInput("listinmargins", label = label, value=defaultvalue),
               checkboxInput("ital", "Italicize Names?", value = defaultvalue)
          )
     	
     })
     
     output$leaveonebox<-renderUI({
          inFile <- input$file1
          
          if (is.null(inFile) | is.null(input$sheets) | input$ssdttc=='TTC')
               return(NULL)
          
          if(substr(inFile$name,nchar(inFile$name)-3,nchar(inFile$name))=="xlsx" | 
                  substr(inFile$name,nchar(inFile$name)-2,nchar(inFile$name))=="xls"){
               dataWB <- inFile$datapath
               wb.data <- loadWorkbook(dataWB)
               foo <- read.xlsx(dataWB,input$sheets)
          }
          
          if(substr(inFile$name,nchar(inFile$name)-2,nchar(inFile$name))=="txt"){
               dataTXT <- inFile$datapath
               foo <- read.delim(file=dataTXT, stringsAsFactors=FALSE)
          }
          
          tablerows <- nrow(foo)
          if (tablerows > 30){
               checkboxInput("leaveone", label= "Include Leave-One-Out Analysis (Not suggested when analyzing > 30 rows)", value=FALSE)
          } else {
               checkboxInput("leaveone", label= "Include Leave-One-Out Analysis", value=TRUE)
          }
          
     })
     
     output$addonebox<-renderUI({
          inFile <- input$file1
          
          if (is.null(inFile) | is.null(input$sheets) | input$ssdttc=='TTC')
               return(NULL)
          
          if(substr(inFile$name,nchar(inFile$name)-3,nchar(inFile$name))=="xlsx" | 
                  substr(inFile$name,nchar(inFile$name)-2,nchar(inFile$name))=="xls"){
               dataWB <- inFile$datapath
               wb.data <- loadWorkbook(dataWB)
               foo <- read.xlsx(dataWB,input$sheets)
          }
          
          if(substr(inFile$name,nchar(inFile$name)-2,nchar(inFile$name))=="txt"){
               dataTXT <- inFile$datapath
               foo <- read.delim(file=dataTXT, stringsAsFactors=FALSE)
          }
          
          tablerows <- nrow(foo)
          if (tablerows > 30){
               checkboxInput("addone", label= "Include Add-One-In Analysis (Not suggested when analyzing > 30 rows)", value=FALSE)
          } else {
               checkboxInput("addone", label= "Include Add-One-In Analysis", value=TRUE)
          }
          
     })
     
     output$nrowmessage <- renderUI({
          req(loadtable())
          p("There are currently ", strong(nrow(loadtable())), " rows in this dataset.")
     })
     
     output$FullTable <- 
       renderTable(width = "auto",spacing = "xs",
                   expr = {
                     df <- loadtable()
                     if(input$ssdttc =="TTC"){
                       names(df)[1] <- "Chemical"
                       print(names(df))
                     }
                     #print("df in FullTable")
                     #print(str(df))
                     #print(df)
                     return(df)
                   })
     
     output$analysisinfo2 <- renderTable({
          #incubator implementation
          withProgress(session, min=0, max=10, {
               setProgress(message="Analysis in progress.")
               return(analysis())
          })
          
     })
          
     output$analysisinfo <- renderImage({
    
          
          #if (input$run==0)
          #     return(NULL)
          input$run
          isolate(withProgress(session, min=0, max=10, {

               setProgress(message="Analysis in progress.")
               analysis()

          }
          ))
          filename <- "temp/displayplot.png"
          
          return(list(src=filename,
          	contentType = "image/png",
          	alt = "plot"))
          
     }, delete=FALSE)
     
     output$QSimage <- renderImage({ 
          filename <- "../images/QSInformaticsInnovationSuite21.png"
          list(src=filename)
     }, deleteFile=FALSE)
     
     output$downloadExcel <- downloadHandler(
          
          filename = function() { 

               fname <- paste0(input$ssdttc, " Results ",input$sheets,".xlsx")
               return(fname) },
          
          content = function(file){
               wb <- loadWorkbook("temp/tempout.xlsx")
               #wb$writeWorksheet(filedata(), sheet=1, startRow=1, startCol=1, header=TRUE)
               randomNumber <- trunc(runif(1, 10000000, 99999999))
               outPath <- paste0("temp/output", randomNumber, ".xls")
               saveWorkbook(wb, outPath)
               
               file.copy(outPath, file)
          }
     )
     
     output$downloadPDF <- downloadHandler(
          
          filename = function() { 
               
               fname <- paste0( input$ssdttc, " Results ",input$sheets,".pdf")
               return(fname) },
          
          content = function(file){
               inFile <- isolate(input$file1)

               #tagString <- substring(inFile$name,first=1,last=nchar(inFile$name)-5)
               tagString <- "tempout"
               randomNumber <- trunc(runif(1, 10000000, 99999999))
               pdfname <- paste0("temp/",tagString,".pdf")
               outPath <- paste0("temp/output", randomNumber, ".pdf")
               file.copy(pdfname, outPath)
               
               file.copy(outPath, file)
          }
     )
     
     output$Excelbutton <- renderUI({
          if (input$run==0)
               return(NULL)
          return(downloadButton('downloadExcel', 'Download Excel File'))
            
     })
     
     output$PDFbutton <- renderUI({
          if (input$run==0)
               return(NULL)
          return(downloadButton('downloadPDF', 'Download PDF File'))
          
     })
     
     output$PlotImage <- renderUI({
          if (input$run==0){
          	return(p("Select a file and click the 'Run Analysis' button to see results."))
          }
               
          return(imageOutput("analysisinfo", width="600px", height="400px"))
          
     })
     
     output$analysisinfo2 <- renderTable({
          #incubator implementation          
          
          if (input$run==0)
               return(NULL)
          input$run
          isolate(withProgress(session, min=0, max=10, {
               
               setProgress(message="Analysis in progress.")
               return(analysis())
               
          }
          ))

          
     })
     
     output$PlotImage2 <- renderUI({
          if (input$run==0)
               return(p("Select a file and click the 'Run Analysis' button to see results."))
          return(tableOutput("analysisinfo2"))
          
     })
     
     output$testing <- renderTable({
          inFile <- input$file1
          
          if (is.null(inFile))
               return(NULL)
          
          if(substr(inFile$name,nchar(inFile$name)-3,nchar(inFile$name))=="xlsx" | 
                  substr(inFile$name,nchar(inFile$name)-2,nchar(inFile$name))=="xls"){
               dataWB <- inFile$datapath
               sheetNames <- getSheetNames(dataWB)
               return (as.data.frame(sheetNames))
          }
          
          if(substr(inFile$name,nchar(inFile$name)-2,nchar(inFile$name))=="txt"){
               txt <- read.delim(file="AllSpeciesC12-C13.txt")
               datatxt <- inFile$datapath
             
               txt.data <- read.delim(file=datatext)
               
               return (as.data.frame(txt.data))
          }
                   
     })
     
     observe({
          if(input$run!=0)
               updateTabsetPanel(session,"inTabset",selected="Output")
     })
     

})