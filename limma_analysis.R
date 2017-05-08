contrasts_matrix <- function(filePath, basePlotPath, outputPath,factors, plotFactors, factorColumnValues, factorColumnValueSet, samples, sampleNames, individuals, contrasts,  factorSetMap, projectId, userId, analysisDetailId, startingColumn, commaOrTabDelimitedFlag, allContrasts) {

  # GIVEN: The counts matrix 
  #        basePath
  #        filePath
  #        basePlotPath
  #        factors
  #        samples
  #        sampleNames
  #        individuals
  #        listPlotFactors
  #        contrasts
  # 
  # RETURNS: 

  require(limma)

  require("reshape2")
  require("ggplot2")
  require("gplots")
  require("plyr")

  library(Biobase)
  library(statmod)
  library(affy)
  library(edgeR)
  
  if (commaOrTabDelimitedFlag)  
  {
    data = read.delim2(filePath,as.is=T,check.names=F,sep=",")
  }
  else{
    data = read.delim2(filePath,as.is=T,check.names=F)
  }
  countTable=data[,-seq(1:startingColumn)]
  genenames=data$Geneid
  design_limma <-model.matrix(~0+factors)
  
  #plotMDS(y,col=ifelse(design$Treatment=="P","red","blue”))
  #nc = cpm(counts)

  colnames(design_limma)<-levels(factor(factors))  
  contrasts_matrix <- makeContrasts(contrasts=allContrasts, levels=design_limma)

  sedIdx = match(samples,colnames(countTable))
  print (sedIdx)

  countTable=countTable[,sedIdx]
  
  colnames(countTable) = sampleNames  
  
  counts = DGEList(counts=countTable,genes=genenames)
  isexpr<-rowSums(cpm(counts)>1)>=3
  counts <- counts[isexpr,keep.lib.sizes=FALSE]
  counts <-calcNormFactors(counts)  
  print (" before normalized")
  normalizedCountsMatrix = cpm(counts, log=TRUE)
  
  rownames(normalizedCountsMatrix) <- counts$genes[,1]
  
  normalizedCountsMatrix <- cbind(genes = rownames(normalizedCountsMatrix), normalizedCountsMatrix)
  
  write.csv(normalizedCountsMatrix,file=paste(outputPath, "/NormalizedData", ".csv", sep=""))

  y <- voom(counts,design_limma)
  
  counter = 1
  for (plotFactor in plotFactors) {  
  
    png(paste(basePlotPath, "/plotMDS_", plotFactor, "_Project_", projectId, "_user_", userId , "_Analysis_", analysisDetailId, ".png", sep=""), width=6, height=4, units = "in", res=300)  
    
     # title  main = paste( "plotMDS_", plotFactor, sep="")    
    
    plotMDS(y,col=ifelse(eval(parse(text=factorColumnValues [plotFactor])) == factorColumnValueSet[counter],"blue","red"), xlab = NA, ylab = NA , main=paste( "plotMDS_", plotFactor, sep=""))

    legend("topleft", legend=eval(parse(text=(factorSetMap[plotFactor]))), col=c("red", "blue"), pch=15)

    dev.off()    
    
    counter = counter + 1

  }
  
  #MA Plot
  png(paste(basePlotPath, "/plotMA", "_Project_", projectId, "_user_", userId , "_Analysis_", analysisDetailId,".png", sep=""), width=4, height=4, units="in", res=300)
  limma::plotMA(y, main=paste( "plotMA", sep=""))
  dev.off()

  if (length(individuals) > 0){
   corfit <- duplicateCorrelation(y,design_limma,block=individuals)
   v <- voom(counts,design_limma,plot=FALSE,block=samples,correlation=corfit$consensus)
   fit <- lmFit(v,design_limma)
  }
  else{
   fit <- lmFit(y,design_limma)
  }

  fit2<-contrasts.fit(fit,contrasts_matrix)
  fit2 <- eBayes(fit2)

  #results <-decideTests(fit2,method="separate",adjust.method="none",p.value = as.numeric(decidetestsPValue),lfc = as.numeric(decidetestsLFC) )
  #write.csv(results,file = paste(outputPath,"/","decideTests.csv", sep=""))

  fitDataOutputPath = paste(outputPath, "/fitCoefficients.csv", sep="")
  
  fitData = data.frame(fit2$genes, fit2$coefficients)  

  write.csv(fitData, file = fitDataOutputPath)
  
  #print (fit2$coefficients)

  for (contrast in contrasts){
  
    res=topTable(fit2,coef=contrast,number=length(which(isexpr==TRUE)),sort="none")
    
    write.csv(res,file=paste(outputPath, "/", contrast, ".csv", sep=""))

    #print (fit2$coefficients[contrast])
    
    #dfContrast = data.frame(res, fit2$coefficients)  
    
    #write.csv(dfContrast,file=paste(outputPath, "/", contrast, ".csv", sep=""))
  
  }    

  return(contrasts_matrix)

}
