import rpy2.robjects as robjects
import pandas as pd
import numpy

import rpy2.robjects.numpy2ri as rpyn

#robjects.r('require(DESeq2)')
#robjects.r('require(reshape2)')
#robjects.r('require(ggplot2)')
#robjects.r('require(gplots)')
#robjects.r('require(plyr)')
#robjects.r('require(genefilter)')
#robjects.r('require(edgeR)')
	
#robjects.r('data <- read.delim2("/Users/mitras/self/teleshova/counts.matrix",as.is=T,check.names=F)')
#robjects.r('col<- sub("_S\\\d+_R1_001", "", colnames(data),perl=TRUE)')
#robjects.r('col.n <- sub("\\\.bam", "",col,perl=T)')
#robjects.r('colnames(data)=col.n')

#robjects.r('genenames<-data$Geneid')
#robjects.r('mart<-useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset="mmulatta_gene_ensembl",host="www.ensembl.org")')
#robjects.r('genemap<-getBM(attributes=c("external_gene_name","ensembl_gene_id"), filters = "ensembl_gene_id",values=genenames,mart=mart)')
#robjects.r('idx <- match(genenames, genemap$ensembl_gene_id )')
#robjects.r('genesymbol <- genemap$external_gene_name[idx]')

robjects.r('design <- read.csv(file = "/Users/mitras/self/teleshova/Teleshova_Phenotype_csv.csv", header=TRUE, sep=",")')
robjects.r('rownames(design) <- design$Sample')

df = pd.DataFrame.from_csv("/Users/mitras/self/teleshova/Teleshova_Phenotype_csv.csv")

groupCol = "Group"
conditionCol = "Condition"

group = df[groupCol]
condition = df[conditionCol]

groupset = set(group)
conditionset = set(condition)

print str(groupset)
print str(conditionset)

x = list(itertools.product(*c))

xl = [".".join(x) for x in combinationStringList]

print xl

#['Placebo.baseline-Placebo.postgel', 'Placebo.baseline-lowgel.postgel', 'lowgel.baseline-Placebo.postgel', 'lowgel.baseline-lowgel.postgel']

#Placebo.postgel-Placebo.baseline
#lowgel.postgel-lowgel.baseline

##makeContrasts(lowgel=lowgel.postgel-lowgel.baseline,
            ##Placebo=Placebo.postgel-Placebo.baseline,diff=(lowgel.postgel-lowgel.baseline)-(Placebo.postgel-Placebo.baseline),
##baseline_diff=lowgel.baseline-Placebo.baseline,
##postgel_diff=lowgel.postgel-Placebo.postgel,levels=design.limma)


#['lowgel.postgel-Placebo.postgel', 'lowgel.postgel-Placebo.baseline', 'lowgel.baseline-Placebo.postgel', 'lowgel.baseline-Placebo.baseline']

#lowgel.baseline-Placebo.baseline
#lowgel.postgel-Placebo.postgel

#makeContrasts(lowgel=lowgel.postgel-lowgel.baseline,
            #Placebo=Placebo.postgel-Placebo.baseline,diff=(lowgel.postgel-lowgel.baseline)-(Placebo.postgel-Placebo.baseline),
#baseline_diff=lowgel.baseline-Placebo.baseline,
#postgel_diff=lowgel.postgel-Placebo.postgel,levels=design.limma)


#robjects.r('countTable <- data[,-seq(1:6)]')

#vec1 = robjects.StrVector(group)
#vec2 = robjects.StrVector(condition)

#robjects.globalenv["vec1"] = vec1
#robjects.globalenv["vec2"] = vec2

##fx11 <- c('Placebo', 'Placebo', 'lowgel','lowgel','lowgel')
##fx12 <- c('postgel', 'baseline', 'baseline','postgel','postgel')

#robjects.r('f1 <-factor(paste(vec1,vec2,sep="."))')

#robjects.r('design.limma <-model.matrix(~0+f1)')

#d = robjects.globalenv["design.limma"]

##vector=rpyn.ri2numpy(d)

#print str ( pd.DataFrame(numpy.array(d )))

#print str ( type(d) )

#colnames(design.limma)<-levels(f1)

#robjects.r('sedIdx <- match(rownames(design),colnames(countTable))')
#robjects.r('countTable <- countTable[,sedIdx]')

#robjects.r('counts <- DGEList(counts=countTable,genes=genesymbol)')
#robjects.r('isexpr <- rowSums(cpm(counts)>1)>=3')
#robjects.r('counts <- counts[isexpr,keep.lib.sizes=FALSE]')
#robjects.r('counts <-calcNormFactors(counts)')

##robjects.r('f <-factor(paste(design$Group,design$Condition,sep="."))')	

#robjects.globalenv["datafileName"] = datafileName	 

#robjects.r('f <-factor(paste(design$Condition,sep="."))')

##design.limma <-model.matrix(~0+f)

#robjects.r('design.limma <-model.matrix(~0+f1)')

## colnames(design.limma)<-levels(f)

#robjects.r('colnames(design.limma)<-levels(f1)')

#robjects.r('y <- voom(counts,design.limma)')
#robjects.r('corfit <- duplicateCorrelation(y,design.limma,block=design$Animal)')
#robjects.r('v <- voom(counts,design.limma,plot=TRUE,block=design$Individual,correlation=corfit$consensus)')
#robjects.r('fit <- lmFit(v,design.limma)')
##plotMDS(y,col=ifelse(design$Group=="Placebo","blue","red"))
#robjects.r('plotMDS(y,col=ifelse(design$Condition=="baseline","blue","red"))')

##cont.matrix<-makeContrasts(lowgel=lowgel.postgel-lowgel.baseline,Placebo=Placebo.postgel-Placebo.baseline,diff=(lowgel.postgel-lowgel.baseline)-(Placebo.postgel-Placebo.baseline),baseline_diff=lowgel.baseline-Placebo.baseline,postgel_diff=lowgel.postgel-Placebo.postgel,levels=design.limma)

#robjects.r('cont.matrix<-makeContrasts(postgel=postgel,baseline=baseline,diff=(postgel-baseline),levels=design.limma)')

#robjects.r('fit2<-contrasts.fit(fit,cont.matrix)')
#robjects.r('fit2 <- eBayes(fit2)')
#robjects.r('results <-decideTests(fit2,method="separate",adjust.method="none",p.value=0.05,lfc=0.585)')
#robjects.r('vennDiagram(results)')
#robjects.r('res=topTable(fit2,coef="lowgel",number=length(which(isexpr==TRUE)),sort="p")')
#robjects.r('write.csv(res,file="voom_lowgel.csv")')
#robjects.r('res=topTable(fit2,coef="Placebo",number=length(which(isexpr==TRUE)),sort="p")')
#robjects.r('write.csv(res,file="voom_Placebo.csv")')
