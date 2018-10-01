## load required packages
library(limma);library(edgeR); library(EMMREML)
library(gridExtra); library(grid); library(cobs); library(scales)
library(parallel); library(doParallel); library(ggplot2)

## load in data
load(url("https://github.com/nsmackler/Dex_status_2018/blob/master/dex_challenge.RData?raw=true"))

## there are 6 R Objects
# ATAC_counts: raw ATAC-seq read counts for all 86 samples after removal of regions with a median CPM ≤ 1 in either the Dex or control conditions
dim(ATAC_counts)
## RNA_counts: raw RNA-seq read counts for all 86 samples after removal of genes with median RPKM ≤ 2 in either the Dex or control conditions
dim(RNA_counts)
# ca_info: metadata for all 86 samples for chromatin accessibility data
dim(ca_info)
# ge_info: metadata for all 86 samples for gene expression data
dim(ge_info)
# kin: kinship matrix for all 45 individuals
dim(kin)
# Z_matrix: Z matrix for EMMA
dim(Z_matrix)

## normalize the count matrix based on the number of reads mapped to the nuclear genome
voom_CA <- voom(calcNormFactors(DGEList(counts=ATAC_counts[,6:91],lib.size = ca_info$nuc_reads_q10)),plot=FALSE)

## Use limma to remove technical effects associated with the group (i.e., cage) and library prepartion (proportion of mtDNA-mapped reads)
design <- model.matrix(~ca_info$p_mtDNA+ca_info$group)
fit <-lmFit(voom_CA,design)
intercepts=data.frame(eBayes(fit))[,1]
resid_CA=apply(residuals.MArrayLM(object=fit, voom_CA),2,function(x){x+intercepts})
rm(design); rm(fit)

#######################################################################################
###                          Main model (Chromatin accessibility)               #######
#######################################################################################
clus <- makeCluster(4)
registerDoParallel(cores=4)  
clusterExport(clus,varlist=c("resid_CA","ca_info","Z_matrix","kin"),envir=environment())
model=t(parApply(clus,resid_CA,1,function(y){
  library(EMMREML)
  design=model.matrix(~ca_info$trt+ca_info$TC1+ca_info$TC2+ca_info$TC3+ca_info$trt:scale(ca_info$elo))
  K=as.matrix(kin)
  Z=as.matrix(Z_matrix)
  emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  p=emma$pvalbeta
  varb=emma$varbetahat
  b=emma$betahat
  return(c(b,varb,p[,"none"]))
}))
model=as.data.frame(model)
colnames(model)=c('beta_intercept','beta_condition','beta_TC1','beta_TC2','beta_TC3','beta_rank_in_control','beta_rank_in_dex','var_beta_intercept','var_beta_condition','var_beta_TC1','var_beta_TC2','var_beta_TC3','var_beta_rank_in_control','var_beta_rank_in_dex','pval_intercept','pval_condition','pval_TC1','pval_TC2','pval_TC3','pval_rank_in_control','pval_rank_in_dex')
EMMA_CA_nested=model
rm(model)


#################################################################################################
###                                PCA of Dex and rank effects
#################################################################################################

## remove random effects:
## generate uhat matrix of random effects
clus <- makeCluster(4)
registerDoParallel(cores=4)  
clusterExport(clus,varlist=c("resid_CA","ca_info","Z_matrix","kin"),envir=environment())
EMMA_CA_random_effects=t(parApply(clus,resid_CA,1,function(y){
  library(EMMREML)
  design=model.matrix(~ca_info$trt+scale(ca_info$elo)+ca_info$TC1+ca_info$TC2+ca_info$TC3)
  K=as.matrix(kin)
  Z=as.matrix(Z_matrix)
  emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  return(t(emma$uhat))
}))

design=model.matrix(~ca_info$trt+scale(ca_info$elo)+ca_info$TC1+ca_info$TC2+ca_info$TC3)
EMMA_CA_random_effects=as.data.frame(EMMA_CA_random_effects)
colnames(EMMA_CA_random_effects)=rownames(emmreml(y=resid_CA[2,],X=design,Z=as.matrix(Z_matrix),K=as.matrix(kin),varbetahat=T,varuhat=T,PEVuhat=T,test=T)$uhat)

## remove tissue composition effects
reduced_CA_matrix=resid_CA-
  EMMA_CA_nested$beta_TC1%*%t(design[,"ca_info$TC1"])-
  EMMA_CA_nested$beta_TC2%*%t(design[,"ca_info$TC2"])-
  EMMA_CA_nested$beta_TC3%*%t(design[,"ca_info$TC3"])
design_mat = model.matrix(~0+ID,data=ca_info)
colnames(design_mat)=gsub("ID","",colnames(design_mat))
design_mat=as.matrix(design_mat[,colnames(EMMA_CA_random_effects)])
design_mat=design_mat[rownames(ca_info),]
## remove random (individual) effects
for (x in 1:ncol(design_mat)) {
  reduced_CA_matrix=reduced_CA_matrix-EMMA_CA_random_effects[,x]%*%t(design_mat[,x])
}
rm(x);rm(design_mat)

ca_pca=prcomp(cor(reduced_CA_matrix))
qplot(ca_pca$x[,1],ca_pca$x[,2],col=ca_info$trt)
qplot(ca_info$elo,ca_pca$x[,3],col=ca_info$trt)+stat_smooth(method="lm")

#######################################################################################
###                                delta model                                 #######
#######################################################################################
## create log2(Dex+/Dex-) matrix:
## based on the principle that log(x/y)=log(x)-log(y)
voom_CA_delta=voom_CA$E[,seq(1,86,2)]-voom_CA$E[,seq(2,86,2)]
design <- model.matrix(~ca_info[seq(1,86,2),"group"])
fit <-lmFit(voom_CA_delta,design)
intercepts=data.frame(eBayes(fit))[,1]
resid_CA_delta=apply(residuals.MArrayLM(object=fit, voom_CA_delta),2,function(x){x+intercepts})
rm(design); rm(fit); rm(intercepts)

clus <- makeCluster(4)
registerDoParallel(cores=4)  
clusterExport(clus,varlist=c("resid_CA_delta","ca_info","Z_matrix","kin"),envir=environment())
EMMA_CA_delta=t(parApply(clus,resid_CA_delta,1,function(y){
  library(EMMREML)
  design=model.matrix(~scale(elo)+TC1+TC2+TC3,data=ca_info[seq(1,86,2),])
  K=as.matrix(kin)
  Z=as.matrix(Z_matrix[seq(1,86,2),])
  emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  p=emma$pvalbeta
  varb=emma$varbetahat
  b=emma$betahat
  return(c(b,varb,p[,"none"]))
}))
EMMA_CA_delta=as.data.frame(EMMA_CA_delta)
colnames(EMMA_CA_delta)=c('beta_intercept','beta_rank','beta_TC1','beta_TC2','beta_TC3','var_beta_intercept','var_beta_rank','var_beta_TC1','var_beta_TC2','var_beta_TC3','pval_intercept','pval_rank','pval_TC1','pval_TC2','pval_TC3')

#######################################################################################
###                                RNA-seq analysis                             #######
#######################################################################################
## Normalize using limma
voom_RNA <- voom(calcNormFactors(DGEList(counts=RNA_counts,lib.size = ge_info$mapped_reads)),plot=FALSE)
## Use limma to remove technical effects associated with the group 
design <- model.matrix(~ge_info$group)
fit <-lmFit(voom_RNA,design)
intercepts=data.frame(eBayes(fit))[,1]
resid_RNA=apply(residuals.MArrayLM(object=fit, voom_RNA),2,function(x){x+intercepts})
rm(design); rm(fit); rm(intercepts)

#######################################################################################
###                          Main model (Gene expression)                       #######
#######################################################################################
clus <- makeCluster(4)
registerDoParallel(cores=4)  
clusterExport(clus,varlist=c("resid_RNA","ge_info","Z_matrix","kin"),envir=environment())
EMMA_RNA_nested=t(parApply(clus,resid_RNA,1,function(y){
  library(EMMREML)
  design=model.matrix(~ge_info$trt+ge_info$TC1+ge_info$TC2+ge_info$TC3+ge_info$trt:scale(ge_info$elo))
  K=as.matrix(kin)
  Z=as.matrix(Z_matrix)
  emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  p=emma$pvalbeta
  varb=emma$varbetahat
  b=emma$betahat
  return(c(b,varb,p[,"none"]))
}))
EMMA_RNA_nested=as.data.frame(EMMA_RNA_nested)
colnames(EMMA_RNA_nested)=c('beta_intercept','beta_condition','beta_TC1','beta_TC2','beta_TC3','beta_rank_control','beta_rank_dex','var_beta_intercept','var_beta_condition','var_beta_TC1','var_beta_TC2','var_beta_TC3','var_beta_rank_control','var_beta_rank_dex','pval_intercept','pval_condition','pval_TC1','pval_TC2','pval_TC3','pval_rank_control','pval_rank_dex')

#################################################################################################
###                                PCA of Dex and rank effects
#################################################################################################
## remove random effects:
## generate uhat matrix of random effects
clus <- makeCluster(4)
registerDoParallel(cores=4)  
clusterExport(clus,varlist=c("resid_RNA","ge_info","Z_matrix","kin"),envir=environment())
EMMA_RNA_random_effects=t(parApply(clus,resid_RNA,1,function(y){
  library(EMMREML)
  design=model.matrix(~ge_info$trt+scale(ge_info$elo)+ge_info$TC1+ge_info$TC2+ge_info$TC3)
  K=as.matrix(kin)
  Z=as.matrix(Z_matrix)
  emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  return(t(emma$uhat))
}))
design=model.matrix(~ge_info$trt+scale(ge_info$elo)+ge_info$TC1+ge_info$TC2+ge_info$TC3)
EMMA_RNA_random_effects=as.data.frame(EMMA_RNA_random_effects)
colnames(EMMA_RNA_random_effects)=rownames(emmreml(y=resid_RNA[2,],X=design,Z=as.matrix(Z_matrix),K=as.matrix(kin),varbetahat=T,varuhat=T,PEVuhat=T,test=T)$uhat)

## remove tissue composition effects
reduced_expression_matrix=resid_RNA-
  EMMA_RNA_nested$beta_TC1%*%t(design[,"ge_info$TC1"])-
  EMMA_RNA_nested$beta_TC2%*%t(design[,"ge_info$TC2"])-
  EMMA_RNA_nested$beta_TC3%*%t(design[,"ge_info$TC3"])
design_mat = model.matrix(~0+ID,data=ge_info)
colnames(design_mat)=gsub("ID","",colnames(design_mat))
design_mat=as.matrix(design_mat[,colnames(EMMA_RNA_random_effects)])
design_mat=design_mat[rownames(ge_info),]
## remove random (individual) effects
for (x in 1:ncol(design_mat)) {
  reduced_expression_matrix=reduced_expression_matrix-EMMA_RNA_random_effects[,x]%*%t(design_mat[,x])
}
rm(x);rm(design_mat)

ge_pca=prcomp(cor(reduced_expression_matrix))
qplot(ge_pca$x[,1],ge_pca$x[,2],col=ge_info$trt)
qplot(ge_info$elo,ge_pca$x[,3],col=ge_info$trt)+stat_smooth(method="lm")

#######################################################################################
###                                delta model                                 #######
#######################################################################################
voom_RNA_delta=voom_RNA$E[,seq(1,86,2)]-voom_RNA$E[,seq(2,86,2)]
design <- model.matrix(~ge_info[seq(1,86,2),"group"])
fit <-lmFit(voom_RNA_delta,design)
intercepts=data.frame(eBayes(fit))[,1]
resid_RNA_delta=apply(residuals.MArrayLM(object=fit, voom_RNA_delta),2,function(x){x+intercepts})
rm(design); rm(fit); rm(intercepts)
clus <- makeCluster(4)
registerDoParallel(cores=4)  
clusterExport(clus,varlist=c("resid_RNA_delta","ge_info","Z_matrix","kin"),envir=environment())
EMMA_RNA_delta=t(parApply(clus,resid_RNA_delta,1,function(y){
  library(EMMREML)
  design=model.matrix(~scale(elo)+TC1+TC2+TC3,data=ge_info[seq(1,86,2),])
  K=as.matrix(kin)
  Z=as.matrix(Z_matrix[seq(1,86,2),])
  emma=emmreml(y=y,X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
  p=emma$pvalbeta
  varb=emma$varbetahat
  b=emma$betahat
  return(c(b,varb,p[,"none"]))
}))
EMMA_RNA_delta=as.data.frame(EMMA_RNA_delta)
colnames(EMMA_RNA_delta)=c('beta_intercept','beta_rank','beta_TC1','beta_TC2','beta_TC3','var_beta_intercept','var_beta_rank','var_beta_TC1','var_beta_TC2','var_beta_TC3','pval_intercept','pval_rank','pval_TC1','pval_TC2','pval_TC3')

#######################################################################################
###                           GO enrichment with topGO                          #######
#######################################################################################
library(topGO)
library(biomaRt)
mmul = useMart(biomart = "ENSEMBL_MART_ENSEMBL",dataset="mmulatta_gene_ensembl")

##ensembl_gene_names is a vector of ensembl gene names from your expression matrix
mmul_GO<-getBM(attributes=c('ensembl_gene_id', 'go_id','name_1006'), filters = 'ensembl_gene_id', values = rownames(EMMA_RNA_nested), mart = mmul)

geneID2GO=lapply(unique(mmul_GO$ensembl_gene_id),function(x){sort(mmul_GO[mmul_GO$ensembl_gene_id==x,"go_id"])})
names(geneID2GO)=unique(mmul_GO$ensembl_gene_id)

rm(mmul);rm(mmul_GO)

## rank_tvals is a vector of standardized rank betas (beta/se(beta))
## the names(rank_tvals) needs to be the names of the ensembl genes (in the same order as the standardized betas)
rank_tvals=EMMA_RNA_nested$beta_rank_control/sqrt(EMMA_RNA_nested$var_beta_rank_control)
names(rank_tvals)=rownames(EMMA_RNA_nested)
## change the number to be whatever the standardized beta value is that corresponds to your FDR threshold:
#FDR=20%; tval=2.26318

## note that the tvals vector is multiplied by -1 here because the KS enrichment tests for the lower end of the distribution (it is assuming p-values are used)
GOdata_highrank<-new("topGOdata",description="Simple session",ontology="BP",allGenes=(-1*rank_tvals),geneSel=function(x) x<(-2.26318),nodeSize=10,annot=annFUN.gene2GO,gene2GO = geneID2GO)
## vector tvals not multiplied by -1 in this example (because low values == more highly expressed in low ranking animals) 
GOdata_lowrank<-new("topGOdata",description="Simple session",ontology="BP",allGenes=rank_tvals,geneSel=function(x) x<(-2.26318),nodeSize=10,annot=annFUN.gene2GO,gene2GO = geneID2GO)

#Calculate enrichment using the KS weight01 algorithm
## GO categories enriched in genes more highly expressed in high-ranking animals
KS_weight01<-runTest(GOdata_highrank,algorithm="weight01",statistic="ks")
topnodes=max(KS_weight01@geneData[4])
GOdata_highrank_results=GenTable(GOdata_highrank,KS_weight01=KS_weight01,orderBy="KS_weight01",ranksOf="KS_weight01",topNodes=topnodes,numChar=200)
rownames(GOdata_highrank_results)=GOdata_highrank_results$GO.ID

## GO categories enriched in genes more highly expressed in low-ranking animals
KS_weight01<-runTest(GOdata_lowrank,algorithm="weight01",statistic="ks")
topnodes=max(KS_weight01@geneData[4])
GOdata_lowrank_results=GenTable(GOdata_lowrank,KS_weight01=KS_weight01,orderBy="KS_weight01",ranksOf="KS_weight01",topNodes=topnodes,numChar=200)
rownames(GOdata_lowrank_results)=GOdata_lowrank_results$GO.ID

## merge the outputs
a=merge(GOdata_lowrank_results[,c(1,2,6)],GOdata_highrank_results[,c(1,6)],by=1)
a$KS_weight01.x=as.numeric(a$KS_weight01.x)
a$KS_weight01.y=as.numeric(a$KS_weight01.y)
colnames(a)=c("GO.ID","Term","enrichment_lowrank","enrichment_highrank")
a=a[order(a$enrichment_lowrank),]
## write the file
write.table(a,file="dex_GO_enrichment.txt",col.names = T,row.names = F,quote = F,sep="\t")
rm(a)
                    
#######################################################################################
###                                WGCNA analysis                                 #######
#######################################################################################                 
library(WGCNA)
### Rename rows and columns.
rownames(ge_info)=paste0(ge_info$ID,"_",ge_info$trt)
colnames(RNA_counts)=rownames(ge_info)

### Subset Dex- samples
ge_info=ge_info[which(ge_info$trt=="control"),]
RNA_counts=RNA_counts[,which(colnames(RNA_counts) %in% rownames(ge_info))]

## Use limma to remove technical effects associated with the group from just the Dex- samples
voom_RNA <- voom(calcNormFactors(DGEList(counts=RNA_counts,lib.size = ge_info$mapped_reads)),plot=TRUE)
design <- model.matrix(~ge_info$group)
fit <-lmFit(voom_RNA,design)
intercepts=data.frame(eBayes(fit))[,1]
resid_RNA=apply(residuals.MArrayLM(object=fit, voom_RNA),2,function(x){x+intercepts})

### Put everything into WGCNA friendly format (samples as rows and genes as columns) 
datExpr = as.data.frame(t(resid_RNA))
### QC routine from WGCNA, says everything's fine.
gsg = goodSamplesGenes(datExpr, verbose = 3)
print(paste("Are samples and genes OK?",gsg$allOK))

if (!gsg$allOK)
{
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
    if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))
    # Remove the offending genes and samples from the data:
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

traitRows=ge_info[which(rownames(ge_info) %in% colnames(resid_RNA)),c(4,5,6,7)]
datExpr=datExpr[order(rownames(datExpr)),]
traitRows=traitRows[order(rownames(traitRows)),]

allowWGCNAThreads()

# Select power to rise the coexpression matrix to, in order to get the adjacency.
powers = c(c(1:11), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,blockSize=ncol(datExpr))
pow=sft$powerEstimate

#Calculate TOM and obtain modules.
net = blockwiseModules(datExpr, blocks=rep(1,ncol(datExpr)),power = pow,TOMType = "unsigned", maxBlockSize=ncol(datExpr),minModuleSize = 30,reassignThreshold = 0, mergeCutHeight = 0.25,numericLabels = TRUE, pamRespectsDendro = FALSE,saveTOMs = TRUE,saveTOMFileBase = "CTL_TOM",verbose = 3)

## convert label names to colors
colnames(net$MEs)=labels2colors(colnames(net$MEs))
net$named_colors=labels2colors(net$colors)

## calculate association between rank and eigengene values controlling for tissue composition (TC1-3)
apply(net$MEs,2,function(y){
  return(summary(lm(y~scale(traitRows$elo)+traitRows$TC1+traitRows$TC2+traitRows$TC3))$coef[2,])
  })
 
## prepare "gene_modules" data frame for GO enrichment analysis
gene_modules=data.frame(wgcna_input_color=net$named_colors,cluster=net$colors); rownames(gene_modules)=colnames(datExpr)

## calculate GO enrichment using topGO
## WGCNA GO enrichment
WGCNA_rank_GO=lapply(unique(gene_modules$wgcna_input_color),function(color){
module_col_vector=as.integer(gene_modules$wgcna_input_color==color)
names(module_col_vector)=rownames(gene_modules)
GOdata<-new("topGOdata",description="Simple session",ontology="BP",allGenes=module_col_vector,geneSel=function(x) x>0.5,nodeSize=10,annot=annFUN.gene2GO,gene2GO = geneID2GO)
FET_weight01<-runTest(GOdata,algorithm="weight01",statistic="fisher")
return(GenTable(GOdata,FET_weight01=FET_weight01,orderBy="FET_weight01",ranksOf="FET_weight01",topNodes=max(FET_weight01@geneData[4]),numChar=200))
})
names(WGCNA_rank_GO)=unique(gene_modules$wgcna_input_color)

## write GO results
WGCNA_rank_GO_table=Reduce(function(...) merge(..., by=c(1,2),all=T),lapply(WGCNA_rank_GO,function(x){x[,c(1,2,6)]}))
colnames(WGCNA_rank_GO_table)[3:19]=names(WGCNA_rank_GO_table)
