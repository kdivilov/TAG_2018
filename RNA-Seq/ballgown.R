library(ballgown)
library(plyr)
library(genefilter)
library(GenomicRanges)

#Read data. Note that you must enter the ballgown directory from the final StringTie command in the RNA-Seq pipeline.
data = ballgown(dataDir='~/InsertDirectory/RNA-Seq/ballgown', samplePattern='BG', meas='all')
#Define covariates according to whether the QTL is homozygous or heterozygous in a genotype
pData(data) = data.frame(id=sampleNames(data), group=c(1,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,1,0,0,0,0,0,0,0,0,1,0,0))
#Keep only transcripts with > 1 FPKM
data_filt = subset(data,"rowVars(texpr(data)) >1",genomesubset=TRUE)

#Test for differentially transcribed genes
results_genes = stattest(data_filt, feature="gene", covariate="group", getFC=TRUE, meas="FPKM",libadjust = F)
#Arrange by q value
results_genes = arrange(results_genes,qval)
#Only keep genes with q value > 0.05
results_genes_sig = subset(results_genes,results_genes$qval<0.05)

#Get genes of transcripts
transtogenes = indexes(data_filt)$t2g

chr = c()
start = c()
end = c()
for(i in 1:nrow(results_genes_sig)){
  #Get chromosome info for gene i
  gene_chr = seqnames(structure(data_filt)$trans[as.character(transtogenes[which(transtogenes$g_id == as.character(results_genes_sig$id[i])),][1,1])])@unlistData@values
  #Get the ranges of the transcripts for gene i
  gene_range = ranges(structure(data_filt)$trans[as.character(transtogenes[which(transtogenes$g_id == as.character(results_genes_sig$id[i])),][1,1])])
  
  chr = c(chr,as.character(gene_chr))
  #We will bound a gene by the minimum and maximum values among all its transcripts
  start = c(start,as.numeric(min(start(gene_range))))
  end = c(end,as.numeric(max(end(gene_range))))
}

results_genes_sig_location = cbind(results_genes_sig,chr,start,end)

#Write the file containing information about the significantly differentially transcribed genes.
#Further filtering was done manually.
write.csv(results_genes_sig_location,"ballgownoutput.csv",row.names = F)



