# Read in vcf file
library(data.table)
library(tidyr)
library(dplyr)

# Read in your vcf file
dlist=fread('BC0309.sorted.rg.dedup.split.realigned.recal.bam.cof20.g.vcf',header=TRUE,skip='#CHROM')

# Extract only autosome SNPs
chr_list=sapply(seq(1,22),function(x)paste0('chr',x),simplify = T)
dlist <- dlist %>% filter(`#CHROM`%in%chr_list)

Format <- dlist%>%select(FORMAT)
Info <- dlist %>% select(`#CHROM`,POS,REF,ALT)
dlist <- dlist %>% select(-`#CHROM`,-POS,-ID,-REF,-ALT,-QUAL,-FILTER,-INFO,-FORMAT)

# By checking the head of the format, we decided the position of GT (genotype info), AD (allele read depth), and DP (total read depth) separated by ":"
# As shown by the example below, GT position is 1, AD position is 2 and DP position is 3.
# Different version of GATK tool may assign different names of the varaibles. Please refer to the GATK manual for more information.
head(Format)
GT_pos=1
AD_pos=2
DP_pos=3

exinfo_x <- function(info){
  return(sapply(strsplit(info,':'),function(x){ifelse(x[GT_pos]=='./.',NA,as.numeric(strsplit(x[AD_pos],',')[[1]][2]))},simplify=T))
}
exinfo_n <- function(info){
  return(sapply(strsplit(info,':'),function(x){ifelse(x[GT_pos]=='./.',NA,as.numeric(x[DP_pos]))},simplify=T))
}
exinfo_z <- function(info){
  return(sapply(strsplit(info,':'),function(x){ifelse(x[GT_pos]=='./.',NA,ifelse(x[GT_pos]!='0/0',1,0))},simplify=T))
}

# Because the vector is too big. I run it on the cluster instead. All you need is to store the object dlist
X=sapply(dlist,exinfo_x)
N=sapply(dlist,exinfo_n)
Z=sapply(dlist,exinfo_z)

# Quick fitering
# If 99% of the cells did not detect any read, we remove such entry.
thres=floor(0.01*ncol(X))
sel=rowSums(!is.na(Z))>thres
X=X[sel,]
N=N[sel,]
Z=Z[sel,]
save(Info,X,N,Z,file='DENDRO_input.rda')
