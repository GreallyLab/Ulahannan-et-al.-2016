#RNA-seq Analysis - #R code
#read in htseq count table
read.table("table",sep="\t",header=T)->exp
exp[,c(2:7)]->hff2
#DEseq Bioconductor Package used to obtain differentially expressed genes and normalized gene counts
library("DESeq")
conds <- factor( c( "I1" ,"I2", "I3", "U1", "U2", "U3" ) )
cds <- newCountDataSet( table, conds )
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
hff3<-( counts( cds, normalized=TRUE ) )
cds <- estimateDispersions( cds )
res <- nbinomTest( cds, "table1", "table2" )

#HELP-tagging/HELP-GT - #R code
#read in reference mspI count table
read.table("MspRef.txt",header=T)->msp
#read in hpaII hg19 annotation table
read.table("anno_hpaii_hg19.txt",skip=1)->hpaii
#read in hpaII count table
read.table("hg19.BGT_U1_HFF.AC263LACXX.lane_7.hcount",skip=1)->U
#merge hpaII count table and annotations table
merge(U1,hpaii,by=c(1),all=TRUE)->table1
#merge with mspI count table 
merge(table1,msp,by=c(1),all=TRUE)->table2
#convert NAs to 0
table2[is.na(table2)]<-0
table2[,c(1,2,3,6,19)]->table3
#calculate angles
colnames(table3)<-c("id","chr","pos","U1","msp")
sum(table3[,c("U1")])
(table3[,c("msp")])/40550229->msp_norm
(table3[,c("U1")])/7034891->U1_norm
cbind(table3,U1_norm,msp_norm)->table4
atan2(table4[,c("U1_norm")],table4[,c("msp_norm")])->U1ANGLE
U1_angle<-U1ANGLE*(2/pi)*100
cbind(hpaii,msp_norm,U1_norm,U1_angle)->BGT_Angle
#filter down to most conident loci 
(BGT_Angle[,c("U1_norm")] < BGT_Angle[,c("msp_norm")])->Filtered_U1
BGT_Angle2[which(BGT_Angle2[,c("Filtered_U1")] == TRUE),]->BGTangle_Final
BGTangle_Final[,c(1:11,14,16,18)]->BGTangle_final2

#linear model to calculate differentially methylated loci
#read in angle file
read.table("HELPtagHFFall_ANGLE.txt",header=T)->hff
#read in covariate file 
read.table("HELPtag_cv.txt",header=T)->CV
for (i in c(1:1656979)) {
as.numeric(hff[i,c(12:14)])->U
as.numeric(hff[i,c(17:19)])->I
cbind(U,I)->Methylation
factor(CV[,c(4)])->I
lmfit = lm(Methylation ~ Infection)
t(coef(summary(lmfit))[,c(4)]) -> PVal
t(coef(summary(lmfit))[,c(1)]) -> MethEstimate
write.table(PVal,file="~/LM_Infection_Pvalue.txt",quote=F,sep="\t",row.names=F,col.names=F,append=T) 
write.table(MethEstimate,file="~/LM_Infection_MethEstimate.txt",quote=F,sep="\t",row.names=F,col.names=F,append=T) 
}

#Permutation Tests - #R code
bed_boot_scoot<-function(bed1, bed2, n.iter=100){
bedTools.2in<-function(functionstring="/apps1/bedtools/current/intersectBed -wao",bed1,bed2,opt.string="")
{
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
 
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=try(read.table(out,header=F),silent=T) 
  unlink(a.file);unlink(b.file);unlink(out)
  if (is(res, "try-error")) return(mat.or.vec(0,3)) else return(res)
}
bedTools2.2in<-function(functionstring="/apps1/bedtools/current/shuffleBed",bed1,bed2,opt.string="")
{
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
 
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
  
  #write bed formatted dataframes to tempfile
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
 
  # create the command string and call the command using system()
  command=paste(functionstring,"-i",a.file,"-g",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=try(read.table(out,header=F),silent=T) 
  unlink(a.file);unlink(b.file);unlink(out)
  if (is(res, "try-error")) return(mat.or.vec(0,3)) else return(res)
}

mainobs<-try(bedTools.2in(bed1=bed1, bed2=bed2, opt.string="-wao"), silent=T)
A<-ifelse(!is.null(ncol(mainobs)),
 { sum(mainobs[,ncol(mainobs)])},

 {1})
mainobs<-A
simulations<-sapply(1:n.iter, function(x){
print(x)
newbed<-bedTools2.2in(bed1=read.table("~/SecondStrandCounts/IMR90minus10"),bed2=read.table("~/fakechrfile"), opt.string="-excl RDIP1/hektracks/gaps.bed")
temp<-try(bedTools.2in(bed1=newbed, bed2=bed2), silent=T)
ifelse(!is.null(ncol(temp)),
 {return(sum(temp[,ncol(temp)]))},
 {1})})
return(list(mainobs, simulations))}

result1<-bed_boot(bed1=bed1, bed2=bed2, bed3=bed3)


#ATAC-seq Analysis - #bash script
#alignment to combined genome (hg19 + Toxoplasma genomes)
qsub bwa mem -p data_from_oxford/Genomes/TgME49_9.0+hg19.fasta data_from_oxford/Human_ATAC/I1_ATAC.fq.gz
samtools view -Sb I1_ATAC_hg19TgME49v9.sam  >  I1_ATAC_hg19TgME49v9.bam
#remove PCR duplicates
qsub -l mem_free=20G -b y -V -cwd -m e -M netha.ulahannan@phd.einstein.yu.edu samtools rmdup I1_ATAC_hg19TgME49v9.bam I1_ATAC_hg19TgME49v9_NoDup.bam
samtools sort I2_ATAC_NoDup.bam I2_ATAC_NoDupSort
samtools index I2_ATAC_NoDupSort.bam
samtools view -H I2_ATAC_NoDupSort.bam | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools reheader - I2_ATAC_NoDupSort.bam > I2_ATAC_NoDupSort_reheader.bam
samtools index I2_ATAC_NoDupSort_reheader.bam
#differentiate between reads aligning to the human and t.gondii genomes
samtools view -b I2_ATAC_NoDupSort_reheader.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY | samtools sort - I2_ATAC_NoDupSort-chrM
samtools view -b I1sort_ATAC_hg19TgME49v9_NoDup.bam TGME49_chrIa TGME49_chrIb TGME49_chrII TGME49_chrIII TGME49_chrIV TGME49_chrV TGME49_chrVI TGME49_chrVIIa TGME49_chrVIIb TGME49_chrVIII TGME49_chrIX TGME49_chrX TGME49_chrXI TGME49_chrXII  | samtools sort - I1_ATACtoxo_hg19TgME49v9
#call peaks using Macs2 
macs2 callpeak -t /oxford/nulahannan/I1_ATACtoxo_hg19TgME49v9.bam -f BAM -g 7e7 -n I1toxomacs -B -q 0.05 --nomodel --nolambda
macs2 callpeak -t /oxford/nulahannan/I1_ATAChuman_hg19TgME49v9_-chrM.bam -f BAM -g hs -n I1human_macs -B -q 0.05 --nomodel --nolambda
bamToBed -i I1_ATAChuman_hg19TgME49v9_-chrM.bam > I1_ATAChuman_hg19TgME49v9.bed

## Versions of software used in unix environment
module list
Currently Loaded Modulefiles:
  1) samtools/1.2/gcc.4.4.7            4) numpy/1.9.0/python.2.7.8
  2) bwa/0.7.13/gcc.4.4.7              5) MACS2/2.1.0-update/python.2.7.8
  3) python/2.7.8/gcc.4.4.7			   6) bedtools2/2.24.0/gcc.4.4.7










