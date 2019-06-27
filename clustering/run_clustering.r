library(methods)
scClustering_dir="scClustering/"
setwd("~/scRNA/proj/ibd/")
source(paste(scClustering_dir,"/clustering3.r",sep=""))


rps=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/ribsomial_genes.csv",sep=""),stringsAsFactors = F)[,1])
igk=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/IGK_genes.csv",sep=""),stringsAsFactors = F)[,1])
igh=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/IGH_genes.csv",sep=""),stringsAsFactors = F)[,1])
igl=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/IGL_genes.csv",sep=""),stringsAsFactors = F)[,1])
ep=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/EP_markers.csv",sep=""),stringsAsFactors = F)[,1])
variable_seg_light_chain=c(igk,igl,igh)
malat=c("MALAT1")
xist="XIST"
jchain="JCHAIN"
#highly_expressed=c("IGKC","IGHA1","JCHAIN")
cell_cycle=strsplit("STMN1,CDC20,ANP32E,CKS1B,NUF2,ASPM,UBE2T,CENPF,RRM2,SGOL2,SGOL1,SMC4,CENPE,CCNA2,CCNB1,PTTG1,HMMR,MXD3,HIST1H4C,CENPW,H2AFV,CKS2,SMC2,ZWINT,CDK1,MKI67,C12orf75,CDKN3,NUSAP1,CCNB2,KIAA0101,PRC1,ARL6IP1,PLK1,CENPN,AURKB,TOP2A,KPNA2,HN1,TK1,BIRC5,TYMS,TPX2,UBE2C,CALM3,UBE2S,GTSE1,CLSPN,CDCA7,CENPU,GMNN,MCM7,CCDC34,CDCA5,HELLS,TMEM106C,IDH2,SNRNP25,CDT1,PCNA,UHRF1,ASF1B",",")[[1]]
hla=strsplit("HLA-F,HLA-G,HLA-A,HLA-E,HLA-C,HLA-B,HLA-DRA,HLA-DRB5,HLA-DRB1,HLA-DQA1,HLA-DQB1,HLA-DQB1-AS1,HLA-DQA2,HLA-DQB2,HLA-DOB,HLA-DMB,HLA-DMA,HLA-DOA,HLA-DPA1,HLA-DPB1",",")[[1]]
mts=c('MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8','MT-ATP6','MT-CO3','MT-ND3','MT-ND4L','MT-ND4','MT-ND5','MT-ND6','MT-CYB')
metallothionein=c('MT1HL1','MT2A','MT1E','MT1M','MT1A','MT1B','MT1F','MT1G','MT1H','MT1X')
insilico_gating=list()

insilico_gating$MITO=list()
insilico_gating$MITO$genes=mts
insilico_gating$MITO$interval=c(0,0.25)

insilico_gating$EP=list()
insilico_gating$EP$genes=ep
insilico_gating$EP$interval=c(0,0.01)

insilico_gating$RBC=list()
insilico_gating$RBC$genes=c("HBB","HBA1","HBA2")
insilico_gating$RBC$interval=c(0,0.1)

#libnames=setdiff(samples,c("122","123","180","181","186","187","197"))

#data_l=read_multiple_mtx("~/GoogleDrive/work/shared/data/IBD_10x_data/human/Grch38/",libnames,cell_interval = c(1000,1e6),noise_interval = c(100,1e6))
#save(data_l,file="/tmp/tmp_data_l.rd")


model_name="combined_081318_50_1e5_kmreg01_wjchain"
pdf(paste(model_name,".pdf",sep=""))
cluster(data_l_path="tmp/tmp_data_l_minumi_1000.rd",model_name = model_name,k=50,load_seed=F,running_mode="LSF_seeding",params=list(train_set_size  = 2000,test_set_size = 2000,insilico_gating = insilico_gating,genes_excluded = c(mts,malat,xist,metallothionein,variable_seg_light_chain,cell_cycle,hla,jchain),seeding_varmean_quantile=.92,min_umis_per_seeding_gene=50,genes_excluded_from_seeding = c(rps,ep),init_min_umis_quantile_range=c(0,.2),fixed_downsampling_min_umis=NA,reg=1e-5,n_init_seeds=1000,max_n_cores=1,init_method="TGL_kmeans",km_reg=.1))
dev.off()
