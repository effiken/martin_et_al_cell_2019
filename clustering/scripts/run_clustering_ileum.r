#install.packages('tglkmeans', repos=c(getOption('repos'), 'https://tanaylab.bitbucket.io/repo'))


library(tglkmeans)
main_clustering=function(clustering_root="martin_et_al_cell_2019/clustering/",download_umitabs=T,raw_umitabs_path=NULL){


scClustering_dir<<-paste(clustering_root,"/scripts",sep="")
source(paste(scClustering_dir,"/clustering3.r",sep=""))


rps=unique(read.csv(file=paste(clustering_root,"/gene_lists/ribsomial_genes.csv",sep=""),stringsAsFactors = F)[,1])
igk=unique(read.csv(file=paste(clustering_root,"/gene_lists/IGK_genes.csv",sep=""),stringsAsFactors = F)[,1])
igh=unique(read.csv(file=paste(clustering_root,"/gene_lists/IGH_genes.csv",sep=""),stringsAsFactors = F)[,1])
igl=unique(read.csv(file=paste(clustering_root,"/gene_lists/IGL_genes.csv",sep=""),stringsAsFactors = F)[,1])
ep=unique(read.csv(file=paste(clustering_root,"/gene_lists/EP_markers.csv",sep=""),stringsAsFactors = F)[,1])
variable_seg_light_chain=c(igk,igl,igh)
malat=c("MALAT1")
xist="XIST"
jchain="JCHAIN"
cell_cycle=strsplit("STMN1,CDC20,ANP32E,CKS1B,NUF2,ASPM,UBE2T,CENPF,RRM2,SGOL2,SGOL1,SMC4,CENPE,CCNA2,CCNB1,PTTG1,HMMR,MXD3,HIST1H4C,CENPW,H2AFV,CKS2,SMC2,ZWINT,CDK1,MKI67,C12orf75,CDKN3,NUSAP1,CCNB2,KIAA0101,PRC1,ARL6IP1,PLK1,CENPN,AURKB,TOP2A,KPNA2,HN1,TK1,BIRC5,TYMS,TPX2,UBE2C,CALM3,UBE2S,GTSE1,CLSPN,CDCA7,CENPU,GMNN,MCM7,CCDC34,CDCA5,HELLS,TMEM106C,IDH2,SNRNP25,CDT1,PCNA,UHRF1,ASF1B",",")[[1]]
hla=strsplit("HLA-F,HLA-G,HLA-A,HLA-E,HLA-C,HLA-B,HLA-DRA,HLA-DRB5,HLA-DRB1,HLA-DQA1,HLA-DQB1,HLA-DQB1-AS1,HLA-DQA2,HLA-DQB2,HLA-DOB,HLA-DMB,HLA-DMA,HLA-DOA,HLA-DPA1,HLA-DPB1",",")[[1]]
mts=c('MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8','MT-ATP6','MT-CO3','MT-ND3','MT-ND4L','MT-ND4','MT-ND5','MT-ND6','MT-CYB')
metallothionein=c('MT1HL1','MT2A','MT1E','MT1M','MT1A','MT1B','MT1F','MT1G','MT1H','MT1X')
stress=strsplit("DUSP10,CSRNP1,IRF1,HSP90AB1,DNAJA1,TUBB4B,DDIT4,RGCC,SOCS1,UBB,PMAIP1,SERTAD1,MYADM,HSPA6,HSPE1,HSPA1A,HSPA1B,BAG3,HSPH1,HSP90AA1,DNAJB1,PIK3R3",",")[[1]]
mtrnr=strsplit("MTRNR2L11,MTRNR2L12,MTRNR2L13,MTRNR2L6,MTRNR2L10,MTRNR2L8,MTRNR2L7,MTRNR2L5,MTRNR2L4,MTRNR2L1,MTRNR2L3",",")[[1]]

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


#


  
  dir.create("tmp",showWarnings = F)
  if (download_umitabs){
    download.file("https://www.dropbox.com/s/8vaezppfzno993l/tmp_data_l_minumi_1000.rd?dl=1",destfile = "tmp/tmp_data_l_minumi_1000.rd")
  }
  else {
    data_l=read_multiple_mtx(raw_umitabs_path,libnames,cell_interval = c(1000,1e6),noise_interval = c(100,1e6))
    save(data_l,file="/tmp/tmp_data_l_minumi_1000.rd")
  }
  
  model_name="combined_081718_50_5e6_kmreg01"
  
  pdf(paste(model_name,".pdf",sep=""))
  cluster(data_l_path="tmp/tmp_data_l_minumi_1000.rd",model_name = model_name,k=50,load_seed=F,
          running_mode="LSF_seeding",
#            running_mode="local",
            params=list(
            train_set_size  = 1000,
            test_set_size = 1000,
            insilico_gating = insilico_gating,
            genes_excluded = c(mts,mtrnr,malat,xist,metallothionein,variable_seg_light_chain,cell_cycle,hla,jchain,stress),
            seeding_varmean_quantile=.92,
            min_umis_per_seeding_gene=50,
            genes_excluded_from_seeding = c(rps,ep),
            init_min_umis_quantile_range=c(.1,.4),
            fixed_downsampling_min_umis=NA,
            reg=5e-6,n_init_seeds=10000,
            max_n_cores=1,
            init_method="TGL_kmeans",
            km_reg=.1,
            init_alpha_noise=0.015))
  dev.off()
}