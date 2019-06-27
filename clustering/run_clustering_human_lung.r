#scDissector_dir="/home/ubuntu/GitHub/scDissector/"
#scClustering_dir="/home/ubuntu/GitHub/scClustering/"
#work_dir="/home/ubuntu/proj/ibd/"

scDissector_dir="~/Documents/GitHub/scDissector/"
scClustering_dir="~/Documents/GitHub/scClustering/"
work_dir="~/GoogleDrive/work/shared/analysis/analysis_iCD_paper_2018/"

setwd(work_dir)
source(paste(scClustering_dir,"clustering3.r",sep=""))

rps=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/ribsomial_genes.csv",sep=""),stringsAsFactors = F)[,1])
igk=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/IGK_genes.csv",sep=""),stringsAsFactors = F)[,1])
igh=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/IGH_genes.csv",sep=""),stringsAsFactors = F)[,1])
igl=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/IGL_genes.csv",sep=""),stringsAsFactors = F)[,1])
ep=unique(read.csv(file=paste(scClustering_dir,"/gene_lists/EP_markers.csv",sep=""),stringsAsFactors = F)[,1])
variable_seg_light_chain=c(igk,igl)
malat=c("MALAT1")
xist="XIST"
#highly_expressed=c("IGKC","IGHA1","JCHAIN")

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

libnames=c("normal_370", "normal_371","normal_403","normal_406","normal_410","normal_408","normal_413",
           "tumor_370","tumor_371","tumor_377","tumor_378_GRCh38","tumor_406","tumor_408","tumor_410","tumor_413")


data_l=read_multiple_mtx("~/GoogleDrive/work/shared/data/lung_human_data/data/",libnames,cell_interval = c(1000,1e6),noise_interval = c(100,1e6))
save(data_l,file="/tmp/tmp_data_l_1000.rd")


load("tmp/tmp_data_l_1000.rd")
#local_multi_core

#cluster(data_l,model_name = "test",k=4,load_seed=F,params=list(train_set_size  = 2000 ,init_test_set_size = 1000,test_set_size = 1000,insilico_gating = insilico_gating,genes_excluded = c(mts,malat,xist,metallothionein,variable_seg_light_chain),seeding_varmean_quantile=.92,min_umis_per_seeding_gene=25,genes_excluded_from_seeding = c(rps,ep),n_init_seeds=100,init_min_umis_quantile_range=c(0,.2),fixed_downsampling_min_umis=NA,reg=1e-5,running_mode="local_multi_core",max_n_cores=4,init_method="TGL_kmeans",km_reg=.1))


cluster(data_l,model_name = "combined_071718_40_2",k=40,load_seed=F,params=list(train_set_size  = 2000,init_test_set_size = 1000,test_set_size = 1000,insilico_gating = insilico_gating,genes_excluded = c(mts,malat,xist,metallothionein,variable_seg_light_chain),seeding_varmean_quantile=.94,min_umis_per_seeding_gene=25,genes_excluded_from_seeding = c(rps,ep),n_init_seeds=200,init_min_umis_quantile_range=c(.1,.3),fixed_downsampling_min_umis=NA,reg=1e-5,running_mode="local_multi_core",max_n_cores=8,init_method="TGL_kmeans",km_reg=.1))


