library(methods)
library(scDissector)
library(Matrix)
library(Matrix.utils)
scClustering_dir<<-"scClustering/"
source("scClustering/clustering3.r")
args = commandArgs(trailingOnly=TRUE)
if (length(args)==3){
  data_l_path=args[1]
  model_name=args[2]
  k=as.numeric(args[3])
}

pdf(paste("saved_clustering/",model_name,"/tmp/figures_EM.pdf",sep=""))
cluster(data_l_path=data_l_path,model_name = model_name,k=k,load_seed=F,running_mode="LSF_EM")
dev.off()
