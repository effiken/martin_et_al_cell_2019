library(methods)
library(scDissector)
library(Matrix)
library(Matrix.utils)
scClustering_dir<<-"scClustering/"
source("scClustering/clustering3.r")
args = commandArgs(trailingOnly=TRUE)
if (length(args)==2){
  model_name=args[1]
  km_i=as.numeric(args[2])
} else { 
  km_i=1
  model_name="lung_190529"
}
pdf(paste("saved_clustering/",model_name,"/tmp/figures_",km_i,".pdf",sep=""))
tmp_seed_rd_filename=paste("saved_clustering/",model_name,"/tmp/init_input_",model_name,".rd",sep="")
load(tmp_seed_rd_filename)
print(paste("res=get_seed(k=",k,", train_umitab_fn, cell_to_batch,noise_models,high_var_genes,params,model_name=",model_name,", km_i=",km_i,")",sep=""))
set.seed(km_i)

cell_to_cluster=run_kmeans(k,train_umitab_fn,ds,high_var_genes,params=params,model_name=model_name,km_i=km_i,trace=trace)
res=get_seed(cell_to_cluster,train_umitab_fn,cell_to_batch,noise_models,params=params,model_name=model_name,km_i=km_i,trace=trace)

save(list=ls(envir = res),file=paste("saved_clustering/",model_name,"/tmp/data_",km_i,".rd",sep=""),envir=res)
writeLines(as.character(res$train_totll),con =paste("saved_clustering/",model_name,"/tmp/totll_",km_i,".txt",sep=""))
dev.off()
