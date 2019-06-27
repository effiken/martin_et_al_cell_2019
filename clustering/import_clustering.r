source("~/Documents/GitHub/scDissector/projector.R")
source("~/GoogleDrive/work/scClustering/umitab_utils.r")


import_data_and_clustering_no_noise=function(data_path,libnames,cell_to_cluster,output_model_fn,reg=1e-5,insilico_gating=NULL){
  data_l=read_multiple_mtx(data_path,libnames,cell_interval = c(200,25000),noise_interval = c(100,25000),type=NA,rename_cell_ids = F)

  for (i in 1:length(data_l[[1]])){
    message(names(data_l[[1]])[i])
    if (i==1){
      umitab=data_l[[1]][[1]]
      cell_to_batch=rep(names(data_l[[1]])[1],ncol(data_l[[1]][[1]]))
    }
    else{
      umitab=cbind(umitab,data_l[[1]][[i]])
      cell_to_batch=c(cell_to_batch,rep(names(data_l[[1]])[i],ncol(data_l[[1]][[i]])))
    }
  }
  print(table(cell_to_cluster))
  print(table(cell_to_batch))
  names(cell_to_batch)=colnames(umitab)
  mask=intersect(colnames(umitab),names(cell_to_cluster))
  models=update_models(umitab[,mask],cell_to_cluster[mask])
  cell_to_batch_train=cell_to_batch
  params=list(reg=reg)
  save(umitab,models,cell_to_cluster,cell_to_batch_train,params,insilico_gating=insilico_gating,file=fn)
}
