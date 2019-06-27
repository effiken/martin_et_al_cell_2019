source("~/Documents/GitHub/scDissector/load_dataset.R")
source("~/Documents/GitHub/scDissector/projector.R")

#69,122,128,138,158,181,187,CD13IN,CD14IN

main=function(){
  path="~/GoogleDrive/work/shared/clustering_data/clustering_data_IBD3/"
  libs=c("69","122","128","138","158","181","187","CD13IN","CD14IN")
  sample_fns=paste(path,"/data_",libs,".rd",sep="")
  names(sample_fns)=libs
  model="model_combined_010318_40"
  model_fn=paste(path,model,".rd",sep="")
  lm<<-load_dataset_and_model(model_fn,sample_fns = sample_fns,min_umis = 800)
}

for (samp in lm$dataset$samples){
  u=
}