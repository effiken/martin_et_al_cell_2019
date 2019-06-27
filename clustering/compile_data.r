source("~/Documents/Github/scClustering/umitab_utils.r")


# compile data
# reads multiple mtx's and convert them to the scDissector format
#
# sample_IDs - indices of sample_annots table
# input_path - the root of the data dir. mtx will be read from: input_path/libnames[i]/filtered_or_raw/matrix.mtx. If sample_annots table has a "path" column input_path is not required. Otherwise, length should be 1 or as the number of unique amplifacations batches assocaites with the provided sample_IDs
# output_path - output clustering_data dir
# filtered_or_raw - String that determines the type of the matrix that will be read. Should be either "filtered" or "raw".
# ds_numis - A vector containing intergers with the number of UMIs that the cells should be downsampled to.
# min_umis - The minimial number of umis per cell. Cells with lower number of UMIs will be excluded.
# sample_annots_path - path of sample_annots table. table has the following columns:
#  Mandatory:
#     sample_ID - specimen sample
#     amp_batch_ID - amplification batch.
#     HTO - comma delimited HTO used for this sample. NA if CITEseq wasn't used
#     patient_ID
#  Optionl :
#     disease
#     tissue
#     treatment
#     disease
#
#  cite_seq_params= If null default parameters are used. List with the following parameters: 
#  hto_UMI_ratio_thresh - minimal ratio between the UMIs of the best and second best HTOs 
#  min_umi_maxhto - minimal number of UMIs for the best HTO

compile_data=function(sample_IDs,input_paths,output_path,filtered_or_raw="raw",prefix="",noise_interval=c(100,25000),cell_interval=c(200,25000),ds_numis=c(200,500,1000,2000),sample_annots_path=NULL,cite_seq_params=NULL){
  
  read_and_compile_umitabs(sample_IDs=sample_IDs,input_paths=input_paths,output_path=output_path,filtered_or_raw=filtered_or_raw,prefix=prefix,noise_interval=noise_interval,cell_interval=cell_interval,ds_numis=ds_numis,sample_annots_path=sample_annots_path,cite_seq_params=cite_seq_params,only_load=F)
  
}


create_clustering_data_dir=function(path,samples=c(),model_names=c()){
  if(dir.exists(path)){
    stop("Error! Cannot create ",path,". It already exists!")
  }
  dir.create(path)
  dir.create(paste(path,"/metadata",sep=""))
  samples_tab=matrix("",length(samples),3,dimnames = list(NULL,c("index","path","title")))
  samples_tab[,1]=samples
  model_versions_tab=matrix("",length(model_names),2,dimnames = list(NULL,c("title","path")))
  model_versions_tab[,1]=model_names
  sample_annots_tab=matrix("",length(samples),15,dimnames=list(NULL,c("sample_ID","amp_batch_ID","old_lib_name","Disease","Origin","tissue","status","biotherapy_status","Inclusion date (M/D/Y)","drug","Patient ID","GRID ID","COMPASS ID","CHEMISTRY","details")))
  sample_annots_tab[,1]=samples
  write.csv(file=paste(path,"/samples.csv",sep=""),samples_tab,row.names = F)
  write.csv(file=paste(path,"/model_versions.csv",sep=""),model_versions_tab,row.names = F)
  write.csv(file=paste(path,"/metadata/sample_annots.csv",sep=""),sample_annots_tab,row.names = F)
}