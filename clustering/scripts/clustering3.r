library(gplots)
library(Matrix.utils)
library(seriation)
library(scDissector)

 #devtools::install_bitbucket('tanaylab/tgstat@default');devtools::install_bitbucket("tanaylab/tglkmeans", ref='default')

source(paste(scClustering_dir,"umitab_utils.r",sep="/"))
source(paste(scClustering_dir,"make_LSF_seeding_job_file.r",sep="/"))
source(paste(scClustering_dir,"make_LSF_EM_job_file.r",sep="/"))
set.seed(42)





colgrad=c(colorRampPalette(c("white",colors()[378]))(50),colorRampPalette(c(colors()[378],"orange", "tomato","mediumorchid4"))(50))



getNoiseModels=function(data_l,low_umi_count_barcode_weight=0.9){
  batch_umi_counts=sapply(data_l$umitabs_by_batch,rowSums)
  noise_models=low_umi_count_barcode_weight*data_l$noise_counts/sum(data_l$noise_counts)+(1-low_umi_count_barcode_weight)*batch_umi_counts/sum(batch_umi_counts)
  noise_models=t(t(noise_models)/colSums(noise_models))
  return(noise_models)
}




update_alpha=function(umitab,models,noise_models,cell_to_batch,reg,max_noise_fraction=.2,cell_to_cluster=NULL,max_n_cores=1){
  
  update_alpha_single_batch_wrapper=function(bi){
    return(update_alpha_single_batch(umitab[,cell_to_batch==bi],models,cell_to_cluster=cell_to_cluster[cell_to_batch==bi],noise_models[,bi],reg,max_noise_fraction))
  }
  
 # alpha_noise=rep(NA,ncol(noise_models))
  #names(alpha_noise)=colnames(noise_models)
  #   message("Updating noise params")
  if (max_n_cores==1){
      alpha_noise=unlist(lapply(colnames(noise_models),update_alpha_single_batch_wrapper))
     names(alpha_noise)=colnames(noise_models)
  }
  else{
    library(parallel)
    numCores=detectCores()
    gc()
    cl <-makeCluster(numCores,type="FORK")
    njobs=min(numCores,max_n_cores)
    alpha_noise=unlist(mclapply(colnames(noise_models),update_alpha_single_batch_wrapper,mc.cores =njobs,mc.preschedule = F))
    stopCluster(cl)
    names(alpha_noise)=colnames(noise_models)
  }
  return(alpha_noise)
}


getBatchCorrectedLikelihood=function(umitab,models,noise_models,cell_to_batch,alpha_noise=NULL,reg){
  
  l=split_sparse(umitab,cell_to_batch)
  ll=matrix(NA,ncol(umitab),ncol(models))
  rownames(ll)=colnames(umitab)
  colnames(ll)=colnames(models)
  
  for (bi in 1:ncol(noise_models)){
    res_l=getOneBatchCorrectedLikelihood(umitab=l[[bi]],models,noise_models[,bi],alpha_noise=alpha_noise[bi],reg=reg)
    ll[rownames(res_l$ll),]=res_l$ll
    ll_noise=res_l$ll_noise[,1]
  }
  
  return(list(ll=ll,ll_noise=ll_noise))
}


get_total_likelihood=function(ll){
  return(mean(apply(ll,1,max)))
}


update_cell_to_batch=function(fn,dwt){
  load(fn)
  cell_to_batch=paste(dwt$exp,dwt$Mouse_ID,dwt$pcr_pool,sep=".")
  names(cell_to_batch)=rownames(dwt)
  save(umitab,ds,models,cell_to_cluster,ll,l_cells_per_cluster,cell_to_batch,chisq_res,params,file=fn)
}





compare_clustering_versions=function(f1,f2){
  e1=new.env()
  e2=new.env()
  load(f1,envir = e1)
  load(f2,envir = e2)
  return(table(e1$cell_to_cluster,e2$cell_to_cluster))
}


run_kmeans=function(k,train_umitab_fn,ds_fixed=NA,high_var_genes_fixed=NA,params=NULL,model_name="",km_i=NA,trace=F){
  
  require(tglkmeans)
    
  load(train_umitab_fn)
  if (params$redownsample){
    min_umis=pmax(params$hard_min_umis,round(quantile(params$numis_seeding[params$numis_seeding>params$hard_min_umis],runif(n=1,min=params$min_umis_quantile_range[1],max=params$min_umis_quantile_range[2]))))
    load(train_umitab_fn)
    ds=downsample(train_umitab[params$seeding_genes,],min_umis,parallel=F)
    rm(train_umitab)
    gc()
    high_var_genes=get_highly_variable_genes(ds,params$seeding_genes,params$min_umis_per_seeding_gene,params$varmean_quantile)
  }
  else{
    ds=ds_fixed
    high_var_genes=high_var_genes_fixed
    
  }
  
  print(paste("selected ",length(high_var_genes)," variable genes",sep=""))
  
  print(paste("KM #",km_i,sep=""))
  cells=colnames(ds)
#  df_to_km=as.data.frame(t(as.matrix(log2(ds[high_var_genes,]+params$km_reg))))
  df_to_km=t(as.matrix(log2(ds[high_var_genes,]+params$km_reg)))
  rm(ds)
  gc()
 
  km1=TGL_kmeans(df_to_km, k, metric = "euclid", max_iter = 100, min_delta = 1e-04,verbose = F, keep_log = F, id_column = F,reorder_func = NULL, hclust_intra_clusters = FALSE, seed = NULL,bootstrap = FALSE)
  cell_to_cluster=as.character(km1$cluster)
  names(cell_to_cluster)=cells
  print("TGLkmeans done")
  rm(df_to_km)
  rm(km1)
  gc()
  res=list()
  res$cell_to_cluster=cell_to_cluster
  res$high_var_genes=high_var_genes
  res$min_umis=min_umis
  return(res)
}

get_seed=function(cell_to_cluster,train_umitab_fn,cell_to_batch,noise_models,params=NULL,model_name="",km_i=NA,trace=F){
  load(train_umitab_fn)
  #  print
  alpha_noise=params$init_alpha_noise
  for (i in 1:10){
    prev_alpha_noise=alpha_noise
    models=update_models_debatched(umis = train_umitab[params$clustering_genes,names(cell_to_cluster)],cell_to_cluster = cell_to_cluster, batch=cell_to_batch[names(cell_to_cluster)],noise_models=noise_models[params$clustering_genes,,drop=F],alpha_noise =alpha_noise)
    alpha_noise=update_alpha(umitab = train_umitab[params$clustering_genes,names(cell_to_cluster)],models = models,noise_models = noise_models,cell_to_batch = cell_to_batch[names(cell_to_cluster)],cell_to_cluster=cell_to_cluster,reg=params$reg)
   
    alpha_dist=sqrt(mean((alpha_noise-prev_alpha_noise)^2))
    print(round(alpha_dist,digits=6))
    if (alpha_dist<1e-3){
      print(alpha_noise)
      break
    }
  }
  
  models=update_models_debatched(train_umitab[params$clustering_genes,names(cell_to_cluster)],cell_to_cluster,cell_to_batch[names(cell_to_cluster)],noise_models=noise_models[params$clustering_genes,,drop=F],alpha_noise =alpha_noise)
  
  train_bcl=getBatchCorrectedLikelihood(train_umitab[params$clustering_genes,],models[params$clustering_genes,],noise_models[params$clustering_genes,,drop=F],cell_to_batch[colnames(train_umitab)],alpha_noise = alpha_noise,reg = params$reg)
  train_totll=get_total_likelihood(train_bcl$ll)
  cell_to_cluster=MAP(train_bcl$ll)
  print(table(cell_to_batch[names(cell_to_cluster)],cell_to_cluster))   
 
  
  print(paste("Done #",km_i," train LL= ",round(train_totll,digits=4)))
  res=new.env()
  res$cell_to_cluster=cell_to_cluster
  res$train_totll=train_totll
  res$alpha_noise=alpha_noise
  return(res)
}


get_seed_launcher=function(train_umitab_fn,ds,high_var_genes,cell_to_batch,noise_models,k,data_l_path=NULL,params,model_name,running_mode,parallel_params=NULL,trace=F){
  
  get_seed_wrapper=function(x){
    res1=run_kmeans(k,train_umitab_fn,ds,high_var_genes,params=params,model_name=model_name,km_i=x,trace=trace)
    res=get_seed(res1$cell_to_cluster,train_umitab_fn,cell_to_batch,noise_models,params=params,model_name=model_name,km_i=x,trace=trace)
    res$min_umis=res1$min_umis
    res$high_var_genes=res1$high_var_genes
    save(list=ls(envir = res),file=paste("saved_clustering/",model_name,"/tmp/data_",x,".rd",sep=""),envir=res)
    writeLines(as.character(res$train_totll),con =paste("saved_clustering/",model_name,"/tmp/totll_",x,".txt",sep=""))
    return(res$train_totll)
  }
  
  default_params=list(
    seeding_genes=rownames(ds),
    clustering_genes=rownames(ds),
    reg=1e-5,
    min_umis_per_seeding_gene=30,
    varmean_quantile=0.95,
    minimum_cluster_size=3,
    min_umis_quantile_range=c(0,.2),
    hard_min_umis=150,
    init_alpha_noise=0.01,
    km_reg=.1
  )
  
  params2=default_params
  if (length(params)>0){
    for (param in names(params)){
      params2[[param]]=params[[param]]
    }
  }
  params=params2
  
  
  
 
  if (running_mode=="local"){
    l=lapply(1:params$n_init_seeds,get_seed_wrapper)
  }
  else if (running_mode=="local_multi_core"){

    library(parallel)
    numCores=detectCores()
    gc()
    
    njobs=min(numCores,params$n_init_seeds,params$max_n_cores)
    print(paste("Running ",params$n_init_seeds," kmeans (",njobs," in parallel)",sep=""))
    cl <-makeCluster(njobs,type="FORK")

    l=mclapply(1:params$n_init_seeds,get_seed_wrapper,mc.cores =njobs,mc.preschedule = F)
    stopCluster(cl)
  }
  else if (running_mode=="LSF_seeding"){
 
    
    tmp_path=paste("saved_clustering/",model_name,"/tmp/",sep="")
    tmp_seed_rd_filename=paste(tmp_path,"/init_input_",model_name,".rd",sep="")
    save(k,train_umitab_fn,ds,cell_to_batch,noise_models,high_var_genes,params,model_name,file = tmp_seed_rd_filename)
    



    make_seeding_job_file(job_name="sc_seeding",account=parallel_params$account,queue=parallel_params$queue,ncores=1,memk=parallel_params$memk,wall_time=parallel_params$wall_time,headnode=parallel_params$headnode,log_prefix="seed_",path=tmp_path,model_name=model_name)
    make_EM_job_file(job_name="sc_EM",account=parallel_params$account,queue=parallel_params$queue,ncores=1,memk=parallel_params$memk,wall_time=parallel_params$wall_time_EM,headnode=parallel_params$headnode,log_prefix="em_",path=tmp_path,data_l_path,model_name,k)


    cmd1=paste("bsub -cwd . -J sc_seeding[1-",params$n_init_seeds,"] < ",paste(tmp_path,"/job_sc_seeding.sh",sep=""),sep="")
    print(cmd1)
  
    s=system(cmd1,intern=T)
    jobid_seeding=strsplit(s,"<|>")[[1]][2]
    cmd2=paste("bsub -cwd . -J sc_EM < ",paste(tmp_path,"/job_sc_EM.sh",sep=""), " -w \'ended(",jobid_seeding,")\'" ,sep="")
    system(cmd2)
    l=NULL
  }
  else{
    stop("No such mode ",running_mode)
  }
  return(l)
}

get_highly_variable_genes=function(ds,seeding_genes,min_umis_per_var_gene,varmean_quantile=.95,log_numis_bin_size=0.2){
  varmin_min_thresh=.1
  s=rowSums(ds)
  m=s/ncol(ds)
  v=apply(ds[s>min_umis_per_var_gene,],1,var)
  m=m[s>min_umis_per_var_gene]
  x=log10(m)
  breaks=seq(min(x),max(x),log_numis_bin_size)
  log_varmean=log2(v/m)
  binned_log_varmean=split(log_varmean,cut(x,breaks))
  mask_binned_log_varmean=sapply(binned_log_varmean,length)>0   
  min_per_bin=sapply(binned_log_varmean[mask_binned_log_varmean],min,na.rm=T)
  b=breaks[-length(breaks)]
  b=b[mask_binned_log_varmean]
  
  lo=loess(min_per_bin~b)
  var_diff=log_varmean-predict(lo,newdata =x)
  diff_thresh=quantile(var_diff,varmean_quantile,na.rm=T)

  high_var_genes=intersect(names(which(var_diff>diff_thresh&var_diff>varmin_min_thresh)),seeding_genes)
  high_var_genes=high_var_genes[!is.na(high_var_genes)]
  if (length(high_var_genes)>ncol(ds)){
    high_var_genes=head(names(sort(var_diff,decreasing = T)),ncol(ds)-1)
  }
  plot(log10(m),log2(v/m),col=ifelse(names(m)%in%high_var_genes,4,1),xlab="Log10(mean_#UMIs",ylab="log2(var/mean)",cex.axis=1.5)
  return(high_var_genes)
}




sample_cells_by_batch=function(cells,cell_to_batch,proportions,max_cells_to_sample){
  sampled_cells=list()
  l=split(cells,cell_to_batch)
  for (i in names(l)){
    n=pmin(length(l[[i]]),max_cells_to_sample)
    intervals=c(1,n*(cumsum(proportions)))
    permuted_cells=sample(l[[i]],n,replace = F)
    sampled_cells$training_set=c(sampled_cells$training_set,permuted_cells[intervals[1]:intervals[2]])
    sampled_cells$test_set=c(sampled_cells$test_set,permuted_cells[(intervals[2]+1):intervals[3]])
  }
  return(sampled_cells)
}


select_best_seeding=function(model_name,params){
  prev_seedtest_total_ll=-Inf
  best_seed=-1 

   l=list()
    for (li in 1:params$n_init_seeds){
      cur_seed=new.env()
      fn=paste("saved_clustering/",model_name,"/tmp/totll_",li,".txt",sep="")
      if (file.exists(fn)){
	l[[li]]=as.numeric(readLines(con=fn))
      } 
      else{
	l[[li]]=NA
      }
    }
  
    for (li in 1:length(l)){
      if (is.na(l[[li]])){
	    next
    }
      if (prev_seedtest_total_ll<l[[li]]){
	      prev_seedtest_total_ll=l[[li]]
	      best_seed=li
      }
    }
    ncrushed=sum(is.na(unlist(l)))
    if (ncrushed>0){
      print (paste(ncrushed, "jobs have crushed"))
    }
    seed=new.env()
    best_seed_file=paste("saved_clustering/",model_name,"/tmp/data_",best_seed,".rd",sep="")
    load(file=best_seed_file,env=seed)
    system(paste("mv",best_seed_file,paste("saved_clustering/",model_name,sep="")))
    system(paste("mv",paste("saved_clustering/",model_name,"/tmp/figures_",best_seed,".pdf",sep=""),paste("saved_clustering/",model_name,sep="")))
    system(paste("rm",paste("saved_clustering/",model_name,"/tmp/data_*.rd",sep="")))
    system(paste("rm",paste("saved_clustering/",model_name,"/tmp/figures_*.pdf",sep="")))
    system(paste("rm",paste("saved_clustering/",model_name,"/tmp/totll_*",sep="")))

    tab=table(seed$cell_to_cluster)
    print(tab[order(as.numeric(names(tab)))])

    gc()
    print("----------")

    print(paste("Best seed is #",best_seed, " (seed down-sampled to ",seed$min_umis,")",sep=""))
    return(list(models=seed$models,cell_to_cluster=seed$cell_to_cluster,alpha_noise=seed$alpha_noise,high_var_genes=seed$high_var_genes,min_umis=seed$min_umis,params=params))

}



init_clusters_kmeans=function(train_umitab,k,cell_to_batch,noise_models,params,model_name,running_mode,parallel_params=NULL,data_l_path,trace=F){
  
  clustering_genes=rownames(train_umitab)[!(rownames(train_umitab)%in%params$genes_excluded)]
  init_alpha_noise=rep(params$init_alpha_noise,ncol(noise_models))
  names(init_alpha_noise)=colnames(noise_models)
  numis_per_gene=rowSums(train_umitab)
  seeding_genes=rownames(train_umitab)[(!(rownames(train_umitab)%in%params$genes_excluded_from_seeding))&(!(rownames(train_umitab)%in%params$genes_excluded))&numis_per_gene>params$min_umis_per_seeding_gene]
  
   cell_to_cluster=c()
  
   numis_seeding=colSums(train_umitab[seeding_genes,])
  alpha_noise=init_alpha_noise
  
  if (!params$fixed_downsampling_min_umis==F){
    ds=downsample(train_umitab[seeding_genes,],params$fixed_downsampling_min_umis)

    high_var_genes=get_highly_variable_genes(ds,seeding_genes,params$min_umis_per_seeding_gene,params$seeding_varmean_quantile)
    print(paste("selected ",length(high_var_genes)," variable genes",sep=""))
    min_umis=params$fixed_downsampling_min_umis
  }
  else{
    ds=NA
    high_var_genes=NA
  }
  
  if (running_mode!="LSF_EM"){
    train_umitab_fn=paste("saved_clustering/",model_name,"/tmp/train_umitab.rd",sep="")
    save(train_umitab,file=train_umitab_fn)
    rm(train_umitab)

    l=get_seed_launcher(train_umitab_fn=train_umitab_fn,
    ds=ds,
    high_var_genes=high_var_genes,
    cell_to_batch=cell_to_batch,
    noise_models=noise_models,
    k=k,
    data_l_path=data_l_path,
    running_mode=running_mode,
    params=list(numis_seeding=numis_seeding,
      seeding_genes=seeding_genes,
      clustering_genes=clustering_genes,
      reg=params$reg,
      min_umis_per_seeding_gene=params$min_umis_per_seeding_gene,
      varmean_quantile=params$seeding_varmean_quantile,
      minimum_cluster_size=params$minimum_cluster_size,
      hard_min_umis = 150,
      min_umis_quantile_range=params$init_min_umis_quantile_range,
      init_alpha_noise=params$init_alpha_noise,
      km_reg=params$km_reg,
      redownsample=(params$fixed_downsampling_min_umis==F),
      n_init_seeds=params$n_init_seeds,
      max_n_cores=params$max_n_cores),
    model_name=model_name,
    parallel_params=parallel_params,
    trace=trace)
  }
   if (running_mode=="LSF_seeding"){
      return()
    }

  seed=select_best_seeding(model_name,params)
  cell_to_cluster=seed$cell_to_cluster
  alpha_noise=seed$alpha_noise
  
  models=update_models_debatched(train_umitab[clustering_genes,names(cell_to_cluster)],cell_to_cluster,cell_to_batch[names(cell_to_cluster)],noise_models=noise_models[clustering_genes,],alpha_noise =alpha_noise)

  return(list(models=models,alpha_noise=seed$alpha_noise,high_var_genes=seed$high_var_genes,min_umis=seed$min_umis,params=params))
}

# cluster
#
# umitab
# model_name - user defined string
# minimum_cluster_size - clusters smaller than this value will be purged
# k - number of clusters in the seeding step
# max_iter - number of likelihood update iterations
# init_expr_lim - expression range of genes participating in seeding step
# genes_excluded_from_seeding
# genes_excluded - genes to ignore
# cell_to_batch - a vector mapping from cell names to batch - currently not used


cluster=function(data_l_path,k,model_name,running_mode="local",params=list(),parallel_params=NULL,max_iter=20,load_seed=F,trace=F){
  default_params=
    list(
      minimum_cluster_size=3,
      disable_cross_validation=F,
      train_set_size=1000,
      seed_test_set_size=1000,
      test_set_size=1000,
      reg=1e-5,
      init_alpha_noise=0.01,
      genes_excluded_from_seeding=c(),
      genes_excluded=c(),
      n_init_seeds=100,
      min_umis_per_seeding_gene=40,
      insilico_gating=NULL,
      seeding_varmean_quantile=.95,
      init_min_umis_quantile_range=c(0,.2),
      fixed_downsampling_min_umis=F,
      km_reg=1,
      max_n_cores=3,
      thresh_delta_ll=2e-4,
      thresh_fraction_cells_moved=.01)
  
  params2=default_params
  if (length(params)>0){
    for (param in names(params)){
      params2[[param]]=params[[param]]
    }
  }
  params=params2

  
  load(data_l_path)
  l_cells_per_cluster=list()
  clustering_genes=rownames(data_l$umitab)[!(rownames(data_l$umitab)%in%params$genes_excluded)]
  orig_noise_models=data_l$noise_models
  noise_models=data_l$noise_models[clustering_genes,,drop=F]
  noise_models=t(t(noise_models)/colSums(noise_models))
  
  if (running_mode=="LSF_EM"||running_mode=="LSF_seeding"){
    if (is.null(parallel_params)){
      parallel_params=list()
    }
    if (is.null(parallel_params$account)){
      parallel_params$account="acc_Chessa01a"
    }
    if (is.null(parallel_params$queue)){
      parallel_params$queue="premium"
    }
    if (is.null(parallel_params$memk)){
      parallel_params$memk=20000
    }
    if (is.null(parallel_params$wall_time)){
      parallel_params$wall_time="01:00"
    }
    if (is.null(parallel_params$wall_time_EM)){
      parallel_params$wall_time_EM="04:00"
    }
    if (is.null(parallel_params$headnode)){
      parallel_params$headnode="mothra"
    }
    
  }
  
  
  if (is.na(params$fixed_downsampling_min_umis)){
    params$fixed_downsampling_min_umis=F
  }
  
  if (!dir.exists("log/")){
    dir.create("log/")
  }
  if (!dir.exists("saved_clustering/")){
    dir.create("saved_clustering/")
  }
  model_dir=paste("saved_clustering/",model_name,sep="")
  if (!dir.exists(model_dir)){
    dir.create(model_dir)
  }
  figures_dir=paste(model_dir,"/figures",sep="")
  if (!dir.exists(figures_dir)){
    dir.create(figures_dir)
  }
  seed_rd_fn=paste(model_dir,"/seed_",model_name,".rd",sep="") 
  model_rd_fn=paste(model_dir,"/model_",model_name,".rd",sep="") 
  params_rd_fn=paste(model_dir,"/tmp/params_",model_name,".rd",sep="") 
  sink(file=paste(model_dir,"/log.out",sep=""),split=TRUE)
  

  if (load_seed){
    if (!file.exists(seed_rd_fn)){
      stop("Error! " ,seed_rd_fn, " wasn't not found")
    }
    input_params=params
    load(file = seed_rd_fn)
    params$reg=input_params$reg
    sorted_umitab=insilico_sorter(data_l$umitab,params$insilico_gating,data_l$cell_to_batch)$umitab
    cell_to_batch=data_l$cell_to_batch[colnames(sorted_umitab)]
    train_umitab=data_l$umitab[,names(init_res$cell_to_cluster)]
    rm(data_l)
    gc()
  }
  else if (running_mode=="LSF_EM"){
    if (!file.exists(params_rd_fn)){
      stop("Error! " ,params_rd_fn, " wasn't not found")
    }
    load(file = params_rd_fn)
    clustering_genes=rownames(data_l$umitab)[!(rownames(data_l$umitab)%in%params$genes_excluded)]
    
    sorted_umitab=insilico_sorter(data_l$umitab,params$insilico_gating,data_l$cell_to_batch)$umitab
    cell_to_batch=data_l$cell_to_batch[colnames(sorted_umitab)]
    init_res=select_best_seeding(model_name,params)
    cell_to_cluster=init_res$cell_to_cluster
    alpha_noise=init_res$alpha_noise
    train_umitab=data_l$umitab[,names(cell_to_cluster)]
    noise_models=data_l$noise_models[clustering_genes,,drop=F]
    noise_models=t(t(noise_models)/colSums(noise_models))
    init_res$models=update_models_debatched(train_umitab[clustering_genes,],cell_to_cluster,cell_to_batch[names(cell_to_cluster)],noise_models=noise_models,alpha_noise =alpha_noise)

    save(file = seed_rd_fn,init_res,params)
    
    rm(data_l)
    gc()
  }
  else{ #running LSF_seeding or running locally
    sorted_umitab=insilico_sorter(data_l$umitab,params$insilico_gating,data_l$cell_to_batch)$umitab
    cell_to_batch=data_l$cell_to_batch[colnames(sorted_umitab)]
    rm(data_l)
    gc()

    
    batch_tab=table(cell_to_batch)
    batches=as.character(unique(cell_to_batch))
    
    if (!params$disable_cross_validation){
      ncells=c(params$train_set_size,params$test_set_size)
      max_cells_to_sample=sum(ncells)
      proportions=ncells/max_cells_to_sample
      sampled_cells=sample_cells_by_batch(colnames(sorted_umitab),cell_to_batch,proportions,max_cells_to_sample)
      training_set=sampled_cells$training_set
      test_set=sampled_cells$test_set
      all_the_rest=setdiff(colnames(sorted_umitab),c(training_set,test_set))
      cell_counts=matrix(0,length(unique(cell_to_batch)),3,dimnames = list(unique(cell_to_batch),c("training_set","test_set","unused")))
      tab=table(cell_to_batch[training_set])
      cell_counts[names(tab),"training_set"]=tab
      tab=table(cell_to_batch[test_set])
      cell_counts[names(tab),"test_set"]=tab
      tab=table(cell_to_batch[all_the_rest])
      cell_counts[names(tab),"unused"]=tab
      print(cell_counts)
    }
    else{
      training_set=colnames(sorted_umitab)
      test_set=colnames(sorted_umitab)
    }
    
    if (!dir.exists(paste("saved_clustering/",model_name,sep=""))){
      dir.create(paste("saved_clustering/",model_name,sep=""))
    }
    
    if (!dir.exists(paste("saved_clustering/",model_name,"/tmp",sep=""))){
      dir.create(paste("saved_clustering/",model_name,"/tmp",sep=""))
    }
    test_umitab=sorted_umitab[,test_set]
    save(file=paste("saved_clustering/",model_name,"/tmp/test_umitab.rd",sep=""),test_umitab)
    train_umitab=sorted_umitab[,training_set]
    rm(sorted_umitab)
    rm(test_umitab)
    gc()
    init_res=init_clusters_kmeans(train_umitab,k=k,cell_to_batch = cell_to_batch,noise_models=noise_models,params=params,model_name=model_name,running_mode=running_mode,parallel_params = parallel_params,data_l_path=data_l_path,trace=trace)
     if (running_mode=="LSF_seeding"){
        save(params,file=params_rd_fn)
	      return()
      }
      else{
        save(file = seed_rd_fn,init_res,params)
      }
  }
  
  
  load(file=paste("saved_clustering/",model_name,"/tmp/test_umitab.rd",sep=""))
  
  
  tab=table(init_res$cell_to_cluster)
  l_cells_per_cluster[[1]]=tab[order(as.numeric(names(tab)))]
  
  
  
  alpha_noise=init_res$alpha_noise
  models=init_res$models
  numis_train=colSums(train_umitab)
  params=init_res$params

  ll_train=getBatchCorrectedLikelihood(train_umitab[clustering_genes,],models,noise_models,cell_to_batch[colnames(train_umitab)],alpha_noise = alpha_noise,reg = params$reg)$ll
  cell_to_cluster=MAP(ll_train)
  prev_cell_to_cluster=cell_to_cluster
  ll_test=getBatchCorrectedLikelihood(test_umitab[clustering_genes,],models,noise_models,cell_to_batch[colnames(test_umitab)],alpha_noise = alpha_noise,reg = params$reg)$ll
  cell_to_cluster_test=MAP(ll_test)
  cur_test_total_ll=get_total_likelihood(ll_test)

  prev_train_total_ll=get_total_likelihood(ll_train)
  #ll_test_rel=ll_test-apply(ll_test,1,max)
  #ll_train_rel=ll_train-apply(ll_train,1,max)
  #ll_test_s=split(as.data.frame(ll_test_rel),c2c_test)
  #ll_train_s=split(as.data.frame(ll_train_rel),c2c_train)
  #clust="15";layout(matrix(1:2,1,2));boxplot(ll_train_s[[clust]],ylim=c(-5,0));boxplot(ll_test_s[[clust]],ylim=c(-5,0))
  
  print(paste("iter #0 LL=", round(prev_train_total_ll,digits=3)," ; test LL=",round(cur_test_total_ll,digits=4),sep=""))
  for (i in 1:max_iter){
    print("")
    print(paste("iter #",i,sep=""))
    #updating model parameters
    if (running_mode!="local"){
      max_n_cores=params$max_n_cores
    }
    else{
      max_n_cores=1
    }
    for (j in 1:10){
      prev_alpha_noise=alpha_noise
      models=update_models_debatched(train_umitab[clustering_genes,names(cell_to_cluster)],cell_to_cluster,batch=cell_to_batch[names(cell_to_cluster)],noise_models,alpha_noise)
      alpha_noise=update_alpha(umitab = train_umitab[clustering_genes,],models = models,noise_models = noise_models,cell_to_batch = cell_to_batch[colnames(train_umitab)],cell_to_cluster = cell_to_cluster,reg=params$reg,max_n_cores = max_n_cores)
      alpha_dist=sqrt(mean((alpha_noise-prev_alpha_noise)^2))
      print(round(alpha_dist,digits=6))
      if (alpha_dist<1e-4){
        break
      }
    }
    if (exists("cell_to_cluster_test")){
      alpha_noise_test=update_alpha(umitab = test_umitab[clustering_genes,],models = models,noise_models = noise_models,cell_to_batch = cell_to_batch[colnames(test_umitab)],cell_to_cluster=cell_to_cluster_test,reg=params$reg,max_n_cores = max_n_cores)
    }
    else{
      alpha_noise_test=update_alpha(umitab = test_umitab[clustering_genes,],models = models,noise_models = noise_models,cell_to_batch = cell_to_batch[colnames(test_umitab)],reg=params$reg,max_n_cores = max_n_cores)
    }
    mat_alpha=round(100*cbind(noise_train=alpha_noise,noise_test=alpha_noise_test),digits=2)
    pdf(paste(figures_dir,"/noise_test_vs_train_iter_",i,".pdf",sep=""))
    plot(mat_alpha,panel.first={grid(lty=1);abline(0,1,lty=2,col="gray")},main=paste("Iter",i," ; % Noise"))
    dev.off()
    plot(mat_alpha,panel.first={grid(lty=1);abline(0,1,lty=2,col="gray")},main=paste("Iter",i," ; % Noise"))
      # calculating log likelihoods
    
    ll_train=getBatchCorrectedLikelihood(train_umitab[clustering_genes,],models,noise_models,cell_to_batch[colnames(train_umitab)],alpha_noise = alpha_noise,reg = params$reg)$ll
    
    # Mapping cells to clusters
    cell_to_cluster=MAP(ll_train)
    
    l_cells_per_cluster[[i+1]]=rep(0,ncol(models))
    names( l_cells_per_cluster[[i+1]])=colnames(models)
    tab=table(cell_to_cluster)
    l_cells_per_cluster[[i+1]][names(tab)]=tab
    
    clusters_to_filter= l_cells_per_cluster[[i+1]]<params$minimum_cluster_size
    # filtering small clusters
    nclusters_filtered=0
    if (sum(clusters_to_filter)>0){
      print(paste("Removing ",sum(clusters_to_filter)," clusters",sep=""))
      ll_train=ll_train[,names(which(!clusters_to_filter))]
      cell_to_cluster=MAP(ll_train)
      models=models[,names(which(!clusters_to_filter))]
      nclusters_filtered=sum(clusters_to_filter)
    }
    n_cells_moved=sum(cell_to_cluster!=prev_cell_to_cluster)
    if (trace){
      table(prev_cell_to_cluster,cell_to_cluster)
    }

    cur_train_total_ll=get_total_likelihood(ll_train)
 
    ll_test=getBatchCorrectedLikelihood(test_umitab[clustering_genes,],models,noise_models,cell_to_batch[colnames(test_umitab)],alpha_noise = alpha_noise,reg = params$reg)$ll
    cell_to_cluster_test=MAP(ll_test)
    models_test=update_models_debatched(test_umitab[clustering_genes,],cell_to_cluster_test,cell_to_batch[colnames(test_umitab)],noise_models=noise_models,alpha_noise =alpha_noise_test)

    train_validation_bcl=getBatchCorrectedLikelihood(train_umitab[clustering_genes,],models_test,noise_models,cell_to_batch[colnames(train_umitab)],alpha_noise = alpha_noise,reg = params$reg)
    cell_to_cluster_validation=MAP(train_validation_bcl$ll)
    tab=table(cell_to_cluster,cell_to_cluster_validation)
    validation_score=tab[diag(nrow(tab))==1]/rowSums(tab)

    test_ll_per_cell=apply(ll_test,1,max)
    test_total_ll_per_batch=sapply(split(test_ll_per_cell,cell_to_batch[names(test_ll_per_cell)]),mean)
    cur_test_total_ll=get_total_likelihood(ll_test)
    
    train_ll_per_cell=apply(ll_train,1,max)
    train_total_ll_per_batch=sapply(split(train_ll_per_cell,cell_to_batch[names(train_ll_per_cell)]),mean)
    train_ll_diff=cur_train_total_ll-prev_train_total_ll
    
    mat_ll_per_batch=round(cbind(ll_train=train_total_ll_per_batch,ll_test=test_total_ll_per_batch),digits=3)
  
    print(cbind(mat_ll_per_batch,mat_alpha))
    print(l_cells_per_cluster[[i+1]])
    
    print(paste("k=",length(l_cells_per_cluster[[i+1]])," training LL= ",round(cur_train_total_ll,digits=4)," ; test LL=",round(cur_test_total_ll,digits=4)," ; diff=",round(train_ll_diff,digits=4)," ; clusters flitered=",sum(clusters_to_filter)," ; cells moved=",n_cells_moved,sep=""))    
    prev_cell_to_cluster=cell_to_cluster
    if (params$thresh_delta_ll>=train_ll_diff&sum(clusters_to_filter)==0&n_cells_moved<length(cell_to_cluster)*params$thresh_fraction_cells_moved){
      print("Converged!")
      break
    }
    prev_train_total_ll=cur_train_total_ll
    
  }
  cell_to_cluster=MAP(ll_train)
  # updating model with all genes (including all filtered by user)
  train_models=update_models_debatched(train_umitab[,names(cell_to_cluster)],cell_to_cluster,batch=cell_to_batch[names(cell_to_cluster)],orig_noise_models,alpha_noise,make_plots=T,figure_prefix=paste(figures_dir,"/UMI_counts_observed_noise_vs_expected",i,".pdf",sep=""))

  cell_to_cluster_test=MAP(ll_test)
  models_test=update_models_debatched(test_umitab[clustering_genes,],cell_to_cluster_test,cell_to_batch[colnames(test_umitab)],noise_models=noise_models,alpha_noise =alpha_noise_test)
  train_validation_bcl=getBatchCorrectedLikelihood(train_umitab[clustering_genes,],models_test,noise_models,cell_to_batch[colnames(train_umitab)],alpha_noise = alpha_noise,reg = params$reg)
  cell_to_cluster_validation=MAP(train_validation_bcl$ll)
  
  tab=table(cell_to_cluster,cell_to_cluster_validation)
  validation_score=round(tab[diag(nrow(tab))==1]/rowSums(tab),digits=2)
  print(validation_score)
  
 
  
  params$clustering_version="v3.0"
  print(paste("%Noise = ",paste(round(100*alpha_noise,digits=3),collapse=","),sep=""))
  pdf(paste(figures_dir,"/noise_test_vs_train_iter_final.pdf",sep=""))
  plot(mat_alpha,panel.first={grid(lty=1);abline(0,1,lty=2,col="gray")},main="% Noise")
  
  if (!params$disable_cross_validation){
    umitab=cbind(train_umitab,test_umitab)
  }
  else {
    umitab=train_umitab
  }
  ll=getBatchCorrectedLikelihood(umitab[clustering_genes,],train_models[clustering_genes,],noise_models,cell_to_batch[colnames(umitab)],alpha_noise = alpha_noise,reg = params$reg)$ll
  cell_to_cluster_test=MAP(ll)
  models=train_models
  print(paste("Saving model to ", model_rd_fn,sep=""))
  save(umitab,models,alpha_noise,noise_models,cell_to_cluster,ll,l_cells_per_cluster,cell_to_batch,params,validation_score,file=model_rd_fn)
  
  cells_per_cluster=table(cell_to_cluster)
  
  print(paste("Finished! (Detected ", length(cells_per_cluster) ," clusters)" ,sep=""))
  print(l_cells_per_cluster)
  sink()
  return(l_cells_per_cluster)
}


