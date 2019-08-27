library(viridis)
library(RColorBrewer)
library(matrixStats)
library(RColorBrewer)
library(gplots)
library(Matrix)
library(Matrix.utils)
library(seriation)
library(scDissector)



make_martin_et_al_figures=function(pipeline_path,download_data=T,make_figures=T){
  
  
  stats_list<<-list()
  pipeline_path<<-pipeline_path
  path_scripts<<-"scripts/"
  source(paste(pipeline_path,"/scripts/graphics.r",sep=""))
  
  if (download_data){
    download_files(pipeline_path)
  }
  
 
  read_sc_data(pipeline_path =pipeline_path)
  create_global_vars(pipeline_path =pipeline_path)
  make_common_stats()
  
  if (make_figures){
    dir.create(paste(pipeline_path,"output",sep="/"))
    dir.create(paste(pipeline_path,"output/main_figures",sep="/"))
    dir.create(paste(pipeline_path,"output/supp_figures",sep="/"))
    dir.create(paste(pipeline_path,"output/tables",sep="/"))
    source(paste(pipeline_path,"scripts/figure1.R",sep="/"))
    source(paste(pipeline_path,"scripts/figure2.R",sep="/"))
    source(paste(pipeline_path,"scripts/figure3.R",sep="/"))
    source(paste(pipeline_path,"scripts/figure4.R",sep="/"))
    source(paste(pipeline_path,"scripts/figure5.R",sep="/"))
    
    message("Making figure 1")
    make_figure1()
    
    message("Making figure 2")
    make_figure2()
    
    message("Making figure 3")
    make_figure3()
    
    message("Making figure 4")
    make_figure4()
    
    message("Making figure 5")
    make_figure5()
    table_s4() 
  }
}










download_files=function(pipeline_path){
  dir.create(paste(pipeline_path,"input","clustered_scRNA_data",sep="/"),showWarnings = F)
  dir.create(paste(pipeline_path,"input","external_cohorts_data",sep="/"),showWarnings = F)
  dir.create(paste(pipeline_path,"input","DE",sep="/"),showWarnings = F)
  
  message("Downloading file 1/4")
  download.file("https://www.dropbox.com/s/v0fr7pejzfynmbr/model_and_samples_ileum.rd?dl=1",destfile = paste(pipeline_path,"input","clustered_scRNA_data","model_and_samples_ileum.rd",sep="/"))
  message("Downloading file 2/4")
  download.file("https://www.dropbox.com/s/79wdka2vtgv7lrr/model_and_samples_blood.rd?dl=1",destfile = paste(pipeline_path,"input","clustered_scRNA_data","model_and_samples_blood.rd",sep="/"))
  message("Downloading file 3/4")
  download.file("https://www.dropbox.com/s/4emrgk7g2xb9o0m/DE_inf_pat1_vs_pat2_total.csv?dl=1",destfile = paste(pipeline_path,"input","DE","DE_inf_pat1_vs_pat2_total.csv",sep="/"))
  message("Downloading file 3/4")
  download.file("https://www.dropbox.com/s/zvuoq8hkxt3cwll/bulk_data.rd?dl=1",destfile = paste(pipeline_path,"input","external_cohorts_data","bulk_data.rd",sep="/"))
}

get_edges=function(l,name=NA){
  if (is.null(l)){
    return()
  }
  else{
    node=c()
    parent=c()
    for (i in 1:length(names(l))){
      if (!is.na(name)&!is.null(names(l))){
        parent[i]=name
        node[i]=names(l)[i]
      }
    }
    m=data.frame(node=node,parent=parent)
    for (i in length(l):1){
      if (!is.null(names(l))){
        m=rbind(get_edges(l[[i]],names(l)[i]),m)
      }
    }
    return(m)
  }
}



randomly_select_cells=function(ldm,nrandom_cells_per_sample_choices){
  randomly_selected_cells=list()
  for (ds_i in 1:length(ldm$dataset$ds_numis)){
    randomly_selected_cells[[ds_i]]<-list()
    for (nrandom_cells in nrandom_cells_per_sample_choices){
      randomly_selected_cells[[ds_i]][[nrandom_cells]]=c()
      for (sampi in ldm$dataset$samples){
        if (nrandom_cells=="All"||pmax(0,as.numeric(nrandom_cells),na.rm=T)>=ncol(ldm$dataset$ds[[ds_i]])){
          randomly_selected_cells[[ds_i]][[nrandom_cells]]<-c(randomly_selected_cells[[ds_i]][[nrandom_cells]],colnames(ldm$ld$ds[[ds_i]]))
        }
        else{
          randomly_selected_cells[[ds_i]][[nrandom_cells]]<-c(randomly_selected_cells[[ds_i]][[nrandom_cells]],sample(colnames(ldm$dataset$ds[[ds_i]]),size=as.numeric(nrandom_cells),replace=F))
        }
      }
    }
  }
  
  return(randomly_selected_cells)
}

read_sc_data=function(pipeline_path){#,ileum_model_version="model_combined_081518_50_5e6_kmreg01"){
  sample_tab<<-read.csv(paste(pipeline_path,"/input/tables/sample_index.csv",sep=""),stringsAsFactors = F)
  rownames(sample_tab)=sample_tab$index
  sample_tab$Exclude_scRNAseq[is.na(sample_tab$Exclude_scRNAseq)]=F
  
  sample_tab<<-sample_tab[!sample_tab$Exclude_scRNAseq,]
  mask_ileum=as.character(sample_tab$index)[(!sample_tab$Exclude_scRNAseq)&sample_tab$tissue=="ILEUM"&(sample_tab$CHEMISTRY=="V1"|sample_tab$CHEMISTRY=="V2")]
  mask_blood=as.character(sample_tab$index)[(!sample_tab$Exclude_scRNAseq)&sample_tab$tissue=="BLOOD"&(sample_tab$CHEMISTRY=="V1"|sample_tab$CHEMISTRY=="V2")]
  
  ileum_rd_filename=paste(pipeline_path,"input","clustered_scRNA_data","model_and_samples_ileum.rd",sep="/")
  blood_rd_filename=paste(pipeline_path,"input","clustered_scRNA_data","model_and_samples_blood.rd",sep="/")
  #if (project){
  #  samples_fn=paste(clustering_data_path,"/data_",sample_tab$index,".rd",sep="")
  #  names(samples_fn)=sample_tab$index
  #  ileum_ldm<<-load_dataset_and_model(model_fn = paste(clustering_data_path,ileum_model_version,".rd",sep=""),sample_fns = samples_fn[mask_ileum],min_umis = 800,excluded_clusters = c("9","23","44"))
  #  blood_ldm<<-load_dataset_and_model(model_fn = paste(clustering_data_path,"model_blood_111718_20_5e6_kmreg01.rd",sep=""),sample_fns = samples_fn[mask_blood],min_umis = 800,excluded_clusters = c("12"))
  #  if (!dir.exists(paste(pipeline_path,"input","clustered_scRNA_data",sep="/"))){
  #    dir.create(paste(pipeline_path,"input","clustered_scRNA_data",sep="/"))
  #  }
  #  save(file=ileum_rd_filename,ileum_ldm)
  #  save(file=blood_rd_filename,blood_ldm)
  #}
  #else{
    env=new.env()
    load(ileum_rd_filename,envir = env)
    ileum_ldm<<-env$ileum_ldm
    load(blood_rd_filename,envir = env)
    blood_ldm<<-env$blood_ldm
  #}
  
  ileum_ldm$dataset$randomly_selected_cells<<-randomly_select_cells(ileum_ldm,c("200","1000","5000"))
  blood_ldm$dataset$randomly_selected_cells<<-randomly_select_cells(blood_ldm,c("200","1000","5000"))
  
}

create_gene_expression_matrices=function(){
  get_one_gene_expression_matrix=function(ldm,samples){
    counts=apply(ldm$dataset$counts[samples,,],2:3,sum)-apply(ldm$dataset$noise_counts[samples,,],2:3,sum)
    s=pmax(aggregate.Matrix(t(counts),groupings = cluster_to_subtype[colnames(counts)],fun = "sum"),0)
    m=s/rowSums(s)
    return(t(m))
  }
  l=list()
  l[["inf_pat1"]]=get_one_gene_expression_matrix(ileum_ldm,inflamed_samples_v2[sample_to_patient[inflamed_samples_v2]%in%pat1])
  l[["inf_pat2"]]=get_one_gene_expression_matrix(ileum_ldm,inflamed_samples_v2[sample_to_patient[inflamed_samples_v2]%in%pat2])
  l[["uninf_pat1"]]=get_one_gene_expression_matrix(ileum_ldm,uninflamed_samples_v2[sample_to_patient[uninflamed_samples_v2]%in%pat1])
  l[["uninf_pat2"]]=get_one_gene_expression_matrix(ileum_ldm,uninflamed_samples_v2[sample_to_patient[uninflamed_samples_v2]%in%pat2])
  return(l)
}

create_global_vars=function(pipeline_path){
  
  main_figures_path<<-paste(pipeline_path,"output/main_figures/",sep="/")
  supp_figures_path<<-paste(pipeline_path,"output/supp_figures/",sep="/")
  tables_path<<-paste(pipeline_path,"output/tables/",sep="/")
  
  samples_by_status<<-split(sample_tab$index,sample_tab$status) 
  sample_to_patient<<-gsub("CD Biop","bp ",sample_tab$Patient.ID)
  names(sample_to_patient)<<-sample_tab$index
  sample_names<<-names(sample_to_patient)
  patientid_to_gridid<<-unique(sample_tab$GRID.ID)
  names(patientid_to_gridid)<<-sample_to_patient[match(unique(sample_tab$GRID.ID),sample_tab$GRID.ID)]
  
  
  condition_to_col<<-c(4,2)
  names(condition_to_col)<<-c("Uninflamed","Inflamed")
  
  inflamed_samples<<-as.character(sample_tab$index[sample_tab$status=="Inflamed"])
  uninflamed_samples<<-as.character(sample_tab$index[sample_tab$status=="Uninflamed"])
  control_samples<<-as.character(sample_tab$index[grep("CONTROL",sample_tab$status)])
  blood_samples<<-sample_tab$index[sample_tab$status=="PBMC"]
  
  inflamed_samples_filtered<<-setdiff(inflamed_samples,inflamed_samples[sample_to_patient[inflamed_samples]%in%c("rp 6","rp 16")])
  uninflamed_samples_filtered<<-setdiff(uninflamed_samples,uninflamed_samples[sample_to_patient[uninflamed_samples]%in%c("rp 6","rp 16")])
  
  inflamed_samples_v2<<-as.character(sample_tab$index[sample_tab$status=="Inflamed"&sample_tab$CHEMISTRY=="V2"])
  uninflamed_samples_v2<<-as.character(sample_tab$index[sample_tab$status=="Uninflamed"&sample_tab$CHEMISTRY=="V2"])
  
  inflamed_samples_v2_filtered<<-intersect(inflamed_samples_v2,inflamed_samples_filtered)
  uninflamed_samples_v2_filtered<<-intersect(uninflamed_samples_v2,uninflamed_samples_filtered)
  
  sample_to_condition<<-sample_tab$status
  names(sample_to_condition)<<-sample_names
  
  cols_inf_uninf<<-brewer.pal(3,"Set1")[1:2]
  
  cluster_to_cluster_set<<-rep(names(ileum_ldm$cluster_sets),sapply(ileum_ldm$cluster_sets,function(x){length(unlist(x))}))
  
  edges=get_edges(ileum_ldm$cluster_sets)
  cluster_to_subtype<<-as.character(edges[match(colnames(ileum_ldm$model$models),edges[,1]),"parent"])
  names(cluster_to_subtype)<<-colnames(ileum_ldm$model$models)
  
  
  cluster_to_subtype1<<-rep(unlist(sapply(ileum_ldm$cluster_sets, names)),unlist(sapply(ileum_ldm$cluster_sets,sapply,length)))
  names(cluster_to_subtype1)<<-unlist(ileum_ldm$cluster_sets)
  cluster_to_subtype1<<-cluster_to_subtype1[!cluster_to_subtype1%in%unlist(ileum_ldm$cluster_sets$`Not good`)]
  
  names(cluster_to_cluster_set)<<-unlist(ileum_ldm$cluster_sets)
  cluster_to_cluster_set<<-cluster_to_cluster_set[names(cluster_to_subtype)]
  
  cluster_to_cluster_set_with_pdc<<-cluster_to_cluster_set
  cluster_to_cluster_set_with_pdc[names(ileum_ldm$clustAnnots)[ileum_ldm$clustAnnots=="pDC"]]<<-"pDC"
  
  
  
  cluster_to_broader_subtype<<-cluster_to_subtype1
  cluster_to_broader_subtype[cluster_to_broader_subtype%in%c("Inf. Macrophages","Resident macrophages")]<<-"Macrophages"
  cluster_to_broader_subtype[cluster_to_broader_subtype%in%c("DC1","DC2","Activated DC","moDC")]<<-"DC"
  cluster_to_broader_subtype[cluster_to_broader_subtype%in%c("Activated fibroblasts")]<<-"Fibroblasts"
  cluster_to_broader_subtype[cluster_to_broader_subtype%in%c("ACKR1+ endothelial cells","CD36+ endothelial cells")]<<-"Endothelial"
  cluster_to_broader_subtype[cluster_to_broader_subtype%in%c("IgA plasma cells","IgG plasma cells","IgM plasma cells")]<<-"Plasma"
  cluster_to_broader_subtype[unlist(ileum_ldm$cluster_sets$T)]<<-"T"
  
  
  
  stats<<-list()
  
  brew_cols<<-rep(brewer.pal(12,"Set3"),10)
  colgrad<<-c(colorRampPalette(c("white",colors()[378],"orange", "tomato","mediumorchid4"))(100))
  sample_cols<<-rep(paste("#",read.table(paste(pipeline_path,"input/colors/sample_colors.txt",sep=""),stringsAsFactors = F)[,1],sep=""),10)
  
  colgrad_rel_file=paste(pipeline_path,"input/colors/colors_brewer_RdBu.txt",sep="")
  colgrad_rel<<-read.table(colgrad_rel_file,stringsAsFactors=F)[,1]
  
  celltypes<<-c("T","MNP","Plasma","B","ILC","Stromal","Mast","pDC")
  celltypes_cols2<<-c(brewer.pal(12,"Paired")[c(1,3,5,7,9)],colors()[40],"lightgray","yellow")
  celltypes_cols1<<-c(brewer.pal(12,"Paired")[c(2,4,6,8,10,12)],"darkgray","gold1")
  names(celltypes_cols1)<<-celltypes
  names(celltypes_cols2)<<-celltypes
  
  
  alt_cols<<-rep(c("#e6194b","#3cb44b","#ffe119","#0082c8","#f58231","#911eb4","#46f0f0","#f032e6","#d2f53c","#fabebe","#008080","#e6beff","#aa6e28","#fffac8","#800000","#aaffc3","#808000","#ffd8b1","#000080","#808080",
                   "#FFFFFF","#000000"),10)
  pat_cols<<-c("#2c7fb8","#31a354")

  pat1<<-paste("rp",c(7,8,5,11,12))
  pat2<<-paste("rp",c(10,13,14,15))
  
  gsc=load_gene_symbol_converters()
  gene_symbol_old2new<<-gsc$old2new
  gene_symbol_new2old<<-gsc$new2old
  rm(gsc)
  
  colgrad_abs<<-c(1,viridis(99))
}

make_common_stats=function(){
  cluster_condition_counts<<-table(sample_to_condition[ileum_ldm$dataset$cell_to_sample],ileum_ldm$dataset$cell_to_cluster)
}




table_s4=function(){
  gene_lists_fn=c(
    fig_1c="input/gene_lists/gene_list_figure_1c.txt",
    fig_1f="input/gene_lists/gene_list_figure_1f.txt",
    fig_2a="input/gene_lists/gene_list_figure_2a.txt",
    fig_2b="input/gene_lists/gene_list_figure_2b.txt",
    fig_2c="input/gene_lists/gene_list_figure_2c.txt",
    fig_2d="input/gene_lists/gene_list_figure_2d.txt",
    fig_2f="input/gene_lists/gene_list_figure_2f.txt",
    fig_2g="input/gene_lists/gene_list_figure_2g.txt",
    fig_2h="input/gene_lists/gene_list_figure_2h.txt",
    fig_s2b="input/gene_lists/gene_list_figure_s2b.txt",
    fig_s3a="input/gene_lists/gene_list_figure_s3a.txt",
    fig_s3e="input/gene_lists/gene_list_figure_s3e.txt",
    fig_s3f="input/gene_lists/gene_list_figure_s3f.txt",
    fig_s6e="input/gene_lists/gene_list_figure_s6e.txt")

  gene_lists=sapply(paste(pipeline_path,gene_lists_fn,sep="/"),function(x){strsplit(readLines(x),",|, ")})
  names(gene_lists)=names(gene_lists_fn)
  
  max_ngenes=max(sapply(gene_lists,length))
  m=sapply(gene_lists,function(x){c(x,rep("",max_ngenes-length(x)))})
  write.csv(file=paste(pipeline_path,"output/tables/table_s4.csv",sep="/"),m,quote = F,row.names = F)
  
 
}


normalize_by_clusterset_frequency2=function(lm,selected_samples){
  freqs=get_freqs(lm,selected_samples)
  
  
  return(sapply(lm$cluster_sets[names(lm$cluster_sets)!="Not good"],norm_one_clusterset))
}
