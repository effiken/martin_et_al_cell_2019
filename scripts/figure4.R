source(paste(pipeline_path,"/scripts/ligand_receptor_analysis.R",sep=""))
plot_pbmc_truth=function(zlim=c(0,3)){
  clusters_to_exclude=c()
  clusters=setdiff(blood_ldm$cluster_order,clusters_to_exclude)
  ds=blood_ldm$dataset$ds[[match("2000",blood_ldm$dataset$ds_numis)]]
  ds=ds[,blood_ldm$dataset$cell_to_cluster[colnames(ds)]%in%clusters]
  l=list()
  for (i in clusters){
    mask=colnames(ds)[blood_ldm$dataset$cell_to_cluster[colnames(ds)]==i]
    if (100<length(mask)){
      l[[i]]=sample(mask,size = 100,replace = F)
    }
    else{
      l[[i]]=mask
    }
  }
  cells=unlist(l)
  ds=ds[,cells]
  gene_listp=readLines(paste(pipeline_path,"input/gene_lists/gene_list_figure_s6e.txt",sep=""))
  gene_listp_adj=adjust_gene_names(gene_listp,row.names(ds))
  gene_listp_adj=gene_listp_adj[gene_listp_adj%in%rownames(ds)]

  open_plot(path=supp_figures_path,fn = "figure_s6e",plot_type = "png",width = 3000,height = 2000)
  plot_truth_heatmap(ds,cell_to_sample =blood_ldm$dataset$cell_to_sample[colnames(ds)],
                     cell_to_cluster = blood_ldm$dataset$cell_to_cluster[colnames(ds)],
                     insamples = blood_samples,
                     ingenes = gene_listp_adj,
                     inclusts =clusters,
                     zlim=zlim,cols=colgrad_abs,plot_batch_bar = F,gene_text_cex = 2, cluster_text_cex = 2,)
  close_plot()
}



# classical monocytes in blood vs inf Macrophages in inflamed ileum
figure_4e=function(){
  
  inf_pat1=inflamed_samples[sample_to_patient[inflamed_samples]%in%pat1]
  inf_pat2=inflamed_samples[sample_to_patient[inflamed_samples]%in%pat2]
  
  tab=table(ileum_ldm$dataset$cell_to_sample,ileum_ldm$dataset$cell_to_cluster)[c(inf_pat1,inf_pat2),]
  rownames(tab)=sample_to_patient[rownames(tab)]
  tab=tab[c(pat1,pat2),]
  infMNP_gut_freq=tab[,unlist(ileum_ldm$cluster_sets$MNP$`Inf. Macrophages`)]/rowSums(tab[,unlist(ileum_ldm$cluster_sets$MNP)])
  
  tab2=table(blood_ldm$dataset$cell_to_sample,blood_ldm$dataset$cell_to_cluster)
  rownames(tab2)=sample_to_patient[rownames(tab2)]
  tab2=tab2[c(pat1,pat2),]
  
  open_plot(path=main_figures_path,fn = "figure_4e",plot_type = "pdf",width = 4,height = 4)
  par(mar=c(4,4,1,1))
  classical_mono_clusters=names(blood_ldm$clustAnnots)[blood_ldm$clustAnnots%in%c('Classical Mono-1','Classical Mono-2')]
  classical_mono_freq=rowSums(tab2[,classical_mono_clusters,drop=F])/rowSums(tab2)
  plot(infMNP_gut_freq,classical_mono_freq,col=ifelse(names(infMNP_gut_freq)%in%pat1,pat_cols[1],pat_cols[2]),panel.first = grid(lty=1),pch=20,cex=2)
  close_plot()
  print(cor.test(infMNP_gut_freq,classical_mono_freq,alternative = "less",method="spe"))
  

  
}


plot_grad=function(ldm,path,prefix,clusts,genes_to_plot=NULL,gene_list=NULL,decreasing_v=NULL,mix_clusters=F){
  u=ldm$dataset$ds[[match("1000",ldm$dataset$ds_numis)]]
  reg=.1
  score_classical_mono=colMeans(log2(reg+u[markers_classical_mono,]))
  score_mature_mac=colMeans(log2(reg+u[c(markers_mature_mac),]))
  score_inf_module=colMeans(log2(reg+u[c(markers_inf),]))
  
  ds=u[,order(score_mature_mac-score_classical_mono)]
  samps=unique(ldm$dataset$cell_to_sample[colnames(ds)])
  cell_to_cluster=ldm$dataset$cell_to_cluster[colnames(ds)]
  mask=colnames(ds)[cell_to_cluster%in%clusts&ldm$dataset$cell_to_sample[colnames(ds)]%in%samps]
  if (mix_clusters){
    cell_to_cluster[mask]=paste(clusts,collapse = ",")
    clusts=paste(clusts,collapse = ",")
  }
  reorder_genes=function(genes,decreasing){genes[order(rowMeans(ds[genes,mask]),decreasing=decreasing)]}
  
  if (is.null(genes_to_plot)){
    genes=c()
    if (is.null(decreasing_v)){
      decreasing_v=rep(T,length(gene_list))
    }
    for (i in 1:length(gene_list)){
      genes=c(genes,reorder_genes(gene_list[[i]],decreasing=decreasing_v[i]))
    }
  }
  else {
    genes=genes_to_plot
  }
  
  sample_cols1=sample_cols[1:length(samps)]
  names(sample_cols1)=samps
  
  open_plot(path=path,paste(prefix,"gradient_truth",sep="_"),"png",width=1500,height=1200)
  plot_truth_heatmap(ds,cell_to_sample = ldm$dataset$cell_to_sample[colnames(ds)],cell_to_cluster =  cell_to_cluster,insamples = samps,inclusts =clusts,ingenes = genes,zlim=c(0,3),sample_cols = sample_cols1,plot_batch_bar =T,gene_text_cex = .7,cols = colgrad_abs)
  close_plot()
  
  open_plot(path=path,fn = paste(prefix,"gradient_truth_genes",sep="_"),plot_type = "pdf",width =10,height = 1)
  par(mar=c(1,1,1,1),bg=0)
  plot.new()
  mtext(text = genes,side=1,line = -2,las=2,at = seq(0,1,l=length(genes)),cex = .7)
  close_plot()
  
  return(genes)
}

plot_macs_grads=function(){

  global_markers_inf_macs=c("CXCL3","IER3","SOD2","C15orf48","G0S2","PLAUR","NFKBIA")
  markers_mature_mac=c("MMP12","CCL18","C1QA","C1QB","C1QC","ACP5","CTSD",'KLF6','FYB','PLA2G7','CTSB','ATP6V1G1','PSAP','CD63','CD68','JUND','NFKBIZ','TXN','LGALS3',"RNASE1","APOE","CAPG","GLUL","CTSL","CNBP","ZFP36L2")
  markers_classical_mono=c("CCR2","S100A8","S100A9","S100A12","VCAN","CD52","OSM",'S100A4','NBEAL1','ATP1B3','PHACTR1','GK','TIMP1','YWHAZ','SLC2A3','CORO1A','C5AR1','IL1B',"ACSL1",'DUSP2',"AREG","EREG","THBS1","FOSL2","IRF8","HIF1A","SUB1","ANXA1","LRRFIP1")
  markers_inf=c("CCL3","CCL3L3","CCL4L2","BCL2A1","PTGS2","CXCL2","CXCL8","CCL20","ZFAND5","IFITM3","OLR1","CHMP1B","CD83","MAP3K8","SAMSN1","CEBPB")
  markers_res_macs=c("LGMN","SLC40A1","SEPP1","CEBPD")
  
  gene_list=list(classical_mono=markers_classical_mono,mature_mac=markers_mature_mac,global_inf=global_markers_inf_macs)

  genes=plot_grad(ileum_ldm,main_figures_path,"4d-top",clusts=unlist(ileum_ldm$cluster_sets$MNP$`Inf. Macrophages`),gene_list=gene_list,decreasing_v = c(T,F,F))
  a=plot_grad(blood_ldm,main_figures_path,"4d-bottom",genes_to_plot = genes,clusts=unlist(blood_ldm$cluster_sets$`Classical Mono-1`))
 
  gene_list=list(classical_mono=markers_classical_mono,mature_mac=markers_mature_mac,global_inf=global_markers_inf_macs,res_macs=markers_res_macs)
  genes=plot_grad(ileum_ldm,supp_figures_path,"6d",clusts=unlist(ileum_ldm$cluster_sets$MNP$`Resident macrophages`),gene_list=gene_list,decreasing_v = c(T,F,F,F),mix_clusters = T)
}

make_figure4=function(){
  plot_pbmc_truth()
  figure_4e()
  main_ligand_receptor()
}
