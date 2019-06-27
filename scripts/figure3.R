

plot_subtype_freqs_per_patient=function(freq_norm,celltype,plot_legend=T,ordc=NULL,cex.names=1,cex.axis=1,ordp=c(pat1,pat2),cols=NULL,samp_labels=NULL){
  if (is.null(ordc)){
    ordc=order(colMeans(freq_norm[[celltype]][sample_to_patient[rownames(freq_norm[[celltype]])]%in%pat1,])/colMeans(freq_norm[[celltype]][sample_to_patient[rownames(freq_norm[[celltype]])]%in%pat2,]),decreasing=T)
  }
  if (!is.null(samp_labels)){
    rownames(freq_norm[[celltype]])=samp_labels
  }
  m=freq_norm[[celltype]][ordp,ordc]

  if (is.null(cols)){
    cols=alt_cols[1:ncol(m)]
  }else{
    cols=cols[ordc]
  }
  if (plot_legend){
    
    barplot(t(m),col=cols,las=2,cex.names = cex.names,cex.axis = cex.axis)
    plot.new()
    legend("bottomleft",legend=rev(colnames(m)),col=rev(cols),pch=15,cex=.9,bty = "n")
  }
  else{
    barplot(t(m),col=cols,las=2,cex.names=cex.names,cex.axis=1)
  }
  return(ordc)
}



plot_all_subtype_freqs_per_patient=function(){
  cols_l=list(
    MNP=c('Inf. Macrophages'=3, 'Resident macrophages'=7, 'pDC'=6, 'moDC'=5, 'Activated DC'=4, 'DC2'=2, 'DC1'=1),
    Plasma=c('IgA plasma cells'=1, 'Plasmablasts'=4, 'IgM plasma cells'=3, 'IgG plasma cells'=2),
    "T"=c('Tregs'=8, 'TFH-like'=7, 'Naive/CM T cells'=5, 'Highly activated T cells'=4, 'Type 1 cytokines Trm'=9, 'Type 3 cytokines Trm'=2, 'Cytotoxic T cells'=3, 'Cytokines low Trm'=1, 'T (gd)'=6),
    Stromal=c('Smooth muscle cells'=8, 'Pericytes'=7, 'Lymphatics'=6, 'Fibroblasts'=5, 'Enteric neurons'=4, 'CD36+ endothelial cells'=3, 'Activated fibroblasts'=2, 'ACKR1+ endothelial cells'=1)
  )
  
  
  cluster_sets=ileum_ldm$cluster_sets
  freq_norm_inf=normalize_by_clusterset_frequency(ileum_ldm$dataset,inflamed_samples_filtered,cluster_sets = cluster_sets,pool_subtype = T,reg = 0)
  freq_norm_uninf=normalize_by_clusterset_frequency(ileum_ldm$dataset,uninflamed_samples_filtered,cluster_sets =cluster_sets,pool_subtype = T,reg = 0)
  cormat=(cor(t(do.call(cbind,freq_norm_inf))))
  ordp=sample_to_patient[rownames(cormat)[get_order(seriate(as.dist(1-cormat)))]]
  celltypes=c("MNP","Stromal","T","Plasma")
  open_plot(path = main_figures_path,fn="figure_3a",plot_type = "pdf",width=6,height=length(celltypes)*2.5)
  layout(matrix(1:(length(celltypes)*3),length(celltypes),3,byrow = T)[,c(3,1,2)])
  par(mar=c(4,3,1,1))
  for (cluster_set_name in celltypes){
    print(cluster_set_name)
    cluster_set=cluster_sets[cluster_set_name]
    if (length(unlist(cluster_set))>1){

      ordc=plot_subtype_freqs_per_patient(freq_norm_inf,cluster_set_name,plot_legend = F,ordp=ordp,cols=alt_cols[cols_l[[cluster_set_name]]],samp_labels=sample_to_patient[inflamed_samples_filtered])
      plot_subtype_freqs_per_patient(freq_norm_uninf,cluster_set_name,plot_legend = T,ordc=ordc,ordp=ordp,cols=alt_cols[cols_l[[cluster_set_name]]],samp_labels=sample_to_patient[uninflamed_samples_filtered])
   
    }
  }
  close_plot()
  return(ordp)
}



geometric_mean_panel=function(ordp=paste("rp",c(5,7,8,12,11,14,15,10,13))){
  module_cell_types=c("IgG plasma cells","Activated DC","Inf. Macrophages","Activated fibroblasts","ACKR1+ endothelial cells","Highly activated T cells")
 # module_cell_types=c("IgG plasma cells","mature DC","Inf. Macrophages","Activated fibroblasts","ACKR1+ endothelial cells")
  reg=0.001
  freqs_inf=do.call(cbind,normalize_by_clusterset_frequency(ileum_ldm$dataset,samples = inflamed_samples_filtered,pool_subtype = T,cluster_sets = ileum_ldm$cluster_sets,reg = reg))
  score_inf=exp(rowMeans(log(reg+freqs_inf[,module_cell_types])))
  names(score_inf)=sample_to_patient[names(score_inf)]
  freqs_uninf=do.call(cbind,normalize_by_clusterset_frequency(ileum_ldm$dataset,samples = uninflamed_samples_filtered,pool_subtype = T,cluster_sets = ileum_ldm$cluster_sets,reg = reg))
  score_uninf=exp(rowMeans(log(reg+freqs_uninf[,module_cell_types])))
  names(score_uninf)=sample_to_patient[names(score_uninf)]
  
  ttest_res=t.test(score_inf[pat1],score_inf[pat2])
  print(ttest_res$statistic)
  print(ttest_res$p.value)
  open_plot(path = main_figures_path,"figure_3b",plot_type = "pdf",5,3)
  layout(matrix(1:2,1,2),widths=c(4.5,3))
  par(mar=c(4,5,2,2))
  barplot(score_inf[ordp],las=2,ylim=c(0,.4),border=F,ylab="Geometric mean(freqs)",col="#3F6D9B",main="Inflamed",cex.axis = .7,cex.names = .8)
  box()
  par(mar=c(4,1,2,1))
  barplot(score_uninf[ordp],las=2,ylim=c(0,.4),border=F,ylab="Geometric mean(freqs)",col="#3F6D9B",main="Uninflamed",axes=T,cex.axis = .7,cex.names = .8)
  box()
  close_plot()
}



gene_per_compartment=function(genes,put_legend=T){
  get_one=function(clusters,gene="TNF"){
    return(apply(counts[,gene,unlist(clusters),drop=F],1,sum)/sum(counts[,,unlist(clusters),drop=F]))
  }
  counts=ileum_ldm$dataset$counts-ileum_ldm$dataset$noise_counts
  for (gene in genes){
    print(gene)
    m=t(sapply(ileum_ldm$cluster_sets,get_one,gene=gene))[intersect(names(ileum_ldm$cluster_sets),names(celltypes_cols1)),inflamed_samples]
    colnames(m)=sample_to_patient[colnames(m)]
    open_plot("output/supp_figures/",paste(gene,"_per_compartment",sep=""),6,6,plot_type = "pdf")
    print(paste(gene,"_per_compartment",sep=""))
    par(mar=c(5,8,1,1))
    barplot(m[,c(pat1,pat2)],col=celltypes_cols1[rownames(m)],las=2)
    mtext(paste("Total ",gene),side = 2,line = 6)
    if (put_legend){
      legend("topright",legend=names(celltypes_cols1[rownames(m)]),col=celltypes_cols1[rownames(m)],pch=15)
    }
    close_plot()
  }
}



make_figure3=function(){
   
  ordp=plot_all_subtype_freqs_per_patient()
  geometric_mean_panel(ordp=ordp)
 
#  gene_per_compartment(c("OSM","IL1B","CXCL8","IL6","CCL3","CCL4","IL11","CSF3","IL1A","CSF2","IL22","IFNG","CCL2","IL17A","CXCL2","CXCL3","CXCL5","CXCL6","CXCL1","CXCL9","CXCL10","CXCL11"))
  }