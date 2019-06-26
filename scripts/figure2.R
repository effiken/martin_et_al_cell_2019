make_figure_correlation_between_subtypes=function(path,prefix,lm,selected_samples,show_celltype_colors=F){
  l=sapply(lm$cluster_sets,unlist)
  v=rep(names(l),sapply(l,length))
  names(v)=unlist(l)
  freqs=do.call(cbind,normalize_by_clusterset_frequency(lm$dataset,selected_samples,cluster_set=lm$cluster_sets,reg=0.001))
  cormat=cor(log2(freqs))

  ord=rev(get_order(seriate(as.dist(1-cormat),method = "GW")))
  large_margin=8
  small_margin=.5
  open_plot(path = path,fn=prefix,plot_type = "pdf",width = 5,height = 5)
  if (show_celltype_colors){
    subtype=sapply(strsplit(cluster_to_subtype," \\| "),head,1)
    celltype=sapply(strsplit(cluster_to_subtype," \\| "),tail,1)
    subtype_to_celltype=celltype[match(unique(subtype),subtype)]
    names(subtype_to_celltype)=unique(subtype)
    celltype_ind=match(subtype_to_celltype[colnames(cormat)[ord]],names(celltypes_cols1))
    layout(matrix(c(1,3,4,2),2,2),heights=c(.7,20),widths = c(20,.7))
    par(mar=c(.2,large_margin,.2,small_margin))
    image(matrix(celltype_ind,,1),col=celltypes_cols1,axes=F,breaks=0.5+c(0:length(celltypes_cols1)))
    par(mar=c(large_margin,.2,small_margin,.2))
    image(matrix(celltype_ind,1,),col=celltypes_cols1,axes=F,breaks=0.5+c(0:length(celltypes_cols1)))
  }
  par(mar=c(large_margin,large_margin,small_margin,small_margin))

  image(cormat[ord,ord],col=colorRampPalette(c("blue","white","red"))(100),axes=F,breaks=seq(-1,1,l=101))
  box(lwd=2)
  mtext(rownames(cormat)[ord],1,at = seq(0,1,l=ncol(cormat)),las=2,cex=.6,line =.5)
  mtext(rownames(cormat)[ord],2,at = seq(0,1,l=ncol(cormat)),las=2,cex=.6,line =.5)
  close_plot()
  
}


freq_sig_test=function(m,design,nperm=1e5){
  mean1=colMeans(m[design==1,])
  mean2=colMeans(m[design==2,])
  obs_fc=abs(log2((1e-5+mean1)/(1e-5+mean2)))
  v=rep(0,ncol(m))
  for (i in 1:nperm){
    samped_design=sample(design,size = length(design),replace = F)
    v=v+obs_fc<abs(log2((1e-5+colMeans(m[samped_design==1,]))/(1e-5+colMeans(m[samped_design==2,]))))
  }
  return(v/nperm)
}

inf_uninf_freq_barplot=function(freqs,compartment){
  m=freqs[[compartment]]
  inf=m[inflamed_samples_filtered,]
  uninf=m[uninflamed_samples_filtered,]
  browser()
  #freq_sig_test(m[c(inflamed_samples_filtered,uninflamed_samples_filtered),],c(inflamed_samples_filtered,uninflamed_samples_filtered)%in%inflamed_samples_filtered+1)
  
  ord=order(log2(colMeans(inf)/colMeans(uninf)),decreasing=T)
  mean_inf=colMeans(inf)
  mean_uninf=colMeans(uninf)
  se_inf=apply(inf,2,sd)/sqrt(nrow(inf))
  se_uninf=apply(uninf,2,sd)/sqrt(nrow(uninf))
  df=data.frame(subtype=factor(c(colnames(inf),colnames(uninf)),levels=colnames(m)[ord]),frequency=c(mean_inf,mean_uninf),se=c(se_inf,se_uninf),status=c(rep("Inflamed",ncol(inf)),rep("Uninflamed",ncol(uninf))))
 # open_plot(main_figures_path,fn=paste("figure2_",compartment,"_inf_uninf_freq",sep=""),plot_type = "pdf",width = 4,height = 4)
  #par(mar=c(2,3,1,1),lwd = .5 )
  #barplot(cbind(uninf,inf)[c("Resident macrophages","moDC","DC2","DC1","pDC","Inf. Macrophages","mature DC"),],col=c("#984EA3","#4DAF4A","#4DAF4A","#E41A1C","#FF7F00","#FFFF33","#377EB8"),cex.axis =.7,space = .5)  
  #open_plot(main_figures_path,fn=paste("figure2_",compartment,"_inf_uninf_freq",sep=""),plot_type = "pdf",width = 4,height = 4)
  p<-ggplot(df, aes(x=subtype, y=frequency, fill=status,width=.7)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=frequency-se, ymax=frequency+se),width=.2,position=position_dodge(.7)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=15),plot.title = element_text(hjust = 0.5),axis.title=element_text(size=14)) +
    scale_fill_manual("legend", values = c("Inflamed" = cols_inf_uninf[1], "Uninflamed" = cols_inf_uninf[2])) +
    labs(title=compartment)
  p
  ggsave(filename = paste("output/main_figures/figure2_",compartment,"_inf_uninf_freq.pdf",sep=""),p)
  pvalues=p.adjust(sapply(1:ncol(inf),function(i){wilcox.test(uninf[,i],inf[,i])$p.value}))
  names(pvalues)=colnames(inf)
  print(pvalues)
  #return(pvalues)
}


make_inf_uninf_freq_barplot=function(){
  freqs=normalize_by_clusterset_frequency(ileum_ldm$dataset,samples = c(inflamed_samples_filtered,uninflamed_samples_filtered),cluster_sets = ileum_ldm$cluster_sets,pool_subtype = T,reg = 0)
  inf_uninf_freq_barplot(freqs,"MNP")
  inf_uninf_freq_barplot(freqs,"Plasma")
  inf_uninf_freq_barplot(freqs,"T")
  inf_uninf_freq_barplot(freqs,"Stromal")
  m=do.call(cbind,freqs)
  print(freq_sig_test(m,rownames(m)%in%inflamed_samples_filtered+1))
}

tregs_de=function(){
  mask=bulk$design$dataset=="RISK"
  genes_to_exclude=c("JCHAIN","IGLC2","IGKC","FABP6","PHGR1","DEFA5","DEFA6","REG1A","FABP1","IGHG3","IGHG2","IGHG1","XIST")
  de_tregs=read.csv("output/tables/DE_patterns/DE_inf_pat1_vs_pat2_Tregs.T.csv",row.names = 1)
  gene_mask=rownames(de_tregs)[de_tregs[,"log2_FC"]<(-1)&de_tregs[,"adj.p.value"]<(0.01)]
  gene_mask=setdiff(gene_mask,genes_to_exclude)
  z=t(bulk$z[intersect(gene_mask,rownames(bulk$z)),mask])
  z=z[rowSums(!is.na(z))>1,colSums(!is.na(z))>1]
  cor(z,bulk$scores[rownames(z),"projected",drop=F],use="comp")
}


expression_profiles_corr=function(){
  max_exprssion=1e-3
  reg=1e-6
  thresh=7
  mask=log2((reg+rowMaxs(ileum_ldm$model$models))/(reg+rowMins(ileum_ldm$model$models)))>thresh&rowMaxs(ileum_ldm$model$models)<max_exprssion
  sum(mask)
  m=ileum_ldm$model$models
  m=log2((reg+m)/(reg+rowMeans(m)))
  cormat=cor(ileum_ldm$model$models[mask,],method="pearson")
  #ord=hclust(as.dist(1-cormat))$order
  ord=ileum_ldm$cluster_order
  open_plot("output/supp_figures/",fn = "figure_s1h_expression_profiles_corr",plot_type = "pdf",nrow(cormat)/6,height =ncol(cormat)/6)
  par(mar=c(7,7,1,1))
  image(cormat[ord,ord],col=colorRampPalette(c("blue","white","red"))(100),breaks=c(-1,seq(-1,1,l=99),1),axes=F)
  mtext(text = ileum_ldm$clustAnnots[ord],side = 1,at = seq(0,1,l=nrow(cormat)),las=2,cex=.7)
  mtext(text = ileum_ldm$clustAnnots[ord],side = 2,at = seq(0,1,l=ncol(cormat)),las=2,cex=.7)
  close_plot()
  
}


make_tf_truth=function(){
  tf_list_s=paste(unique(read.csv("input/Gene lists/Human_TFs_Lambert_et_al_Cell_2018.csv",stringsAsFactors = F)[,1]),collapse=",")
  tf_list_adj=adjust_gene_names(tf_list_s,row.names(ileum_ldm$model$models))
  tf_list_adj=unique(tf_list_adj[tf_list_adj%in%rownames(ileum_ldm$model$models)])
  reg=1e-5
  tf_list=tf_list_adj[rowMaxs(ileum_ldm$model$models[tf_list_adj,])>1e-4&log2((reg+rowMaxs(ileum_ldm$model$models[tf_list_adj,]))/(reg+apply(ileum_ldm$model$models[tf_list_adj,],1,quantile,.7)))>1]
  
  make_truth_plot(clusters = ileum_ldm$cluster_order[ileum_ldm$cluster_order%in%colnames(ileum_ldm$model$models)],gene_list =tf_list,path="output/supp_figures/",figure_fn ="truth_tfs",reorder_genes = T,gene_text_cex = 1.5,seperatorBars_lwd=1,zlim = c(0,3),ncell_per_cluster = 100)
}


# correlation_between_subtypes_axcross inflamed samples
figure_2j=function(){
  make_figure_correlation_between_subtypes(main_figures_path,"figure_2j",ileum_ldm,inflamed_samples_filtered)
 
  open_plot(main_figures_path,fn="figure_2j_colorscale",plot_type = "png",width = 200,height =50)
  par(mar=c(2,1,0,1))
  image(matrix(seq(-1,1,l=100),,1),col=colorRampPalette(c("blue","white","red"))(100),axes=F,breaks=seq(-1,1,l=101))
  axis(1,at=seq(0,1,l=5),labels=F,lwd = 2)
  close_plot()
  
}

freq_barplot=function(fn,cell_type,subtypes=NA){
  inf=colMeans(normalize_by_clusterset_frequency(ileum_ldm$dataset,samples = inflamed_samples_filtered,pool_subtype = T,cluster_sets = ileum_ldm$cluster_sets,reg = 0)[[cell_type]])
  uninf=colMeans(normalize_by_clusterset_frequency(ileum_ldm$dataset,samples = uninflamed_samples_filtered,pool_subtype = T,cluster_sets = ileum_ldm$cluster_sets,reg = 0)[[cell_type]])
  if (is.na(subtypes)){
    subtypes=names(inf)[order(inf/uninf)]
  }
  open_plot(main_figures_path,fn=fn,plot_type = "pdf",width = 3,height = 3)
  par(mar=c(2,3,1,1),lwd = .5 )
  barplot(cbind(uninf,inf)[subtypes,],col=brewer.pal(length(subtypes),"Set3"),cex.axis =.7,space = .5,ylim=c(0,1))  
  close_plot()
  open_plot(main_figures_path,fn=paste(fn,"_legend",sep=""),plot_type = "pdf",width = 3,height = 3)
  par(mar=c(2,3,1,1),lwd = .5 )
  plot.new()
  legend("bottomright",legend=rev(subtypes),pch=15,col=rev(brewer.pal(length(subtypes),"Set3")))  
  close_plot()
}




make_avg_heatmaps=function(zlim=c(-2,3)){
    read_genes=function(gene_list_fn){  
      gene_list2=readLines(paste(pipeline_path,"/input/gene_lists/",gene_list_fn,sep=""))
      gene_list2_adj=adjust_gene_names(gene_list2,row.names(ileum_ldm$dataset$umitab))
      gene_list2_adj=gene_list2_adj[gene_list2_adj%in%rownames(ileum_ldm$dataset$umitab)]
      return(gene_list2_adj)
    }
    
    
    
    #avg_heatmap_MNP
    genes=read_genes("gene_list_figure_2a.txt")
    genes=intersect(genes,rownames(ileum_ldm$model$models))
    clusters=unlist(c(rev(ileum_ldm$cluster_sets$MNP$`Resident macrophages`),ileum_ldm$cluster_sets$MNP$`Inf. Macrophages`,ileum_ldm$cluster_sets$MNP$moDC,ileum_ldm$cluster_sets$MNP$DC2,ileum_ldm$cluster_sets$MNP$DC1,ileum_ldm$cluster_sets$MNP$`Activated DC`))
    open_plot(main_figures_path,fn='figure_2a',plot_type = "pdf",width = 5,height = 2)
    par(mar=c(4,5,1,1))
    plot_avg_heatmap(ileum_ldm$model$models[genes,clusters],genes=genes,gene.cols = 1,clusters = clusters,clusters_text = "",annots = cluster_to_subtype1[clusters],zlim=c(-2,2),main_title = "",Relative_or_Absolute = "Relative",colgrad = colgrad_rel,reg=1e-6,cex.genes = .6,cex.clusters = .5,line.genes = .1)
    close_plot()
    
    #avg_heatmap_MNP cytokines
    genes=read_genes("gene_list_figure_2d.txt")
    genes=intersect(genes,rownames(ileum_ldm$model$models))
    clusters=unlist(c(ileum_ldm$cluster_sets$MNP$DC1,ileum_ldm$cluster_sets$MNP$moDC,ileum_ldm$cluster_sets$MNP$DC2,ileum_ldm$cluster_sets$MNP$`Resident macrophages`,ileum_ldm$cluster_sets$MNP$`Inf. Macrophages`,ileum_ldm$cluster_sets$MNP$`Activated DC`))
    open_plot(main_figures_path,fn='figure_2d',plot_type = "pdf",width = 4.2,height = 1.5)
    par(mar=c(4,5,1,1))
    plot_avg_heatmap(ileum_ldm$model$models[genes,clusters],genes=genes,gene.cols = 1,clusters = clusters,clusters_text = "",annots = cluster_to_subtype1[clusters],zlim=zlim,main_title = "",Relative_or_Absolute = "Relative",colgrad = colgrad_rel,reg=1e-6,cex.genes = .6,cex.clusters = .5,line.genes = .1)
    close_plot()
    
    # avg_heatmap_t_cells
    genes=read_genes("gene_list_figure_2f.txt")
    genes=intersect(genes,rownames(ileum_ldm$model$models))
    clusters=c("19","29","22","17","43","34","47","36","42","48","46","30","26","4","21","38","28","35")
    open_plot(main_figures_path,fn='figure_2f',plot_type = "pdf",width = 3.5,height = 2.5)
    par(mar=c(4,5,1,1))
    plot_avg_heatmap(ileum_ldm$model$models[genes,clusters],genes=genes,gene.cols = 1,clusters = clusters,clusters_text = "",annots = cluster_to_subtype1[clusters],zlim=zlim,main_title = "",Relative_or_Absolute = "Relative",colgrad = colgrad_rel,reg=1e-6,cex.genes = .7,cex.clusters = .5,line.genes = .1)
    close_plot()
    
    # avg_heatmap_stromal_cells
    genes=read_genes("gene_list_figure_2h.txt")
    genes=intersect(genes,rownames(ileum_ldm$model$models))
    clusters=unlist(c(ileum_ldm$cluster_sets$Stromal$`CD36+ endothelial cells`,ileum_ldm$cluster_sets$Stromal$`ACKR1+ endothelial cells`,ileum_ldm$cluster_sets$Stromal$Lymphatics,ileum_ldm$cluster_sets$Stromal$`Smooth muscle cells`,ileum_ldm$cluster_sets$Stromal$Fibroblasts,ileum_ldm$cluster_sets$Stromal$`Activated fibroblasts`))
    open_plot(main_figures_path,fn='figure_2h',plot_type = "pdf",width = 3.5,height = 1.5)
    par(mar=c(4,5,1,1))
    plot_avg_heatmap(ileum_ldm$model$models[genes,clusters],genes=genes,gene.cols = 1,clusters = clusters,clusters_text = "",annots = cluster_to_subtype1[clusters],zlim=zlim,main_title = "",Relative_or_Absolute = "Relative",colgrad = colgrad_rel,reg=1e-6,cex.genes = .7,cex.clusters = .5,line.genes = .1)
    close_plot()

}



numbers_figures2=function(){
  l=list()
  
  ### Inf.MAcs inf uninf fold change
  mask_inf=ileum_ldm$dataset$cell_to_sample%in%inflamed_samples_filtered
  mask_uninf=ileum_ldm$dataset$cell_to_sample%in%uninflamed_samples_filtered
  cluster_sets=ileum_ldm$cluster_sets
  freq_norm_inf=(normalize_by_clusterset_frequency(ileum_ldm$dataset,inflamed_samples_filtered,pool_subtype = T,cluster_sets = ileum_ldm$cluster_sets,reg = 1e-3)$MNP)
  freq_norm_uninf=(normalize_by_clusterset_frequency(ileum_ldm$dataset,uninflamed_samples_filtered,pool_subtype = T,cluster_sets = ileum_ldm$cluster_sets,reg =1e-3)$MNP)
  
  mean_inf=colMeans(freq_norm_inf)
  mean_uninf=colMeans(freq_norm_uninf)
  fc=mean_inf/mean_uninf
  
  l[["Inf. Macrophages Freq Inflamed/Uninflamed"]]=fc["Inf. Macrophages"]
  l[["Activated DC Freq Inflamed/Uninflamed"]]=fc["Activated DC"]
  return(l)
}



# make_average_gene_expression_table
table_s9=function(){
  counts=pmax(apply(ileum_ldm$dataset$counts[c(uninflamed_samples_v2_filtered,inflamed_samples_v2_filtered),,],2:3,sum)-apply(ileum_ldm$dataset$noise_counts[c(uninflamed_samples_v2_filtered,inflamed_samples_v2_filtered),,],2:3,sum),0)
#  counts=apply(ileum_ldm$dataset$counts[samples,,],2:3,sum)
  counts=counts[,names(cluster_to_subtype)]
  s=aggregate.Matrix(t(counts),groupings = cluster_to_subtype[colnames(counts)],fun = "sum")
  m=t(s/rowSums(s))
  m=m[apply(m,1,max)>1e-6,]
  lm=log2(2^-20+(m))
  # table s9 - gene_expression_per_subtype_log2_average
  write.csv(file=paste(pipeline_path,"output/tables/table_s9.csv",sep="/"),as.matrix(lm),quote = F,row.names = T)
}



make_figure2=function(){
  make_avg_heatmaps()
  
  freq_barplot("figure_2i_MNP_inf_uninf_freq",cell_type = "MNP")
  freq_barplot("figure_2i_T_inf_uninf_freq",cell_type = "T")
  freq_barplot("figure_2i_Plasma_inf_uninf_freq",cell_type = "Plasma")
  freq_barplot("figure_2i_Stromal_inf_uninf_freq",cell_type = "Stromal")
  
  #ILC_inf_uninf_freq
  freq_barplot("figure_S4g",cell_type = "ILC")
  
  gene_list_path=paste(pipeline_path,"input/gene_lists/",sep="/")
  #T-cells 
  make_truth_plot(clusters = unlist(ileum_ldm$cluster_sets$T),gene_list_fn ="gene_list_figure_s4e.txt" ,figure_fn ="figure_s4e",gene_list_path=gene_list_path,path=supp_figures_path)
  #DCs
  make_truth_plot(clusters = c(unlist(ileum_ldm$cluster_sets$MNP$moDC),unlist(ileum_ldm$cluster_sets$MNP$DC2),unlist(ileum_ldm$cluster_sets$MNP$DC1),unlist(ileum_ldm$cluster_sets$MNP$`Activated DC`),unlist(ileum_ldm$cluster_sets$MNP$pDC)),gene_list_fn ="gene_list_figure_2c.txt" ,figure_fn ="figure_2b",use_default_cluster_order=F,ncell_per_cluster = 500,,gene_list_path=gene_list_path,path=main_figures_path)
  #Plasma and B cells
  make_truth_plot(clusters = c(unlist(ileum_ldm$cluster_sets$B),unlist(ileum_ldm$cluster_sets$Plasma)),gene_list_fn ="gene_list_figure_s4a.txt" ,figure_fn ="figure_s4a",gene_list_path=gene_list_path,path=supp_figures_path)
  #ILC
  make_truth_plot(clusters = unlist(ileum_ldm$cluster_sets$ILC),gene_list_fn ="gene_list_figure_s4f.txt" ,figure_fn ="figure_s4f",ncell_per_cluster = 300,gene_text_cex=2,gene_list_path=gene_list_path,path=supp_figures_path)
  #Macrophages
  make_truth_plot(clusters = c(unlist(ileum_ldm$cluster_sets$MNP$`Resident macrophages`),unlist(ileum_ldm$cluster_sets$MNP$`Inf. Macrophages`)),gene_list_fn ="gene_list_figure_2b.txt" ,figure_fn ="figure_2b",ncell_per_cluster = 300,gene_list_path=gene_list_path,path=main_figures_path)
  #Stroma
  make_truth_plot(clusters = unlist(ileum_ldm$cluster_sets$Stromal),gene_list_fn ="gene_list_figure_2g.txt" ,figure_fn ="figure_2g",gene_list_path=gene_list_path,path=main_figures_path)
  print(numbers_figures2())
  table_s9()
  figure_2j()
}


