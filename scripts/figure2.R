make_figure_correlation_between_subtypes=function(path,prefix,lm,selected_samples,show_celltype_colors=F){
  l=sapply(lm$cluster_sets,unlist)
  v=rep(names(l),sapply(l,length))
  names(v)=unlist(l)
  freqs=do.call(cbind,normalize_by_clusterset_frequency(lm$dataset$cell_to_cluster,lm$dataset$cell_to_sample,selected_samples,cluster_set=lm$cluster_sets,reg=0.001))
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




inf_uninf_freq_barplot=function(freqs,compartment,fig,path=supp_figures_path){
  m=freqs[[compartment]]
  inf=m[inflamed_samples_filtered,]
  uninf=m[uninflamed_samples_filtered,]
  
  ord=order(log2(colMeans(inf)/colMeans(uninf)),decreasing=T)
  mean_inf=colMeans(inf)
  mean_uninf=colMeans(uninf)
  se_inf=apply(inf,2,sd)/sqrt(nrow(inf))
  se_uninf=apply(uninf,2,sd)/sqrt(nrow(uninf))
  df=data.frame(subtype=factor(c(colnames(inf),colnames(uninf)),levels=colnames(m)[ord]),frequency=c(mean_inf,mean_uninf),se=c(se_inf,se_uninf),status=c(rep("Inflamed",ncol(inf)),rep("Uninflamed",ncol(uninf))))
  p<-ggplot(df, aes(x=subtype, y=frequency, fill=status,width=.7)) + 
    geom_bar(position=position_dodge(), stat="identity") +
    geom_errorbar(aes(ymin=frequency-se, ymax=frequency+se),width=.2,position=position_dodge(.7)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5,size=15),plot.title = element_text(hjust = 0.5),axis.title=element_text(size=14)) +
    scale_fill_manual("legend", values = c("Inflamed" = cols_inf_uninf[1], "Uninflamed" = cols_inf_uninf[2])) +
    labs(title=compartment)
  p
  ggsave(filename = paste(path,"figure_",fig,".pdf",sep=""),p)

}


make_inf_uninf_freq_barplot=function(){
  freqs=normalize_by_clusterset_frequency(ileum_ldm$dataset$cell_to_cluster,ileum_ldm$dataset$cell_to_sample,samples = c(inflamed_samples_filtered,uninflamed_samples_filtered),cluster_sets = ileum_ldm$cluster_sets,pool_subtype = T,reg = 0)
  inf_uninf_freq_barplot(freqs,"MNP",fig="s4a")
  inf_uninf_freq_barplot(freqs,"Plasma",fig="s4d")
  inf_uninf_freq_barplot(freqs,"T",fig="s4b")
  inf_uninf_freq_barplot(freqs,"Stromal",fig="s4c")
}




figure_2e=function(ncell_per_cluster=300){
  
  tcell_list=read.table(paste(pipeline_path,"input/gene_lists/gene_list_figure_2e.txt",sep="/"),row.names = 1,header=T,stringsAsFactors = F,sep="\t")
  genes=strsplit(tcell_list[,1],",| ,|, ")
  names(genes)=rownames(tcell_list)
  ds=ileum_ldm$dataset$ds[[3]]
  byvec=rep(0,nrow(ds))
  names(byvec)=rownames(ds)
  for (i in 1:length(genes)){
    byvec[genes[[i]]]=names(genes)[i]
  }
  clusters=unlist(ileum_ldm$cluster_sets$T)
  tcells=intersect(colnames(ds),names(ileum_ldm$dataset$cell_to_cluster)[ileum_ldm$dataset$cell_to_cluster%in%clusters])
  l=list()
  
  for (i in clusters){
    mask=tcells[ileum_ldm$dataset$cell_to_cluster[tcells]==i]
    if (ncell_per_cluster<length(mask)){
      l[[i]]=sample(mask,size = ncell_per_cluster,replace = F)
    }
    else{
      l[[i]]=mask
    }
  }
  tcells=unlist(l)
  
  
  scores=t(aggregate.Matrix(ds[,tcells],groupings = byvec,fun="mean"))[,-1]/1000
  scores=as.matrix(scores[order(match(ileum_ldm$dataset$cell_to_cluster[tcells],ileum_ldm$cluster_order)),])
  normed_scores=scores/colSums(scores)
  reg=1e-6
  normed_scores=t(log2((t(normed_scores)+reg)/(apply(normed_scores,2,function(x){quantile(x[x>0],.75)})+reg)))
  normed_scores=normed_scores[,get_order(seriate(as.dist(1-cor(normed_scores)),method="GW"))]
  #normed_scores=normed_scores[,c("Cycling","CD8_Cytotoxic","Tissue_resident_memory","Regulatory","Naive_CM")]]
  open_plot(path=main_figures_path,fn = "figure_2e",plot_type = "png",width = 2400,height = 500)
  par(mar=c(1,1,3,20))
  zlim=c(-1,1)
  image(normed_scores,col=colgrad_abs,,axes=F,breaks=c(-100,seq(zlim[1],zlim[2],l=99),100))
  mtext(colnames(normed_scores),4,at = seq(0,1,l=ncol(normed_scores)),las=2,cex=2,line=0.5)
  mtext(intersect(ileum_ldm$cluster_order,clusters),3,at = seq(par()$usr[1],par()$usr[2],l=length(clusters)+1)[-1]-.5/(length(clusters)+1),cex=2,line=.5)
  abline(v=seq(par()$usr[1],par()$usr[2],l=length(clusters)+1),col="gray",lwd=4)
  close_plot()
  
  open_plot(path = main_figures_path,fn="figure_2e_colorscale",plot_type = "png",width = 200,height =50)
  par(mar=c(2,1,0,1))
  image(matrix(seq(zlim[1],zlim[2],l=100),,1),col=colgrad_abs,axes=F,breaks=seq(zlim[1],zlim[2],l=101))
  axis(1,at=seq(0,1,l=5),labels = seq(zlim[1],zlim[2],l=5),lwd = 2)
  close_plot()
  
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

freq_barplot=function(fn,cell_type,subtypes=NA,path=main_figures_path){
  inf=colMeans(normalize_by_clusterset_frequency(ileum_ldm$dataset$cell_to_cluster,ileum_ldm$dataset$cell_to_sample,samples = inflamed_samples_filtered,pool_subtype = T,cluster_sets = ileum_ldm$cluster_sets,reg = 0)[[cell_type]])
  uninf=colMeans(normalize_by_clusterset_frequency(ileum_ldm$dataset$cell_to_cluster,ileum_ldm$dataset$cell_to_sample,samples = uninflamed_samples_filtered,pool_subtype = T,cluster_sets = ileum_ldm$cluster_sets,reg = 0)[[cell_type]])
  if (is.na(subtypes)){
    subtypes=names(inf)[order(inf/uninf)]
  }
  open_plot(path,fn=fn,plot_type = "pdf",width = 3,height = 3)
  par(mar=c(2,3,1,1),lwd = .5 )
  barplot(cbind(uninf,inf)[subtypes,],col=brewer.pal(length(subtypes),"Set3"),cex.axis =.7,space = .5,ylim=c(0,1))  
  close_plot()
  open_plot(path,fn=paste(fn,"_legend",sep=""),plot_type = "pdf",width = 3,height = 3)
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
    
    # avg_heatmap_stromal_cells
    genes=read_genes("gene_list_figure_s2b.txt")
    genes=intersect(genes,rownames(ileum_ldm$model$models))
    clusters=c("41","15","2","8","31","24","32","40")
    open_plot(supp_figures_path,fn='figure_s2b',plot_type = "pdf",width = 3.5,height = 1.5)
    par(mar=c(4,5,1,1))
    plot_avg_heatmap(ileum_ldm$model$models[genes,clusters],genes=genes,gene.cols = 1,clusters = clusters,clusters_text = "",annots = cluster_to_subtype1[clusters],zlim=zlim,main_title = "",Relative_or_Absolute = "Relative",colgrad = colgrad_rel,reg=1e-6,cex.genes = .7,cex.clusters = .5,line.genes = .1)
    close_plot()

}

# Gene modules could be slightly differenet between runs due to random seeding
run_gene_modules_analysis=function(){
  genes_to_exclude=c(grep("RP",rownames(ileum_ldm$dataset$umitab),val=T),grep("MT-",rownames(ileum_ldm$dataset$umitab),val=T))
  cormat=gene_cor_analysis(ileum_ldm,"2000",.1,30,genes_to_exclude=genes_to_exclude,clusters=unlist(ileum_ldm$cluster_sets$T),samples = setdiff(c(inflamed_samples,uninflamed_samples),c("68","69")))
  parse_modules(ileum_ldm,cormat,"2000","ileum_T",nmods =50,figures_path=supp_figures_path,tables_path = tables_path,zlim = c(-.5,.5))
  
  genes_to_exclude=c(grep("RP",rownames(ileum_ldm$dataset$umitab),val=T),grep("MT-",rownames(ileum_ldm$dataset$umitab),val=T))
  cormat=gene_cor_analysis(ileum_ldm,"2000",.1,30,genes_to_exclude=genes_to_exclude,clusters=unlist(ileum_ldm$cluster_sets$MNP),samples = setdiff(c(inflamed_samples,uninflamed_samples),c("68","69")))
  parse_modules(ileum_ldm,cormat,"2000","ileum_MNP",nmods =50,figures_path=supp_figures_path,tables_path = tables_path,zlim = c(-.5,.5))
  
}


numbers_figures2=function(){

  
  ### Inf.MAcs inf uninf fold change
  mask_inf=ileum_ldm$dataset$cell_to_sample%in%inflamed_samples_filtered
  mask_uninf=ileum_ldm$dataset$cell_to_sample%in%uninflamed_samples_filtered
  cluster_sets=ileum_ldm$cluster_sets
  freq_norm_inf=(normalize_by_clusterset_frequency(ileum_ldm$dataset$cell_to_cluster,ileum_ldm$dataset$cell_to_sample,inflamed_samples_filtered,pool_subtype = T,cluster_sets = ileum_ldm$cluster_sets,reg = 1e-3)$MNP)
  freq_norm_uninf=(normalize_by_clusterset_frequency(ileum_ldm$dataset$cell_to_cluster,ileum_ldm$dataset$cell_to_sample,uninflamed_samples_filtered,pool_subtype = T,cluster_sets = ileum_ldm$cluster_sets,reg =1e-3)$MNP)
  
  mean_inf=colMeans(freq_norm_inf)
  mean_uninf=colMeans(freq_norm_uninf)
  fc=mean_inf/mean_uninf
  
  stats$"Inf. Macrophages Freq Inflamed/Uninflamed"<<-fc["Inf. Macrophages"]
  stats$"Activated DC Freq Inflamed/Uninflamed"<<-fc["Activated DC"]

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
  
  freq_barplot("figure_2i_MNP",cell_type = "MNP")
  freq_barplot("figure_2i_T",cell_type = "T")
  freq_barplot("figure_2i_Plasma",cell_type = "Plasma")
  freq_barplot("figure_2i_Stromal",cell_type = "Stromal")
  
  #ILC_inf_uninf_freq
  freq_barplot("figure_s3g",cell_type = "ILC",path=supp_figures_path)
  
  gene_list_path=paste(pipeline_path,"input/gene_lists/",sep="/")
  #T-cells 
  make_truth_plot(clusters = unlist(ileum_ldm$cluster_sets$T),gene_list_fn ="gene_list_figure_s3e.txt" ,figure_fn ="figure_s3e",gene_list_path=gene_list_path,path=supp_figures_path)
  #DCs
  make_truth_plot(clusters = c(unlist(ileum_ldm$cluster_sets$MNP$moDC),unlist(ileum_ldm$cluster_sets$MNP$DC2),unlist(ileum_ldm$cluster_sets$MNP$DC1),unlist(ileum_ldm$cluster_sets$MNP$`Activated DC`),unlist(ileum_ldm$cluster_sets$MNP$pDC)),gene_list_fn ="gene_list_figure_2c.txt" ,figure_fn ="figure_2b",use_default_cluster_order=F,ncell_per_cluster = 500,,gene_list_path=gene_list_path,path=main_figures_path)
  #Plasma and B cells
  make_truth_plot(clusters = c(unlist(ileum_ldm$cluster_sets$B),unlist(ileum_ldm$cluster_sets$Plasma)),gene_list_fn ="gene_list_figure_s3a.txt" ,figure_fn ="figure_s3a",gene_list_path=gene_list_path,path=supp_figures_path)
  #ILC
  make_truth_plot(clusters = unlist(ileum_ldm$cluster_sets$ILC),gene_list_fn ="gene_list_figure_s3f.txt" ,figure_fn ="figure_s3f",ncell_per_cluster = 300,gene_text_cex=2,gene_list_path=gene_list_path,path=supp_figures_path)
  #Macrophages
  make_truth_plot(clusters = c(unlist(ileum_ldm$cluster_sets$MNP$`Resident macrophages`),unlist(ileum_ldm$cluster_sets$MNP$`Inf. Macrophages`)),gene_list_fn ="gene_list_figure_2b.txt" ,figure_fn ="figure_2b",ncell_per_cluster = 300,gene_list_path=gene_list_path,path=main_figures_path)
  #Stroma
  make_truth_plot(clusters = unlist(ileum_ldm$cluster_sets$Stromal),gene_list_fn ="gene_list_figure_2g.txt" ,figure_fn ="figure_2g",gene_list_path=gene_list_path,path=main_figures_path)

  table_s9()
  figure_2e()
  figure_2j()

  #s5a-d
  make_inf_uninf_freq_barplot()
  run_gene_modules_analysis()
}


