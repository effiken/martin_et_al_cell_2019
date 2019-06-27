open_plot=function(path="figures/",fn,plot_type="png",width=NA,height=NA,bg="white"){
  if (plot_type=="png"){
    if (is.na(width)){
      width=400
    }
    if (is.na(height)){
      height=400
    }
    png(paste(path,"/",fn,".png",sep=""),width,height,bg=bg)
  }
  else if(plot_type=="pdf"){
    if (is.na(width)){
      width=4
    }
    if (is.na(height)){
      height=4
    }
    pdf(paste(path,"/",fn,".pdf",sep=""),width,height,bg=bg)
  }
  par(cex.axis=1.5)
}

close_plot=function(){
  dev.off()
}


write_data=function(path="tables/",fn,tab,round_nums=T){
  if (round_nums){
    tab=round(tab,digits=2)
  }
  write.csv(tab,paste(path,fn,".csv",sep=""),quote=F)
}



plot_truth=function(){
  ds=ileum_ldm$dataset$ds[[match("2000",ileum_ldm$dataset$ds_numis)]]
  plot_truth_heatmap(ds,cell_to_sample = ileum_ldm$dataset$cell_to_sample[colnames(ds)],cell_to_cluster = ileum_ldm$dataset$cell_to_sample[colnames(ds)],cols = c("black",colgrad_abs[-1]),insamples = insamples,ingenes = ingenes,inclusts = inclusts,zlim = zlim,sample_cols = sample_cols,showSeparatorBars = T)
}


make_truth_plot=function(clusters=NULL,gene_list_fn="",path="output/main_figures/",figure_fn="truth",zlim=c(0,3),ncell_per_cluster=100,downsampling_version="2000",use_default_cluster_order=T,gene_list=NULL,gene_text_cex = 3,cluster_text_cex=2,seperatorBars_lwd=3,reorder_genes=F,samples=NULL,gene_list_path="")
{
  if (is.null(clusters)){
    clusters=unique(ileum_ldm$dataset$cell_to_cluster)
  }
  if (use_default_cluster_order){
    clusters=intersect(ileum_ldm$cluster_order,clusters)
  }
  ds=ileum_ldm$dataset$ds[[match(downsampling_version,ileum_ldm$dataset$ds_numis)]]
  ds=ds[,ileum_ldm$dataset$cell_to_cluster[colnames(ds)]%in%clusters]
  if (!is.null(samples)){
    ds=ds[,ileum_ldm$dataset$cell_to_sample[colnames(ds)]%in%samples]
  }
  
  l=list()
  for (i in clusters){
    mask=colnames(ds)[ileum_ldm$dataset$cell_to_cluster[colnames(ds)]==i]
    if (ncell_per_cluster<length(mask)){
      l[[i]]=sample(mask,size = ncell_per_cluster,replace = F)
    }
    else{
      l[[i]]=mask
    }
  }
  cells=unlist(l)
  ds=ds[,cells]
  if (is.null(gene_list)){
    gene_list2=readLines(paste(gene_list_path,gene_list_fn,sep="/"))
    gene_list2_adj=adjust_gene_names(gene_list2,row.names(ds))
    gene_list2_adj=gene_list2_adj[gene_list2_adj%in%rownames(ds)]
  }
  else{
    gene_list2_adj=gene_list
  }
  
  if (reorder_genes){
    gene_list2_adj=gene_list2_adj[order(apply(ileum_ldm$model$models[gene_list2_adj,ileum_ldm$cluster_order],1,which.max))]
  }
  
  # should be automatic!
  open_plot(path=path,fn = figure_fn,plot_type = "png",width = 3000,height = 2000)
  plot_truth_heatmap(ds,cell_to_sample =ileum_ldm$dataset$cell_to_sample[colnames(ds)],
                     cell_to_cluster = ileum_ldm$dataset$cell_to_cluster[colnames(ds)],
                     insamples = c(inflamed_samples,uninflamed_samples),
                     ingenes = gene_list2_adj,
                     inclusts =clusters,
                     zlim=zlim,cols=colgrad_abs,plot_batch_bar = F,lower_mar=15,seperatorBars_lwd=seperatorBars_lwd,gene_text_cex = gene_text_cex, cluster_text_cex = cluster_text_cex)
  close_plot()
  open_plot(path=path,fn = paste(figure_fn,"genes",sep="_"),plot_type = "pdf",width = 8*length(gene_list2_adj)/100,height = 1)
  par(mar=c(1,1,1,1))
  plot.new()
  mtext(text = gene_list2_adj,side=1,line = -2,las=2,at = seq(0,1,l=length(gene_list2_adj)),cex = .4)
  close_plot()
}


error.bar <- function(x,y, upper, lower=upper, length=0.1,horiz=T,...){
  if( length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  if (horiz){
    arrows(y+upper,x, y-lower,x, angle=90, code=3, length=length, ...)
  }
  else{
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
  }
}



plot_subtype_freqs=function(freq_norm,celltype,plot_legend=T,ordc=NULL,cex.names=1,cex.axis=1,ordp=c(pat1,pat2),cols=NULL,samp_labels=NULL){
  if (is.null(ordc)){
    ordc=order(colMeans(freq_norm[[celltype]][sample_to_patient[rownames(freq_norm[[celltype]])]%in%pat1,])/colMeans(freq_norm[[celltype]][sample_to_patient[rownames(freq_norm[[celltype]])]%in%pat2,]),decreasing=T)
  }
  if (!is.null(samp_labels)){
    rownames(freq_norm[[celltype]])=samp_labels
  }
  m=freq_norm[[celltype]][ordp,ordc]
  #cols=brewer.pal(length(ordc),"Set3")
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

