

# Generating figures 1c and 1f
make_truth_plots=function(clusters_to_exclude=c("23","44"),zlim=c(0,3))
{
  clusters=setdiff(ileum_ldm$cluster_order,clusters_to_exclude)
  ds=ileum_ldm$dataset$ds[[match("2000",ileum_ldm$dataset$ds_numis)]]
  ds=ds[,ileum_ldm$dataset$cell_to_cluster[colnames(ds)]%in%clusters]

  l=list()
  for (i in clusters){
    mask=colnames(ds)[ileum_ldm$dataset$cell_to_cluster[colnames(ds)]==i]
    if (100<length(mask)){
      l[[i]]=sample(mask,size = 100,replace = F)
    }
    else{
      l[[i]]=mask
    }
  }
  cells=unlist(l)
  ds=ds[,cells]
  gene_list2=readLines(paste(pipeline_path,"input/gene_lists/gene_list_figure_1f.txt",sep=""))
  gene_list2_adj=adjust_gene_names(gene_list2,row.names(ds))
  gene_list2_adj=gene_list2_adj[gene_list2_adj%in%rownames(ds)]

  main_figures_path=paste(pipeline_path,"output/main_figures/",sep="/")
  open_plot(path=main_figures_path,fn = "truth_color_legend",plot_type = "png",width =500,height = 100)
  par(mar=c(1,1,1,1))
  image(matrix(1:100,,1),col=colgrad_abs,pch=20,cex=3,axes=F)
  close_plot()
  
  
  open_plot(path=main_figures_path,fn = "figure_1f",plot_type = "png",width = 3000,height = 2000)
  plot_truth_heatmap(ds,cell_to_sample =ileum_ldm$dataset$cell_to_sample[colnames(ds)],
                   cell_to_cluster = ileum_ldm$dataset$cell_to_cluster[colnames(ds)],
                   insamples = c(inflamed_samples,uninflamed_samples),
                   ingenes = gene_list2_adj,
                   inclusts =clusters,
                   zlim=zlim,cols=colgrad_abs,plot_batch_bar = F,gene_text_cex = 1.5, cluster_text_cex = 2)
  close_plot()
  
  gene_list1=readLines(paste(pipeline_path,"/input/gene_lists/gene_list_figure_1c.txt",sep=""))
  gene_list1_adj=adjust_gene_names(gene_list1,row.names(ds))
  gene_list1_adj=gene_list1_adj[gene_list1_adj%in%rownames(ds)]
  open_plot(path=main_figures_path,fn = "figure_1c",plot_type = "png",width = 1000,height = 2000)
  plot_truth_heatmap(ds,cell_to_sample =ileum_ldm$dataset$cell_to_sample[colnames(ds)],
                     cell_to_cluster = ileum_ldm$dataset$cell_to_cluster[colnames(ds)],
                     insamples = c(inflamed_samples,uninflamed_samples),
                     ingenes = gene_list1_adj,
                     inclusts =clusters,
                     zlim=zlim,cols=colgrad_abs,plot_batch_bar = F,gene_text_cex = 2, cluster_text_cex =2)
  close_plot()
  
 
  open_plot(path=main_figures_path,fn = "figure_1f",plot_type = "pdf",width = 15,height = 10)
  plot_truth_heatmap(ds,cell_to_sample =ileum_ldm$dataset$cell_to_sample[colnames(ds)],
                     cell_to_cluster = ileum_ldm$dataset$cell_to_cluster[colnames(ds)],
                     insamples = c(inflamed_samples,uninflamed_samples),
                     ingenes = gene_list2_adj,
                     inclusts =clusters,
                     zlim=zlim,cols=colgrad_abs,plot_batch_bar = F,gene_text_cex = .5, cluster_text_cex = 1)
  close_plot()
  
  open_plot(path=main_figures_path,fn = "figure_1c",plot_type = "pdf",width = 5,height = 10)
  plot_truth_heatmap(ds,cell_to_sample =ileum_ldm$dataset$cell_to_sample[colnames(ds)],
                     cell_to_cluster = ileum_ldm$dataset$cell_to_cluster[colnames(ds)],
                     insamples = c(inflamed_samples,uninflamed_samples),
                     ingenes = gene_list1_adj,
                     inclusts =clusters,
                     zlim=zlim,cols=colgrad_abs,plot_batch_bar = F,gene_text_cex = .8, cluster_text_cex =1)
  close_plot()
  
  
  
}

get_pooled_freqs=function(samples,celltypes_to_include=c()){
  scrna=as.data.frame.matrix(table(ileum_ldm$dataset$cell_to_cluster,ileum_ldm$dataset$cell_to_sample))
  scrna_subset=scrna[,samples]

  dcs=setdiff(grep("DC",ileum_ldm$clustAnnots,val=T),"pDC")
  pool_by=cluster_to_cluster_set_with_pdc[rownames(scrna_subset)]
  pool_by[ileum_ldm$clustAnnots[names(pool_by)]%in%dcs]="cDC"
  ag=aggregate(scrna_subset,by=list(pool_by),sum)
  rownames(ag)=ag[,1]
  ag=ag[,-1]
  scrna_counts=ag[celltypes_to_include,]
  scrna_perc=t(100*t(scrna_counts)/colSums(scrna_counts))
  colnames(scrna_perc)=sample_to_patient[colnames(scrna_perc)]
  return(scrna_perc)
}


# Generating figure 1e and s2B
cytof_comparison=function(){
  cytof=read.csv(paste(pipeline_path,"input/tables/CyTOF_cell_perecentage.csv",sep=""))
  colnames(cytof)=gsub(pattern = "\\."," ",colnames(cytof))
  cytof_inf=cytof[cytof[,1]=="Inflamed",-1]
  cytof_uninf=cytof[cytof[,1]=="Uninflamed",-1]
  rownames(cytof_inf)=cytof_inf[,1]
  cytof_inf=cytof_inf[,-1]
  rownames(cytof_uninf)=cytof_uninf[,1]
  cytof_uninf=cytof_uninf[,-1]
  cytof_to_scRNA_annots=c("T","T","T","B","Plasma","ILC","ILC","pDC","cDC","MNP","Neutrophils","Eosinophils","Basophils","Mast")
  names(cytof_to_scRNA_annots)=c('CD4 T cells','CD8 T cells','DN T cells','B cells','Plasma cells','NK cells','ILCs','pDC','cDC','Macrophages','Neutrophils','Eosinophils','Basophils','Mast cells')
  
  ag=aggregate(cytof_inf,list(cytof_to_scRNA_annots[rownames(cytof_inf)]),FUN = sum);rownames(ag)=ag[,1];ag=ag[,-1]
  cytof_inf_pooled=ag
  ag=aggregate(cytof_uninf,list(cytof_to_scRNA_annots[rownames(cytof_uninf)]),FUN = sum);rownames(ag)=ag[,1];ag=ag[,-1]
  cytof_uninf_pooled=ag
  
  
  scrna_inf_perc=get_pooled_freqs(inflamed_samples,c('B','cDC','ILC','Mast','MNP','pDC','Plasma','T'))
  scrna_uninf_perc=get_pooled_freqs(uninflamed_samples,c('B','cDC','ILC','Mast','MNP','pDC','Plasma','T'))
  
  
  patients=intersect(colnames(scrna_inf_perc),colnames(cytof_inf_pooled))
  celltypes=unique(c(rownames(scrna_inf_perc),rownames(cytof_inf_pooled)))
  
  add_missing_celltypes=function(m,celltypes){
    m2=matrix(0,length(celltypes),ncol(m),dimnames = list(celltypes,colnames(m)))
    m2[rownames(m),]=as.matrix(m)
    return(m2)
  }
  scrna_inf_perc=add_missing_celltypes(scrna_inf_perc,celltypes)
  scrna_uninf_perc=add_missing_celltypes(scrna_uninf_perc,celltypes)
  cytof_inf_pooled=add_missing_celltypes(cytof_inf_pooled,celltypes)
  cytof_uninf_pooled=add_missing_celltypes(cytof_uninf_pooled,celltypes)
  

  options(scipen=5)
  col=c(celltypes_cols1[c("B","MNP","ILC","Mast","MNP","pDC","Plasma","T")],"gray","gray","gray")
  pch=                c(20,  3,20,20,20,20,20,20,4,3,18)
  legend_text=c('B cells','cDC','ILC','Mast cells','Macrophages','pDC','Plasma cells','T cells','Basophils','Eosinophils','Neutrophils')
  open_plot(path=main_figures_path,fn = "figure_1e",plot_type = "pdf",width = 4,height = 4)
  par(mar=c(4,4,1,1))
  plot(1e-2+scrna_inf_perc[celltypes,patients],1e-2+cytof_inf_pooled[celltypes,patients],log="xy",col=col,pch=pch,panel.first=grid(lty=1),xlab="scRNA",ylab="CyTOF",xlim=c(0.01,100),ylim=c(0.01,100),cex=1,cex.axis=1)
  legend("bottomright",legend = legend_text,col = col,pch=pch,cex=.7)
  close_plot()
  open_plot(path=supp_figures_path,fn = "figure_s2b",plot_type = "pdf",width = 4,height = 4)
  par(mar=c(4,4,1,1))
  plot(1e-2+scrna_uninf_perc[celltypes,patients],1e-2+cytof_uninf_pooled[celltypes,patients],log="xy",col=col,pch=pch,panel.first=grid(lty=1),xlab="scRNA",ylab="CyTOF",xlim=c(0.01,100),ylim=c(0.01,100),cex=1,cex.axis=1)
  legend("bottomright",legend = legend_text,col =col,pch=pch ,cex=.7)
  close_plot()
  cell_types_fitered=setdiff(celltypes,c("Basophils","Eosinophils","Neutrophils"))
  res_cor_inf=cor.test(log(1e-2+scrna_inf_perc[cell_types_fitered,patients]),log(1e-2+cytof_inf_pooled[cell_types_fitered,patients]))
  res_cor_uninf=cor.test(log(1e-2+scrna_uninf_perc[cell_types_fitered,patients]),log(1e-2+cytof_uninf_pooled[cell_types_fitered,patients]))
  stats_list[["figure_1e"]]<<-list(inflamed_pvalue=res_cor_inf$p.value,inflamed_r=as.vector(res_cor_inf$estimate),uninflamed_pvalue=res_cor_uninf$p.value,uninflamed_r=as.vector(res_cor_uninf$estimate))

}
# Cell-lineage distribution barplot
figure_1d=function(){
  scrna_inf_perc=get_pooled_freqs(inflamed_samples,c("T","ILC","Plasma","B","MNP","pDC","Mast","Stromal"))
  scrna_uninf_perc=get_pooled_freqs(uninflamed_samples,c("T","ILC","Plasma","B","MNP","pDC","Mast","Stromal"))
  scrna_uninf_perc=scrna_uninf_perc[,order(scrna_inf_perc[1,])]
  scrna_inf_perc=scrna_inf_perc[,order(scrna_inf_perc[1,])]
  m=matrix(0,nrow(scrna_inf_perc),2*ncol(scrna_inf_perc),dimnames = list(rownames(scrna_inf_perc),paste(rep(colnames(scrna_inf_perc),each=2),c("Inf","Uninf"))))
  m[,(1:ncol(scrna_inf_perc))*2-1]=scrna_inf_perc
  m[,(1:ncol(scrna_inf_perc))*2]=scrna_uninf_perc
  
  open_plot(path=main_figures_path,fn = "figure_1d",plot_type = "pdf",width = 4,height = 4)
  par(mar=c(6,3,1,1),lwd = .3)
  barplot(m,col=celltypes_cols1[rownames(m)],las=2,cex.axis=.8,cex.names=.6,space = c(.8,.2))
  grid(nx=NA,ny=NULL,lty=1)
  barplot(m,col=celltypes_cols1[rownames(m)],las=2,cex.axis=.8,cex.names=.6,space = c(.8,.2),add=T)
  close_plot()
}


get_cell_counts=function(lm,selected_samples){
  scRNA_tab=table(lm$dataset$cell_to_cluster,(lm$dataset$cell_to_sample))
  not_good_clusters=names(lm$clustAnnots[rownames(scRNA_tab)])[grep("Not good",lm$clustAnnots[rownames(scRNA_tab)])]
  scRNA_tab=scRNA_tab[setdiff(rownames(scRNA_tab),not_good_clusters),as.character(selected_samples)]
  
  scRNA_tab=t(scRNA_tab)
  rownames(scRNA_tab)=sample_to_patient[rownames(scRNA_tab)]
  return(scRNA_tab)
}

# Clusters_fractions_barplots_inflamded_vs_uninflamed
figure_1g=function(){
  
  
  get_freqs=function(lm,selected_samples){
    scRNA_tab=get_cell_counts(lm,selected_samples)
    freqs=(scRNA_tab)/rowSums(scRNA_tab)
  
    return(freqs)
  }
  freqs_inf=get_freqs(ileum_ldm,inflamed_samples_filtered[sample_to_patient[inflamed_samples_filtered]%in%c(pat1,pat2)])
  freqs_uninf=get_freqs(ileum_ldm,uninflamed_samples_filtered[sample_to_patient[uninflamed_samples_filtered]%in%c(pat1,pat2)])
  border_col=1

  y_inf=colMeans(freqs_inf)
  y_uninf=colMeans(freqs_uninf)
  reg=0.005
  ranks_l=sapply(split(data.frame(names(y_inf),(reg+y_uninf)/(reg+y_inf)),cluster_to_cluster_set_with_pdc[colnames(freqs_inf)]),function(x){rankv=rank(x[,2]);names(rankv)=x[,1];rankv})
  ranks_l=ranks_l[rev(c("T","ILC","Plasma","B","MNP","pDC","Mast","Stromal"))]
  v=rep((1:length(ranks_l))*1000,sapply(ranks_l,length))+unlist(ranks_l)
  names(v)=unlist(sapply(ranks_l,names))
  ord=names(sort(v))
  open_plot(path=main_figures_path,fn="figure_1g",plot_type = "pdf",width=4,height=4)
  layout(matrix(1:2,1,2))
  par(mar=c(5,5,2,0),lwd = .3 )
  infy=colMeans(freqs_inf)[ord]
  se=apply(freqs_inf,2,sd)[ord]/sqrt(nrow(freqs_inf))
  barpx<-barplot(infy,beside=T,las=2,col=celltypes_cols1[cluster_to_cluster_set_with_pdc[names(infy)]],xlim=c(.14,0),horiz=T,cex.names = .3,xlab="Inflamed",cex.axis = .5,border = border_col)
  grid(nx=NULL, ny=NA)
  barpx<-barplot(infy,beside=T,las=2,col=celltypes_cols1[cluster_to_cluster_set_with_pdc[names(infy)]],xlim=c(.14,0),horiz=T,cex.names = .3,xlab="Inflamed",cex.axis = .5,border = border_col,add=T)
  error.bar(barpx,infy,se,col=border_col,length = .01)
  par(mar=c(5,0,2,5))
  uninfy=colMeans(freqs_uninf)[ord]
  se=apply(freqs_uninf,2,sd)[ord]/sqrt(nrow(freqs_uninf))
  barpx<-barplot(uninfy,beside=T,las=2,col=celltypes_cols2[cluster_to_cluster_set_with_pdc[names(uninfy)]],xlim=c(0,.14),horiz=T,names.arg ="",xlab= "Uninflamed",cex.axis=.5,border = border_col)
  grid(nx=NULL, ny=NA)
  barpx<-barplot(uninfy,beside=T,las=2,col=celltypes_cols2[cluster_to_cluster_set_with_pdc[names(uninfy)]],xlim=c(0,.14),horiz=T,names.arg ="",xlab= "Uninflamed",cex.axis=.5,border = border_col,add=T)
  error.bar(barpx,uninfy,se,col=border_col,length = .01)

  close_plot()
}

# distance_between_inf_uninf
figure_s1o=function(){
  pool_subtypes_frequencies=function (lm, samples, cluster_sets, pool_subtype = T) {
    cluster_sets=cluster_sets[!names(cluster_sets)%in%"Not good"]
    get_freqs=function(lm,selected_samples){
      scRNA_tab=get_cell_counts(lm,selected_samples)
      freqs=(scRNA_tab)/rowSums(scRNA_tab)
      return(freqs)
    }
    freqs = get_freqs(lm, samples)
    pool_subtype_freqs = function(one_subtype) {
      return(rowSums(freqs[, unlist(one_subtype), drop = F]))
    }
    pool_one_clusterset = function(one_clusterset) {
      subtypes_freqs = sapply(one_clusterset,  pool_subtype_freqs)
      colnames(subtypes_freqs) = names(one_clusterset)
      return(subtypes_freqs)
    }
    return(sapply(cluster_sets, pool_one_clusterset, simplify = F))
  }
  freqs_inf= 1e-2+t(do.call(cbind,pool_subtypes_frequencies(ileum_ldm,inflamed_samples,cluster_sets = ileum_ldm$cluster_sets ,pool_subtype=T)))[,-2]
  freqs_uninf= 1e-2+t(do.call(cbind,pool_subtypes_frequencies(ileum_ldm,uninflamed_samples,cluster_sets = ileum_ldm$cluster_sets ,pool_subtype=T)))[,-2]
  open_plot(path =supp_figures_path,fn = "figure_s1o",plot_type = "pdf",6,6)
  par(mar=c(5,5,1,1))
  barplot(sqrt(colSums((freqs_inf-freqs_uninf)^2)),ylab="dissimilarity (inf vs uninf",border=F,cex.names = .7)
  close_plot()
}

# pca analysis

figure_1h_s2f_s2g=function(){
  
  reg=1e-3
  mat_inf_freq= t(do.call(cbind,normalize_by_clusterset_frequency(ileum_ldm$dataset$cell_to_cluster,ileum_ldm$dataset$cell_to_sample,inflamed_samples_filtered,cluster_sets = ileum_ldm$cluster_sets[setdiff(names(ileum_ldm$cluster_sets),"Mast")] ,pool_subtype=T,reg=reg)))
  colnames(mat_inf_freq)=paste(colnames(mat_inf_freq),"inf")
  mat_uninf_freq=t(do.call(cbind,normalize_by_clusterset_frequency(ileum_ldm$dataset$cell_to_cluster,ileum_ldm$dataset$cell_to_sample,uninflamed_samples_filtered,cluster_sets =  ileum_ldm$cluster_sets[setdiff(names(ileum_ldm$cluster_sets),"Mast")] ,pool_subtype=T,reg=reg)))
  colnames(mat_uninf_freq)=paste(colnames(mat_uninf_freq),"uninf")
  
  z=log2(cbind(mat_inf_freq,mat_uninf_freq));z2=apply(z,1,quantile,c(.15,.85))
  zmask=z2[2,]-z2[1,]>2
  pca_res=princomp(t(z[zmask,]))
  pca_scores=pca_res$scores
  pca_loadings=pca_res$loadings

  plot_one_pca=function(xi,yi,with_labels=T){
    col='grey'
    
    
    par(mar=c(5,5,3,3))
    ylim=c(min(unclass(pca_scores[,yi])),0.03+max(unclass(pca_scores[,yi])))
    xlim=range(unclass(pca_scores[,xi]))
    inf_mask=1:ncol(mat_inf_freq)
    plot(unclass(pca_scores[inf_mask,c(xi,yi)]),col=0,pch=16,ylim=ylim,xlim=xlim,panel.first=grid(lty=1))
     arrows(unclass(pca_scores[-1*inf_mask,xi]),unclass(pca_scores[-1*inf_mask,yi]),unclass(pca_scores[inf_mask,xi]),unclass(pca_scores[inf_mask,yi]),col =col,lty=1,lwd=0.5,length = 0)
    if (with_labels){
      text(unclass(pca_scores[inf_mask,xi]),unclass(pca_scores[inf_mask,yi])+.3,labels = gsub("rp ","",sample_to_patient[gsub(pattern = "rp | inf| uninf",replacement = "",rownames(pca_scores)[inf_mask])]),cex=.8)
    }
    points(unclass(pca_scores[-1*inf_mask,xi]),unclass(pca_scores[-1*inf_mask,yi]),pch=20,cex=1.5,col="cornflowerblue")
    points(unclass(pca_scores[1*inf_mask,xi]),unclass(pca_scores[1*inf_mask,yi]),pch=20,cex=1.5,col=2)
    legend("topleft",legend = c("Inflamed","Uninflamed"),col = c(2,"cornflowerblue"),pch=20,cex=1,pt.cex = 1.5)
  }
  
  plot_one_pca_inf_uninf=function(xi,with_labels=T){
    col='grey'
    inf_mask=1:ncol(mat_inf_freq)
    scores_inf=unclass(pca_scores[inf_mask,xi])
    scores_uninf=unclass(pca_scores[-inf_mask,xi])
    
    deltas=scores_inf-scores_uninf
    par(mar=c(5,5,3,3))

    ylim=c(0,.5+max(deltas))
    xlim=range(scores_inf)
    inf_mask=1:ncol(mat_inf_freq)
    plot(scores_inf,deltas,col=0,pch=16,ylim=ylim,xlim=xlim,panel.first=grid(lty=1))
    if (with_labels){
      text(scores_inf,deltas+.5,labels = gsub("rp ","",sample_to_patient[gsub(pattern = "rp | inf| uninf",replacement = "",rownames(pca_scores)[inf_mask])]),cex=.8)
    }
    points(scores_inf,deltas,pch=20,cex=1.5,col=1)
  }
  
  
  
  open_plot(path=main_figures_path,fn="figure_1h_with_labels",plot_type = "pdf",width = 5,height = 5)
  plot_one_pca_inf_uninf(1)
  title(paste(round(100*sum(pca_res$sdev[1]^2)/sum(pca_res$sdev^2),digits=1),"% of variance explained by PC1"))
  close_plot()
  
  open_plot(path=main_figures_path,fn="figure_1h",plot_type = "pdf",width = 5,height = 5)
  plot_one_pca_inf_uninf(1,with_labels = F)
  title(paste(round(100*sum(pca_res$sdev[1]^2)/sum(pca_res$sdev^2),digits=1),"% of variance explained by PC1"),cex=.7)
  close_plot()
  
  open_plot(path=supp_figures_path,fn="figure_s1p_with_labels",plot_type = "pdf",width = 5,height = 5)
  #  plot_one_pca(1,2)
  plot_one_pca(1,2)
  title(paste(round(100*sum(pca_res$sdev[1:2]^2)/sum(pca_res$sdev^2),digits=1),"% of variance explained by PC1+PC2"))
  close_plot()
  
  open_plot(path=supp_figures_path,fn="figure_s1p",plot_type = "pdf",width = 5,height = 5)
  plot_one_pca(1,2,with_labels = F)
  title(paste(round(100*sum(pca_res$sdev[1:2]^2)/sum(pca_res$sdev^2),digits=1),"% of variance explained by PC1+PC2"),cex=.7)
  close_plot()
  
  
  v=cluster_to_cluster_set[names(cluster_to_subtype1)]
  names(v)=cluster_to_subtype1
  subtype_to_clusterset=v[unique(names(v))]
  plot_comp=function(i,ylim=c(-.5,.5)){
    par(mar=c(7,3,3,1),lwd=.3)
    m=sort(pca_loadings[,i])
    barplot(m,col=celltypes_cols1[subtype_to_clusterset[names(m)]],names.arg = names(m),ylim=ylim,las=2,cex.names = .4)
    grid(nx=NA, ny=NULL,lwd=2)
    barplot(m,col=celltypes_cols1[subtype_to_clusterset[names(m)]],names.arg = names(m),ylim=ylim,add=T,las=2,cex.names = .4)
  }
  
  open_plot(path=supp_figures_path,fn="figure_s1q",plot_type = "pdf",width = 4,height = 4)
  plot_comp(1)
  title(main="loadings PC1")
  close_plot()
  

  stats$pca_t_test<<-t.test(pca_scores[grep(" inf",rownames(pca_scores),val=T),1]-pca_scores[grep(" uninf",rownames(pca_scores),val=T),1])


}






numis_threshold_analysis=function(){
  mask=ileum_ldm_400umis$dataset$cell_to_sample%in%inflamed_samples_v2[sample_to_patient[inflamed_samples_v2]%in%pat1][2];tab=apply(table(floor((1+colSums(ileum_ldm_400umis$dataset$umitab[,mask])/100)),cluster_to_subtype[ileum_ldm_400umis$dataset$cell_to_cluster[mask]]),2,function(x){rev(cumsum(rev(x)))});tab=round(100*tab/rowSums(tab),digits=1)
  mask=ileum_ldm_400umis$dataset$cell_to_sample%in%inflamed_samples_v2[sample_to_patient[inflamed_samples_v2]%in%pat1];tab=apply(table(floor((1+colSums(ileum_ldm_400umis$dataset$umitab[,mask])/100)),cluster_to_subtype[ileum_ldm_400umis$dataset$cell_to_cluster[mask]]),2,function(x){rev(cumsum(rev(x)))});tab=round(100*tab/rowSums(tab),digits=1)

  samps=c(inflamed_samples_v2[sample_to_patient[inflamed_samples_v2]%in%pat1][1:3])
  tabl=list()
  for (samp in samps){
    mask=ileum_ldm_400umis$dataset$cell_to_sample%in%samp
    tab=apply(table(floor((1+colSums(ileum_ldm_400umis$dataset$umitab[,mask])/100)),cluster_to_cluster_set[ileum_ldm_400umis$dataset$cell_to_cluster[mask]]),2,function(x){rev(cumsum(rev(x)))})
    tabl[[samp]]=t(round(100*tab/rowSums(tab),digits=1))[,c(2,7,12,17)]
  }
  open_plot(path="output/supp_figures/","celltype_piecharts",plot_type = "pdf",4,3)
  par(mar=c(1,1,1,1))
  layout(matrix(1:(length(samps)*ncol(tab1)),length(samps)))
    for (i in 1:ncol(tab1)){
      for (j in 1:length(samps)){
        pie(tabl[[j]][,i],col = celltypes_cols1[rownames(tabl[[j]])],labels = labels,border="gray")
      }
    }
  close_plot()
}



# Cumulative distributions of UMI counts/cells
figure_s1a=function(){
  open_plot(supp_figures_path,"figure_s1a",plot_type = "pdf",width = 5,height = 5)
  samp_fact=as.factor(sample_to_patient[names(ileum_ldm$dataset$numis_before_filtering)])
  matplot(sapply(ileum_ldm$dataset$numis_before_filtering,quantile,0:100/100),0:100/100,log="x",type='l',lty=ifelse(names(ileum_ldm$dataset$numis_before_filtering)%in%inflamed_samples,1,2),lwd=2,xlab="#UMIs",ylab="Cumulative Fraction",col=sample_cols[samp_fact])
  grid(lty=1)
  matplot(sapply(ileum_ldm$dataset$numis_before_filtering,quantile,0:100/100),0:100/100,log="x",type='l',lty=ifelse(names(ileum_ldm$dataset$numis_before_filtering)%in%inflamed_samples,1,2),lwd=2,xlab="#UMIs",ylab="Cumulative Fraction",col=sample_cols[samp_fact],add=T)
  legend("bottomright",legend=paste(samp_fact,c("I","U")),col=sample_cols[samp_fact],lty=1:2,ncol=2,lwd=2,cex=.7)
  close_plot()
}

ncells_per_celltype=function(){
  tab=table(cluster_to_subtype[ileum_ldm$dataset$cell_to_cluster])
  clusterset=sapply(strsplit(names(tab)," | "),function(x){tail(x,1)})
  names(tab)=sapply(strsplit(names(tab)," | "),function(x){head(x,1)})
  ord=order(clusterset)
  
  open_plot("output/supp_figures/","figure_sx_ncells_per_celltype",plot_type = "pdf",width = 12,height = 8)
  par(mar=c(10,5,4,1))
  barplot(tab[ord],log="y",las=2,ylim=c(100,10000),col=celltypes_cols1[clusterset[ord]])
  close_plot()
}

# Number of cells per cluster
figure_s1j=function(){
  tab=table(ileum_ldm$dataset$cell_to_cluster)[ileum_ldm$cluster_order]
  clusterset=sapply(strsplit(cluster_to_cluster_set[names(tab)]," | "),function(x){tail(x,1)})
  ord=order(clusterset)
  open_plot(supp_figures_path,"figure_s1j",plot_type = "pdf",width = 15,height = 8)
  par(mar=c(5,7,4,1))
  barplot(tab[ord],log="y",las=2,ylim=c(1,10000),col=celltypes_cols1[clusterset[ord]],ylab="", xlab="")
  mtext("Cell Counts",side = 2,line=4,cex=2)
  mtext("Cluster",side=1,line=4,cex=2)
  close_plot()
}

# Median UMI counts per cluster
figure_s1k=function(){
  numis=colSums(ileum_ldm$dataset$umitab)
  med=sapply(split(numis,ileum_ldm$dataset$cell_to_cluster),median)[ileum_ldm$cluster_order]
  clusterset=sapply(strsplit(cluster_to_cluster_set[names(med)]," | "),function(x){tail(x,1)})
  ord=order(clusterset)
  open_plot(supp_figures_path,"figure_s1k",plot_type = "pdf",width = 15,height = 8)
  par(mar=c(5,7,4,1))
  barplot(med[ord],log="y",las=2,col=celltypes_cols1[clusterset[ord]],ylab="", xlab="",ylim=c(100,10000))
  mtext("Median (#UMIs)",side = 2,line=4,cex=2)
  mtext("Cluster",side=1,line=4,cex=2)
  close_plot()
}

figure_s1b_c=function(){
  # in-silico_sorting by mitochondrial genes
  open_plot(supp_figures_path,"figure_s1b",plot_type = "pdf",width = 6,height = 6)
  samp=inflamed_samples_v2[2]
  prof="MITO"
  numis=ileum_ldm$dataset$numis_before_filtering[[samp]]
  plot(numis,ileum_ldm$dataset$insilico_gating_scores[[prof]][names(numis)],log="x",ylab=prof,xlab="#UMIs",cex=.8,panel.first = grid(lty=1),col=ifelse(numis<1000,"gray",1))
  abline(h=ileum_ldm$model$params$insilico_gating[[prof]]$interval[2],col=2,lty=2)
  close_plot()
  
  # in-silico_sorting by epithelial genes
  open_plot(supp_figures_path,"figure_s1c",plot_type = "pdf",width = 6,height = 6)
  samp=inflamed_samples_v2[2]
  prof="EP"
  numis=ileum_ldm$dataset$numis_before_filtering[[samp]]
  plot(numis,ileum_ldm$dataset$insilico_gating_scores[[prof]][names(numis)],log="x",ylab=prof,xlab="#UMIs",cex=.8,panel.first = grid(lty=1),col=ifelse(numis<1000,"gray",1))
  abline(h=ileum_ldm$model$params$insilico_gating[[prof]]$interval[2],col=2,lty=2)
  close_plot()
  
}

# Number of cells per sample 
figure_s1d=function(){
  tab_inf=table(sample_to_patient[ileum_ldm$dataset$cell_to_sample[ileum_ldm$dataset$cell_to_sample%in%inflamed_samples]])
  tab_uninf=table(sample_to_patient[ileum_ldm$dataset$cell_to_sample[ileum_ldm$dataset$cell_to_sample%in%uninflamed_samples]])
  
  open_plot(supp_figures_path,"figure_s1d",plot_type = "pdf",width = 6,height = 6)
  par(mar=c(5,5,1,1))
  tab=rbind(tab_inf,tab_uninf)[,order(as.numeric(sapply(strsplit(names(tab_inf)," "),tail,1)))]
  barplot(tab,beside=T,las=2,space = c(.5,2),col=cols_inf_uninf,ylim=c(0,10000))
  grid(nx=NA, ny=NULL,lty=1)
  barplot(tab,beside=T,las=2,space = c(.5,2),col=cols_inf_uninf,add=T)
  legend("topleft",c("Inflamed","Uninflamed"),col=cols_inf_uninf,pch=15,cex=1)
  close_plot()
}


# Estimated noise percentage barplots
figure_s1h=function(){
  open_plot(supp_figures_path,"figure_s1h",plot_type = "pdf",width = 6,height = 6)
  par(mar=c(5,5,1,1))
  v=100*ileum_ldm$dataset$alpha_noise
  v_inf=v[names(v)%in%inflamed_samples]
  v_uninf=v[names(v)%in%uninflamed_samples]
  names(v_inf)=sample_to_patient[names(v_inf)]
  names(v_uninf)=sample_to_patient[names(v_uninf)]
  tab=rbind(v_inf,v_uninf)[,order(as.numeric(sapply(strsplit(names(v_inf)," "),tail,1)))]
  barplot(tab,beside=T,ylim=c(0,10),las=2,space=c(0.5,2),ylab="% Noise",col=cols_inf_uninf)
  grid(nx=NA, ny=NULL,lty=1)
  barplot(tab,beside=T,ylim=c(0,10),las=2,space=c(0.5,2),add=T,ylab="% Noise",col=cols_inf_uninf)
  legend("topright",c("Inflamed","Uninflamed"),col=cols_inf_uninf,pch=15)
  close_plot()
}


doublets_analysis=function(){
  all_samples=ileum_ldm$dataset$samples
  a=sapply(split(as.data.frame(t(ileum_ldm$dataset$ll)),cluster_to_cluster_set_with_pdc[colnames(ileum_ldm$dataset$ll)],drop=F),apply,2,max)
  a=a[,setdiff(colnames(a),c("ILC","pDC"))]
  
  cluster_sets_models=update_models_debatched(umis=ileum_ldm$dataset$umitab,cell_to_cluster = cluster_to_cluster_set_with_pdc[ileum_ldm$dataset$cell_to_cluster],batch = ileum_ldm$dataset$cell_to_sample,noise_models = ileum_ldm$dataset$noise_models,alpha_noise = ileum_ldm$dataset$alpha_noise)
  numis=colSums(ileum_ldm$dataset$umitab)
  l=list()
  for (i in 1:ncol(a)){
    numis_i=numis[cluster_to_cluster_set_with_pdc[ileum_ldm$dataset$cell_to_cluster]==colnames(a)[i]]
    for (j in 1:ncol(a))
    {
      if (j>i){
        numis_j=numis[cluster_to_cluster_set_with_pdc[ileum_ldm$dataset$cell_to_cluster]==colnames(a)[j]]
        n=1000
        s=mean(numis_i)*cluster_sets_models[,colnames(a)[i]]+mean(numis_j)*cluster_sets_models[,colnames(a)[j]]
        l[[paste(colnames(a)[i],colnames(a)[j],sep=" X ")]]=s/sum(s)
      }
    }
  }
  doublet_models=do.call(cbind,l)
  dll=getBatchCorrectedLikelihood(ileum_ldm$dataset$umitab,doublet_models,noise_models = ileum_ldm$dataset$noise_models,ileum_ldm$dataset$cell_to_sample,alpha_noise = ileum_ldm$dataset$alpha_noise,reg=ileum_ldm$model$params$reg)[[1]]
  for (samp in all_samples){
    open_plot(path="output/supp_figures/",paste("doublet_analayis_",samp,sep=""),"pdf",10,10)
    layout(matrix(1:(ncol(a)^2),ncol(a),ncol(a),byrow=T))
    par(mar=c(4,4,0.1,0.1))
    for (i in 1:ncol(a)){
      for (j in 1:ncol(a))
      {
        if (i<=j){
          plot.new()
        }
        else{
          pair=c(i,j)
          mask=rownames(a)[ileum_ldm$dataset$cell_to_sample==samp&cluster_to_cluster_set[ileum_ldm$dataset$cell_to_cluster]%in%colnames(a)[pair]]
          ax=a[mask,pair]
          plot(dll[mask,paste(rev(colnames(ax)),collapse=" X ")],ax[,2]-ax[,1],pch=".",xlab="Doublet LL",ylab=paste(rev(colnames(ax)),collapse = " - "))
          abline(h=0,lty=1)
          abline(h=c(-.5,.5),lty=3)
          print(sum(abs(ax[,2]-ax[,1])<.5)/nrow(ax))
        }
      }
    }
    close_plot()
  }
}

modules_supps=function(){
  modsa=read.table("output/tables/ileum_all_0809_modules.txt",header=F,stringsAsFactors = F,row.names = 1)
  modsl=strsplit(modsa[,1],",")
  names(modsl)=rownames(modsa)
  s=t(sapply(modsl,function(x){colSums(ileum_ldm$dataset$ds[[match("2000",ileum_ldm$dataset$ds_numis)]][x,,drop=F])}))
  l=list()
  for (i in 1:length(modsl)){
    l[[i]]=split(t(log2(1+s[i,])),cluster_to_cluster_set_with_pdc[ileum_ldm$dataset$cell_to_cluster[colnames(s)]])
  }
  ord=order(sapply(l,function(x){median(x$Plasma)-median(unlist(x))}))
  open_plot(path="output/supp_figures/","module_per_clusterset","pdf",14,14)
  layout(matrix(1:100,10,10,byrow=T))
  par(mar=c(3,3,1,1))
  for (i in ord){
    boxplot(l[[i]],main=paste("module",i),lwd=.5,las=2,cex.axis=.7)
  }
  close_plot()
  
  
  
  doublet_plot=function(module_x="12",module_y="74",celltype_x="T",celltype_y="Plasma",b){
    maskb=intersect(colnames(s),names(ileum_ldm$dataset$cell_to_sample)[ileum_ldm$dataset$cell_to_sample==b])
    mask=maskb[cluster_to_cluster_set[ileum_ldm$dataset$cell_to_cluster[maskb]]%in%c(celltype_x,celltype_y)]
    plot(s[module_x,maskb],s[module_y,maskb],log="xy",pch=20,col="gray",main=paste(sample_to_patient[b],ifelse(b%in%inflamed_samples,"I","U")),xlab=paste("Module",module_x),ylab=paste("Module",module_y))
    points(s[module_x,mask],s[module_y,mask],col=celltypes_cols1[cluster_to_cluster_set[ileum_ldm$dataset$cell_to_cluster[mask]]],pch=20)
    
  }
  
  open_plot(path="output/supp_figures/","Plasma_T_doublets","pdf",10,10)
  for (b in c(inflamed_samples,uninflamed_samples)){
    doublet_plot(b=b)
  }
  close_plot()
}

# ncells per sample per cluster
table_s3=function(){
  write.csv(file=paste(pipeline_path,"output/tables/table_S3",sep="/"),table(ileum_ldm$dataset$cell_to_sample,ileum_ldm$dataset$cell_to_cluster),quote=F)
}

# correlation_between_subtypes_by_gene_expression
figure_1b=function(){
  
  gene_mask=apply(ileum_ldm$model$models,1,max)>5e-5
#  print(sum(gene_mask))
  m=log2(1e-5+ileum_ldm$model$models[gene_mask,unlist(ileum_ldm$cluster_sets[!names(ileum_ldm$cluster_sets)%in%"Not good"])])
  m=m-rowMeans(m)
  cormat=cor(m)
  d=as.dist(1-cormat)
  order=seriate(d,method = "GW")
  ord=get_order(order)
  
  
  large_margin=4
  small_margin=.5
  open_plot(path = main_figures_path,fn="figure_1b",plot_type = "pdf",width = 5,height = 5)
  celltype_ind=match(cluster_to_cluster_set_with_pdc[colnames(cormat)[ord]],names(celltypes_cols1))
  layout(matrix(c(1,3,4,2),2,2),heights=c(.7,20),widths = c(20,.7))
  par(mar=c(.2,large_margin,.2,small_margin))
  image(matrix(celltype_ind,,1),col=celltypes_cols1,axes=F,breaks=0.5+c(0:length(celltypes_cols1)))
  par(mar=c(large_margin,.2,small_margin,.2))
  image(matrix(celltype_ind,1,),col=celltypes_cols1,axes=F,breaks=0.5+c(0:length(celltypes_cols1)))
  
  par(mar=c(large_margin,large_margin,small_margin,small_margin))
  
  image(cormat[ord,ord],col=colorRampPalette(c("blue","white","red"))(100),axes=F,breaks=seq(-1,1,l=101))
  box(lwd=2)
  mtext(rownames(cormat)[ord],1,at = seq(0,1,l=ncol(cormat)),las=2,cex=.6,line =.5)
  mtext(rownames(cormat)[ord],2,at = seq(0,1,l=ncol(cormat)),las=2,cex=.6,line =.5)
  close_plot()
  
  open_plot(path = main_figures_path,fn="figure_1b_colorscale",plot_type = "pdf",width = 4,height = 1)
  par(mar=c(2,1,0,1))
  image(matrix(seq(-1,1,l=100),,1),col=colorRampPalette(c("blue","white","red"))(100),axes=F,breaks=seq(-1,1,l=101))
  axis(1)
  close_plot()
}




# stats per sample
table_s2=function(){
  annots=read.csv(paste(pipeline_path,"input/tables/Cohort_technical_summary.csv",sep="/"),stringsAsFactors = F)
  tenx=read.csv(paste(pipeline_path,"input/tables/tenx_stats_per_sample.csv",sep="/"),stringsAsFactors = F)
  all=merge(annots,tenx,by.y="sample_ID",by.x="index")

  selected_columns=c('index', 'Disease', 'Origin', 'tissue', 'status', 'Patient.ID', 'scRNAseq', 'X10.CHEMISTRY', 'CyTOF', 'MICSSS', 
                   'Number.of.Reads', 'Valid.Barcodes', 'Reads.Mapped.Confidently.to.Transcriptome', 'Reads.Mapped.Confidently.to.Exonic.Regions', 'Reads.Mapped.Confidently.to.Intronic.Regions', 'Reads.Mapped.Confidently.to.Intergenic.Regions', 'Q30.Bases.in.Barcode', 'Q30.Bases.in.RNA.Read', 'Q30.Bases.in.Sample.Index', 'Q30.Bases.in.UMI')

  tab=all[,selected_columns]

  ncells_per_sample=table(c(ileum_ldm$dataset$cell_to_sample,blood_ldm$dataset$cell_to_sample))
 
  nfiltered_cells_per_sample_ileum=cbind(sapply(ileum_ldm$dataset$gated_out_umitabs,function(l){sapply(l,ncol)}))
  nfiltered_cells_per_sample_blood=cbind(sapply(blood_ldm$dataset$gated_out_umitabs,function(l){sapply(l,ncol)}),EP=NA)
  colnames(nfiltered_cells_per_sample_blood)[colnames(nfiltered_cells_per_sample_blood)=="HB"]="RBC"
  nfiltered_cells_per_sample_blood=nfiltered_cells_per_sample_blood[,c("MITO","EP","RBC")]
  nfiltered_cells_per_sample=rbind(nfiltered_cells_per_sample_ileum,nfiltered_cells_per_sample_blood)
  colnames(nfiltered_cells_per_sample)=paste("N_cells_filtered_",colnames(nfiltered_cells_per_sample),sep="")
  numis_per_cell=round(t(sapply(split(c(colSums(ileum_ldm$dataset$umitab),colSums(blood_ldm$dataset$umitab)),c(ileum_ldm$dataset$cell_to_sample,blood_ldm$dataset$cell_to_sample)),quantile,c(.25,.5,.75))))

  colnames(numis_per_cell)=c("N_UMIs_percentile_25","N_UMIs_percentile_50","N_UMIs_percentile_75")

  tab=cbind(tab,nfiltered_cells_per_sample[as.character(tab$index),],N_lamina_propria_cells=as.numeric(ncells_per_sample[as.character(tab$index)]),numis_per_cell[as.character(tab$index),])
  write.csv(file=paste(pipeline_path,"output/tables/table_s2.csv",sep="/"),tab,row.names = F)
 
}

make_gene_reference_table=function(){
  l=list()
  gd=c()
  fig_ind=read.csv("input/Gene lists/Figures_Gene_lists.csv",stringsAsFactors = F)
  fns=paste("input/Gene lists/",fig_ind[,2],".txt",sep="")
  for (i in 1:nrow(fig_ind)){
    l[[fig_ind[i,1]]]=gsub(" ","",read.csv(fns[i],stringsAsFactors = F,header = F)[1,])
  #  write.xlsx(l[[fig_ind[i,1]]], file="output/tables/gene_index_table.xlsx",
  #             sheetName=fig_ind[i,1], append=i>1)
    
    #for (g in l[[fig_ind[i,1]]]){
    #  if (is.null(gd[g])){
    #    gd[g]=fig_ind[i,1]
    #  }
    #  else if (is.na(gd[g])){
    #    gd[g]=fig_ind[i,1]
    #  }
    #  else{
    #    gd[g]=paste(gd[g],fig_ind[i,1],sep=", ")
    #  }
    #}
  }
  maxl=max(sapply(l,length))
  tab=matrix("",maxl,length(l))
  colnames(tab)=names(l)
  for (i in 1:length(l)){
    tab[1:length(l[[i]]),i]=l[[i]]
  }
  write.csv(file = "output/tables/gene_index_table.csv",tab,row.names = F,quote=F)
  
}


numbers_figures1=function(){
  
  stats$total_number_of_lamina_propria_cells<<-sum(table(ileum_ldm$dataset$cell_to_sample)[c(inflamed_samples,uninflamed_samples)])
  stats$cell_per_cluster_range<<-range(table(ileum_ldm$dataset$cell_to_cluster[ileum_ldm$dataset$cell_to_sample%in%c(inflamed_samples,uninflamed_samples)]))
  stats$number_of_clusters<<-length(table(ileum_ldm$dataset$cell_to_cluster[ileum_ldm$dataset$cell_to_sample%in%c(inflamed_samples,uninflamed_samples)]))
  
}

make_figure1=function(){
  #also making supp figures 1-2 and supp tables 2-3
  message("Making Fig 1, Fig S1-2 and tables 2-3")
  
  figure_1b()
  # figure 1c, 1f:
  make_truth_plots()
  figure_1d()
  #figure 1e s1l:
  cytof_comparison()
  figure_1g()
  figure_1h_s2f_s2g()
  
  figure_s1a()
  figure_s1b_c()
  figure_s1d()
  # figure s1e-g - clustering
  figure_s1h() 
  # figure s1i   - clustering
  figure_s1j()
  figure_s1k()
  
  figure_s1o()
  table_s2()
  table_s3()
  
  numbers_figures1()
}
