
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
  gene_listp=readLines("input/Gene lists/Genes_list_for_PBMC_truth_plot v2.txt")
  gene_listp_adj=adjust_gene_names(gene_listp,row.names(ds))
  gene_listp_adj=gene_listp_adj[gene_listp_adj%in%rownames(ds)]
  # should be automatic!
  open_plot(path="output/main_figures/",fn = "truth_pbmc",plot_type = "png",width = 3000,height = 2000)
  plot_truth_heatmap(ds,cell_to_sample =blood_ldm$dataset$cell_to_sample[colnames(ds)],
                     cell_to_cluster = blood_ldm$dataset$cell_to_cluster[colnames(ds)],
                     insamples = blood_samples,
                     ingenes = gene_listp_adj,
                     inclusts =clusters,
                     zlim=zlim,cols=colgrad_abs,plot_batch_bar = F,gene_text_cex = 2, cluster_text_cex = 2,)
  close_plot()
}




plot_pbmc_MNP_scatter2=function(){
  
  inf_pat1=inflamed_samples[sample_to_patient[inflamed_samples]%in%pat1]
  inf_pat2=inflamed_samples[sample_to_patient[inflamed_samples]%in%pat2]
  
  tab=table(ileum_ldm$dataset$cell_to_sample,ileum_ldm$dataset$cell_to_cluster)[c(inf_pat1,inf_pat2),]
  rownames(tab)=sample_to_patient[rownames(tab)]
  tab=tab[c(pat1,pat2),]
  infMNP_gut_freq=tab[,unlist(ileum_ldm$cluster_sets$MNP$`Inf. Macrophages`)]/rowSums(tab[,unlist(ileum_ldm$cluster_sets$MNP)])
  
  tab2=table(blood_ldm$dataset$cell_to_sample,blood_ldm$dataset$cell_to_cluster)
  rownames(tab2)=sample_to_patient[rownames(tab2)]
  tab2=tab2[c(pat1,pat2),]
  
  open_plot(path="output/main_figures/",fn = "figure3_scatter_inf_MNP_classical_mono_pooled",plot_type = "pdf",width = 4,height = 4)
  par(mar=c(4,4,1,1))
  classical_mono_clusters=names(blood_ldm$clustAnnots)[blood_ldm$clustAnnots%in%c('Classical Mono-1','Classical Mono-2')]
  classical_mono_freq=rowSums(tab2[,classical_mono_clusters,drop=F])/rowSums(tab2)
  plot(infMNP_gut_freq,classical_mono_freq,col=ifelse(names(infMNP_gut_freq)%in%pat1,pat_cols[1],pat_cols[2]),panel.first = grid(lty=1),pch=20,cex=2)
  close_plot()
  print(cor.test(infMNP_gut_freq,classical_mono_freq,alternative = "less",method="spe"))
  
  open_plot(path="output/main_figures/",fn = "figure3_scatter_inf_MNP_classical_mono",plot_type = "pdf",width = 8,height = 4)
  par(mar=c(4,4,1,1))
  layout(matrix(1:2,1,2))
  clusts='Classical Mono-1'
  classical_mono_clusters=names(blood_ldm$clustAnnots)[blood_ldm$clustAnnots%in%c(clusts)]
  classical_mono_freq=rowSums(tab2[,classical_mono_clusters,drop=F])/rowSums(tab2)
  plot(infMNP_gut_freq,classical_mono_freq,col=ifelse(names(infMNP_gut_freq)%in%pat1,2,1),panel.first = grid(lty=1),pch=20,main=clusts)
  clusts='Classical Mono-2'
  classical_mono_clusters=names(blood_ldm$clustAnnots)[blood_ldm$clustAnnots%in%c(clusts)]
  classical_mono_freq=rowSums(tab2[,classical_mono_clusters,drop=F])/rowSums(tab2)
  plot(infMNP_gut_freq,classical_mono_freq,col=ifelse(names(infMNP_gut_freq)%in%pat1,2,1),panel.first = grid(lty=1),pch=20,main=clusts)
  close_plot()
  
  
  open_plot(path="output/main_figures/",fn = "figure3_scatter_inf_MNP_intermediate_mono",plot_type = "pdf",width = 4,height = 4)
  par(mar=c(4,4,1,1))
  c('Classical Mono-1','Classical Mono-2','Intermediate monocytes','Non classical monocytes')
  classical_mono_clusters=names(blood_ldm$clustAnnots)[blood_ldm$clustAnnots%in%c('Intermediate monocytes')]
  classical_mono_freq=rowSums(tab2[,classical_mono_clusters,drop=F])/rowSums(tab2)
  plot(infMNP_gut_freq,classical_mono_freq,col=ifelse(names(infMNP_gut_freq)%in%pat1,2,1),panel.first = grid(lty=1),pch=20)
  close_plot()
  
  open_plot(path="output/main_figures/",fn = "figure3_scatter_inf_MNP_non_classical_mono",plot_type = "pdf",width = 4,height = 4)
  par(mar=c(4,4,1,1))
  c('Classical Mono-1','Classical Mono-2','Intermediate monocytes','Non classical monocytes')
  classical_mono_clusters=names(blood_ldm$clustAnnots)[blood_ldm$clustAnnots%in%c('Non classical monocytes')]
  classical_mono_freq=rowSums(tab2[,classical_mono_clusters,drop=F])/rowSums(tab2)
  plot(infMNP_gut_freq,classical_mono_freq,col=ifelse(names(infMNP_gut_freq)%in%pat1,2,1),panel.first = grid(lty=1),pch=20)
  close_plot()
  
  
}


pbmc_cluster_freqs=function(){
  inf_pat1=inflamed_samples[sample_to_patient[inflamed_samples]%in%pat1]
  inf_pat2=inflamed_samples[sample_to_patient[inflamed_samples]%in%pat2]
  
  tab=table(ileum_ldm$dataset$cell_to_sample,ileum_ldm$dataset$cell_to_cluster)[c(inf_pat1,inf_pat2),]
  rownames(tab)=sample_to_patient[rownames(tab)]
  tab=tab[c(pat1,pat2),]
  infMNP_gut_freq=tab[,unlist(ileum_ldm$cluster_sets$MNP$`Inf. Macrophages`)]/rowSums(tab[,unlist(ileum_ldm$cluster_sets$MNP)])
  
  tab2=table(blood_ldm$dataset$cell_to_sample,blood_ldm$dataset$cell_to_cluster)
  rownames(tab2)=sample_to_patient[rownames(tab2)]
  tab2=tab2[c(pat1,pat2),]
  m2=tab2/rowSums(tab2)
  reg=1e-2
  m2=t(log2(t(reg+m2)/rowMeans(t(reg+m2))))
  ord=get_order(seriate(as.dist(1-cor(m2))))
  open_plot(path = "output/supp_figures/","blood_freq_heatmap",plot_type = "pdf",6,8)
  par(mar=c(7,12,1,1))
  image(m2[,ord],col=colorRampPalette(c("blue","white","red"))(100),axes=F,breaks=c(-100,seq(-2,2,l=99),100))
  mtext(text = rownames(m2),side = 1,at = seq(0,1,l=nrow(m2)),las=2,cex=1,col=1,line=.5)
  mtext(text = blood_ldm$clustAnnots[colnames(m2)[ord]],side = 2,at = seq(0,1,l=ncol(m2)),las=2,cex=1,col=1,line=.5)
  close_plot()
}



correlation_blood_gut=function(method="spe"){
  b=do.call(cbind,normalize_by_clusterset_frequency(blood_ldm,pool_subtype = T));b=b[intersect(c(pat1,pat2),rownames(b)),]
  g_inf=do.call(cbind,normalize_by_clusterset_frequency(ileum_ldm,pool_subtype = T,selected_samples = inflamed_samples));g_inf=g_inf[intersect(c(pat1,pat2),rownames(g_inf)),]
  g_uninf=do.call(cbind,normalize_by_clusterset_frequency(ileum_ldm,pool_subtype = T,selected_samples = uninflamed_samples));g_uninf=g_uninf[intersect(c(pat1,pat2),rownames(g_uninf)),]
  
  
  cormat_inf=cor(g_inf,b,method=method)
  cormat_uninf=cor(g_uninf,b,method=method)
  ord1=hclust(dist(cormat_inf))$order
  ord2=hclust(dist(t(cormat_inf)))$order
  
  
  reg=0.01
  b2=t(log2((reg+t(b))/(reg+apply(b,2,max))))
  cormat=cor(b2,method=method)
  ord=order(colMeans(b2[pat1,])/colMeans(b2[pat2,]))
  par(mar=c(5,10,2,2))
  image(b2[,ord],col=colorRampPalette(c("white","red"))(100),breaks=c(-100,seq(-2,0,l=99),100),axes=F);box()
  mtext(text = colnames(cormat)[ord],side = 2,at = seq(0,1,l=ncol(cormat)),las=2,cex=1,col=1,line=.5)
}


main_figure_4=function(){
  plot_pbmc_truth()
  plot_pbmc_MNP_scatter2()
}