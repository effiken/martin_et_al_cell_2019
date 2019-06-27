
DE_between_two_sets=function(ldm,mask_bg,mask_fg,nmin_umi_thresh=0,nmin_cells_with_min_umi=20,reg=1e-5,nchunks=100,n_per_chunk=1000,noise_correction=F){
  
  u=ldm$dataset$umitab[,c(mask_bg,mask_fg)]
  n_cells_with_numis_above_thresh=rowSums(u>nmin_umi_thresh)
  names(n_cells_with_numis_above_thresh)=rownames(u)
  gene_mask=n_cells_with_numis_above_thresh>nmin_cells_with_min_umi
  u=u[gene_mask,c(mask_bg,mask_fg)]
  message("Testing ",sum(gene_mask)," genes")
  if (noise_correction){
    u=as.matrix(u)
    samples=unique(ldm$dataset$cell_to_sample[colnames(u)])
    for (sampi in 1:length(samples)){
      mask_samp=colnames(u)[colnames(u)%in%names(which(ldm$dataset$cell_to_sample==samples[sampi]))]
      ncells=table(ldm$dataset$cell_to_cluster[mask_samp])
      mat_noise=t(t(ldm$dataset$noise_counts[sampi,gene_mask,ldm$dataset$cell_to_cluster[mask_samp]])/as.vector(ncells[ldm$dataset$cell_to_cluster[mask_samp]]))
#    noise=matrix(ldm$dataset$beta_noise[ldm$dataset$cell_to_sample[colnames(u)]],nrow(u),ncol(u),byrow=T)*ldm$dataset$noise_models[rownames(u),ldm$dataset$cell_to_sample[colnames(u)]]
    #u[,mask_samp]=pmax(u[,mask_samp]-mat_noise,0)
      u[,mask_samp]=u[,mask_samp]-mat_noise
    }
  }
  obs_s=pmax(rowSums(u),0)
  obs_s_bg=pmax(rowSums(u[,mask_bg,drop=F]),0)
  obs_s_fg=pmax(obs_s-obs_s_bg,0)
  obs_m_bg=obs_s_bg/sum(obs_s_bg)
  obs_m_fg=obs_s_fg/sum(obs_s_fg)
  obs_log2_fc=log2((reg+obs_m_fg)/(reg+obs_m_bg))
  ncounts_bigger=rep(0,nrow(u))
  n1=length(mask_bg)
  mat=matrix(c(rep(T,n1),rep(F,ncol(u)-n1)),n_per_chunk,ncol(u),byrow =T)
 
  ntot=ncol(u)
  for (i in 1:nchunks){
    message(i)

    print(system.time({
    s_bg=pmax(u%*%apply(mat,1,sample,ntot),0)
    }))
    s_fg=pmax(obs_s-s_bg,0)
    m_bg=t(t(s_bg)/colSums(s_bg))
    m_fg=t(t(s_fg)/colSums(s_fg))
    log2_fc=log2((reg+m_fg)/(reg+m_bg))
    ncounts_bigger=ncounts_bigger+rowSums(abs(log2_fc)>=abs(obs_log2_fc))

#    save(file="de_current_iter.rd",i)
  }
  p.value=ncounts_bigger/(i*n_per_chunk)
  adj.p.value=p.adjust(p.value,method = "BH")
  de_res=data.frame(counts_bg=obs_s_bg,counts_fg=obs_s_fg,freq_bg=obs_m_bg,freq_fg=obs_m_fg,n_cells_with_min_umis_above_thresh=n_cells_with_numis_above_thresh[gene_mask],log2_FC=obs_log2_fc,p.value=p.value,adj.p.value=adj.p.value)

  return(de_res)
}


plot_de_volcano=function(de,fn,xsize=8,ysize=8,xlim=c(-1,1),main=""){
  pdf(fn,xsize,ysize)
  plot(de$log2_FC,-log10(de$adj.p.value),xlim=c(-1,1),panel.first = grid(lty=1),main=main,xlab="Log2(TET/WT)",ylab="-Log10(adj.pvalue)")
  dev.off()
  
}


#DE_pat1_vs_pat2=function(){
#  inf_pat1=inflamed_samples[sample_to_patient[inflamed_samples]%in%pat1]
#  inf_pat2=inflamed_samples[sample_to_patient[inflamed_samples]%in%pat2]
#  mask1=colnames(ileum_ldm$dataset$umitab)[ileum_ldm$dataset$cell_to_sample%in%inf_pat1]
#  mask2=colnames(ileum_ldm$dataset$umitab)[ileum_ldm$dataset$cell_to_sample%in%inf_pat2]
#  de_res=DE_between_two_sets(ileum_ldm,mask1,mask2,nchunks=1000,n_per_chunk=100,reg=1e-7)
#}

DE_tregs_pat1_vs_pat2=function(){
  inf_pat1=inflamed_samples[sample_to_patient[inflamed_samples]%in%pat1]
  inf_pat2=inflamed_samples[sample_to_patient[inflamed_samples]%in%pat2]
  clusts=cluster_sets$`T cells`$Tregs
  mask1=colnames(ileum_ldm$dataset$umitab)[ileum_ldm$dataset$cell_to_sample%in%inf_pat1&ileum_ldm$dataset$cell_to_cluster%in%clusts]
  mask2=colnames(ileum_ldm$dataset$umitab)[ileum_ldm$dataset$cell_to_sample%in%inf_pat2&ileum_ldm$dataset$cell_to_cluster%in%clusts]
  de_res=DE_between_two_sets(ileum_ldm,mask1,mask2,nchunks=1000,n_per_chunk=100,reg=1e-7)
  write.csv(de_res[order(de_res$log2_FC),],file="tables/DE_pat1_vs_pat2_tregs.csv")
  }
