
get_ligand_receptor_data=function(ldm,samps1,samps2,reg=1e-6,exprs_thresh=5e-6,ligand_receptor_table_fn){
    samps=c(samps1,samps2)
    a=read.csv(ligand_receptor_table_fn,stringsAsFactors = F)[,1:5]
    a=rbind(a,"-")
    a=apply(a,2,function(x){adjust_gene_names(paste(x,collapse = ","),rownames(ldm$dataset$umitab))})
    a=a[-nrow(a),]
    get_one=function(i){
      mask=a[,i]!=""
      return(cbind(ligand=a[mask,1],receptor=a[mask,i]))
    }
    validated_pairs<<-do.call(rbind,sapply(2:ncol(a),get_one))


  genes=rownames(ldm$model$models)
  validated_pairs<<-validated_pairs[validated_pairs[,1]%in%genes&validated_pairs[,2]%in%genes,]
  subtype_models=sapply(split(as.data.frame(t(ldm$model$models)),cluster_to_subtype1[colnames(ldm$model$models)],drop=F),colMeans)
  
  cc2<<-unique(validated_pairs[,1])
  corrected_counts<<-pmax(ldm$dataset$counts-ldm$dataset$noise_counts,0)
  m=apply(corrected_counts[samps,,],1:2,sum)
  m_bulk=t(m/rowSums(m))
  cc2<<-cc2[rowMaxs(subtype_models[cc2,])>exprs_thresh]
  rec<<-unique(as.vector(validated_pairs[match(cc2,validated_pairs[,1]),2]))
  validated_pairs<<-validated_pairs[validated_pairs[,1]%in%cc2&validated_pairs[,2]%in%rec,]
  
  counts_inf_pat1=apply(corrected_counts[samps1,,],2:3,sum)
  counts_inf_pat2=apply(corrected_counts[samps2,,],2:3,sum)

  a_rec=aperm(corrected_counts[samps,rec,],c(1,3,2))
  a_tot=100+array(apply(corrected_counts[samps,,],c(1,3),sum),dimnames = dimnames(a_rec),dim = dim(a_rec))
  a_rec2=(a_rec/a_tot)

 
  total_counts_per_patient<<-apply(corrected_counts,1,sum)
  total_counts_per_gene_per_patient<<-apply(corrected_counts,1:2,sum)
  total_counts_per_cluster_per_patient=apply(corrected_counts,c(1,3),sum)
  total_counts_per_subtype_per_patient<<-aggregate(t(total_counts_per_cluster_per_patient),list(cluster_to_subtype1[colnames(total_counts_per_cluster_per_patient)]),sum)
  rownames(total_counts_per_subtype_per_patient)<<-total_counts_per_subtype_per_patient[,1]


  total_counts_per_subtype_per_patient<<-total_counts_per_subtype_per_patient[,-1]
  l=apply(corrected_counts[,cc2,],1,function(x){x2=aggregate(t(x),list(cluster_to_subtype1[colnames(x)]),sum);rownames(x2)=x2[,1];x2=x2[,-1];return(x2)})
  counts_per_gene_per_subtype_per_patient<<-array(unlist(l),dim=c(dim(l[[1]]),length(l)),dimnames=list(rownames(l[[1]]),colnames(l[[2]]),names(l)))
  freqs=normalize_by_clusterset_frequency(cell_to_cluster = ldm$dataset$cell_to_cluster,cell_to_sample = ldm$dataset$cell_to_sample,samples = samps,cluster_sets=ldm$cluster_sets,pool_subtype = T,reg = 0)
  freqs<<-do.call(cbind,freqs)

  m_rec2=ldm$model$models[validated_pairs[,2],]
  m_rec2=m_rec2[match(unique(rownames(m_rec2)),rownames(m_rec2)),]
  m_rec2=m_rec2[,intersect(names(cluster_to_subtype1),colnames(m_rec2))]
  m_rec2=sapply(split(as.data.frame(t(m_rec2)),cluster_to_subtype1[colnames(m_rec2)],drop=F),colMeans)
  m_rec2=m_rec2[rowMaxs(m_rec2)>exprs_thresh,]

  m_rec2_normed<<-(log2((reg+(m_rec2))/(reg+rowMeans(m_rec2))))
  
  m_rec2_blood=blood_ldm$model$models[validated_pairs[,2],]
  m_rec2_blood=m_rec2_blood[match(unique(rownames(m_rec2_blood)),rownames(m_rec2_blood)),]
  m_rec2_blood=sapply(split(as.data.frame(t(m_rec2_blood)),blood_ldm$clustAnnots[colnames(m_rec2_blood)],drop=F),colMeans)
  m_rec2_blood=m_rec2_blood[rowMaxs(m_rec2_blood)>exprs_thresh,]
  m_rec2_blood_normed<<-(log2((reg+(m_rec2_blood))/(reg+rowMeans(m_rec2_blood))))
  
  
  rec_to_cc<<-sapply(split(validated_pairs[,1],validated_pairs[,2]),paste,collapse=",")
  cc_to_rec<<-split(validated_pairs[,2],validated_pairs[,1])
  
  
  m_cc2=subtype_models[cc2,]

  m_cc2_normed<<-(log2((reg+(m_cc2))/(reg+rowMeans(m_cc2))))
  
  cormatall=cor(m_cc2_normed)
  cormatall=cormatall[rowSums(!is.finite(cormatall))==0,]
  celltype_ord<<-colnames(cormatall)[get_order(seriate((as.dist(1-cor(cormatall)))))]
  
  m_rec2_normed<<-m_rec2_normed[intersect(unlist(cc_to_rec[rownames(m_cc2_normed)]),rownames(m_rec2_normed)),]
  
  
  ligand_ord=rev(rownames(m_cc2_normed)[order(apply(m_cc2_normed[,celltype_ord],1,which.max))])
  m_cc2_normed<<-m_cc2_normed[ligand_ord,]
  
  cormat_genes_rec=cor(t(m_rec2_normed))
  receptor_ord=get_order(seriate(as.dist(1-cormat_genes_rec),method = "GW"))
  m_rec2_normed=m_rec2_normed[receptor_ord,]
  m_rec2_normed=m_rec2_normed[intersect(rownames(m_rec2_normed),unlist(cc_to_rec[rownames(m_cc2_normed)])),]
  m_rec2_normed=m_rec2_normed[rowMaxs(abs(m_rec2_normed))>.25,]
  m_rec2_normed<<-m_rec2_normed
  
  
  
  corrected_tot_mean=corrected_counts/apply(corrected_counts,1,sum)
  n_expressed=(corrected_counts>0)[,c(validated_pairs[,1],validated_pairs[,2]),]
  
  
  total_umi_frac1=apply(corrected_tot_mean[samps1,,],2:3,mean)
  total_umi_frac2=apply(corrected_tot_mean[samps2,,],2:3,mean)
  
  
  pool_clusters=function(x){sapply(split(as.data.frame(t(x)),cluster_to_broader_subtype[colnames(x)],drop=F),colSums)}
  
  
  total_umi_frac_per_celltype1=pool_clusters(total_umi_frac1)
  total_umi_frac_per_celltype2=pool_clusters(total_umi_frac2)
  
  n_expressed_pooled=array(0,dim = c(dim(n_expressed)[1:2],dim(total_umi_frac_per_celltype1)[2]),dimnames = list(dimnames(n_expressed)[[1]],dimnames(n_expressed)[[2]],dimnames(total_umi_frac_per_celltype1)[[2]]))
  for (samp in dimnames(n_expressed)[[1]]){
    n_expressed_pooled[samp,,]=pool_clusters(n_expressed[samp,,])
  }
  
  n_patient_expressing1<<-colSums(n_expressed_pooled[samps1,,]>0)
  n_patient_expressing2<<-colSums(n_expressed_pooled[samps2,,]>0)
  
  
  
  m1=t(t(total_umi_frac_per_celltype1)/colSums(total_umi_frac_per_celltype1))
  m2=t(t(total_umi_frac_per_celltype2)/colSums(total_umi_frac_per_celltype2))
  
  total_umi_frac_per_celltype_ligands1<<-total_umi_frac_per_celltype1[cc2,]
  m1_rec<<-m1[rec,]
  total_umi_frac_per_celltype_ligands2<<-total_umi_frac_per_celltype2[cc2,]
  m2_rec<<-m2[rec,]
  
  
  
  
  validated_pairs_matrix=matrix(NA,length(rec),length(cc2))
  rownames(validated_pairs_matrix)=rec
  colnames(validated_pairs_matrix)=cc2
  for (i in 1:nrow(validated_pairs)){
    if (validated_pairs[i,2]%in%rec&(validated_pairs[i,1]%in%cc2)){
      validated_pairs_matrix[validated_pairs[i,2],validated_pairs[i,1]]=1
    }
  }
  validated_pairs_matrix<<-validated_pairs_matrix
  
}


figure_4b_c=function(){
  cormat_genes_rec=cor(t(m_rec2_normed))
  receptor_ord=get_order(seriate(as.dist(1-cormat_genes_rec),method = "GW"))
  m_rec2_normed<-m_rec2_normed[receptor_ord,]
  open_plot(main_figures_path,fn = "figure_4c",plot_type = "pdf",nrow(m_cc2_normed)/6,height =2+length(celltype_ord)/8)
  par(mar=c(1,1,4,15),bg=NA)
  image(m_rec2_normed[,celltype_ord],col=colgrad_rel,axes=F,breaks=c(-200,seq(-4,4,l=99),200))
  mtext(text = rownames(m_rec2_normed),side = 3,at = seq(0,1,l=nrow(m_rec2_normed)),las=2,cex=.7,line=.5)
  mtext(text = sapply(strsplit(celltype_ord," \\| "),head,1),side = 4,at = seq(0,1,l=ncol(m_rec2_normed)),las=2,cex=.7,line=.5)
  close_plot()


  open_plot(main_figures_path,fn = "figure_4b",plot_type = "pdf",nrow(m_cc2_normed)/6,height =2+length(celltype_ord)/8)
  par(mar=c(6,1,1,15),bg=NA)
  image(m_cc2_normed[rownames(m_cc2_normed),celltype_ord],col=colgrad_rel,axes=F,breaks=c(-200,seq(-4,4,l=99),200))
  mtext(text = rownames(m_cc2_normed),side = 1,at = seq(0,1,l=nrow(m_cc2_normed)),las=2,cex=.7,line=.5)
  mtext(text = sapply(strsplit(celltype_ord," \\| "),head,1),side = 4,at = seq(0,1,l=ncol(m_cc2_normed)),las=2,cex=.7,line=.5)
  close_plot()


  lrank=sapply(strsplit(rec_to_cc[rownames(m_rec2_normed)],","),function(x){y=match(x,rownames(m_cc2_normed));y=y[!is.na(y)]})

  a=cbind(rep(1:length(lrank),sapply(lrank,length))/length(lrank),unlist(lrank)/nrow(m_cc2_normed))

  l=list()

  ligands_to_highlight_tab=read.table(paste(pipeline_path,"/input/tables/figure_4bc_ligands_to_highlight.txt",sep=""),sep="\t",header=T,stringsAsFactors = F)
  ligands_to_highlight=strsplit(ligands_to_highlight_tab[,1],",")
  ligands_to_highlight_cols=ligands_to_highlight_tab[,2]
  open_plot(main_figures_path,fn = "figure_4bc_connectors",plot_type = "pdf",width =nrow(m_cc2_normed)/6,height =(2+length(celltype_ord)/8)/3)
  par(mar=c(0,1,0,15),bg=NA)
  par(xaxs="i")
  plot.new()
  cols=rep(colors()[350],length(unlist(lrank)))
  lwds=rep(1,length(unlist(lrank)))
  for (i in 1:length(ligands_to_highlight)){
    ranks=which(rownames(m_cc2_normed)%in%ligands_to_highlight[[i]])
    mask=unlist(lrank)%in%ranks
    cols[mask]=ligands_to_highlight_cols[i]
  }
  lwds[mask]=2
  arrows(a[,1],0,a[,2],1,code=0,col=cols,lwd=lwds)
  mask2=cols!="grey89"
  arrows(a[mask2,1],0,a[mask2,2],1,code=0,col=cols[mask2],lwd=lwds)

  close_plot()
}


DE_ligand_receptor_interactions=function(ldm,mask_pat2_ligands,mask_pat1_ligands,clusters_ligands,mask_pat2_clusters_receptors,mask_pat1_clusters_receptors,ligands,receptors,nmin_umi_thresh=0,nmin_cells_with_min_umi=20,reg=1e-14,nchunks=100,n_per_chunk=1000,noise_correction=F){
  numis_ligands=colSums(ldm$dataset$umitab[,c(mask_pat2_ligands,mask_pat1_ligands)])
  mask_pat2_clusters_ligands=mask_pat2_ligands[ldm$dataset$cell_to_cluster[mask_pat2_ligands]%in%clusters_ligands]
  mask_pat1_clusters_ligands=mask_pat1_ligands[ldm$dataset$cell_to_cluster[mask_pat1_ligands]%in%clusters_ligands]
  
  ntot_umis_ligands=colSums(ldm$dataset$umitab[,c(mask_pat2_ligands,mask_pat1_ligands)])
  u_ligands=ldm$dataset$umitab[,c(mask_pat2_ligands,mask_pat1_ligands)]
  u_ligands_clusters=ldm$dataset$umitab[,c(mask_pat2_clusters_ligands,mask_pat1_clusters_ligands)]
  n_cells_with_numis_above_thresh_ligands=rowSums(u_ligands_clusters>nmin_umi_thresh)
  names(n_cells_with_numis_above_thresh_ligands)=rownames(u_ligands_clusters)
  gene_mask_ligands=n_cells_with_numis_above_thresh_ligands>nmin_cells_with_min_umi&rownames(u_ligands_clusters)%in%c(ligands)
  if (sum(gene_mask_ligands)>0){
    u_ligands_clusters=u_ligands_clusters[gene_mask_ligands,,drop=F]
    u_ligands=u_ligands[gene_mask_ligands,,drop=F]
    if (noise_correction){
      u_ligands_clusters=as.matrix(u_ligands_clusters)
      samples=unique(ldm$dataset$cell_to_sample[colnames(u_ligands_clusters)])
      for (sampi in 1:length(samples)){
        mask_samp=colnames(u_ligands_clusters)[colnames(u_ligands_clusters)%in%names(which(ldm$dataset$cell_to_sample==samples[sampi]))]
        ncells=table(ldm$dataset$cell_to_cluster[mask_samp])
        mat_noise=t(t(ldm$dataset$noise_counts[sampi,,][gene_mask_ligands,ldm$dataset$cell_to_cluster[mask_samp],drop=F])/as.vector(ncells[ldm$dataset$cell_to_cluster[mask_samp]]))
        u_ligands_clusters[,mask_samp]=u_ligands_clusters[,mask_samp]-mat_noise
      }
    }
    
    
    u_receptors=ldm$dataset$umitab[,c(mask_pat2_clusters_receptors,mask_pat1_clusters_receptors)]
    n_cells_with_numis_above_thresh_receptors=rowSums(u_receptors>nmin_umi_thresh)
    names(n_cells_with_numis_above_thresh_receptors)=rownames(u_receptors)
    gene_mask_receptors=n_cells_with_numis_above_thresh_receptors>nmin_cells_with_min_umi&rownames(u_receptors)%in%c(receptors)
    
    if (sum(gene_mask_receptors>0)){
      u_receptors=u_receptors[gene_mask_receptors,c(mask_pat2_clusters_receptors,mask_pat1_clusters_receptors)]
      
      if (noise_correction){
        u_receptors=as.matrix(u_receptors)
        samples=unique(ldm$dataset$cell_to_sample[colnames(u_receptors)])
        for (sampi in 1:length(samples)){
          mask_samp=colnames(u_receptors)[colnames(u_receptors)%in%names(which(ldm$dataset$cell_to_sample==samples[sampi]))]
          ncells=table(ldm$dataset$cell_to_cluster[mask_samp])
          mat_noise=t(t(ldm$dataset$noise_counts[sampi,gene_mask_receptors,ldm$dataset$cell_to_cluster[mask_samp]])/as.vector(ncells[ldm$dataset$cell_to_cluster[mask_samp]]))
          u_receptors[,mask_samp]=u_receptors[,mask_samp]-mat_noise
        }
      }
      
      message("Testing ",sum(gene_mask_ligands)," ligands vs. ",sum(gene_mask_receptors)," receptors")
      
      obs_s_ligands=pmax(rowSums(u_ligands_clusters[,c(mask_pat2_clusters_ligands,mask_pat1_clusters_ligands)]),0)
      obs_s_pat2_ligands=pmax(rowSums(u_ligands_clusters[,mask_pat2_clusters_ligands,drop=F]),0)
      obs_s_pat1_ligands=pmax(obs_s_ligands-obs_s_pat2_ligands,0)
      
      obs_m_pat2_ligands=obs_s_pat2_ligands/sum(ntot_umis_ligands[mask_pat2_ligands])
      obs_m_pat1_ligands=obs_s_pat1_ligands/sum(ntot_umis_ligands[mask_pat1_ligands])
      
      obs_s_receptors=pmax(rowSums(u_receptors),0)
      obs_s_pat2_receptors=pmax(rowSums(u_receptors[,mask_pat2_clusters_receptors,drop=F]),0)
      obs_s_pat1_receptors=pmax(obs_s_receptors-obs_s_pat2_receptors,0)
      obs_m_pat2_receptors=obs_s_pat2_receptors/sum(obs_s_pat2_receptors)
      obs_m_pat1_receptors=obs_s_pat1_receptors/sum(obs_s_pat1_receptors)
      
      
      obs_interaction_intensity_pat1=matrix(obs_m_pat1_ligands,length(obs_m_pat1_ligands),length(obs_m_pat1_receptors))*matrix(obs_m_pat1_receptors,length(obs_m_pat1_ligands),length(obs_m_pat1_receptors),byrow = T,dimnames=list(names(obs_m_pat1_ligands),names(obs_m_pat1_receptors)))
      obs_interaction_intensity_pat2=matrix(obs_m_pat2_ligands,length(obs_m_pat2_ligands),length(obs_m_pat2_receptors))*matrix(obs_m_pat2_receptors,length(obs_m_pat2_ligands),length(obs_m_pat2_receptors),byrow = T,dimnames=list(names(obs_m_pat2_ligands),names(obs_m_pat2_receptors)))
      
      obs_log2_fc=log2((reg+obs_interaction_intensity_pat1)/(reg+obs_interaction_intensity_pat2))
      obs_log2_fc_arr=array(obs_log2_fc,dim=c(dim(obs_log2_fc)[1],dim(obs_log2_fc)[2],n_per_chunk))
      
      ncounts_bigger=obs_interaction_intensity_pat1*0
      n1_ligands=length(mask_pat2_ligands)
      n1_receptors=length(mask_pat2_clusters_receptors)
      mat_ligands=matrix(c(rep(T,n1_ligands),rep(F,length(ntot_umis_ligands)-n1_ligands)),n_per_chunk,length(ntot_umis_ligands),byrow =T)
      mat_receptors=matrix(c(rep(T,n1_receptors),rep(F,ncol(u_receptors)-n1_receptors)),n_per_chunk,ncol(u_receptors),byrow =T)
      
      ntot_receptors=ncol(u_receptors)
      ntot_ligands=ncol(u_ligands)
      
      pb=txtProgressBar(min = 1,max=nchunks)
      for (i in 1:nchunks){
        setTxtProgressBar(pb,value = i)
        
          mat_resampled_ligands=apply(mat_ligands,1,sample,ntot_ligands)
          
          s_bg_ligands=pmax(u_ligands%*%(mat_resampled_ligands*colnames(u_ligands)%in%c(mask_pat2_clusters_ligands,mask_pat1_clusters_ligands)),0)
          s_bg_ligands_ntot=ntot_umis_ligands%*%mat_resampled_ligands
          s_bg_receptors=pmax(u_receptors%*%apply(mat_receptors,1,sample,ntot_receptors),0)
        
        s_fg_ligands=pmax(obs_s_ligands-s_bg_ligands,0)
        s_fg_ligands_ntot=sum(ntot_umis_ligands)-s_bg_ligands_ntot
        m_bg_ligands=t(t(s_bg_ligands)/t(s_bg_ligands_ntot))
        m_fg_ligands=t(t(s_fg_ligands)/t(s_fg_ligands_ntot))
        
        s_fg_receptors=pmax(obs_s_receptors-s_bg_receptors,0)
        m_bg_receptors=t(t(s_bg_receptors)/colSums(s_bg_receptors))
        m_fg_receptors=t(t(s_fg_receptors)/colSums(s_fg_receptors))
        
        arr_fg_ligands=array(m_fg_ligands,dim=c(dim(m_fg_ligands)[1],n_per_chunk,dim(m_fg_receptors)[1]))
        arr_fg_receptors=aperm(array(m_fg_receptors,dim=c(dim(m_fg_receptors)[1],n_per_chunk,dim(m_fg_ligands)[1])),c(3,2,1))
        arr_bg_ligands=array(m_bg_ligands,dim=c(dim(m_bg_ligands)[1],n_per_chunk,dim(m_bg_receptors)[1]))
        arr_bg_receptors=aperm(array(m_bg_receptors,dim=c(dim(m_bg_receptors)[1],n_per_chunk,dim(m_bg_ligands)[1])),c(3,2,1))
        
        log2_fc=aperm(log2((reg+arr_fg_ligands*arr_fg_receptors)/(reg+arr_bg_ligands*arr_bg_receptors)),c(1,3,2))
        ncounts_bigger=ncounts_bigger+apply(abs(log2_fc)>abs(obs_log2_fc_arr),1:2,sum)
      }
      
      p.value=ncounts_bigger/(i*n_per_chunk)
      adj.p.value=matrix(p.adjust(p.value,method = "BH"),nrow(p.value),ncol(p.value),dimnames = dimnames(p.value))
      de_res=data.frame(ligand= rownames(p.value),receptor=rep(colnames(p.value),each=nrow(p.value)),counts_pat2_ligands=obs_s_pat2_ligands,counts_pat1_ligands=obs_s_pat1_ligands,freq_pat2_ligands=obs_m_pat2_ligands,freq_pat1_ligands=obs_m_pat1_ligands,n_cells_with_min_umis_above_thresh=n_cells_with_numis_above_thresh_ligands[gene_mask_ligands],log2_FC=as.vector(obs_log2_fc),p.value=as.vector(p.value),adj.p.value=as.vector(adj.p.value))
      return(de_res)
    }
    else{
      de_res=c()
    }
  }
  else{
    de_res=c()
  }
  
  
  
}





DE_ligand_receptor_interactions_patterns=function(clusters_ligand,clusters_receptors,samples,ligands,receptors,nmin_umi_thresh=0,nmin_cells_with_min_umi=20,reg=1e-14,nchunks=100,n_per_chunk=1000,ncells_per_sample=1000,noise_correction=T,ldm=ileum_ldm,min_n_cells=100){
  samp_by_sample=function(mask,cell_to_sample,ncells_per_sample){
    batches=unique(cell_to_sample[mask])
    mask2=c()
    for (b in batches){
      maskb=mask[cell_to_sample[mask]==b]
      if (length(maskb)>ncells_per_sample){
        maskb=sample(maskb,size = ncells_per_sample,replace = F)
      }
    #  message(b," ",length(maskb))
      mask2=c(mask2,maskb)
    }
    return(mask2)
  }
  
  pat1_samples=names(sample_to_patient)[sample_to_patient%in%pat1]
  pat2_samples=names(sample_to_patient)[sample_to_patient%in%pat2]
  pat1_samples_v2=intersect(samples,pat1_samples)
  pat2_samples_v2=intersect(samples,pat2_samples)
  
  # mask_pat1_clusters_ligands=colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample%in%pat1_samples_v2&ldm$dataset$cell_to_cluster%in%clusters_ligand]
  #  mask_pat2_clusters_ligands=colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample%in%pat2_samples_v2&ldm$dataset$cell_to_cluster%in%clusters_ligand]
  mask_pat1_ligands=colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample%in%pat1_samples_v2]
  mask_pat2_ligands=colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample%in%pat2_samples_v2]
  
  mask_pat1_clusters_receptors=colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample%in%pat1_samples_v2&ldm$dataset$cell_to_cluster%in%clusters_receptors]
  mask_pat2_clusters_receptors=colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample%in%pat2_samples_v2&ldm$dataset$cell_to_cluster%in%clusters_receptors]
  
  mask_pat1_ligands=samp_by_sample(mask_pat1_ligands,ldm$dataset$cell_to_sample,ncells_per_sample)
  mask_pat2_ligands=samp_by_sample(mask_pat2_ligands,ldm$dataset$cell_to_sample,ncells_per_sample)
  mask_pat1_clusters_receptors=samp_by_sample(mask_pat1_clusters_receptors,ldm$dataset$cell_to_sample,ncells_per_sample)
  mask_pat2_clusters_receptors=samp_by_sample(mask_pat2_clusters_receptors,ldm$dataset$cell_to_sample,ncells_per_sample)
  
  
  if (min(c(length(mask_pat1_ligands),length(mask_pat2_ligands)))>=min_n_cells&min(c(length(mask_pat1_clusters_receptors),length(mask_pat2_clusters_receptors)))>=min_n_cells){
    return(DE_ligand_receptor_interactions(ldm,mask_pat1_ligands,mask_pat2_ligands,clusters_ligand,mask_pat1_clusters_receptors,mask_pat2_clusters_receptors,ligands = ligands,receptors = receptors,nmin_umi_thresh=nmin_umi_thresh,nmin_cells_with_min_umi=nmin_cells_with_min_umi,reg=reg,nchunks=nchunks,n_per_chunk=n_per_chunk,noise_correction=noise_correction))
  }
}



get_two_subtypes_interactions=function(subtype1,subtype2,rec_thresh=1e-5,lig_thresh=1e-6,n_expressed_thresh=1,tables_output_path){
  reg=1e-14
  rec_crit1=m1_rec[,subtype2]>rec_thresh&(n_patient_expressing1[rownames(m1_rec),subtype2]>n_expressed_thresh)&rowSums(!is.na(validated_pairs_matrix))>0
  rec_crit2=m2_rec[,subtype2]>rec_thresh&(n_patient_expressing2[rownames(m2_rec),subtype2]>n_expressed_thresh)&rowSums(!is.na(validated_pairs_matrix))>0
 
  ligmat1=matrix(total_umi_frac_per_celltype_ligands1[,subtype1],nrow(m1_rec),nrow(total_umi_frac_per_celltype_ligands1),byrow = T)
  ligmat2=matrix(total_umi_frac_per_celltype_ligands2[,subtype1],nrow(m2_rec),nrow(total_umi_frac_per_celltype_ligands2),byrow = T)
  
  recmat1=matrix(m1_rec[,subtype2],nrow(m1_rec),nrow(total_umi_frac_per_celltype_ligands1))
  recmat2=matrix(m2_rec[,subtype2],nrow(m2_rec),nrow(total_umi_frac_per_celltype_ligands2))
  
  intensity_score1=log10(reg+recmat1*ligmat1*validated_pairs_matrix)
  intensity_score2=log10(reg+recmat2*ligmat2*validated_pairs_matrix)
  lig_crit1=(ligmat1>lig_thresh)&matrix(n_patient_expressing1[colnames(intensity_score1),subtype1]>n_expressed_thresh&colSums(!is.na(validated_pairs_matrix[,colnames(intensity_score1)]))>0,nrow(intensity_score1),ncol(intensity_score1),byrow=T)
  lig_crit2=(ligmat2>lig_thresh)&matrix(n_patient_expressing2[colnames(intensity_score2),subtype1]>n_expressed_thresh&colSums(!is.na(validated_pairs_matrix[,colnames(intensity_score2)]))>0,nrow(intensity_score2),ncol(intensity_score2),byrow=T)
  
  
  pair_mask_mat=(lig_crit1|lig_crit2)&(rec_crit1|rec_crit2)&ifelse(is.na(validated_pairs_matrix),F,T)
  
  z1=intensity_score1*ifelse(pair_mask_mat,pair_mask_mat,NA)
  z2=intensity_score2*ifelse(pair_mask_mat,pair_mask_mat,NA)
  row_mask=pmax(rowSums(pair_mask_mat,na.rm=T),0)>0
  column_mask=pmax(colSums(pair_mask_mat,na.rm=T),0)>0
  if (sum(row_mask)==0|sum(column_mask)==0){
    return(c())
  }
  z1=z1[which(row_mask),which(column_mask),drop=F]
  z2=z2[which(row_mask),which(column_mask),drop=F]
  

  if (!dir.exists(tables_output_path)){
    dir.create(tables_output_path)
  }
  
  ligands=colnames(z1)
  receptors=rownames(z1)
  message("Running differential intensity analysis ", subtype1, " vs ", subtype2)
  de_res=DE_ligand_receptor_interactions_patterns(names(cluster_to_broader_subtype)[cluster_to_broader_subtype==subtype1],names(cluster_to_broader_subtype)[cluster_to_broader_subtype==subtype2],inflamed_samples_v2_filtered,ligands = ligands,receptors = receptors,nchunks = 100,nmin_cells_with_min_umi=25,ncells_per_sample=2000,nmin_umi_thresh = 0,n_per_chunk = 1e3)
  rownames(de_res)=paste(de_res$ligand,de_res$receptor,sep="_")
  if (!is.null(de_res)){
    write.csv(de_res[order(de_res$log2_FC),],file=paste(tables_output_path,"/DE_inf_pat1_vs_pat2_",subtype1,"_",subtype2,".csv",sep=""))
  }
  
  
  
  mat=z1-z2
  ord1=order(rowSums(mat,na.rm=T),decreasing = T)
  ord2=order(colSums(mat,na.rm=T))
  
  
  mat_conv=as.data.frame(as.table(mat))
  mat_conv=mat_conv[,c(2:1,3)]
  colnames(mat_conv)=c("Ligand","Receptor","Log10_ratio_pat1_pat2")
  stats_mat=cbind(subtype1,subtype2,mat_conv[,1:2],log2_ratio_score1_score2=log2(10)*mat_conv[,3],log10_exprs_ligand_pat1=log10(reg+as.vector(total_umi_frac_per_celltype_ligands1[mat_conv[,1],subtype1])),log10_exprs_ligand_pat2=log10(reg+as.vector(total_umi_frac_per_celltype_ligands2[mat_conv[,1],subtype1])),n_expressed_ligand_pat1=n_patient_expressing1[as.character(mat_conv[,1]),subtype1],n_expressed_ligand_pat1=n_patient_expressing2[as.character(mat_conv[,1]),subtype1],log10_exprs_receptor_pat1=log10(reg+m1_rec[as.character(mat_conv[,2]),subtype2]),log10_exprs_receptor_pat2=log10(reg+m2_rec[as.character(mat_conv[,2]),subtype2]),n_expressed_receptor_pat1=n_patient_expressing1[as.character(mat_conv[,2]),subtype2],n_expressed_receptor_pat2=n_patient_expressing2[as.character(mat_conv[,2]),subtype2],log10_score1=as.vector(z1),log10_score2=as.vector(z2))
  stats_mat=stats_mat[!is.na(stats_mat[,5]),]
  stats_mat$p.value=de_res[paste(stats_mat$Ligand,stats_mat$Receptor,sep="_"),"p.value"]
  stats_mat$adj.p.value=p.adjust(stats_mat$p.value)
  rownames(stats_mat)=paste(stats_mat[,1],stats_mat[,2],stats_mat[,3],stats_mat[,4],sep="_")
  return(list(stats_mat=stats_mat,intensity_score1=intensity_score1,intensity_score2=intensity_score2,pair_mask=pair_mask_mat,plot_list=list(mat=mat,ord1=ord1,ord2=ord2)))
}


bin_plot2=function(subtype1,subtype2,mat,ord1,ord2,stats,fdr_thresh=1e-2,figure_path){
  xlab=paste(subtype2,"Receptors")
  ylab=paste(subtype1,"Ligands")
  open_plot(figure_path,paste("bin_lr",subtype1,subtype2,sep="_"),"pdf",.9+.2*nrow(mat),.9+.2*ncol(mat),2)
  par(mar=c(6,6,1,1))
  image(mat[ord1,ord2,drop=F],axes=F,breaks=c(-100,seq(-2,2,l=99),100),col=colorRampPalette(c("black","gray","red"))(100))
  
  
  stats2=stats[stats$subtype1==subtype1&stats$subtype2==subtype2,]
  pmask=stats2$adj.p.value<fdr_thresh
  t1=cbind((match(stats2[pmask,]$Receptor,rownames(mat)[ord1])-1)/(length(ord1)-1),(match(stats2[pmask,]$Ligand,colnames(mat)[ord2])-1)/(length(ord2)-1))
  text(t1,labels = rep("*",nrow(t1)),cex=1.5)
  mtext(1,text = xlab,line= 5)
  mtext(2,text = ylab ,line =5)
  box()
  mtext(text = rownames(mat)[ord1],side = 1,at = seq(0,1,l=nrow(mat)),las=2,cex=.8,line=.5)
  mtext(text = colnames(mat)[ord2],side = 2,at = seq(0,1,l=ncol(mat)),las=2,cex=.8,line=.5)
  close_plot()
}


pairwise_intensity_maps=function(st,tables_output_path,load_stats=F){
  intensity_score1=c()
  intensity_score2=c()
  pair_mask=c()
  interaction_stats=c()
  if (load_stats){
    load(file=paste(pipeline_path,"/intermediates/ligand_receptor.rd",sep=""))
  }
  else{
    res_l=list()
    for (i in 1:nrow(st)){
      subtype1=as.character(st[i,1])
      subtype2=as.character(st[i,2])
      key_s=paste(subtype1,subtype2,sep="_")
      
      res_l[[key_s]]=get_two_subtypes_interactions(subtype1,subtype2,tables_output_path = tables_output_path)
      interaction_stats=rbind(interaction_stats,res_l[[key_s]][["stats_mat"]])
      intensity_score1=c(intensity_score1,res_l[[key_s]][["intensity_score1"]])
      intensity_score2=c(intensity_score2,res_l[[key_s]][["intensity_score2"]])
      pair_mask=c(pair_mask,res_l[[key_s]][["pair_mask"]])
    }
    interaction_stats$adj.p.value=p.adjust(interaction_stats$p.value)
    save(list = c("res_l","interaction_stats","intensity_score1","intensity_score2","pair_mask"),file=paste(pipeline_path,"/intermediates/ligand_receptor.rd",sep=""))
  }
  for (i in 1:nrow(st)){
    subtype1=as.character(st[i,1])
    subtype2=as.character(st[i,2])
      key_s=paste(subtype1,subtype2,sep="_")
      print(key_s)
      bin_plot2(subtype1,subtype2,res_l[[key_s]]$plot_list[["mat"]],res_l[[key_s]]$plot_list[["ord1"]],res_l[[key_s]]$plot_list[["ord2"]],interaction_stats,figure_path = main_figures_path)
  }
  
  
  breaks=seq(-14,-6,l=17)-.25
  tab=rbind(log10(1+hist(pmax(intensity_score1,intensity_score2)[!pair_mask],plot=F,breaks=breaks)$counts),
            log10(1+hist(pmax(intensity_score1,intensity_score2)[pair_mask],plot=F,breaks=breaks)$counts))
  return(tab)
  
}


main_ligand_receptor=function(){
  samps_pat1=inflamed_samples_v2[sample_to_patient[inflamed_samples_v2]%in%pat1]
  samps_pat2=inflamed_samples_v2[sample_to_patient[inflamed_samples_v2]%in%pat2]
  
  get_ligand_receptor_data(ileum_ldm,samps_pat1,samps_pat2,ligand_receptor_table_fn = paste(pipeline_path,"/input/tables/Ligand_receptor_pairs.csv",sep=""))
  figure_4b_c()
  
  st=expand.grid(c("Macrophages","Fibroblasts","DC","T","Endothelial"),c("Macrophages","Fibroblasts","DC","T","Endothelial"));st=st[st[,1]!=st[,2],];st=st[st[,1]!="Endothelial",]
  score_distrib=pairwise_intensity_maps(st=st,tables_output_path = paste(pipeline_path,"/output/tables/",sep=""),load=T)
  
  open_plot(supp_figures_path,"figure_s6a",plot_type = "pdf",width = 5,height = 4)
  par(mar=c(6,5,1,1))
  barplot(score_distrib,beside=T,xlab="",ylab="log10(#pairs)",names.arg=seq(-14,-6.5,.5),ylim=c(0,3),las=2,col=c("gray",2),border=F,cex.axis = 1.2,cex.names = 1.2,space = c(0,.5))
  #legend("topright",col=c("gray",2),c("Lowly expressed Receptor/Ligand","Ligand & Receptor expression detected"),pch=15,cex=.5)
  close_plot()

}

