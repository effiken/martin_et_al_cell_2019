
       
cluster_sets_for_heatmaps=c("Inf. Macrophages","mature cDCs","Resident Macrophages","Plasma cells_IgG","Plasma cells_IgA","Endothelial cells","Fibroblasts","Resident Memory")
# screen for marker genes diff. expressed between pat1 and pat2.
prediction_thresh<<-0.25


thresh_fc1=1
thresh_fc2=1
freq_thresh=1e-5
thresh=1



malat=c("MALAT1")
xist="XIST"
jchain="JCHAIN"
hla=strsplit("HLA-F,HLA-G,HLA-A,HLA-E,HLA-C,HLA-B,HLA-DRA,HLA-DRB5,HLA-DRB1,HLA-DQA1,HLA-DQB1,HLA-DQB1-AS1,HLA-DQA2,HLA-DQB2,HLA-DOB,HLA-DMB,HLA-DMA,HLA-DOA,HLA-DPA1,HLA-DPB1",",")[[1]]
mts=c('MT-ND1','MT-ND2','MT-CO1','MT-CO2','MT-ATP8','MT-ATP6','MT-CO3','MT-ND3','MT-ND4L','MT-ND4','MT-ND5','MT-ND6','MT-CYB','MTRNR2L11','MTRNR2L12','MTRNR2L13','MTRNR2L6','MTRNR2L10','MTRNR2L8','MTRNR2L7','MTRNR2L5','MTRNR2L4','MTRNR2L1','MTRNR2L3')
igv=grep("IGLV|IGKV|IGHV",rownames(ileum_ldm$dataset$umitab),val=T)
metallothionein=grep("^MT1",rownames(ileum_ldm$dataset$umitab),val=T)
excluded_genes=c(malat,xist,jchain,hla,mts,igv,metallothionein)




sample_set_list=list(
  list(datasets="CERTIFI",diagnosises="CD",tissues=c("Ileum")),
#  list(datasets="UNITI-1",diagnosises="CD",tissues=c("Ileum")),
  list(datasets="UNITI-2",diagnosises="CD",tissues=c("Ileum")),
  list(datasets="RISK",diagnosises="CD",tissues=c("Ileum"))
)



read_RISK=function(genes){
  path=paste(pipeline_path,"input/external_cohorts_data/",sep="")
  fn=paste(path,"/risk.rd",sep="")
  env=new.env()
  load(file=fn,envir=env)
  m=env$m
  sample_tab=env$sample_tab
  clin2=env$clin2
  pcdai=env$pcdai
  mask=env$mask
  
  return(list(raw=m[,mask],design=sample_tab[mask,],clin=clin2,pcdai=pcdai))
}


read_one_microarray_dataset=function(path,fn,genes,reg=10,remove_first_column=F,log_transform=T){
  if (tail(strsplit(fn,"\\.")[[1]],1)=="csv"){
    exprs_tab=read.csv(paste(path,fn,sep=""),row.names = 1)
  }
  else{
    exprs_tab=read.table(paste(path,fn,sep=""))
  }
  if (remove_first_column){
    exprs_tab=exprs_tab[,-1]
  }
  plot(reg+rowMeans(exprs_tab[,sample(colnames(exprs_tab),10)]),reg+rowMeans(exprs_tab[,sample(colnames(exprs_tab),10)]),log="xy",main="Before log transform")
  if (log_transform){
    exprs_tab=log10(reg+exprs_tab)
  }
  
  colnames(exprs_tab)=gsub(".CEL","",colnames(exprs_tab))
  return(exprs_tab[genes,])
}






read_certifi=function(genes){

  path=paste(pipeline_path,"input/external_cohorts_data/",sep="")
  fn=paste(path,"/certifi.rd",sep="")
  
  env=new.env()
  load(file=fn,envir=env)
  exprs_tab=env$exprs_tab
  annots=env$annots
  
  return(list(raw=exprs_tab,design=annots))
}



read_uniti1=function(load_text=F,genes){
  path=paste(pipeline_path,"input/external_cohorts_data/",sep="")
  fn=paste(path,"/uniti-1.rd",sep="")
  if (load_text){
    
    exprs_tab1=read_one_microarray_dataset(path="input/microarray/output/",fn="expr_UNITI-1.csv",genes,remove_first_column = F)
    
    annots1=read.csv(paste(path,"UNITI-1 Biopsy/UNITI1 biopsy_wk0wk8_347s_CDAI_CRP_FCAL_SESCD.csv",sep=""),stringsAsFactors = F,row.names = 1)
    rownames(annots1)=gsub("\\.|-","_",rownames(annots1))
    colnames(exprs_tab1)=gsub("\\.|-","_",colnames(exprs_tab1))
    exprs_tab1=exprs_tab1[,rownames(annots1)]
    
    save(file=fn,exprs_tab1,annots1)
  }
  else {
    env=new.env()
    load(file=fn,envir=env)
    exprs_tab1=env$exprs_tab1
    annots1=env$annots1
  }
  
  return(list(raw=exprs_tab1,design=annots1))
}


read_uniti2=function(load_text=F,genes){
  path=paste(pipeline_path,"input/external_cohorts_data/",sep="")
  fn=paste(path,"/uniti-2.rd",sep="")
  
  if (load_text){
    exprs_tab2=read_one_microarray_dataset(path="input/microarray/output/",fn="expr_UNITI-2.csv",genes,remove_first_column = F)
    
    annots2=read.csv(paste(path,"UNITI-2/UNITI2 biopsy_wk0wk8_582s_CDAI_CRP_FCAL_SESCD.csv",sep=""),stringsAsFactors = F,row.names = 1)[,-2]
    rownames(annots2)=gsub("\\.|-","_",rownames(annots2))
    colnames(exprs_tab2)=gsub("\\.|-","_",colnames(exprs_tab2))
    
    exprs_tab2=exprs_tab2[,rownames(annots2)]
    
    
    save(file=fn,exprs_tab2,annots2)
  }
  else {
    env=new.env()
    load(file=fn,envir=env)
    exprs_tab2=env$exprs_tab2
    annots2=env$annots2
  }
  
  return(list(raw=exprs_tab2,design=annots2))
}

##################################################################################################

zscore_trasformation=function(m){
  return((m-rowMeans(m,na.rm=T))/apply(m,1,sd,na.rm=T))
}



get_scores=function(g1,g2,m,log_scores=F){

  g1=unlist(g1)
  g2=unlist(g2)
  x1=t(m[intersect(g1,rownames(m)),])
  x2=t(m[intersect(g2,rownames(m)),])
 
  score1=rowMeans(x1,na.rm=T)
  score2=rowMeans(x2,na.rm=T)
  if (log_scores){
    score1=log10(1e-3+score1)
    score2=log10(1e-3+score2)
  }
  m2=cbind(score1=score1,score2=score2,projected=NA)
  rownames(m2)=colnames(m)
  return(m2)
}





get_all_zscores=function(m,design,sample_set_list=list()){
  z=m*NA
  for (i in 1:length(sample_set_list)){
    mask= design$dataset%in%sample_set_list[[i]]$datasets&
          design$diagnosis%in%sample_set_list[[i]]$diagnosises&
          design$tissue%in%sample_set_list[[i]]$tissues
    
    if (sum(mask)>=10){
      z[,mask]=zscore_trasformation(m[,mask])
    }

  }
  return(z)
}




plot_dataset=function(bulk,dataset_name="",tissue="Ileum",ngenes=100,genes_to_show=NULL,path=NULL,diagnosis="CD"){
  mask=as.character(bulk$design$dataset)==dataset_name
  statusv=as.character(bulk$design$status)
  statusv[is.na(statusv)]=""
  
  mask_diagnosis=mask&(as.character(bulk$design$diagnosis)==diagnosis)
  mask_tissue=mask_diagnosis&(as.character(bulk$design$tissue)==tissue)
  mask_tissue[is.na(mask_tissue)]=F
  statuses=setdiff(unique(statusv[mask_tissue]),"")
  if (length(statuses)==2){
        
    suffix=paste(dataset_name,diagnosis,tissue,sep="_")
    open_plot(path=path,fn=paste("scatter_",suffix,sep=""),plot_type = "pdf",width=8,height=4)
    layout(matrix(1:2,1,2))
    maski=mask_tissue&(statusv%in%statuses)
    xlim=range(bulk$scores[maski,"score1"],na.rm=T)
    ylim=range(bulk$scores[maski,"score2"],na.rm=T)
  }
  for (status in sort(statuses)){
    maski=mask_tissue&(statusv==status)&(!is.na(bulk$scores[,"score1"]+bulk$scores[,"score2"]))
    col=ifelse(is.na(bulk$design[maski,"response"]),1,ifelse(bulk$design[maski,"response"]=="R",responder_col,non_responder_col))
    if (sum(!is.na(bulk$scores[maski,"score1"]))==0){
      next;
    }
    cort_res=cor.test(bulk$scores[maski,"score1"],bulk$scores[maski,"score2"],method="spe")
    suffix=paste(dataset_name,diagnosis,tissue,status,sep="_")
    if (length(statuses)!=2){
      open_plot(path=,fn=paste("scatter_",suffix,sep=""),plot_type = "pdf",width=4,height=4)
      xlim=range(bulk$scores[maski,"score1"],na.rm=T)
      ylim=range(bulk$scores[maski,"score2"],na.rm=T)
    }
    plot(bulk$scores[maski,"score1"],bulk$scores[maski,"score2"],panel.first=grid(),col=col,pch=20,xlab="score1",ylab="score2",cex.axis=.7,main=status,xlim=xlim,ylim=ylim)
    mtext(paste("rho =",format(cort_res$estimate,digits=2),"; p =",format(cort_res$p.value,digits=2)),3,line = -1,adj = 1,cex=.7)
      

    if (length(statuses)!=2){
      close_plot()
    }
  }
  if (length(statuses)==2){
    close_plot()
  }
  for (status in statuses){
    suffix=paste(dataset_name,diagnosis,tissue,status,sep="_")
    maski=mask_tissue&(statusv==status)
    if (sum(!is.na(bulk$scores[maski,"score1"]))==0){
      next;
    }
    annot_profiles=list()

    gene_groups_heatmaps(bulk$z[,maski],genes_to_show$markers1,genes_to_show$markers2,by=bulk$scores[maski,"score1"],paste("heatmap_",suffix,"_selected_genes",sep=""),annot_profiles = annot_profiles,path=path)
  }
}



gene_groups_heatmaps=function(m,markers1,markers2,by=NULL,s="",reg=1e-5,annot_profiles=NULL,sample_labels=F,show_gene_names=T,zlim=c(-2,2),log_trans=F,show_celltype_bar=T,height=4,path=""){
  markersl1=markers1
  markersl2=markers2
  markers1=unlist(markers1)
  markers2=unlist(markers2)

  genes_to_show=c(markers1,markers2)
  

  if (length(genes_to_show)==0){
    stop("Error! no genes to show!")
  }
  if (!is.null(by)){
    sample_ord= order(by)
  }else{
    sample_ord=1:ncol(m)
  }

  laft_mar=2
  bottom_mar=ifelse(show_gene_names,5,1)
  
  m=as.matrix(m[match(genes_to_show,rownames(m)),sample_ord]) 
  h1=floor(ncol(m)/2)
  if (sum(m[match(markers1,rownames(m)),1:h1],na.rm=T)+sum(m[match(markers2,rownames(m)),-1:-h1],na.rm=T)>sum(m[match(markers1,rownames(m)),-1:-h1],na.rm=T)+sum(m[match(markers2,rownames(m)),1:h1],na.rm=T)){
    m=m[,ncol(m):1]
    sample_ord=rev(sample_ord)
  }
  if (log_trans==T){
    m=log2((reg+m)/(reg+rowMeans(m,na.rm=T)))
  }
  m=m[rowSums(!is.na(m))>0,]
  if (s!=""){
    open_plot(path=path,fn=paste("",s,sep=""),plot_type = "pdf",width=length(genes_to_show)/10,height=height)
  }
  if ((!is.null(annot_profiles))&length(annot_profiles)>0){
    layout_mat=matrix(1:(length(annot_profiles)+1),1,length(annot_profiles)+1)
    par(mar=c(bottom_mar,laft_mar,1,1))
  }
  else{
    layout_mat=matrix(1,1,1)
    par(mar=c(bottom_mar,laft_mar,1,1))
  }
  
  if (show_celltype_bar){
    layout_mat=rbind(c(1,rep(max(layout_mat)+2,ncol(layout_mat)-1)),layout_mat+1)
    layout(layout_mat,widths = c(30,rep(1,length(annot_profiles))),heights=c(1,20))
    par(mar=c(0.2,laft_mar,0.2,1),xpd=F)
    markersl1=sapply(markersl1,function(x){x=x[x%in%rownames(m)]})
    markersl2=sapply(markersl2,function(x){x=x[x%in%rownames(m)]})
    sizes=c(sapply(markersl1,length),sapply(markersl2,length))

    cols=celltypes_cols1[cluster_to_cluster_set_with_pdc[names(cluster_to_subtype1[match(names(sizes),cluster_to_subtype1)])]]
    ctc=celltypes_cols1[names(sizes)]
    barplot(as.matrix(sizes/sum(sizes)), col = cols, axes =F ,beside = F,horiz=T,xaxs='i',yaxs='i',lwd=.5,xlim=c(0,1))
    par(mar=c(bottom_mar,laft_mar,.5,1))
  }
  else if ((!is.null(annot_profiles))&length(annot_profiles)>0){
    layout(layout_mat,widths = c(30,rep(1,length(annot_profiles))))
  }
  
 
  if (nrow(m)==0){
    close_plot()
    return()
  }
  image((as.matrix(m)),axes=F,col=colorRampPalette(c("blue","white","red"))(100),breaks=c(-100,seq(zlim[1],zlim[2],l=99),100))

  abline(v=(par()$usr[2]-par()$usr[1])*(sum(rownames(m)%in%markers1)-.5)/(nrow(m)),lwd=2)
  box()
  if (show_gene_names){
    mtext(text = rownames(m),side = 1,at = seq(0,1,l=nrow(m)),las=2,cex=.5)
  }
  if (sample_labels){
    mtext(text = (colnames(m)),side = 2,at = seq(0,1,l=ncol(m)),las=2,cex=.5,line=0.5)
  }
  if ((!is.null(annot_profiles))&length(annot_profiles)>0){
    for (i in 1:length(annot_profiles)){
      par(mar=c(5,0,1,1))
      u=unique(annot_profiles[[i]])

      image(matrix(match(annot_profiles[[i]][sample_ord],u),1,length(annot_profiles[[i]])),col=u,breaks=0:length(u)+.5,axes=F)
    }
  }
  if (s!=""){
    close_plot()
  }
}









######################################################################################################

load_bulk_datasets=function(load_text=F){
  clin=list()
  scRNAseq_genes=rownames(ileum_ldm$dataset$umitab)
  l=list()
  message("Loading CERTIFI")
  l$certifi=read_certifi(genes=scRNAseq_genes)
#  message("Loading UNITI-1")
#  l$uniti1=read_uniti1(genes=scRNAseq_genes)
  message("Loading UNITI-2")
  l$uniti2=read_uniti2(genes=scRNAseq_genes)
  message("Loading RISK")
  l$risk<-read_RISK(genes=scRNAseq_genes)

  for (i in 1:length(l)){
    if (i==1){
      genes=rownames(l[[i]]$raw)
    }
    else{
      genes=unique(c(genes,rownames(l[[i]]$raw)))
    }
  }  
 
  raw_exprs=c()
  z_exprs=c()
  design=c()
  design_columns=c("subject","geo_accesion","diagnosis","age","tissue")
  
  diagnosis_conversion_table=cbind(c("inflammatory bowel disease","cd","CDc","CDi","iCD","CD","uc","UC","ctrl","IL","IBD-U","Not IBD"),c("CD","CD","CD","CD","CD","CD","UC","UC","CTRL","IL","IBD-U","Not IBD"))
  diagnosis_conversion=diagnosis_conversion_table[,2]
  names(diagnosis_conversion)=diagnosis_conversion_table[,1]
  treatment_conversion=c("Pbo","Pbo","Ust 130 mg","Ust 130 mg","Ust 6 mg/kg","Ust 6 mg/kg")
  names(treatment_conversion)=c("Placebo IV","Pbo","Ustekinumab 130 mg IV","Ust 130 mg","Ustekinumab 6 mg/kg IV","Ust 6 mg/kg")
  tissue_conversion_table=cbind(c('Ascending colon','Rectum',"RECTUM",'Terminal Ileum',"T. ILEUM","ILEUM",'Descending colon','Sigmoid colon','Transverse colon','DELETED','Not Collected','Blood','ileum','colon','sigmoid','CDc','CDi',"iCD",'uc','UC','ctrl','IL','SP. FLEX'),c('Colon','Rectum','Rectum','Ileum','Ileum','Ileum','Colon','Colon','Colon','Exclude','Exclude','Blood','Ileum','Colon','Colon','Colon',"Ileum","Ileum","Colon","Colon","?","?",'SP. FLEX'))
  
   tissue_conversion=tissue_conversion_table[,2]
  names(tissue_conversion)=tissue_conversion_table[,1]
  status_conversion=c("Inflamed","Uninflamed","Inflamed","Inflamed","Uninflamed","Unknown")
  names(status_conversion)=c("Inflamed area","Non-involved area","Endoscopically involved ","Histologically involved","Normal","")
  sex_conversion=c("M","F","M","F","M","F")
  names(sex_conversion)=c(" M"," F","M","F","Male","Female")
  ba_converstion=c("beforeT","beforeT","afterT","afterT")
  names(ba_converstion)=c("I-WEEK 0","I-WK0","I-WEEK 8","I-WK8")
  
  
  ##############
  ##
  ##  certifi
  ##
  ##############
  
  certifi_raw=matrix(NA,length(genes),ncol(l$certifi$raw),dimnames = list(genes,colnames(l$certifi$raw)))
  certifi_raw[rownames(l$certifi$raw),]=as.matrix(l$certifi$raw)
  
  certifi_subject=sapply(strsplit(as.character(l$certifi$design$title),"-"),function(x){x[1]})
  certifi_diagnosis=diagnosis_conversion[l$certifi$design$`diagnosis:ch1`]
  certifi_age=as.numeric(l$certifi$design$`age (yr):ch1`)
  certifi_tissue=tissue_conversion[l$certifi$design$`tissue:ch1`]
  certifi_status=status_conversion[l$certifi$design$`inflam:ch1`]
  certifi_sex=sex_conversion[l$certifi$design$`Sex:ch1`]
  certifi_design2=data.frame(dataset="CERTIFI",subject=certifi_subject,geo_accesion=l$certifi$design$geo_accession,sex=certifi_sex,diagnosis=certifi_diagnosis,age=certifi_age,tissue=certifi_tissue,status=certifi_status,deep_ulcer=NA,BA=NA,treatment=NA,response=NA)
  rownames(certifi_design2)=certifi_design2$geo_accesion
  raw_exprs=cbind(raw_exprs,certifi_raw)
  design=rbind(design,certifi_design2)



  ####################################
  ##
  ## UNITI-1
  ##

  if (!is.null(l$uniti1)){
    uniti1_raw[rownames(l$uniti1$raw),]=as.matrix(l$uniti1$raw)
  
    uniti_subject=sapply(strsplit(as.character(l$uniti1$design$USUBJID),"-"),function(x){x[2]})
    uniti_tissue=tissue_conversion[l$uniti1$design$ANALOC]
    uniti_induction_treatment=treatment_conversion[l$uniti1$design$TR01AG1]
    uniti_treatment=uniti_induction_treatment
    uniti_BA=ba_converstion[l$uniti1$design$visit]
  
    clin$uniti1=l$uniti1$design[,-1:-5]
  
    uniti1_design=data.frame(dataset="UNITI-1",subject=uniti_subject,geo_accesion=NA,sex=NA,diagnosis="CD",age=NA,tissue=uniti_tissue,status="Inflamed",deep_ulcer=NA,BA=uniti_BA,treatment=uniti_treatment,response=NA)
    rownames(uniti1_design)=rownames(l$uniti1$design)
    raw_exprs=cbind(raw_exprs,uniti1_raw)
    design=rbind(design,uniti1_design)
  }
  ####################################
  ##
  ## UNITI-2
  ##
  #################################

  uniti2_raw=matrix(NA,length(genes),ncol(l$uniti2$raw),dimnames = list(genes,colnames(l$uniti2$raw)))
  uniti2_raw[rownames(l$uniti2$raw),]=as.matrix(l$uniti2$raw)
  uniti_subject=sapply(strsplit(as.character(l$uniti2$design$USUBJID),"-"),function(x){x[2]})
  uniti_tissue=tissue_conversion[l$uniti2$design$location]
  uniti_induction_treatment=treatment_conversion[l$uniti2$design$TR01AG1]
  
  uniti_treatment=uniti_induction_treatment
  uniti_BA=ba_converstion[l$uniti2$design$VISIT]
  clin$uniti2=l$uniti2$design[,-1:-5]
  uniti2_design=data.frame(dataset="UNITI-2",subject=uniti_subject,geo_accesion=NA,sex=NA,diagnosis="CD",age=NA,tissue=uniti_tissue,status="Inflamed",deep_ulcer=NA,BA=uniti_BA,treatment=uniti_treatment,response=NA)
 
   rownames(uniti2_design)=rownames(l$uniti2$design)
  raw_exprs=cbind(raw_exprs,uniti2_raw)
  design=rbind(design,uniti2_design)
  

  ################################
  ##
  ##  RISK
  ##  
  ##################################
  
  risk_raw=matrix(NA,length(genes),ncol(l$risk$raw),dimnames = list(genes,colnames(l$risk$raw)))
  risk_raw[rownames(l$risk$raw),]=as.matrix(log(1+l$risk$raw))
  raw_exprs=cbind(raw_exprs,risk_raw)
  risk_design=data.frame(dataset="RISK",subject=l$risk$design$RISK.DNA.ID,geo_accesion=l$risk$design$GEO_ID,sex= sex_conversion[l$risk$design$Gender],diagnosis=diagnosis_conversion[as.character(l$risk$design$Disease.status)],age=l$risk$design$AgeDxYrs.y,tissue="Ileum",status=status_conversion[l$risk$design$Inflammation],deep_ulcer=l$risk$design$Ulceration,BA=NA,treatment=NA,response=l$risk$design$response,B_Dx=l$risk$design$B_Dx,B_Cur=l$risk$design$B_Cur)
  rownames(risk_design)=rownames(l$risk$design)
  design=rbind(cbind(design,B_Dx=NA,B_Cur=NA),risk_design)
  clin$risk=transform(merge(l$risk$clin,l$risk$pcdai,by=0,all=TRUE), row.names=Row.names, Row.names=NULL)
  z_exprs=get_all_zscores(raw_exprs,design,sample_set_list)

  return(list(raw=raw_exprs,z=z_exprs,design=design,clin=clin))
}

################################################################################################


risk_clinical_data=function(bulk){
  columns_res=grep("RES",colnames(bulk$clin$risk),val=T)
  profs=gsub("RES","",columns_res)
  for (i in 1:length(profs)){
    open_plot(path=supp_figures_path,fn = paste("risk_",profs[i],"_vs_scores",sep=""),plot_type = "pdf",width = 8,height = 4)
    layout(matrix(1:2,1,2))
    plot(bulk$scores[rownames(bulk$clin$risk),"score1"],bulk$clin$risk[,columns_res[i]],xlab="Score 1",ylab=profs[i])
    cor_res=cor.test(bulk$scores[rownames(bulk$clin$risk),"score1"],bulk$clin$risk[,columns_res[i]])
    mtext(side=3,paste("r=",round(cor_res$estimate,digits=2)," ; p=",round(cor_res$p.value,digits=4),sep=""))
    plot(bulk$scores[rownames(bulk$clin$risk),"score2"],bulk$clin$risk[,columns_res[i]],xlab="Score 2",ylab=profs[i])
    cor_res=cor.test(bulk$scores[rownames(bulk$clin$risk),"score2"],bulk$clin$risk[,columns_res[i]])
    mtext(side=3,paste("r=",round(cor_res$estimate,digits=2)," ; p=",round(cor_res$p.value,digits=4),sep=""))
    dev.off()  
  }
}

select_genes=function(max_per_clusterset=30,thresh_cluster=0.5,thresh_pattern1=1,thresh_pattern2=0.5,thresh_freq=0,cluster_to_subtype,fdr_thresh=0.01,
                           included_cluster_sets1,included_cluster_sets2,samples=inflamed_samples,restrict_to_included_cluster_sets=T,trace=F)
{

  counts=pmax(apply(ileum_ldm$dataset$counts[samples,,],2:3,sum)-apply(ileum_ldm$dataset$noise_counts[samples,,],2:3,sum),0)
  s=aggregate.Matrix(t(counts[,intersect(colnames(counts),names(cluster_to_subtype))]),groupings = cluster_to_subtype[intersect(colnames(counts),names(cluster_to_subtype))],fun = "sum")+.1
  s_ep=rowSums(sapply(ileum_ldm$dataset$gated_out_umitabs$EP,rowSums))
  m=cbind(t(s/rowSums(s)),EP=s_ep/sum(s_ep))
  ml=(log2(1e-6+m))
  
  de_markers=setdiff(rownames(de_res)[(de_res[,"log2_FC"]>thresh_pattern1|abs(de_res[,"log2_FC"])>thresh_pattern2)&de_res[,"adj.p.value"]<fdr_thresh&(de_res[,"freq_bg"]>thresh_freq|de_res[,"freq_fg"]>thresh_freq)],excluded_genes)

  if (restrict_to_included_cluster_sets){
    ml.max1=apply(ml[,included_cluster_sets1,drop=F],1,max)
    ml.max2=apply(ml[,included_cluster_sets2,drop=F],1,max)
    ml.max_other=apply(ml[,setdiff(colnames(ml),c(included_cluster_sets1,included_cluster_sets2)),drop=F],1,max)
  
    markers=names(which(abs(ml.max1[de_markers]-ml.max2[de_markers])>thresh_cluster&(abs(pmax(ml.max1[de_markers],ml.max2[de_markers])-ml.max_other[de_markers])>0)))
  }
  else{
    markers=de_markers
  }
  
  
  marker_to_cluster_set=colnames(ml)[apply(ml[markers,],1,which.max)]
  names(marker_to_cluster_set)=markers
  if (trace){
    print(table(colnames(ml)[apply(ml[markers[de_res[markers,"log2_FC"]<0],,drop=F],1,which.max)]))
    print(table(colnames(ml)[apply(ml[markers[de_res[markers,"log2_FC"]>0],,drop=F],1,which.max)]))
  }
  

    marker_to_cluster_set1=marker_to_cluster_set[markers[de_res[markers,"log2_FC"]>0]]
    marker_to_cluster_set2=marker_to_cluster_set[markers[de_res[markers,"log2_FC"]<0]]
  

  short_markers_list1=sapply(split(names(marker_to_cluster_set1),marker_to_cluster_set1),function(x){x[order(abs(de_res[x,"log2_FC"]),decreasing = T)][1:pmin(length(x),max_per_clusterset)]},simplify = F)
  short_markers_list2=sapply(split(names(marker_to_cluster_set2),marker_to_cluster_set2),function(x){x[order(abs(de_res[x,"log2_FC"]),decreasing = T)][1:pmin(length(x),max_per_clusterset)]},simplify = F)
  short_markers2=sapply(short_markers_list2,function(x){x[de_res[x,"log2_FC"]<0]})
  short_markers1=sapply(short_markers_list1,function(x){x[de_res[x,"log2_FC"]>0]})
  short_markers1=short_markers1[sapply(short_markers1,length)>0]
  short_markers2=short_markers2[sapply(short_markers2,length)>0]

  return(list(markers1=short_markers1[included_cluster_sets1],markers2=short_markers2[included_cluster_sets2]))
}



#############################################################################################


plot_figure_5b=function(bulk,status="Inflamed"){
  before_mask=bulk$design$BA=="beforeT"
  before_mask[is.na(before_mask)]=T
  mask=bulk$design$tissue=="Ileum"&bulk$design$status==status&bulk$design$diagnosis=="CD"&before_mask

  #datasets=c('RISK',"CERTIFI",'UNITI-1', 'UNITI-2')
  datasets=c('RISK',"CERTIFI", 'UNITI-2')
  datasets=datasets[table(bulk$design[mask,"dataset"])[datasets]>0]

  cols=alt_cols[c(-3,-5)][1:length(datasets)]
  par(mar=c(4,4,4,4))
  plot(bulk$scores[mask,"score2"],bulk$scores[mask,"score1"],col=cols[match(bulk$design$dataset[mask],datasets)],pch=20,ylim=c(-1.2,2.3),xlim=c(-2,1),panel.first=grid(lty=1),xlab="Score 2",ylab="Score 1",cex=.8)
  legend("bottomleft",paste(datasets," (",table(bulk$design$dataset[mask])[datasets],")",sep=""),col=cols[1:length(datasets)],pch=20,cex=1)
  res=list()
  
  for (dataset in datasets){
    maski=mask&bulk$design$dataset==dataset
    res[[dataset]]=unlist(cor.test(bulk$scores[maski,"score1"],bulk$scores[maski,"score2"],method = "pearson")[c("p.value","estimate")])
  }
  cor_stats=as.data.frame.list(res)
  print(cor_stats)
  print(cor.test(bulk$scores[mask,"score1"],bulk$scores[mask,"score2"],method = "spe")[c("p.value","estimate")])
}







plot_risk_figures=function(bulk){
  mask_risk=bulk$design$dataset=="RISK"&bulk$design$tissue=="Ileum"&bulk$design$status=="Inflamed"&bulk$design$diagnosis=="CD"
  mask_risk[is.na(mask_risk)]=F
  by= bulk$scores[mask_risk,"score1"]
  patient_groups=cut(by,quantile(by,0:2/2),labels = c("low","high"))
  l_score=split(by,bulk$design$response[mask_risk])

  simple_roc <- function(labels, scores){
    labels <- labels[order(scores, decreasing=TRUE)]
    data.frame(TPR=cumsum(labels)/sum(labels), FPR=cumsum(!labels)/sum(!labels), labels)
  }
  
  simple_auc <- function(m){
    TPR=m[,1]
    FPR=m[,2]
    # inputs already sorted, best scores first 
    dFPR <- c(diff(FPR), 0)
    dTPR <- c(diff(TPR), 0)
    sum(TPR * dFPR) + sum(dTPR * dFPR)/2
  }
  
  response1=bulk$design$response[mask_risk]
  score1=by[!is.na(response1)]
  response1=response1[!is.na(response1)]
  open_plot(path=supp_figures_path,fn ="risk_roc",plot_type = "pdf",width = 6,height = 6)
  par(mar=c(5,5,1,1))
  auc=simple_auc(simple_roc(response1=="R",-1*score1))
  auc_resamp=replicate(1e4,simple_auc(simple_roc(sample(response1=="R",length(response1),replace = F),-1*score1)))
  auc_pvalue=mean(auc_resamp>auc)
  plot(rbind(c(0,0),simple_roc(response1=="R",-1*score1)[,c(2,1)]),type='l',panel.first=grid(lty=3,col="gray"),main=paste(round(auc,digits=1),"p=",round(auc_pvalue,digits=3)))
 
  abline(0,1,lty=3,col=1)
  close_plot()

  open_plot(path=main_figures_path,fn ="risk_response",plot_type = "pdf",width = 4,height = 4)
#  layout(matrix(1:2,1,2))
  par(mar=c(4,4,4,4))
  plot(1,col=0,xlim=c(-1,1),ylim=c(0,1),panel.first=grid(lty=1),ylab="Cumulative Fraction",xlab="Projection score")
  quants=sapply(l_score,quantile,0:100/100)
  matplot(quants,0:100/100,type='l',lty=1,lwd=2,col=c(non_responder_col,responder_col),xlim=c(-1,1),cex=.7,add=T)

  ks_res=ks.test(l_score$NR,l_score$R,alternative = "less")
  text(-1,0.95,labels = paste("KS D=",round(ks_res$statistic,digits=2),"; p=",round(ks_res$p.value,digits = 4)),pos=4,cex=0.7)
  legend("bottomright",legend=c("Non Responders","Responders"),col=c(non_responder_col,responder_col),lty=1,lwd=2,cex=.7)
  close_plot()
  
  
  n.x <- length(l_score$R)
  n.y <- length(l_score$NR)
  n <- n.x * n.y/(n.x + n.y)
  w <- c(l_score$R, l_score$NR)
  z <- cumsum(ifelse(order(w) <= n.x, 1/n.x, -1/n.y))
 # max(abs(z)) # Dmax
  max.at <- sort(w)[which(abs(z) == max(abs(z)))]
#  message("D is max at",max.at)
  open_plot(path=main_figures_path,fn ="risk_response_threshold",plot_type = "pdf",width = 4,height = 4)
  par(mar=c(4,4,2,2))
  plot(1,col=0,xlim=c(-1,1),ylim=c(0,1),panel.first=grid(lty=1),ylab="Cumulative Fraction",xlab="Projection score")
  quants=sapply(l_score,quantile,0:100/100)
  matplot(quants,0:100/100,type='l',lty=1,lwd=2,col=c(non_responder_col,responder_col),xlim=c(-1,1),cex=.7,add=T)
  ks_res=ks.test(l_score$NR,l_score$R,alternative = "less")
  abline(v=max.at,col=1,lty=3)
  text(min(unlist(l_score))+.1,.95,labels = paste("KS D=",round(ks_res$statistic,digits=2),"; p=",round(ks_res$p.value,digits = 4)),pos=4,cex=1)
  legend("bottomright",legend=c("Non Responders","Responders"),col=c(non_responder_col,responder_col),lty=1,lwd=2,cex=.8)
   close_plot()
  tab=table(bulk$scores[mask_risk,"score1"]>max.at,bulk$design$response[mask_risk])
  
  open_plot(path=main_figures_path,fn ="risk_response_piecharts",plot_type = "pdf",width = 6,height = 6)
  par(mar=c(1,1,4,1))
  layout(matrix(1:2,2,1))
  pie(tab[1,],c("Non responders","Responders"),col=c(non_responder_col,responder_col),main=paste("Signature score <",round(max.at,digits=2)))
  pie(tab[2,],c("Non responders","Responders"),col=c(non_responder_col,responder_col),main=paste("Signature score >",round(max.at,digits=2)))
  close_plot()
  
 

}


response_cumulative=function(vr,vnr){
    quants_R=quantile(vr,0:100/100)
    quants_NR=quantile(vnr,0:100/100)
    quants=cbind(quants_NR,quants_R)
    plot(1,col=0,xlim=c(-1,1),ylim=c(0,1),panel.first=grid(lty=1),ylab="Cumulative Fraction",xlab="Projection score")
    matplot(quants,0:100/100,type='l',lty=1,lwd=2,col=c(non_responder_col,responder_col),cex=.7,add=T)
    if (length(vnr)>=5&length(vr)>=5){
      ks_res=ks.test(vnr,vr,alternative = "two.sided")
    }else{
      ks_res=list(statistic=NA,p.value=NA)
      
    }
    text(-.95,1,labels = paste("KS D=",round(ks_res$statistic,digits=2),"; p=",round(ks_res$p.value,digits = 4)),pos=4,cex=1)
    legend("bottomright",legend=c("Non Responders","Responders"),col=c(non_responder_col,responder_col),lty=1,lwd=2,cex=.8)
  return(ks_res)
}



#######################################################################################################################################


make_figure5=function(load_normalized_data=F){
  
  
  lm=ileum_ldm
  score_thresh<<-list()  
  responder_col<<-brew_cols[5]
  non_responder_col<<-brew_cols[4]
  if (load_normalized_data){
    lbd_res<-load_bulk_datasets(load_text)
    mask=rownames(lbd_res$design[lbd_res$design$diagnosis=="CD"&lbd_res$design$tissue=="Ileum"&(lbd_res$design$BA=="beforeT"|is.na(lbd_res$design$BA)),])
    lbd_res2=list()
    lbd_res2$design=lbd_res$design[mask,]
    lbd_res2$raw=lbd_res$raw[,mask]
    lbd_res2$z=lbd_res$z[,mask]
    lbd_res2$design$deep_ulcer=NULL
    lbd_res2$design$B_Cur=NULL
    lbd_res2$design$B_Dx=NULL
    lbd_res2$design$age=NULL
    bulk<<-lbd_res2
    
    
     
    save(file=paste(pipeline_path,"input","external_cohorts_data","bulk_data.rd",sep="/"),bulk)
  }
  else {
    load(paste(pipeline_path,"input","external_cohorts_data","bulk_data.rd",sep="/"))
  }
  de_res<<-read.csv(paste(pipeline_path,"input/DE/DE_inf_pat1_vs_pat2_total.csv",sep="/"),row.names = 1)
  
  
  inf_pat1<<-inflamed_samples_v2[sample_to_patient[inflamed_samples_v2]%in%pat1]
  inf_pat2<<-inflamed_samples_v2[sample_to_patient[inflamed_samples_v2]%in%pat2]
  uninf_pat1<<-uninflamed_samples_v2[sample_to_patient[uninflamed_samples_v2]%in%pat1]
  uninf_pat2<<-uninflamed_samples_v2[sample_to_patient[uninflamed_samples_v2]%in%pat2]
  ncells_per_sample=table(lm$dataset$cell_to_sample)
  
  s1=apply(lm$dataset$counts[inf_pat1,,],2,sum)
  s2=apply(lm$dataset$counts[inf_pat2,,],2,sum)
  ms=apply(lm$dataset$counts[c(inf_pat1,inf_pat2),,],1:2,sum)
  mm=ms/rowSums(ms)
  mm1=1e-6+t(mm[inf_pat1,])
  mm2=1e-6+t(mm[inf_pat2,])
  
  afc=log2(rowMeans(mm2)/rowMeans(mm1))
  
  s=apply(lm$dataset$counts,2,sum)
  ms=t(apply(lm$dataset$counts,1:2,sum))
  ms=t(t(ms)/colSums(ms))
  
  mz=(ms-rowMeans(ms))/apply(ms,1,sd)
  mz=mz[,inflamed_samples]
  colnames(mz)=sample_to_patient[colnames(mz)]
  p_thresh=1e-4
  
  de_res2=de_res[de_res$adj.p.value<p_thresh,]
  s_ep=rowSums(sapply(ileum_ldm$dataset$gated_out_umitabs$EP,rowSums))
  
  models=cbind(ileum_ldm$model$models,EP=s_ep/sum(s_ep))
  cluster_to_cluster_set=c(cluster_to_cluster_set,"EP")
  names(cluster_to_cluster_set)[cluster_to_cluster_set=="EP"]="EP"

  scores_genes<-select_genes(max_per_clusterset =20,thresh_cluster=3,thresh_pattern1=1,thresh_pattern2 =1,thresh_freq=1e-6,cluster_to_subtype=cluster_to_subtype1,fdr_thresh = 0.001,
                              included_cluster_sets1=c("Inf. Macrophages","Activated DC","pDC","IgG plasma cells","Naive/CM T cells","Tregs","ACKR1+ endothelial cells","Activated fibroblasts"),
                              included_cluster_sets2=c("Resident macrophages","moDC","IgA plasma cells","ILC3","Cytokines low Trm","Type 1 cytokines Trm","Type 3 cytokines Trm","Enteric neurons","CD36+ endothelial cells", "Fibroblasts"),
                              trace=T,restrict_to_included_cluster_sets = T)
  
  
  bulk$scores<-get_scores(scores_genes$markers1,scores_genes$markers2,bulk$z) 
  
  scRNA_scores<-get_scores(scores_genes$markers1,scores_genes$markers2,mz)

  write.csv(file=paste(pipeline_path,"/output//tables/score_genes_figure_5d.csv",sep=""),rbind(cbind(geneSymbol=unlist(scores_genes$markers1[sapply(scores_genes$markers1,length)>0]),group="with GIMATS"),cbind(geneSymbol=unlist(scores_genes$markers2[sapply(scores_genes$markers2,length)>0]),group="without GIMATS")),quote=F,row.names=T)
  
  open_plot(main_figures_path,fn="figure_5a",plot_type="pdf",width = 5,height = 5)
  par(mar=c(4,4,4,4))
  plot(de_res$log2_FC,-log2(de_res$adj.p.value),xlim=c(-1,1),xlab="Log2(freq(pat2)/freq(pat1))",ylab="FDR adjusted p.value (BH)",cex=.7)
  close_plot()
  
  open_plot(main_figures_path,fn="figure_5b",plot_type="pdf",width = 5,height = 5)
  plot_figure_5b(bulk,"Inflamed")
  close_plot()
  
  open_plot(main_figures_path,fn="figure_5c",plot_type="pdf",width = 5,height = 5)
  plot_figure_5b(bulk,"Uninflamed")
  close_plot()
  
  
 
  plot_dataset(bulk,"RISK",tissue="Ileum",genes_to_show = scores_genes,path=main_figures_path)
  
  if (!is.null(bulk$clin)){
    risk_clinical_data(bulk)
  }
  plot_risk_figures(bulk)
  
  plot_dataset(bulk,"CERTIFI",tissue="Ileum",genes_to_show = scores_genes,,path=supp_figures_path)

}

