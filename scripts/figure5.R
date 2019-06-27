
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


de_res=read.csv("~/GoogleDrive/work/shared/analysis/analysis_iCD_paper_2018/output/tables/DE_patterns/DE_inf_pat1_vs_pat2_total.csv",row.names = 1)


sample_set_list=list(
  list(datasets="CERTIFI",diagnosises="CD",tissues=c("Ileum")),
  list(datasets="UNITI-1",diagnosises="CD",tissues=c("Ileum")),
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
  
  return(list(raw=m[,mask],design=sample_tab[mask,],clin=clin2[mask,],pcdai=pcdai))
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
  path="~/GoogleDrive/work/shared/data/public/human/Peters_et_al_NG_2017/GSE100833_RAW/"
  
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




quantile_normalization=function(raw_data,reg=10){
  m=t(t(raw_data)/colSums(raw_data))
  ranked=apply(m,2,rank)
  normed_data=-1*log2(.9999-(reg+ranked)/(reg+nrow(ranked)))
  return(normed_data)
}

zscore_trasformation=function(m){
  return((m-rowMeans(m,na.rm=T))/apply(m,1,sd,na.rm=T))
}


get_score_per_gl=function(gl,m){
  get_avg=function(g){
    return(colMeans(m[intersect(g,rownames(m)),,drop=F],na.rm=T))
  }
  x=sapply(gl,get_avg)
  colnames(x)=names(gl)
  return(x)
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



read_modules=function(){
  tmp_mod_tab=read.table("tables/modules_all_data.txt",stringsAsFactors = F)
  mods_l=strsplit(tmp_mod_tab[,2],",")
  names(mods_l)=tmp_mod_tab[,1]
  mods=rep(names(mods_l),sapply(mods_l,length))
  names(mods)=unlist(mods_l)
  return(mods)
}


get_all_zscores=function(m,design,sample_set_list=list()){
  z=m*NA
  for (i in 1:length(sample_set_list)){
    mask= design$dataset%in%sample_set_list[[i]]$datasets&
          design$diagnosis%in%sample_set_list[[i]]$diagnosises&
          design$tissue%in%sample_set_list[[i]]$tissues
    
    if (sum(mask)>=10){
      message("Zscoring sample set ",i)
      z[,mask]=zscore_trasformation(m[,mask])
    }
    else{
    
      message("Sample set ",i,"was not zscored")
    }
  }
  return(z)
}



project_scores=function(m){
  a=as.matrix(m)
  rotation=prcomp(a)$rotation[1,]
  proj<-(a%*%rotation)[,1]
  if (rotation[1]<0){
    proj<-proj*-1
  }

  return(proj)
}



plot_dataset=function(dataset_name="",ngenes=100,plot_ulceration=F,plot_bcur=F,genes_to_show=NULL){
  mask=as.character(bulk$design$dataset)==dataset_name
  statusv=as.character(bulk$design$status)
  statusv[is.na(statusv)]=""
  
  for (diagnosis in unique(as.character(bulk$design[mask,]$diagnosis))){
    mask_diagnosis=mask&(as.character(bulk$design$diagnosis)==diagnosis)
    tissues=unique(as.character(bulk$design[mask_diagnosis,]$tissue))
    tissues=tissues[!is.na(tissues)]
    for (tissue in tissues){
      mask_tissue=mask_diagnosis&(as.character(bulk$design$tissue)==tissue)
      mask_tissue[is.na(mask_tissue)]=F
      statuses=setdiff(unique(statusv[mask_tissue]),"")
      if (length(statuses)==2){
        
        suffix=paste(dataset_name,diagnosis,tissue,sep="_")
        open_plot(path="output/bulk/",fn=paste("scatter_",suffix,sep=""),plot_type = "pdf",width=8,height=4)
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
          open_plot(path="output/bulk/",fn=paste("scatter_",suffix,sep=""),plot_type = "pdf",width=4,height=4)
          xlim=range(bulk$scores[maski,"score1"],na.rm=T)
          ylim=range(bulk$scores[maski,"score2"],na.rm=T)
        }
        plot(bulk$scores[maski,"score1"],bulk$scores[maski,"score2"],panel.first=grid(),col=col,pch=20,xlab="score1",ylab="score2",cex.axis=.7,main=status,xlim=xlim,ylim=ylim)
        mtext(paste("rho =",format(cort_res$estimate,digits=2),"; p =",format(cort_res$p.value,digits=2)),3,line = -1,adj = 1,cex=.7)
      
        print(suffix)
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
        if (plot_ulceration){
          annot_profiles$ulceration=ifelse(bulk$design$deep_ulcer[maski]=="no deep ulcer",responder_col,non_responder_col)
        }
        if (plot_bcur){
          annot_profiles$B_cur=c("white","blue","red","black")[match(bulk$design$B_Cur[maski],c("B1","B2","B3","B2+B3"))]
        }


        gene_groups_heatmaps(bulk$z[,maski],genes_to_show$markers1,genes_to_show$markers2,by=bulk$scores[maski,"score1"],paste("heatmap_",suffix,"_selected_genes",sep=""),annot_profiles = annot_profiles)
      }
     
    }
  }
}



gene_groups_heatmaps=function(m,markers1,markers2,by=NULL,s="",reg=1e-5,annot_profiles=NULL,sample_labels=F,show_gene_names=T,zlim=c(-2,2),log_trans=F,show_celltype_bar=T,height=4){
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
    open_plot(path="output/bulk/",fn=paste("",s,sep=""),plot_type = "pdf",width=length(genes_to_show)/10,height=height)
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
  #  browser()
    cols=celltypes_cols1[cluster_to_cluster_set_with_pdc[names(cluster_to_subtype1[match(names(sizes),cluster_to_subtype1)])]]
    ctc=celltypes_cols1[names(sizes)]
  #  if (any(!is.na(ctc))){
  #    cols[!is.na(ctc)]=ctc[!is.na(ctc)]
  #  }
    rep(1:length(sizes),each=sizes)
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


plot_arijs_response=function(){
  
  mask_before=(bulk$design$dataset=="Arijs")&bulk$design$diagnosis=="CD"&bulk$design$tissue%in%c("Ileum","Colon")&bulk$design$BA=="beforeT"
  mask_after=(bulk$design$dataset=="Arijs")&bulk$design$diagnosis=="CD"&bulk$design$tissue%in%c("Ileum","Colon")&bulk$design$BA=="afterT"
  mask_after=rownames(bulk$design[mask_after,])[match(bulk$design[mask_before,"subject"],bulk$design[mask_after,"subject"])]
  mask_before=rownames(bulk$design[mask_before,]);mask2=!(is.na(mask_before)|is.na(mask_after));mask_before=mask_before[mask2]
  mask_after=mask_after[mask2]
  open_plot("output/main_figures/",fn="figure_4g_arijs_cd_before_delta_one_scatter",plot_type="pdf",width = 4,height = 4)
  par(mar=c(4,4,4,4))
  plot(bulk$scores[mask_before,"score1"],bulk$scores[mask_after,"score1"]-bulk$scores[mask_before,"score1"],xlab="before=T",ylab="afterT minus beforeT",pch=c(20,18)[match(bulk$design[mask_before,"tissue"],c("Ileum","Colon"))],col=ifelse(bulk$design[mask_before,"response"]=="R",responder_col,non_responder_col),panel.first = {grid(lty=1);abline(0,0,lty=2)})
  close_plot()
  open_plot("output/main_figures/",fn="figure_4g_arijs_cd_before_after_one_scatter",plot_type="pdf",width = 4,height = 4)
  par(mar=c(4,4,4,4))
  lim=range(c(bulk$scores[mask_before,"score1"],bulk$scores[mask_after,"score1"]))
  plot(bulk$scores[mask_before,"score1"],bulk$scores[mask_after,"score1"],xlab="before=T",ylab="afterT",pch=c(20,18)[match(bulk$design[mask_before,"tissue"],c("Ileum","Colon"))],xlim=lim,ylim=lim,col=ifelse(bulk$design[mask_before,"response"]=="R",responder_col,non_responder_col),panel.first = {grid(lty=1);abline(0,1,lty=2)})
  close_plot()
  
  
  
  mask_before=(bulk$design$dataset=="Arijs")&bulk$design$diagnosis=="UC"&bulk$design$tissue%in%c("Ileum","Colon")&bulk$design$BA=="beforeT"
  mask_after=(bulk$design$dataset=="Arijs")&bulk$design$diagnosis=="UC"&bulk$design$tissue%in%c("Ileum","Colon")&bulk$design$BA=="afterT"
  mask_after=rownames(bulk$design[mask_after,])[match(bulk$design[mask_before,"subject"],bulk$design[mask_after,"subject"])]
  mask_before=rownames(bulk$design[mask_before,]);mask2=!(is.na(mask_before)|is.na(mask_after));mask_before=mask_before[mask2]
  mask_after=mask_after[mask2]
  open_plot("output/supp_figures/",fn="figure_sx_arijs_uc_before_delta_one_scatter",plot_type="pdf",width = 4,height = 4)
  par(mar=c(4,4,4,4))
 
  plot(bulk$scores[mask_before,"score1"],bulk$scores[mask_after,"score1"]-bulk$scores[mask_before,"score1"],xlab="before=T",ylab="afterT minus beforeT",pch=c(20,18)[match(bulk$design[mask_before,"tissue"],c("Ileum","Colon"))],col=ifelse(bulk$design[mask_before,"response"]=="R",responder_col,non_responder_col),panel.first = {grid(lty=1);abline(0,0,lty=2)})
  close_plot()
  open_plot("output/supp_figures/",fn="figure_sx_arijs_uc_before_after_one_scatter",plot_type="pdf",width = 4,height = 4)
  par(mar=c(4,4,4,4))
  lim=range(c(bulk$scores[mask_before,"score1"],bulk$scores[mask_after,"score1"]))
  plot(bulk$scores[mask_before,"score1"],bulk$scores[mask_after,"score1"],xlab="before=T",ylab="afterT",pch=c(20,18)[match(bulk$design[mask_before,"tissue"],c("Ileum","Colon"))],xlim=lim,ylim=lim,col=ifelse(bulk$design[mask_before,"response"]=="R",responder_col,non_responder_col),panel.first = {grid(lty=1);abline(0,1,lty=2)})
  close_plot()
}





plot_dataset_with_stratification=function(dataset_name="Arijs",tissues,pool_tissues=F,ngenes=100,genes_to_show=c(),plot_cor=T,stratify_by="ba"){
  plot_tissue=function(mask_tissue,statuses,file_prefix){
    for (status in statuses){
      if (status=="Uninflamed"){
        next()
      }
      suffix=paste(file_prefix,status,sep="_")
      print(suffix)
      
      
      maski=mask_tissue&(statusv==status)
      
      
      mask_stratify_list=list()
      if (stratify_by=="ba"){
        bav=bulk$design$BA
        bav[is.na(bav)]=""
        
        mask_stratify_list$beforeT=maski&bav=="beforeT"
        mask_stratify_list$afterT=maski&bav=="afterT"
      }
      
      if (stratify_by=="treatment"){
        treatmentv=bulk$design$treatment
        treatmentv[is.na(treatmentv)]=""

              mask_stratify_list$Pbo=maski&treatmentv=="Pbo"
        mask_stratify_list$Ust_6_mg_kg=maski&treatmentv=="Ust 6 mg/kg"
        mask_stratify_list$Ust_130_mg=maski&treatmentv=="Ust 130 mg"
      }
      
      open_plot(path="output/bulk/",fn=paste("scatter_",suffix,sep=""),plot_type = "pdf",width=4*length(mask_stratify_list),height=4)
      par(mar=c(4,4,4,4))
      layout(matrix(1:length(mask_stratify_list),1,length(mask_stratify_list)))
      
      xlim=range(bulk$scores[mask_tissue&(statusv==status),"score1"],na.rm=T)
      ylim=range(bulk$scores[mask_tissue&(statusv==status),"score2"],na.rm=T)
      
      for (i in 1:length(mask_stratify_list)){
        mask_stratify=mask_stratify_list[[i]]
        mask_stratify=mask_stratify&!(is.na(bulk$scores[,"score1"])|is.na(bulk$scores[,"score2"]))
        if (sum(mask_stratify)<10){
          plot_cor=F
        }
        if (plot_cor){
          cort_res=cor.test(bulk$scores[mask_stratify,"score1"],bulk$scores[mask_stratify,"score2"],method="spe")
        }
        response=ifelse(bulk$design$response[mask_stratify]=="R",responder_col,non_responder_col)
        response[is.na(response)]=1
        if (sum(mask_stratify)>0){
          tissues=unique(as.character(bulk$design$tissue[mask_stratify]))
          if ("Ileum"%in%tissues){
            tissues=c("Ileum",tissues[tissues!="Ileum"])
          }
          pchs=c(20,18,15,16,17)[match(bulk$design$tissue[mask_stratify],tissues)]
          plot(bulk$scores[mask_stratify,"score1"],bulk$scores[mask_stratify,"score2"],panel.first=grid(),xlab="score1",ylab="score2",cex.axis=.7,col=response,xlim=xlim,ylim=ylim,main=names(mask_stratify_list)[i],pch=pchs)
          if (length(tissues)>1){
            legend("bottomleft",legend = tissues,pch=c(20,18,15,16,17)[1:length(tissues)],col=1,cex=.7)
          }
          if (plot_cor){
            mtext(paste("rho =",format(cort_res$estimate,digits=2),"; p =",format(cort_res$p.value,digits=2)),3,line = -1,adj = 1,cex=.7)
          }
        }
        else {
          plot.new()
        }
      }
      close_plot()
      for (i in 1:length(mask_stratify_list)){
        mask_stratify=mask_stratify_list[[i]]
        mask_stratify=mask_stratify&!(is.na(bulk$scores[,"score1"])|is.na(bulk$scores[,"score2"]))
        if (sum(mask_stratify)==0){
          next
        }
        suffix=paste(dataset_name,diagnosis,status,names(mask_stratify_list)[i],sep="_")
        annot_profiles=list()
        annot_profiles$response=ifelse(bulk$design$response[mask_stratify]=="R",responder_col,non_responder_col)
        
        gene_groups_heatmaps(bulk$z[,mask_stratify],genes_to_show$markers1,genes_to_show$markers2,by=bulk$scores[mask_stratify,"score1"],paste("heatmap_",suffix,"_selected_genes",sep=""),annot_profiles=annot_profiles)
      }
    }
  }
  
  
  mask=bulk$design$dataset==dataset_name
  statusv=as.character(bulk$design$status)
  statusv[is.na(statusv)]=""
  for (diagnosis in unique(bulk$design[mask,]$diagnosis)){
    mask_diagnosis=mask&(bulk$design$diagnosis==diagnosis)
  
      mask_tissue=mask_diagnosis&(bulk$design$tissue%in%tissues)
      mask_tissue[is.na(mask_tissue)]=F
      statuses=setdiff(unique(statusv[mask_tissue]),"")
      
      if (pool_tissues){
        plot_tissue(mask_tissue,statuses,file_prefix=paste(dataset_name,diagnosis,sep="_"))
      }
      else {
        for (tissue in tissues){
          mask_tissue=mask_diagnosis&(bulk$design$tissue==tissue)
          mask_tissue[is.na(mask_tissue)]=F
          file_prefix=paste(dataset_name,diagnosis,sep="_")
          plot_tissue(mask_tissue,statuses,file_prefix=paste(dataset_name,diagnosis,tissue,sep="_"))
        }
      }
  }
}


get_module_scores=function(mat,mods){
  mask=intersect(rownames(mat),names(mods))
  modsums=t(sapply(split(as.data.frame(as.matrix(mat[mask,])),mods[mask],drop=F),colSums,na.rm=T))
  mask_tot=intersect(rownames(mat),rownames(ileum_ldm$dataset$umitab))
  modfreqs=modsums/colSums(mat[mask_tot,],na.rm=T)
  return(modfreqs)
}






get_module_scores=function(mat,mods){
  mask=intersect(rownames(mat),names(mods))
  modsums=t(sapply(split(as.data.frame(as.matrix(mat[mask,])),mods[mask],drop=F),colSums,na.rm=T))
  mask_tot=intersect(rownames(mat),rownames(ileum_ldm$dataset$umitab))
  modfreqs=modsums/colSums(mat[mask_tot,],na.rm=T)
  return(modfreqs)
}

plot_modules_heatmap=function(normed,order_patiens_by,annot_cols=NULL,reg=1e-5){
  mods_scores=get_module_scores(normed,mods)
  if (!is.null(annot_cols)){
    layout(matrix(1:2,1,2),widths = c(20,1))
    par(mar=c(5,7,1,1))
  }
  else{
    par(mar=c(5,7,1,1))
  }
  patient_ord=order(order_patiens_by)
#  mods_ord=order(as.numeric(rownames(mods_scores)))

  mods_scores_normed=log2((reg+mods_scores)/(reg+rowMeans(mods_scores)))
  mod_ord=order(as.numeric(rownames(mods_scores_normed)))
  image(mods_scores_normed[mod_ord,patient_ord],col=greenred(100),axes=F,breaks=c(-1e8,seq(-1.5,1.5,l=99),1e8))
  mtext(text = paste("Module",rownames(mods_scores_normed))[mod_ord],side = 1,at = seq(0,1,l=nrow(mods_scores_normed)),las=2,cex=.5)
  mtext(text = colnames(mods_scores_normed)[patient_ord],side = 2,at = seq(0,1,l=ncol(mods_scores_normed)),las=2,cex=.5)
  if (!is.null(annot_cols)){
    par(mar=c(5,0,1,1))
    u=unique(annot_cols)
    image(matrix(match(annot_cols[patient_ord],u),1,length(annot_cols)),col=u,breaks=0:length(u)+.5,axes=F)
  }
}



cor_fun=function(x,y){
  return(round(cor(x,y,method="pearson"),digits=2))
}

get_m_one_sample=function(samp){
  genes=intersect(rownames(lm$dataset$umitab),rownames(lm$dataset$noise_model))
  s=rowSums(lm$dataset$umitab[genes,lm$dataset$cell_to_sample==samp])
  s_noise=lm$dataset$beta_noise[samp]*lm$dataset$noise_model[genes,samp]*sum(lm$dataset$cell_to_sample==samp)
  s_net=pmax(s-s_noise,0)
  return(s_net/sum(s_net))
}


######################################################################################################
c('IGLV3-10', 'IGHG1', 'IGLJ3', 'WT1-AS', 'DNAAF1', 'KIF3C', 'INHBA', 'CLEC5A', 'MARCO', 'IL6', 'CXCL10', 'AQP9', 'IL1RN', 'GPR84', 'LILRA1', 'SOD2', 'SLAMF9', 'G0S2', 'CD300E', 'IL1A', 'CCL3', 'CXCL2', 'SERPINB2', 'CXCL8', 'S100A8', 'CXCL3', 'S100A9', 'DRAM1', 'CCL3L3', 'NINJ1', 'KMO', 'MIR3945HG', 'FPR2', 'IRAK2', 'IL1B', 'PILRA', 'CCL19', 'CCL22', 'FSCN1', 'CCL17', 'LAMP3', 'SLC1A2', 'HMSD', 'TVP23A', 'NRP2', 'NCCRP1', 'CD274', 'EBI3', 'CD80', 'CCR7', 'GPR157', 'PTGIR', 'TRAFD1', 'TRIP10', 'MSLN', 'B3GALT5', 'HOXA5', 'SPINK2', 'NEO1', 'LINC00299', 'IL22', 'PRMT9', 'KRT86', 'NRIP1', 'SLC4A10', 'NUDT7', 'USP53', 'SEPP1', 'SLC40A1', 'VSIG4', 'CD209', 'IL2', 'ODF2L')

load_bulk_datasets=function(load_text=F){
  clin=list()
  scRNAseq_genes=rownames(ileum_ldm$dataset$umitab)
  l=list()
  message("Loading CERTIFI")
  l$certifi=read_certifi(genes=scRNAseq_genes)
  message("Loading UNITI-1")
  l$uniti1=read_uniti1(genes=scRNAseq_genes)
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
#  tissue_conversion_table=cbind(c('Ascending colon','Rectum',"RECTUM",'Terminal Ileum',"T. ILEUM","ILEUM",'Descending colon','Sigmoid colon','Transverse colon','DELETED','Not Collected','Blood','ileum','colon','sigmoid','CDc','CDi',"iCD",'uc','UC','ctrl','IL','SP. FLEX'),c('Ascending Colon','Rectum','Rectum','Ileum','Ileum','Ileum','Descending Colon','Sigmoid Colon','Transverse Colon','Exclude','Exclude','Blood','Ileum','Colon','Sigmoid Colon','Colon',"Ileum","Ileum","Colon","Colon","?","?",'SP. FLEX'))
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


  print(nrow(design))

  ####################################
  ##
  ## UNITI-1
  ##
  #################################
#  To increase sample size for maintenance related analysis, we would normally combined groups of the same dose from both randomized and non-randomized populations.  For instance:
#  placebo= “NR_CRP_P” + “RD_CRU_P”
#  U12W=” NR_NCRP_U130W0_CRU12W” + “RD_CRU_U12W”
#  U8W=” NR_NCRU_U90W0_CRU8W” + “RD_CRU_U8W”
  
  
  uniti1_raw=matrix(NA,length(genes),ncol(l$uniti1$raw),dimnames = list(genes,colnames(l$uniti1$raw)))
  uniti1_raw[rownames(l$uniti1$raw),]=as.matrix(l$uniti1$raw)
  maintainance_trt_conversion=c("Pbo","Pbo","Ust q8wk","Ust q8wk","Ust q12wk","Ust q12wk")
  names(maintainance_trt_conversion)=c("NR_CRP_P","RD_CRU_P","NR_NCRU_U90W0_CRU8W","RD_CRU_U8W","NR_NCRP_U130W0_CRU12W","RD_CRU_U12W")
 
  
  uniti_subject=sapply(strsplit(as.character(l$uniti1$design$USUBJID),"-"),function(x){x[2]})
  uniti_tissue=tissue_conversion[l$uniti1$design$ANALOC]
  uniti_induction_treatment=treatment_conversion[l$uniti1$design$TR01AG1]
  uniti_maintainance_treatment=maintainance_trt_conversion[l$uniti1$design$TRTGRPAM]
  uniti_treatment=uniti_induction_treatment
  uniti_BA=ba_converstion[l$uniti1$design$visit]
  
  clinical_response=(l$uniti1$design$"CDAI_I_WK8"-l$uniti1$design$"CDAI_I_WK0")<(-100)|l$uniti1$design$"CDAI_I_WK8"<150
  clinical_remission=l$uniti1$design$"CDAI_M_WK44"<(150)
  
  clin$uniti1=l$uniti1$design[,-1:-5]
  
  ses_pooled_wk0=rowMeans(apply(l$uniti1$design[,c('SES_CD_Ileum_WK0','SES_CD_LeftColon_WK0','SES_CD_Rectum_WK0','SES_CD_RightColon_WK0','SES_CD_TransverseColon_WK0')],2,as.numeric),na.rm=T)
  ses_pooled_wk44=rowMeans(apply(l$uniti1$design[,c('SES_CD_Ileum_WK44','SES_CD_LeftColon_WK44','SES_CD_Rectum_WK44','SES_CD_RightColon_WK44','SES_CD_TransverseColon_WK44')],2,as.numeric),na.rm=T)
  uniti_endoscopic_response=ses_pooled_wk44/ses_pooled_wk0<0.5
  
  uniti_endoscopic_remission=ses_pooled_wk44<4

  uniti_response=ifelse(clinical_response,"R","NR")
#  uniti_response=ifelse(clinical_remission,"R","NR")
  #  uniti_endoscopic_response=as.numeric(l$uniti1$design[,c('SES_CD_Ileum_WK44')])/as.numeric(l$uniti1$design[,c('SES_CD_Ileum_WK0')])<.5
  
  #&pmax(-1,l$uniti1$design$SES_CD_Ileum_WK8,na.rm=T)<=2,"R","NR")
  uniti1_design=data.frame(dataset="UNITI-1",subject=uniti_subject,geo_accesion=NA,sex=NA,diagnosis="CD",age=NA,tissue=uniti_tissue,status="Inflamed",deep_ulcer=NA,BA=uniti_BA,treatment=uniti_treatment,response=uniti_response)
  rownames(uniti1_design)=rownames(l$uniti1$design)
  raw_exprs=cbind(raw_exprs,uniti1_raw)
  design=rbind(design,uniti1_design)
  
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
  
  uniti_maintainance_treatment=maintainance_trt_conversion[l$uniti2$design$TRTGRPAM]
  uniti_treatment=uniti_induction_treatment
  uniti_BA=ba_converstion[l$uniti2$design$VISIT]
  clinical_response=(l$uniti2$design$"CDAI_I_WK8"-l$uniti2$design$"CDAI_I_WK0")<(-100)|l$uniti2$design$"CDAI_I_WK8"<150
  clinical_remission=l$uniti2$design$"CDAI_M_WK44"<(150)
  
  ses_pooled_wk0=rowMeans(apply(l$uniti2$design[,c('SES_CD_Ileum_WK0','SES_CD_LeftColon_WK0','SES_CD_Rectum_WK0','SES_CD_RightColon_WK0','SES_CD_TransverseColon_WK0')],2,as.numeric),na.rm=T)
  ses_pooled_wk44=rowMeans(apply(l$uniti2$design[,c('SES_CD_Ileum_WK44','SES_CD_LeftColon_WK44','SES_CD_Rectum_WK44','SES_CD_RightColon_WK44','SES_CD_TransverseColon_WK44')],2,as.numeric),na.rm=T)
  uniti_endoscopic_response=ses_pooled_wk44/ses_pooled_wk0<0.5
  uniti_endoscopic_remission=ses_pooled_wk44<3
 
  uniti_response=ifelse(clinical_response,"R","NR")
  clin$uniti2=l$uniti2$design[,-1:-5]
  uniti2_design=data.frame(dataset="UNITI-2",subject=uniti_subject,geo_accesion=NA,sex=NA,diagnosis="CD",age=NA,tissue=uniti_tissue,status="Inflamed",deep_ulcer=NA,BA=uniti_BA,treatment=uniti_treatment,response=uniti_response)
 
  #  uniti2_design=data.frame(dataset="UNITI-2",subject=uniti_subject,geo_accesion=NA,sex=NA,diagnosis="CD",age=NA,tissue=uniti_tissue,status="Inflamed",deep_ulcer=NA,BA=uniti_BA,treatment=uniti_treatment,response=uniti_response)
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
  clin$risk=cbind(l$risk$clin,l$risk$pcdai)
  
  
  z_exprs=get_all_zscores(raw_exprs,design,sample_set_list)

  return(list(raw=raw_exprs,z=z_exprs,design=design,clin=clin))
}



DE_pat1_vs_pat2=function(){
  source('~/GoogleDrive/work/scClustering/DE.r')
 
  de_res=DE_between_two_sets(ileum_ldm,mask1,mask2,nchunks=1000,n_per_chunk=1000,reg=1e-7)
  write.csv(de_res[order(de_res$log2_FC),],file="output/tables/DE_inflamed_pat1_vs_pat2.csv")
  
  uninf_pat1=inflamed_samples[sample_to_patient[uninflamed_samples]%in%pat1]
  uninf_pat2=inflamed_samples[sample_to_patient[uninflamed_samples]%in%pat2]
  mask1=colnames(ileum_ldm$dataset$umitab)[ileum_ldm$dataset$cell_to_sample%in%uninf_pat1]
  mask2=colnames(ileum_ldm$dataset$umitab)[ileum_ldm$dataset$cell_to_sample%in%uninf_pat2]
  de_res=DE_between_two_sets(ileum_ldm,mask1,mask2,nchunks=1000,n_per_chunk=1000,reg=1e-7)
  write.csv(de_res[order(de_res$log2_FC),],file="output/tables/DE_uninflamed_pat1_vs_pat2.csv")
  
  inf_pat1=inflamed_samples[sample_to_patient[inflamed_samples]%in%pat1]
  inf_pat2=inflamed_samples[sample_to_patient[inflamed_samples]%in%pat2]
  for (cluster_set in setdiff(names(cluster_sets),"")){
    message(cluster_set)
    clusts=unlist(cluster_sets[[cluster_set]])
    mask1=colnames(ileum_ldm$dataset$umitab)[ileum_ldm$dataset$cell_to_sample%in%inf_pat1&ileum_ldm$dataset$cell_to_cluster%in%clusts]
    mask2=colnames(ileum_ldm$dataset$umitab)[ileum_ldm$dataset$cell_to_sample%in%inf_pat2&ileum_ldm$dataset$cell_to_cluster%in%clusts]
    de_res=DE_between_two_sets(ileum_ldm,mask1,mask2,nchunks=100,n_per_chunk=1000,reg=1e-7,noise_correction = T)
    write.csv(de_res[order(de_res$log2_FC),],file=paste("output/tables/DE_inflamed_pat1_vs_pat2_",cluster_set,".csv",sep=""))
  }
  
}


DE_pat1_vs_pat2_with_EP=function(){
  nchunks=1000
  n_per_chunk=100
  reg=1e-7
  nmin_umi_thresh=1
  nmin_cells_with_min_umi=20
  ldm=ileum_ldm
  
  # without patients 6 and 11
  pat1<<-paste("rp",c(7,8,5,11,12))
  pat2<<-paste("rp",c(6,10,13,14,15))
  
  inf_pat1=inflamed_samples[sample_to_patient[inflamed_samples]%in%pat1]
  inf_pat2=inflamed_samples[sample_to_patient[inflamed_samples]%in%pat2]
  mask1=colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample%in%inf_pat1]
  mask2=colnames(ldm$dataset$umitab)[ldm$dataset$cell_to_sample%in%inf_pat2]
  u=ldm$dataset$umitab[,c(mask1,mask2)]
  
  for (si in inf_pat1){
    u=cBind(u,ldm$dataset$gated_out_umitabs$EP[[si]])
    mask1=c(mask1,colnames(ldm$dataset$gated_out_umitabs$EP[[si]]))
  }
  for (si in inf_pat2){
    u=cBind(u,ldm$dataset$gated_out_umitabs$EP[[si]])
    mask2=c(mask2,colnames(ldm$dataset$gated_out_umitabs$EP[[si]]))
  }
  
  gene_mask=rowSums(u>nmin_umi_thresh)>nmin_cells_with_min_umi
  u=u[gene_mask,c(mask1,mask2)]
  
  obs_s=pmax(rowSums(u),0)
  obs_s1=pmax(rowSums(u[,mask1]),0)
  obs_s2=obs_s-obs_s1
  obs_m1=obs_s1/sum(obs_s1)
  obs_m2=obs_s2/sum(obs_s2)
  obs_log2_fc=log2((reg+obs_m2)/(reg+obs_m1))
  ncounts_bigger=rep(0,nrow(u))

  nA=min(length(mask1),length(mask2))
  mat=matrix(c(rep(T,nA),rep(F,ncol(u)-nA)),n_per_chunk,ncol(u),byrow =T)

  ntot=ncol(u)
  for (i in 1:nchunks){
    message(i)
  
    print(system.time({
      sA=u%*%apply(mat,1,sample,ntot)
    }))
    sB=obs_s-sA
    mA=t(t(sA)/colSums(sA))
    mB=t(t(sB)/colSums(sB))
    log2_fc=log2((reg+mB)/(reg+mA))
    ncounts_bigger=ncounts_bigger+rowSums(abs(log2_fc)>=abs(obs_log2_fc))
  
  #    
  #    save(file="de_current_iter.rd",i)
  }
  p.value=ncounts_bigger/(i*n_per_chunk)
  adj.p.value=p.adjust(p.value,method = "BH")
  de_res=data.frame(counts1=obs_s1,counts2=obs_s2,freq1=obs_m1,freq2=obs_m2,log2_FC=obs_log2_fc,p.value=p.value,adj.p.value=adj.p.value)
  write.csv(de_res[order(de_res$log2_FC),],file="output/tables/DE_inflamed_pat1_vs_pat2.csv")
}
  

risk_clinical_data=function(bulk){
  columns_res=grep("RES",colnames(bulk$clin$risk),val=T)
  profs=gsub("RES","",columns_res)
  for (i in 1:length(profs)){
    open_plot(path="output/bulk/",fn = paste("risk_",profs[i],"_vs_scores",sep=""),plot_type = "pdf",width = 8,height = 4)
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
  

  #marker_to_cluster_set1=marker_to_cluster_set[marker_to_cluster_set!="EP"]#&!marker_to_cluster_set%in%grep("Stromal",marker_to_cluster_set)] 
  #marker_to_cluster_set2=marker_to_cluster_set[marker_to_cluster_set!="EP"]#&!marker_to_cluster_set%in%grep("Stromal",marker_to_cluster_set)] 
  short_markers_list1=sapply(split(names(marker_to_cluster_set1),marker_to_cluster_set1),function(x){x[order(abs(de_res[x,"log2_FC"]),decreasing = T)][1:pmin(length(x),max_per_clusterset)]},simplify = F)
  short_markers_list2=sapply(split(names(marker_to_cluster_set2),marker_to_cluster_set2),function(x){x[order(abs(de_res[x,"log2_FC"]),decreasing = T)][1:pmin(length(x),max_per_clusterset)]},simplify = F)
  short_markers2=sapply(short_markers_list2,function(x){x[de_res[x,"log2_FC"]<0]})
  short_markers1=sapply(short_markers_list1,function(x){x[de_res[x,"log2_FC"]>0]})
  short_markers1=short_markers1[sapply(short_markers1,length)>0]
  short_markers2=short_markers2[sapply(short_markers2,length)>0]

  return(list(markers1=short_markers1[included_cluster_sets1],markers2=short_markers2[included_cluster_sets2]))
}



#######################################################################################################################################


main_bulk=function(load_data=F){

  
  lm=ileum_ldm
  score_thresh<<-list()  
  responder_col<<-brew_cols[5]
  non_responder_col<<-brew_cols[4]
  if (load_data){
    lbd_res<-load_bulk_datasets(load_text)
    bulk<<-lbd_res
    #mods<<-read_modules()
    write.csv(file="output/tables/full_bulk_design.csv",bulk$design)
  }
  
  #message("pattern1:",paste(pat1,collaspse=","),"pattern2:",paste(pat2,collaspse=","))
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
  
  gene_tab=data.frame(ifelse(de_res2$log2_FC>thresh,2,1),colnames(models)[apply(models[rownames(de_res2),],1,which.max)],de_res2[,c("freq_bg","freq_fg","log2_FC")])
  rownames(gene_tab)=rownames(de_res2)
  
  gene_tab=gene_tab[de_res2$adj.p.value<1e-4&abs(de_res2$log2_FC)>thresh,]
  
  scores_EP_genes<<-select_genes(max_per_clusterset =20,thresh_cluster=1,thresh_pattern1=1,thresh_pattern2 =.5 ,thresh_freq=1e-6,cluster_to_subtype=cluster_to_subtype,included_cluster_sets1=c("EP"),included_cluster_sets2=c("EP"),restrict_to_included_cluster_sets=F,trace=F)
  
  scores_genes<<-select_genes(max_per_clusterset =20,thresh_cluster=3,thresh_pattern1=1,thresh_pattern2 =1,thresh_freq=1e-6,cluster_to_subtype=cluster_to_subtype1,fdr_thresh = 0.001,
                                included_cluster_sets1=c("Inf. Macrophages","Activated DC","pDC","IgG plasma cells","Naive/CM T cells","Tregs","ACKR1+ endothelial cells","Activated fibroblasts"),
                                included_cluster_sets2=c("Resident macrophages","moDC","IgA plasma cells","ILC3","Cytokines low Trm","Type 1 cytokines Trm","Type 3 cytokines Trm","Enteric neurons","CD36+ endothelial cells", "Fibroblasts"),
                              trace=T,restrict_to_included_cluster_sets = T)
  
  
  bulk$scores<<-get_scores(scores_genes$markers1,scores_genes$markers2,bulk$z) 

  scRNA_scores<<-get_scores(scores_genes$markers1,scores_genes$markers2,mz)
  ep_score<<-get_scores(scores_genes$markers1,scores_EP_genes,bulk$z) 
  
  project_all_scores(sample_set_list,design = bulk$design)
 
  maski=bulk$design$dataset%in%c("UNITI-1","UNITI-2")
  mask_early=maski&bulk$design$BA=="beforeT"&maski&bulk$design$tissue=="Ileum"
  mapping_mask= which(mask_early)[match(paste(bulk$design[maski,"dataset"],bulk$design[maski,"subject"]),paste(bulk$design[mask_early,"dataset"],bulk$design[mask_early,"subject"]))]
  bulk$design[maski,"prediction"]<<-ifelse(bulk$scores[mapping_mask,"score1"]>prediction_thresh,"NR","R")
  
  
  markers<<-scores_genes

  message("pattern1")
  print(markers$markers1[sapply(markers$markers1,length)>0])
  message("pattern2")
  print(markers$markers2[sapply(markers$markers2,length)>0])
  write.csv(file="output/tables/score_genes.csv",rbind(cbind(geneSymbol=unlist(markers$markers1[sapply(markers$markers1,length)>0]),group="with module"),cbind(geneSymbol=unlist(markers$markers2[sapply(markers$markers2,length)>0]),group="No module")),quote=F,row.names=T)
  open_plot("output/main_figures/",fn="figure_4a_volcano_pat1_vs_pat2",plot_type="pdf",width = 5,height = 5)
  par(mar=c(4,4,4,4))
  plot(de_res$log2_FC,-log2(de_res$adj.p.value),xlim=c(-1,1),xlab="Log2(freq(pat2)/freq(pat1))",ylab="FDR adjusted p.value (BH)",cex=.7)
  close_plot()
  
  open_plot("output/main_figures/",fn="figure_4b_scatter_scores_inflamed",plot_type="pdf",width = 5,height = 5)
  plot_figure_4b()
  close_plot()
 
  open_plot("output/main_figures/",fn="figure_4b_scatter_scores_uninflamed",plot_type="pdf",width = 5,height = 5)
  plot_figure_4b("Uninflamed")
  close_plot()
  

  plot_dataset("CERTIFI",genes_to_show = markers)
  plot_dataset_with_stratification(tissues=c("Ileum"),dataset_name="UNITI-1",genes_to_show = markers,stratify_by = "ba")
  plot_dataset_with_two_stratifications(tissues=c("Ileum"),dataset_name="UNITI-2",genes_to_show = markers)
  plot_dataset_with_two_stratifications(tissues=c("Ileum"),dataset_name="UNITI-1",genes_to_show = markers)
  plot_dataset("RISK",plot_ulceration = F,plot_bcur = F,genes_to_show = markers)
  

  
  rownames(mm)=sample_to_patient[rownames(mm)]
  ord=match(rownames(mm),c(pat1,pat2))
  gene_groups_heatmaps(t(mm),markers$markers1,markers$markers2,by=ord,paste("heatmap_scRNA_patients",sep=""),sample_labels=T,zlim=c(-2,2),log_trans = T,reg=1e-6,height = 1.3,show_gene_names=F,show_celltype_bar=F)

  risk_clinical_data(bulk)

  score_thresh$RISK<<-plot_risk_figures()
  

  
  write.csv(file="output/tables/scores.csv",bulk$scores,row.names = T)
}


plot_figure_4b=function(status="Inflamed"){
  before_mask=bulk$design$BA=="beforeT"
  before_mask[is.na(before_mask)]=T
  mask=bulk$design$tissue=="Ileum"&bulk$design$status==status&bulk$design$diagnosis=="CD"&before_mask

  datasets=c('RISK',"CERTIFI",'UNITI-1', 'UNITI-2')
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






compare <- function(x, y) {
  n <- length(x); m <- length(y)
  w <- c(x, y)
  o <- order(w)
  z <- cumsum(ifelse(o <= n, m, -n))
  i <- which.max(abs(z))
  w[o[i]]
}


plot_risk_figures=function(){
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
  open_plot(path="output/bulk/",fn ="risk_roc",plot_type = "pdf",width = 6,height = 6)
  par(mar=c(5,5,1,1))
  auc=simple_auc(simple_roc(response1=="R",-1*score1))
  auc_resamp=replicate(1e4,simple_auc(simple_roc(sample(response1=="R",length(response1),replace = F),-1*score1)))
  auc_pvalue=mean(auc_resamp>auc)
  plot(rbind(c(0,0),simple_roc(response1=="R",-1*score1)[,c(2,1)]),type='l',panel.first=grid(lty=3,col="gray"),main=paste(round(auc,digits=1),"p=",round(auc_pvalue,digits=3)))
 
  abline(0,1,lty=3,col=1)
  close_plot()

  open_plot(path="output/bulk/",fn ="risk_response",plot_type = "pdf",width = 4,height = 4)
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
  max(abs(z)) # Dmax
  max.at <- sort(w)[which(abs(z) == max(abs(z)))]
  message("D is max at",max.at)
  open_plot(path="output/bulk/",fn ="risk_response_threshold",plot_type = "pdf",width = 4,height = 4)
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
  
  open_plot(path="output/bulk/",fn ="risk_response_piecharts",plot_type = "pdf",width = 6,height = 6)
  par(mar=c(1,1,4,1))
  layout(matrix(1:2,2,1))
  pie(tab[1,],c("Non responders","Responders"),col=c(non_responder_col,responder_col),main=paste("Signature score <",round(max.at,digits=2)))
  pie(tab[2,],c("Non responders","Responders"),col=c(non_responder_col,responder_col),main=paste("Signature score >",round(max.at,digits=2)))
  close_plot()
  
  l_score=split(-1*by,bulk$design$deep_ulcer[mask_risk])[c("deep ulcer","no deep ulcer")]
  open_plot(path="output/bulk/",fn ="risk_ulceration",plot_type = "pdf",width = 8,height = 4)
  layout(matrix(1:2,1,2))
  plot(1,col=0,xlim=c(-1,1),ylim=c(0,1),panel.first=grid(lty=1),ylab="Cumulative Fraction",xlab="Projection score")
  matplot(sapply(l_score,quantile,0:100/100),0:100/100,type='l',lty=1,lwd=2,col=c(non_responder_col,responder_col),xlim=c(-1,1),cex=.7,add=T)
  ks_res=ks.test(l_score$`deep ulcer`,l_score$`no deep ulcer`,alternative = "less")
  text(-1,0.95,labels = paste("KS D=",round(ks_res$statistic,digits=2),"; p=",round(ks_res$p.value,digits = 4)),pos=4,cex=1)
  legend("bottomright",legend=c("Deep Ulcer","No Deep Ulcer"),col=c(non_responder_col,responder_col),lty=1,lwd=2,cex=.8)
  data=t(table(patient_groups,bulk$design$deep_ulcer[mask_risk]))[c("no deep ulcer","deep ulcer"),c("low","high")]
  barplot(data,beside=T,col=c(responder_col,non_responder_col),names.arg = c("Low score", "High score"),yaxp=c(0, 10*ceiling(max(data)/10), ceiling(max(data)/10)))
  close_plot()
  
  
  l_score=split(-1*by,bulk$design$B_Cur[mask_risk]=="B1")

  names(l_score)=c("Other","B1")
  open_plot(path="output/bulk/",fn ="B1_",plot_type = "pdf",width = 8,height = 4)
  layout(matrix(1:2,1,2))
  plot(1,col=0,xlim=c(-1,1),ylim=c(0,1),panel.first=grid(lty=1),ylab="Cumulative Fraction",xlab="Projection score")
  matplot(sapply(l_score,quantile,0:100/100),0:100/100,type='l',lty=1,lwd=2,col=c(non_responder_col,responder_col),xlim=c(-1,1),cex=.7,add=T)
  ks_res=ks.test(l_score$Other,l_score$B1,alternative = "less")
  text(-1,0.95,labels = paste("KS D=",round(ks_res$statistic,digits=2),"; p=",round(ks_res$p.value,digits = 4)),pos=4,cex=1)
  legend("bottomright",legend=c("other","B1"),col=c(non_responder_col,responder_col),lty=1,lwd=2,cex=.8)
  data=t(table(patient_groups,ifelse(bulk$design$B_Cur[mask_risk]=="B1","B1","Other")))[c("B1","Other"),c("low","high")]
  barplot(data,beside=T,col=c(responder_col,non_responder_col),names.arg = c("Low score", "High score"),yaxp=c(0, 10*ceiling(max(data)/10), ceiling(max(data)/10)))
  close_plot()
  
  columns_res=grep("RES",colnames(bulk$clin$risk),val=T)
  profs=gsub("RES","",columns_res)
  columns_res=c(columns_res,"ENCRF.18.PCDAI.SCORE.PCDAI")
  profs=c(profs,"pcdai")
  browser()
  for (i in 1:length(profs)){
    open_plot(path="output/bulk/",fn = paste("risk_",profs[i],"_vs_scores",sep=""),plot_type = "pdf",width = 6.5,height = 4)
    layout(matrix(1:2,1,2))
    plot(-1*by,bulk$clin$risk[rownames(bulk$design)[mask_risk],columns_res[i]],xlab="Score 1",ylab=profs[i],panel.first = grid(lty=1))
    cor_res=cor.test(-1*by,bulk$clin$risk[rownames(bulk$design)[mask_risk],columns_res[i]],method = "spe")
    mtext(side=3,paste("r=",round(cor_res$estimate,digits=2)," ; p=",round(cor_res$p.value,digits=4),sep=""),line = 1)
    boxplot(split(bulk$clin$risk[rownames(bulk$design)[mask_risk],columns_res[i]],patient_groups))
    dev.off()
  }

  ##genetics
  apply(bulk$clin$risk[,c("IL23R_imm_1_67475114_rs11465804","ATG16L1_imm_2_233845149_rs3828309","PRDM1_imm_6_106541962_rs7746082","JAK2_imm_9_4971602_rs10758669","NOD2_imm_16_49303427_rs2066844","NOD2_imm_16_49314041_rs2066845","NOD2_rs5743293_rs2066847")],2,function(x){tab=table(x,1+as.numeric(bulk$scores[rownames(bulk$clin$risk),"score1"]<.25));tab=t(t(tab)/colSums(tab));return(tab)})
  apply(bulk$clin$risk[,c("IL23R_imm_1_67475114_rs11465804","ATG16L1_imm_2_233845149_rs3828309","PRDM1_imm_6_106541962_rs7746082","JAK2_imm_9_4971602_rs10758669","NOD2_imm_16_49303427_rs2066844","NOD2_imm_16_49314041_rs2066845","NOD2_rs5743293_rs2066847")],2,table, 1+as.numeric(bulk$scores[rownames(bulk$clin$risk),"score1"]<.25))
 
  return(max.at)
}
plot_cor=function(dataset="Arijs-frma"){
  gs=c(intersect(rownames(bulk$raw),unlist(markers$markers1)),intersect(rownames(bulk$raw),unlist(markers$markers2)))
  x=bulk$z[gs,bulk$design$dataset==dataset]
  x=x[rowSums(is.na(x))<ncol(x),]
  x=x[,colSums(is.na(x))<nrow(x)]
  cormat=cor(t(x),use="comp")
  image(cormat,col=colorRampPalette(c("blue","white","red"))(100),breaks=c(-1,seq(-.7,.7,l=99),1))
}


plot_dataset_with_two_stratifications=function(dataset_name="Arijs",tissues,pool_tissues=F,ngenes=100,genes_to_show=c(),plot_cor=T){
  plot_tissue=function(mask_tissue,statuses,file_prefix){
    for (status in statuses){
      if (status=="Uninflamed"){
        next()
      }
      suffix=paste(file_prefix,status,sep="_")
      print(suffix)
      
      maski=mask_tissue&(statusv==status)
      
     

      mask_stratify_list1=list()
      mask_stratify_list2=list()
  
        bav=bulk$design$BA
        bav[is.na(bav)]=""
  
        mask_stratify_list1$beforeT=maski&bav=="beforeT"
        mask_stratify_list1$afterT=maski&bav=="afterT"
      
      
   
        treatmentv=bulk$design$treatment
        treatmentv[is.na(treatmentv)]=""
        mask_stratify_list2$Pbo=maski&treatmentv=="Pbo"
        mask_stratify_list2$Ust_130_mg=maski&treatmentv=="Ust 130 mg"
        mask_stratify_list2$Ust_6_mg_kg=maski&treatmentv=="Ust 6 mg/kg"
   
      nx=length(mask_stratify_list1)
      ny=length(mask_stratify_list2)
      open_plot(path="output/bulk/",fn=paste("scatter_",suffix,sep=""),plot_type = "pdf",width=4*ny,height=4*nx)
      layout(matrix(1:(nx*ny),nx,ny,byrow = T))
      
      xlim=range(bulk$scores[mask_tissue&(statusv==status),"score1"],na.rm=T)
      ylim=range(bulk$scores[mask_tissue&(statusv==status),"score2"],na.rm=T)
    
      for (i in 1:length(mask_stratify_list1)){
        for (j in 1:length(mask_stratify_list2)){
          
        mask_stratify=mask_stratify_list1[[i]]&mask_stratify_list2[[j]]
        mask_stratify=mask_stratify&!(is.na(bulk$scores[,"score1"])|is.na(bulk$scores[,"score2"]))
        if (sum(mask_stratify)<10){
          plot_cor=F
        }
        if (plot_cor){
          cort_res=cor.test(bulk$scores[mask_stratify,"score1"],bulk$scores[mask_stratify,"score2"],method="spe")
        }
        response=ifelse(bulk$design$response[mask_stratify]=="R",responder_col,non_responder_col)
      #  response=ifelse(bulk$design$prediction=="R",responder_col,non_responder_col)
        response[is.na(response)]=1
        if (sum(mask_stratify)>0){
          plot(bulk$scores[mask_stratify,"score1"],bulk$scores[mask_stratify,"score2"],panel.first=grid(),pch=20,xlab="score1",ylab="score2",cex.axis=.7,col=response,xlim=xlim,ylim=ylim,main=paste(names(mask_stratify_list1)[i],names(mask_stratify_list2)[j]))
          if (plot_cor){
            mtext(paste("rho =",format(cort_res$estimate,digits=2),"; p =",format(cort_res$p.value,digits=2)),3,line = -1,adj = 1,cex=.7)
          }
        }
        else {
          plot.new()
        }
        }
      }
      close_plot()
      for (i in 1:length(mask_stratify_list1)){
     
          mask_stratify=mask_stratify_list1[[i]]
          mask_stratify=mask_stratify&!(is.na(bulk$scores[,"score1"])|is.na(bulk$scores[,"score2"]))
          if (sum(mask_stratify)==0){
            next
          }
          suffix=paste(file_prefix,names(mask_stratify_list1)[i],status,sep="_")
         
          if (names(mask_stratify_list1)[i]=='beforeT'){
            annot_profiles=list()
            annot_profiles$response=ifelse(bulk$design$response[mask_stratify]=="R",responder_col,non_responder_col)
        #    annot_profiles$response=ifelse(bulk$design[mask_stratify,"prediction"]=="NR",non_responder_col,responder_col)
            
            write.csv(cbind(bulk$design[mask_stratify,1:2],bulk$design[mask_stratify,"prediction"]),file = paste("output/tables/prediction_",suffix,".csv",sep=""))
            gene_groups_heatmaps(bulk$z[,mask_stratify],genes_to_show$markers1,genes_to_show$markers2,by=bulk$scores[mask_stratify,"score1"],paste("heatmap_",suffix,"_combined_selected_genes",sep=""),annot_profiles=annot_profiles)
            open_plot(path="output/bulk/",fn=paste("scatter_with_prefictions_",suffix,sep=""),plot_type = "pdf",width=4,height=4)
            plot(bulk$scores[mask_stratify,"score1"],bulk$scores[mask_stratify,"score2"],col=annot_profiles$response,pch=20,xlab="score1",ylab="score2")
            close_plot()
          }
          
          for (j in 1:length(mask_stratify_list2)){
            mask_stratify=mask_stratify_list1[[i]]&mask_stratify_list2[[j]]
            mask_stratify=mask_stratify&!(is.na(bulk$scores[,"score1"])|is.na(bulk$scores[,"score2"]))
            if (sum(mask_stratify)==0){
              next
            }
            suffix=paste(file_prefix,names(mask_stratify_list1)[i],names(mask_stratify_list2)[j],status,sep="_")
            annot_profiles=list()
       #     annot_profiles$response=ifelse(bulk$design[mask_stratify,"prediction"]=="NR",non_responder_col,responder_col)
           annot_profiles$response=ifelse(bulk$design$response[mask_stratify]=="R",responder_col,non_responder_col)
            gene_groups_heatmaps(bulk$z[,mask_stratify],genes_to_show$markers1,genes_to_show$markers2,by=bulk$scores[mask_stratify,"score1"],paste("heatmap_",suffix,"_selected_genes",sep=""),annot_profiles=annot_profiles)
          }
        
      }
    }
  }
  
  
  mask=bulk$design$dataset==dataset_name
  statusv=as.character(bulk$design$status)
  statusv[is.na(statusv)]=""
  for (diagnosis in unique(bulk$design[mask,]$diagnosis)){
    mask_diagnosis=mask&(bulk$design$diagnosis==diagnosis)
    
    if (pool_tissues){
      mask_tissue=mask_diagnosis&(bulk$design$tissue%in%tissues)
      mask_tissue[is.na(mask_tissue)]=F
      statuses=setdiff(unique(statusv[mask_tissue]),"")
      plot_tissue(mask_tissue,statuses,file_prefix=paste(dataset_name,diagnosis,sep="_"))
    }
    else {
      for (tissue in tissues){
        mask_tissue=mask_diagnosis&(bulk$design$tissue==tissue)
        mask_tissue[is.na(mask_tissue)]=F
        statuses=setdiff(unique(statusv[mask_tissue]),"")
        file_prefix=paste(dataset_name,diagnosis,sep="_")
        plot_tissue(mask_tissue,statuses,file_prefix=paste(dataset_name,diagnosis,tissue,sep="_"))
      }
    }
  }
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


  

plot_roc=function(mask){
  get_xy=function(samp){
    tmp_tab=table(bulk$design[mask,]$response=="NR",bulk$scores[mask,"score1"]>=bulk$scores[samp,"score1"])
    tab=matrix(0,2,2,dimnames=list(c(T,F),c(T,F)))
    tab[rownames(tmp_tab),colnames(tmp_tab)]=tmp_tab
    
    x=tab[2,1]/(tab[2,1]+tab[2,2])
    y=tab[1,1]/(tab[1,1]+tab[1,2])
    return(c(x,y))
  }
  z=t(sapply(mask,get_xy))
  z=z[order(z[,1]),]
  plot(z,panel.first=abline(0,1),type='l',panel.first=abline(0,1))
}


