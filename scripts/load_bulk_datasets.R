
read_RISK=function(load_text=F,genes){
  path="~/GoogleDrive/work/shared/data/RISK data files/"
   fn="output/bulk/risk.rd"
  if (load_text){
    sample_tab=read.csv(paste(path,"Updated_sheets/Merged_RISKID_FPKM_SIden.csv",sep=""),stringsAsFactors = F)
    rownames(sample_tab)=gsub("-","",sample_tab$RNA.ID)
    risk_to_rna_id=gsub("-","",sample_tab$RNA.ID)
    names(risk_to_rna_id)=gsub("RS-","",sample_tab$RISK.DNA.ID)
    rna_to_risk_id=names(risk_to_rna_id)
    names(rna_to_risk_id)=risk_to_rna_id
    
    exprs=read.csv(paste(path,"FPKM/FPKMdataSummary_fpkm_genes.csv",sep=""),stringsAsFactors = F)
    converted_symbol=adjust_gene_names(paste(exprs[,1],collapse = ","),genes)
    agg=aggregate(as.matrix(exprs[,-1:-3]),by=list(as.character(converted_symbol)),FUN=mean)
    m=as.matrix(agg[,-1])
    rownames(m)=as.character(unlist(agg[,1]))
    
    clin=read.csv(paste(path,"Updated_sheets/KM_B2B3_FinalChart_DiagEdit_Aug2017.csv",sep=""),stringsAsFactors = F)
    rownames(clin)=gsub("-","",clin[,"Subject.Identifier"])
    #tnf_response=ifelse(!is.na(ryan$Progressor),ifelse(ryan$Progressor==1,"NR",""),ifelse(ryan$B_Dx=="B1"&ryan$B_curr=="B1","R",""))
    #names(tnf_response)=gsub("-","",ryan$SUBJECTIDENTIFIER)
    #  sample_tab$tnf_response=as.vector(tnf_response[names(risk_to_rna_id)])
    sample_tab$B_Dx=clin[names(risk_to_rna_id),]$B_Dx
    sample_tab$B_Cur=clin[names(risk_to_rna_id),]$B_curr
    
    pcdai=read.csv(paste(path,"PCDAI_BLANDFU_0_36_2017-11-19_20-44-05.csv",sep=""))
    pcdai=pcdai[!is.na(risk_to_rna_id[gsub("-","",pcdai[,"SUBJECTIDENTIFIER"])]),]
    rownames(pcdai)=risk_to_rna_id[gsub("-","",pcdai[,"SUBJECTIDENTIFIER"])]
    
    clin2=read.csv(paste(path,"RISK_BASELINE_DSHIFTED_2017-11-19_20-44-56.csv",sep=""),stringsAsFactors = F)
    # rownames(clin2)=risk_to_rna_id[gsub("-","",clin2$SUBJECTIDENTIFIER)]
    rna_ids2=risk_to_rna_id[gsub("-","",clin2$SUBJECTIDENTIFIER)]
    clin2=clin2[!is.na(rna_ids2),]
    rownames(clin2)=rna_ids2[!is.na(rna_ids2)]
    mask=intersect(intersect(rownames(sample_tab),colnames(m)),rownames(clin2))
    pcdai=pcdai[mask,-1]
    #names(tissue)=risk_to_rna_id
    #read.csv("~/GoogleDrive/work/shared/data/RISK data files/Updated_sheets/RNASeq_linking_Ileal.csv",stringsAsFactors = T)[,1]
    #tissue[gsub("-","",read.csv("~/GoogleDrive/work/shared/data/RISK data files/Updated_sheets/RNASeq_linking_Ileal.csv",stringsAsFactors = F)[,1])]="Ileal"
    #read.csv("~/GoogleDrive/work/shared/data/RISK data files/Updated_sheets/RNASeq_linking_Rectal.csv.csv",stringsAsFactors = T)[,1]
    #tissue[gsub("-","",read.csv("~/GoogleDrive/work/shared/data/RISK data files/Updated_sheets/RNASeq_linking_Ileal.csv",stringsAsFactors = F)[,1])]="Rectal"
    #sample_tab$tissue=tissue
    
    jeromes_tnf_response=read.delim(paste(path,"Jerome/TNF Response.txt",sep=""),stringsAsFactors = F)
    jeromes_tnf_response_subjects=risk_to_rna_id[gsub("-","",jeromes_tnf_response$SUBJECTIDENTIFIER)]
    jeromes_tnf_response=jeromes_tnf_response[!is.na(jeromes_tnf_response_subjects),]
    rownames(jeromes_tnf_response)=jeromes_tnf_response_subjects[!is.na(jeromes_tnf_response_subjects)]
    jeromes_tnf_treatment=read.delim(paste(path,"Jerome/Patients with TNF within 1st year.txt",sep=""),stringsAsFactors = F)
    jeromes_tnf_treatment_subjects=risk_to_rna_id[gsub("-","",jeromes_tnf_treatment$FUCD.SUBJECTIDENTIFIER)]
    jeromes_tnf_treatment=jeromes_tnf_treatment[!is.na(jeromes_tnf_treatment_subjects),]
    rownames(jeromes_tnf_treatment)=jeromes_tnf_treatment_subjects[!is.na(jeromes_tnf_treatment_subjects)]
    durable_remission=jeromes_tnf_response$Durable.Remission..no.active.disease.btw.18.24.months.
    tnf_treatment=jeromes_tnf_treatment[rownames(jeromes_tnf_response),"TNF.1st.year"]
    tnf_treatment[is.na(tnf_treatment)]="no"
    response=ifelse(ifelse(durable_remission=="yes",1,ifelse(durable_remission=="no",0,NA))*(ifelse(tnf_treatment=="yes",1,NA))==1,"R","NR")
    names(response)=rownames(jeromes_tnf_response)
    sample_tab[intersect(rownames(jeromes_tnf_response),rownames(sample_tab)),"response"]=response[intersect(rownames(jeromes_tnf_response),rownames(sample_tab))]
    save(file=fn,sample_tab,m,clin2,pcdai,mask)
  }
  else{
    env=new.env()
    load(file=fn,envir=env)
    m=env$m
    sample_tab=env$sample_tab
    clin2=env$clin2
    pcdai=env$pcdai
    mask=env$mask
  }
  return(list(raw=m[,mask],design=sample_tab[mask,],clin=clin2[mask,],pcdai=pcdai))
}


read_habermann=function(genes){
  patient_tab=read.csv("~/GoogleDrive/work/shared/data/public/human/habermann_et_al_JCI_2015/Habermann_iCD_patients.csv",stringsAsFactors = F,row.names = 1)
  gsm_tab=read.csv("~/GoogleDrive/work/shared/data/public/human/habermann_et_al_JCI_2015/Habermann_GSM_list.csv",stringsAsFactors = F,header=F,row.names = 1)
  exprs=read.delim("~/GoogleDrive/work/shared/data/public/human/habermann_et_al_JCI_2015/GSE57945_all_samples_RPKM.txt",stringsAsFactors = F)
  patient_mask=sapply(strsplit(gsm_tab[rownames(patient_tab),],"\\("),function(x){substr(x[2],1,nchar(x[2])-1)})
  names(patient_mask)=row.names(patient_tab)
  patient_mask=patient_mask[patient_mask%in%colnames(exprs)]
  design=patient_tab[names(patient_mask),]
  m=exprs[,patient_mask]
  
  converted_symbol=adjust_gene_names(exprs[,"Gene.Symbol"],ref_genes = genes)
  
  m2=aggregate(m,by=list(converted_symbol),mean)
  rownames(m2)=m2[,1]
  m2=m2[,-1]
  m=m2
  #  symbol_hist=table(converted_symbol)
  #  symbol_unique=names(symbol_hist)[as.vector(unlist(symbol_hist))==1]
  #  m1=m[converted_symbol%in%symbol_unique,]
  
  #  m2=aggregate(m[!(converted_symbol%in%symbol_unique),],by=list(converted_symbol[!(converted_symbol%in%symbol_unique)]),sum)
  
  #  rownames(m2)=m2[,1]
  #  m=rbind(m1,m2[,-1])
  #  m=m[rowSums(m>0)>1,]
  return(list(raw=m,design=design))
}


read_msh=function(load_text=F,genes){
  path="~/GoogleDrive/work/shared/data/public/human/Peters_et_al_NG_2017/"
  out_fn="output/bulk/risk.rd"
  fn=paste(path,"GSE83687.rd",sep="")
  if (load_text){
    msh_raw=read.csv(paste(path,"GSE83687_all_raw.csv",sep=""),row.names = 1)
    exprs_tab=aggregate(msh_raw,by=list(adjust_gene_names(rownames(msh_raw),ref_genes = genes)),mean)
    rownames(exprs_tab)=exprs_tab[,1]
    exprs_tab=exprs_tab[,-1]
    annots=read.csv(paste(path,"GSE83687_RAW/annotations/design.csv",sep=""),row.names = 1)
    save(file=fn,exprs_tab,annots)
  }
  else{
    env=new.env()
    load(file=fn,envir=env)
    exprs_tab=env$exprs_tab
    annots=env$annots
  }
  return(list(raw=exprs_tab,design=annots))
}


###### Not used
#read_one_hgu133plus2=function(path,fn,genes){
#  exprs_tab=read.table(paste(path,fn,sep=""))
#  library(hgu133plus2.db)
#  probe2gene=select(hgu133plus2.db,keys=gsub("PM_","",rownames(exprs_tab)),columns="SYMBOL")
#  gene_convertor1=adjust_gene_names(probe2gene[,2],ref_genes = genes)
#  names(gene_convertor1)=probe2gene[,2]
#  probe2gene[,2]=gene_convertor1
#  mask=!(probe2gene[,2]%in%genes)
#  tmp_tab=aggregate(exprs_tab[mask,],list(probe2gene[mask,2]),mean)
#  exprs_tab=tmp_tab[,-1]
#  rownames(exprs_tab)=tmp_tab[,1]
#  colnames(exprs_tab)=gsub(".CEL","",colnames(exprs_tab))
#  return(exprs_tab)
#}
########### ^^^^^^^


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



read_certifi=function(load_text=F,genes){
  path="~/GoogleDrive/work/shared/data/public/human/Peters_et_al_NG_2017/GSE100833_RAW/"
  
   fn="output/bulk/certifi.rd"
  if (load_text){
    #  library(Biobase)
    #  library(GEOquery)
    #    gset <- getGEO("GSE100833", GSEMatrix =TRUE, getGPL=T)
    #    if (length(gset) > 1) idx <- grep("GPL13158", attr(gset, "names")) else idx <- 1
    #    gset <- gset[[idx]]
    #    exprs_tab1=read_one_hgu133plus2pm(path="microarray/output/",fn="expr_mas5_CERTIFI-biopsy.csv",genes)
    #    exprs_tab2=read_one_hgu133plus2pm(path="microarray/output/",fn="expr_mas5_CERTIFI-blood.csv",genes)
    exprs_tab1=read_one_microarray_dataset(path="input/microarray/output/",fn="expr_CERTIFI-biopsy.csv",genes)
    exprs_tab2=read_one_microarray_dataset(path="input/microarray/output/",fn="expr_CERTIFI-blood.csv",genes)
    
    exprs_tab=cbind(exprs_tab1,exprs_tab2)
    colnames(exprs_tab)=sapply(strsplit(colnames(exprs_tab),"_"),head,1)
    #    annots=pData(phenoData(gset))
    load(file="~/GoogleDrive/work/shared/data/public/human/Peters_et_al_NG_2017/GSE100833_RAW/annotations.rd")
    
    mask=!(annots$`tissue:ch1`%in%c("DELETED","Not Collected"))
    exprs_tab=exprs_tab[,rownames(annots)[mask]]
    annots=annots[mask,]
    save(file=fn,exprs_tab,annots)
  }
  else {
    env=new.env()
    load(file=fn,envir=env)
    exprs_tab=env$exprs_tab
    annots=env$annots
  }
  
  return(list(raw=exprs_tab,design=annots))
}



read_uniti1=function(load_text=F,genes){
  path="~/GoogleDrive/work/shared/data/Janssen/"
  
   fn="output/bulk/uniti-1.rd"
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
  path="~/GoogleDrive/work/shared/data/Janssen/"
  
   fn="output/bulk/uniti-2.rd"
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


read_arijs_plos1=function(load_text=F,genes){
  
  path="~/GoogleDrive/work/shared/data/public/human/Arijs_et_al_PLoS1_2009/"
   fn="output/bulk/arijs_plos1.rd"
  if (load_text){
    exprs_tab=read_one_microarray_dataset(path="input/microarray/output/",fn="expr_Arijs_Plos1.csv",genes,remove_first_column =F)
    des=read.csv(paste(path,"arijs_design.csv",sep=""),row.names=1)
    annots=des
    save(file=fn,exprs_tab,annots)
  }
  else{
    
    env=new.env()
    load(file=fn,envir=env)
    exprs_tab=env$exprs_tab
    annots=env$annots
  }
  return(list(raw=exprs_tab,design=annots))
}

read_arijs_gut_2009=function(load_text=F,genes){
  
  path="~/GoogleDrive/work/shared/data/public/human/Arijs_et_al_Gut_2009/"
   fn="output/bulk/arijs_gut_2009.rd"
  if (load_text){
    exprs_tabA=read_one_microarray_dataset(path="input/microarray/output/",fn="expr_Arijs_Gut_cohortA.csv",genes,remove_first_column =F)
    exprs_tabB=read_one_microarray_dataset(path="input/microarray/output/",fn="expr_Arijs_Gut_cohortB.csv",genes,remove_first_column =F)
    desA=read.csv(paste(path,"Arijs_et_al_CohortA_GSE14580_Gut_Metadata.csv",sep=""),row.names=1)
    desB=read.csv(paste(path,"Arijs_et_al_CohortB_GSE12251_Gut_Metadata.csv",sep=""),row.names=1)
    desA$diagnosis=desA$Subject
    desA$Subject=sapply(strsplit(as.character(desA$Title),"_"),head,1)
    desB[,"diagnosis"]="UC"
    colnames(desA)[7]="Wk4.6_or_8.response"
    colnames(desB)[7]="Wk4.6_or_8.response"
    #exprs_tab=cbind(exprs_tabA,exprs_tabB)
    #annots=rbind(cbind(dataset="Arijs_Gut_A",desA),cbind(dataset="Arijs_Gut_B",desB))
    exprs_tab=exprs_tabB
    annots=cbind(dataset="Arijs_Gut_B",desB)
    mask=intersect(rownames(annots),colnames(exprs_tab))
    exprs_tab=exprs_tab[,mask]
    annots=annots[mask,]
    save(file=fn,exprs_tab,annots)
  }
  else{
    env=new.env()
    load(file=fn,envir=env)
    exprs_tab=env$exprs_tab
    annots=env$annots
  }
  return(list(raw=exprs_tab,design=annots))
}




read_toedter_ajg_2011=function(load_text=F,genes){
  
  path="~/GoogleDrive/work/shared/data/public/human/Toedter_et_al_AJG_2011/"
   fn="output/bulk/toedter_ajg_2011.rd"
  if (load_text){
    exprs_tab=read_one_microarray_dataset(path="input/microarray/output/",fn="expr_Toedter.csv",genes,remove_first_column =F)
    annots=read.csv(paste(path,"/Toedter_et_al_Metadata.csv",sep=""),row.names=1)
    save(file=fn,exprs_tab,annots)
  }
  else{
    env=new.env()
    load(file=fn,envir=env)
    exprs_tab=env$exprs_tab
    annots=env$annots
  }
  return(list(raw=exprs_tab,design=annots))
}



cap <- function(x) {
  paste(toupper(substring(x, 1,1)), tolower(substring(x, 2)),sep="", collapse=" ")
}

#adjust_gene_names=function(in_genes,ref_genes){
#  in_genes2=sapply(strsplit(in_genes," /// "),head,1)
#  names(in_genes)=in_genes
#  mask1=toupper(in_genes2)%in%ref_genes
#  in_genes2[mask1]=toupper(in_genes2[mask1])
#  mask2=unlist(sapply(in_genes2,cap))%in%ref_genes
#  in_genes2[mask2]=unlist(sapply(in_genes2[mask2],cap))
#  mask3=gene_symbol_old2new[in_genes2]%in%ref_genes
#  in_genes2[mask3]=gene_symbol_old2new[in_genes2[mask3]]
#  mask4=gene_symbol_new2old[in_genes2]%in%ref_genes
#  in_genes2[mask4]=gene_symbol_new2old[in_genes2[mask4]]
#  
#  print(sum(mask1))
#  print(sum(mask2))
#  print(sum(mask3))
#  print(sum(mask4))
#  
#  return(in_genes2)
#  
#}
