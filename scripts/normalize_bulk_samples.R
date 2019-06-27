
library(xps) 




init_probe_gene_converters=function(genes){
  probe_to_geneSymbol<<-list()
  tmp=read.csv("~/GoogleDrive/work/affymetrix/HT_HG-U133_Plus_PM-na36-annot-csv/HT_HG-U133_Plus_PM.na36.annot.csv",skip=25)
  tmp_probe_to_geneSymbol<-as.character(tmp[,"Gene.Symbol"])
  names(tmp_probe_to_geneSymbol)<-tmp[,"Probe.Set.ID"]
  probe_to_geneSymbol[["HT_HG-U133_Plus_PM"]]<<-adjust_gene_names(tmp_probe_to_geneSymbol,ref_genes = genes)
  
  tmp=read.csv("~/GoogleDrive/work/affymetrix/HG-U133_Plus_2-na36-annot-csv/HG-U133_Plus_2.na36.annot.csv",skip=25)
  tmp_probe_to_geneSymbol<-as.character(tmp[,"Gene.Symbol"])
  names(tmp_probe_to_geneSymbol)<-tmp[,"Probe.Set.ID"]
  probe_to_geneSymbol[["HG-U133_Plus"]]<<-adjust_gene_names(tmp_probe_to_geneSymbol,ref_genes = genes)
}

norm_frma=function(cels_path,output_csv){
  browser()
  cel_files1=list.files(path="~/GoogleDrive/work/shared/data/Janssen/UNITI-1 Biopsy/cels/",full.names = T,pattern = ".CEL")
  cel_files2=list.files(path="~/GoogleDrive/work/shared/data/Janssen/UNITI-2/PBOwk8RNR subset_wk0wk8samples/",full.names = T,pattern = ".CEL")
  cel_files3=list.files(path="~/GoogleDrive/work/shared/data/Janssen/UNITI-2/USTwk8RNR subset_wk0wk8samples/",full.names = T,pattern = ".CEL")
  
  cel_files=list.files(path=cels_path,full.names = T,pattern = ".CEL")
  Data <- ReadAffy(widget=F,filenames = cel_files[1:20],cdfname="hthgu133pluspmcdf")
  vecs=makeVectorsFeatureSet(files = c(cel_files1[1:10],cel_files2[1:10],cel_files3[1:10]),batch.id = c(rep(1,10),rep(2,10),rep(3,10)),pkgname = "pd.ht.hg.u133.plus.pm")
  frma_output <- frma(Data,input.vecs = vecs)
#  eset <- expresso(Data, normalize.method=normalize.method,bg.correct = T,bgcorrect.method="rma",normalize=F,bgcorrect.param = list(na.rm=T),pmcorrect.param = list(na.rm=T))
  
}


norm_expresso=function(cels_path,output_csv,normalize.method="qspline"){
  Data <- ReadAffy(widget=F,filenames = list.files(path=cels_path,full.names = T,pattern = ".CEL"))
  eset <- expresso(Data, normalize.method=normalize.method,bg.correct = T,bgcorrect.method="rma",normalize=T,pmcorrect.method = "pmonly")
  browser()
  returm(exprs(eset))
}



norm_dataset=function(cels_path,output_txt){
  message("Normalizing ", cels_path)
 # norm_xps(cels_path,output_txt)
  norm_frma(cels_path,output_txt)
#  message("Normalizing ", cels_path)
#  Data <- ReadAffy(widget=F,filenames = list.files(path=cels_path,full.names = T,pattern = ".CEL"))
#  browser()
#  frma_output <- frma(Data)
#  write.csv(exprs(frma_output), file=output_csv)
  
  
  
}

import_xps_dataset=function(dataset_name,cels_paths,root_path="~/GoogleDrive/work/affymetrix/Schemes/HT_HG-U133_Plus_PM.root"){
  message(dataset_name)
  scheme.mine <- root.scheme(root_path)

  data.test3 <- import.data(scheme.mine, paste(dataset_name,"_dt",sep=""),filedir = "~/GoogleDrive/work/affymetrix/data/", celdir=cels_paths[1], verbose=T)
  if (length(cels_paths)>1){
    for (i in 2:length(cels_paths))
      cel_files=list.files(path=cels_paths[i],full.names = F)
    data.test3 <- addData(data.test3, celfiles=cel_files, celdir=cels_paths[i], verbose=T)
  }
}

import_xps_all=function(){
  import_xps_dataset("UNITI-1","~/GoogleDrive/work/shared/data/Janssen/UNITI-1 Biopsy/cels/")
  import_xps_dataset("UNITI-2","~/GoogleDrive/work/shared/data/Janssen/UNITI-2/cels/")
  import_xps_dataset("CERTIFI-biopsy","~/GoogleDrive/work/shared/data/public/human/Peters_et_al_NG_2017/GSE100833_RAW/biopsy/cels")
  import_xps_dataset("CERTIFI-blood","~/GoogleDrive/work/shared/data/public/human/Peters_et_al_NG_2017/GSE100833_RAW/blood/cels")
  import_xps_dataset("Arijs_Plos","~/GoogleDrive/work/shared/data/public/human/Arijs_et_al_PLoS1_2009/GSE16879_RAW/","~/GoogleDrive/work/affymetrix/Schemes/HG-U133_Plus.root")
  import_xps_dataset("Arijs_Gut_cohortB","~/GoogleDrive/work/shared/data/public/human/Arijs_et_al_Gut_2009/GSE12251_RAW/","~/GoogleDrive/work/affymetrix/Schemes/HG-U133_Plus.root")
  import_xps_dataset("Arijs_Gut_cohortA","~/GoogleDrive/work/shared/data/public/human/Arijs_et_al_Gut_2009/GSE14580_RAW/","~/GoogleDrive/work/affymetrix/Schemes/HG-U133_Plus.root")
  import_xps_dataset("Toedter","~/GoogleDrive/work/shared/data/public/human/Toedter_et_al_AMJG_2011/GSE23597_RAW/","~/GoogleDrive/work/affymetrix/Schemes/HG-U133_Plus.root")

}

normalize_xps_one_dataset=function(root_scheme,dataset,array_type="HG-U133_Plus"){
  data.test3 <- root.data(root_scheme,rootfile = paste("~/GoogleDrive/work/affymetrix/data/",dataset,"_dt_cel.root",sep=""))
  
  bgrd.rma <- express(data.test3, paste("tmpdt",dataset,"RMA",sep="_"), filedir="/tmp/", update=FALSE, bgcorrect.method="rma", bgcorrect.select="none", bgcorrect.option="pmonly:epanechnikov", bgcorrect.params=c(16384))
#  norm.low <- express(bgrd.rma, paste("tmpdt",dataset,"normLowess",sep="_"), filedir="/tmp/", update=FALSE,normalize.method="lowess", normalize.select="pmonly", normalize.option="transcript:all",normalize.logbase="log2", normalize.params=c(0.67, 3.0, 0.0, 0.0))
  norm.low <- express(bgrd.rma, paste("tmpdt",dataset,"normLowess",sep="_"), filedir="/tmp/", update=FALSE, normalize.method="quantile", normalize.select="all", normalize.option="transcript:together:none", normalize.logbase="0", normalize.params=c(0.0, 1.0))      
#  expr.mp <- express(norm.low, paste("tmpdt",dataset,"ExprMedpol",sep="_"), filedir="/tmp/", update=FALSE, summarize.method="medianpolish", summarize.select="pmonly", summarize.option="transcript", summarize.logbase="log2", summarize.params=c(10, 0.01, 1.0))
  expr.mp <- express(norm.low, paste("tmpdt",dataset,"ExprFarms",sep="_"), filedir="/tmp/", update=FALSE,summarize.method="farms", summarize.select="pmonly", summarize.option="transcript",summarize.logbase="log2", summarize.params=c(131, 0.5, 0.001 , 1.0, 0.1, 100, 1))
  # expr.mp  <- express(norm.low, paste("tmpdt",dataset,"ExprDFW",sep="_"), filedir="/tmp/", update=FALSE,summarize.method="dfw", summarize.select="pmonly", summarize.option="transcript",summarize.logbase="log2", summarize.params=c(3.0, 1.0, 0.01))
  m_probes <- validData(expr.mp)
  colnames(m_probes)=sapply(strsplit(colnames(m_probes),"\\."),head,1)
  
  probe2gene=probe_to_geneSymbol[[array_type]]
  ag=aggregate(m_probes[names(probe2gene),],list(probe2gene),mean)
  m=ag[,-1]
  rownames(m)=ag[,1]
  write.csv(file=paste("~/GoogleDrive/work/shared/analysis/analysis_iCD_paper_2018/microarray/output/expr_",dataset,".csv",sep=""),x = m,quote = F)
}



normalize_xps=function(){
  init_probe_gene_converters(rownames(ileum_ldm$dataset$umitab))
  scheme_hthgu133ppm <<- root.scheme("~/GoogleDrive/work/affymetrix/Schemes/HT_HG-U133_Plus_PM.root")
  datasets=c("UNITI-1","UNITI-2","CERTIFI-biopsy","CERTIFI-blood")
  for (dataset in datasets){
    message(dataset)
    normalize_xps_one_dataset(root_scheme = scheme_hthgu133ppm,dataset=dataset,array_type="HT_HG-U133_Plus_PM")
  }
  
  scheme_hgu133p <<- root.scheme("~/GoogleDrive/work/affymetrix/Schemes/HG-U133_Plus.root")
 # datasets=c("Arijs","Arijs_Gut1","Arijs_Gut2","Toedter")
  datasets=c("Arijs_Gut1","Arijs_Gut2","Toedter")
  for (dataset in datasets){
    message(dataset)
    normalize_xps_one_dataset(root_scheme = scheme_hgu133p,dataset=dataset,array_type = "HG-U133_Plus")
  }
}


make_scheme=function(){
  scmdir <- "~/GoogleDrive/work/affymetrix/Schemes"
  
  libdir <- "~/GoogleDrive/work/affymetrix/HT_HG-U133_Plus_PM_rev04"
  anndir <- "~/GoogleDrive/work/affymetrix/HT_HG-U133_Plus_PM-na36-annot-csv"
  scheme.test3 <- import.expr.scheme("HT_HG-U133_Plus_PM.scheme",filedir = scmdir, schemefile = file.path(libdir, "HT_HG-U133_Plus_PM.CDF"),probefile = file.path(libdir, "HT_HG-U133_Plus_PM.probe.tab"),annotfile = file.path(anndir, "HT_HG-U133_Plus_PM.na36.annot.csv"))

  libdir <- "~/GoogleDrive/work/affymetrix/CD_HG-U133_Plus_2/Full/HG-U133_Plus_2/LibFiles/"
  anndir <- "~/GoogleDrive/work/affymetrix/HG-U133_Plus_2-na36-annot-csv/"
  
  scheme.test3 <- import.expr.scheme("HG-U133_Plus.scheme",filedir = scmdir, schemefile = file.path(libdir, "HG-U133_Plus_2.CDF"),probefile = file.path(libdir, "HG-U133_Plus_2.probe.tab"),annotfile = file.path(anndir, "HG-U133_Plus_2.na36.annot.csv"))
  
}

#normalize.method="qspline"
#summary.method="liwond"

#norm_dataset("~/GoogleDrive/work/shared/data/Janssen/UNITI-1 Biopsy/cels/",
#              "~/GoogleDrive/work/shared/data/Janssen/UNITI-1 Biopsy/UNITI-1_biopsy.csv")

#norm_dataset("~/GoogleDrive/work/shared/data/Janssen/UNITI-2/PBOwk8RNR subset_wk0wk8samples/",
#              "~/GoogleDrive/work/shared/data/Janssen/UNITI-2/PBOwk8RNR subset_wk0wk8.csv")

#norm_dataset("~/GoogleDrive/work/shared/data/Janssen/UNITI-2/USTwk8RNR subset_wk0wk8samples/",
#              "~/GoogleDrive/work/shared/data/Janssen/UNITI-2/USTwk8RNR_subset_wk0wk8.csv")

#norm_dataset("~/GoogleDrive/work/shared/data/public/human/Peters_et_al_NG_2017/GSE100833_RAW/biopsy/cels",
#              "~/GoogleDrive/work/shared/data/public/human/Peters_et_al_NG_2017/GSE100833_RAW/biopsy/certifi_biopsy.csv")

#norm_dataset("~/GoogleDrive/work/shared/data/public/human/Peters_et_al_NG_2017/GSE100833_RAW/blood/cels",
#              "~/GoogleDrive/work/shared/data/public/human/Peters_et_al_NG_2017/GSE100833_RAW/blood/certifi_blood.csv")

#norm_dataset("~/GoogleDrive/work/shared/data/public/human/Arijs_et_al_PLoS1_2009/GSE16879_RAW/",
#              "~/GoogleDrive/work/shared/data/public/human/Arijs_et_al_PLoS1_2009/Arijs_probes_norm.csv")



norm_arijs_affy=function(reg=10){
  library(affy)
  Data <- ReadAffy(widget=F,filenames = list.files("~/GoogleDrive/work/shared/data/public/human/Arijs_et_al_PLoS1_2009/GSE16879_RAW/",full.names = T))
 
  eset <- expresso(Data, normalize.method="qspline",bgcorrect.method="rma",pmcorrect.method="pmonly",summary.method="liwong")
#  eset <- expresso(Data, normalize.method="quantiles",bgcorrect.method="rma",pmcorrect.method="pmonly",summary.method="liwong")
#  eset <- expresso(Data,bgcorrect.method="rma",normalize.method="quantiles",pmcorrect.method="pmonly",summary.method="medianpolish")
 
  m_probes=log10(reg+exprs(eset))
  probe2gene=probe_to_geneSymbol[["HG-U133_Plus"]]
  ag=aggregate(m_probes[names(probe2gene),],list(probe2gene),mean)
  m=log(ag[,-1])
  rownames(m)=ag[,1]
  colnames(m)=unlist(strsplit(colnames(m),".CEL"))
  write.csv(m,file="~/GoogleDrive/work/shared/analysis/analysis_iCD_paper_2018/microarray/output/expr_Arijs_affy.csv",row.names = T)
}

norm_arijs_frma=function(){
  library(frma)
  frma_output <- frma(Data,summarize = "median_polish")
  m_probes=exprs(frma_output)
  probe2gene=probe_to_geneSymbol[["HG-U133_Plus"]]
  ag=aggregate(m_probes[names(probe2gene),],list(probe2gene),mean)
  m=ag[,-1]
  rownames(m)=ag[,1]
  colnames(m)=unlist(strsplit(colnames(m),".CEL"))
  write.csv(m,file="~/GoogleDrive/work/shared/analysis/analysis_iCD_paper_2018/microarray/output/expr_Arijs_frma.csv",row.names = T)
  
  
}

norm_arijs_xps_test=function(){
  scheme.test3 <- root.scheme("~/GoogleDrive/work/affymetrix/Schemes/HG-U133_Plus.root")
  dataset="Arijs"
  message(dataset)
  normalize_xps_one_dataset(scheme.test3,dataset)
}

norm_certifi_test=function(){
  dataset="CERTIFI-biopsy"
  root_scheme=root.scheme("~/GoogleDrive/work/affymetrix/Schemes/HT_HG-U133_Plus_PM.root")
  array_type="HT_HG-U133_Plus_PM"
  data.test3 <- root.data(root_scheme,rootfile = paste("~/GoogleDrive/work/affymetrix/data/",dataset,"_dt_cel.root",sep=""))
  
  bgrd.rma <- express(data.test3, paste("tmpdt",dataset,"RMA",sep="_"), filedir="/tmp/", update=FALSE, bgcorrect.method="rma", bgcorrect.select="none", bgcorrect.option="pmonly:epanechnikov", bgcorrect.params=c(16384))
  norm.quant <- express(bgrd.rma, paste("tmpdt",dataset,"normQuant",sep="_"), filedir="/tmp/", update=FALSE, normalize.method="mean", normalize.select="all", normalize.option="probeset:all", normalize.logbase="0", normalize.params=c(0.0, -1))     
  expr.mp2 <- express(norm.quant, paste("tmpdt",dataset,"ExprFarms",sep="_"), filedir="/tmp/", update=FALSE,summarize.method="farms", summarize.select="pmonly", summarize.option="transcript",summarize.logbase="log2", summarize.params=c(131, 0.5, 0.001, 1.0, 0.00001, 100, 1))
  m_probes <- validData(expr.mp2)
  colnames(m_probes)=sapply(strsplit(colnames(m_probes),"\\."),head,1)
  
  expr.adf <- express(norm.quant, paste("tmpdt",dataset,"ExprAvgDif",sep="_"), filedir="/tmp/", update=FALSE,summarize.method="avgdiff", summarize.select="none", summarize.option="transcript", summarize.logbase="0", summarize.params=c(3.0))
  m_probes3 <- validData(expr.adf)
  
  probe2gene=probe_to_geneSymbol[[array_type]]
  ag=aggregate(m_probes[names(probe2gene),],list(probe2gene),mean)
  m=ag[,-1]
  rownames(m)=ag[,1]
  write.csv(file=paste("~/GoogleDrive/work/shared/analysis/analysis_iCD_paper_2018/microarray/output/expr_",dataset,".csv",sep=""),x = m,quote = F)
  
  
  data.test3 <- root.data(root_scheme,rootfile = paste("~/GoogleDrive/work/affymetrix/data/",dataset,"_dt_cel.root",sep=""))
  
}
