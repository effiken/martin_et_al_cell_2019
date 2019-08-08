library(Matrix)
library(Matrix.utils)

adt_list=c('CD11b.IH', 'CD11c.IH', 'CD123.IH', 'CD14.IH', 'CD141.IH', 'CD16.IH', 'CD19.IH', 'CD1c.IH', 'CD206.IH', 'CD24.IH', 'CD27.IH', 'CD3.IH', 'CD33.IH', 'CD38.IH', 'CD4.IH', 'CD56.IH', 'CD62L.IH', 'CD64.IH', 'CD66b.IH', 'CD8.IH', 'HLADR.IH','CD10', 'CD103', 'CD117', 'CD11b', 'CD11c', 'CD123', 'CD127', 'CD138', 'CD14', 'CD141', 'CD16', 'CD161', 'CD163', 'CD169', 'CD183..CXCR3.', 'CD185..CXCR5.', 'CD19', 'CD194..CCR4.', 'CD195..CCR5.', 'CD196..CCR6.', 'CD197..CCR7.', 'CD1c', 'CD20', 'CD206', 'CD24', 'CD25', 'CD26', 'CD27', 'CD272..BTLA.', 'CD273..PD.L2.', 'CD274..PD.L1.', 'CD278..ICOS.', 'CD279..PD.1.', 'CD28', 'CD3', 'CD314..NKG2D.', 'CD33', 'CD34', 'CD370', 'CD38', 'CD4', 'CD40', 'CD45', 'CD45RA', 'CD5', 'CD56', 'CD57', 'CD58', 'CD62L', 'CD66b', 'CD69', 'CD8', 'CD86', 'HLA.DR', 'IgD', 'TCR_gamma_delta', 'TIGIT..VSTM3.', 'XCR1', "APC.IH", "CD11c.IHM", "I.A.I.E.IHM", "PE.IH")
#CD11b-IH	CD11c-IH	CD123-IH	CD14-IH	CD141-IH	CD16-IH	CD19-IH	CD1c-IH	CD206-IH	CD24-IH	CD27-IH	CD3-IH	CD33-IH	CD38-IH	CD4-IH	CD56-IH	CD62L-IH	CD64-IH	CD66b-IH	CD8-IH	HLADR-IH
# read_mtx should get the path of the 10x sampled directory. 
#It expects the following directory sturcture:
#/sample_x
#         /filtered/
#                  /matrix.mtx
#				   /barcodes.tsv
#				   /genes.tsv
#         /raw
#                  /matrix.mtx
#				   /barcodes.tsv
#				   /genes.tsv
#
# type can be either "filtered" or "raw"

read_mtx=function(dir,type="raw",min_umis=100,gene_ID_converter=NULL){
  if (!is.na(type)){
    if (!dir.exists(paste(dir,type,sep="/"))){
      orig_type=type
      type_u= paste(toupper(substring(type, 1,1)), substring(type, 2),sep="")
      type_l= paste(tolower(substring(type, 1,1)), substring(type, 2),sep="")
      if (type==type_u){
        type=type_l
      }else{
        type=type_u
      }
      if(!dir.exists(paste(dir,type,sep="/"))){
        stop("Error! Directory does not exist! ",paste(dir,type,sep="/"))
      }
    }
    fn=paste(dir,type,"matrix.mtx",sep="/")
    
    genes_fn=paste(dir,type,"genes.tsv",sep="/")
    ######
    # AL edit 4/5/19
    # Cellranger V3 (necessary for the V3 chemistry kit) outputs the gene names as features.tsv
    if(!file.exists(genes_fn)){
      genes_fn=paste(dir,type,"features.tsv",sep="/")
    }
    ######
    
    barcodes_fn=paste(dir,type,"barcodes.tsv",sep="/")
  }
  else{
    fn=paste(dir,"matrix.mtx",sep="/")
    genes_fn=paste(dir,"genes.tsv",sep="/")
    barcodes_fn=paste(dir,"barcodes.tsv",sep="/")
  }
  if (!file.exists(fn)){
    fn=paste(fn,".txt",sep="")
  }
  message("Reading ",fn)
  if (!file.exists(fn)){
    message("ERROR! ",fn," not found!")
    return(NULL)
  }
  tab=readMM(fn)
  #  message(paste(dim(tab),collapse=","))
  gene_tab=read.delim(file=genes_fn,stringsAsFactors=F,header=F)
  ######
  # AL edit 4/5/19
  # Accommodating Cellranger V3 output
  if(genes_fn==paste(dir,type,"features.tsv",sep="/")){
    
    rownames(tab)=gene_tab$V2
    gene_tab=gene_tab$V2
  }
  else{
    rownames(tab)=gene_tab[,ncol(gene_tab)]
  }
  ######
  
  if (!is.null(gene_ID_converter)){
    rownames(tab)=gene_ID_converter(rownames(tab))
  }

  barcodes=read.delim(file=barcodes_fn,stringsAsFactors=F,header=F)[,1]
  if (length(barcodes)!=ncol(tab)){
    warning("Number of barcodes does not equal to the number of columns!")
    colnames(tab)=paste("B",1:ncol(tab),sep="")
  }
  else{
    colnames(tab)=barcodes
  }
  
  barcode_hist=table(colnames(tab))
  unique_barcodes=names(barcode_hist)[barcode_hist==1]
  tab=tab[,unique_barcodes]
  tab=tab[,Matrix::colSums(tab)>min_umis]
  return(pool_rows_by_id(tab))
}

pool_rows_by_id=function(tab){
  genehist=table(rownames(tab))[rownames(tab)]
  tab2=tab[which(genehist==1),]
  tab3=tab[which(genehist>1),,drop=F]
  if (nrow(tab3)>0){
    tab3u=aggregate(tab3,rownames(tab3))
    return(rbind(tab2,tab3u))
  }
  else{
    return(tab)
  }
}

read_mtx_and_split_by_10x_batch=function(path,libname,type="raw",noise_interval=c(40,100),cell_interval=c(201,20000),gene_ID_converter=NULL,genes=NULL){
  l=list()
  min_umis=min(c(noise_interval,cell_interval))
  umitab=read_mtx(paste(path,libname,sep="/"),type,min_umis,gene_ID_converter=gene_ID_converter)
  
  if (!is.null(genes)){
    umitab=umitab[genes,]
  }
  lib_batch=sapply(strsplit((colnames(umitab)),"-"),function(x){x[2]})
  if (all(is.na(lib_batch))){
    lib_batch[is.na(lib_batch)]=1
  }
  lib_batches=unique(lib_batch)
  if (length(lib_batches)>1){
    libnames=lib_batches
  }
  else{
    libnames=as.character(libname)
  }
  
  for (i in 1:length(lib_batches)){
    tab=umitab[,lib_batch==lib_batches[i]]
    numis=Matrix::colSums(tab)
    noise_counts=Matrix::rowSums(tab[,numis>noise_interval[1]&numis<noise_interval[2]])
    tab=tab[,numis>cell_interval[1]&numis<cell_interval[2]]

    
    ###########
    # A.L. fix 3/28/19
    # read adt files for CITEseq samples that have ADT matrix but not HTOs
    
    txt_files=grep(list.files(paste(path,"/",libnames[i],sep="")),pattern = "\\.txt$",val=T)
    if(length(grep("adt",txt_files)>=1)){
    adt_fn=paste(path,"/",libname,"/",grep(txt_files,pattern = "adt",val=T),sep="")
    adt=read.table(adt_fn,header=T,sep="\t",row.names = 1)
    adt_columns=unique(c(grep("^adt_",colnames(adt),val=T),intersect(colnames(adt),adt_list)))
    l[[libnames[i]]]=list(umitab=tab,noise_counts=noise_counts,adttab=t(adt[colnames(tab),adt_columns]))
    }
    else{
      l[[libnames[i]]]=list(umitab=tab,noise_counts=noise_counts)
    }
    ###########
    
  }
  
  return(l)
}

read_mtx_and_split_by_hto=function(path,libname,type="raw",noise_interval=c(40,100),cell_interval=c(201,20000),hto_UMI_ratio_thresh=5,min_umi_maxhto=20,hto_to_sample=NULL,gene_ID_converter=NULL,genes=NULL){
   min_umi=min(c(noise_interval,cell_interval))
  
  #####
  ## Edit AL 6/21/19 - make function compatible with cellranger feature barcoding
   
   if(dir.exists(paste(path,libname,"features",sep="/"))){
     message("Reading ",paste(path,libname,"features",sep="/"))
     
     alldat=readMM(gzfile(paste(path,libname,"features","matrix.mtx.gz",sep="/")))
     features=read.delim(gzfile(paste(path,libname,"features","features.tsv.gz",sep="/")),stringsAsFactors=F,header=F)
     barcodes = read.table(gzfile(paste(path,libname,"features","barcodes.tsv.gz",sep="/")),stringsAsFactors=F)[,1]
     
     if (!is.null(gene_ID_converter)){
       features$V2[features$V3=="Gene Expression"]=gene_ID_converter(features$V2[features$V3=="Gene Expression"])
     }
     
     
     colnames(alldat) <- barcodes
     rownames(alldat) <- features$V2
     
     barcode_hist=table(colnames(alldat))
     unique_barcodes=names(barcode_hist)[barcode_hist==1]
     alldat=alldat[,unique_barcodes]
     
     ####treat similarly to read_mtx function:
     alldat <- alldat[,Matrix::colSums(alldat[features$V3=="Gene Expression",])>min_umi]
     
     u <- alldat[features$V3=="Gene Expression",]
     
     genehist=table(rownames(u))[rownames(u)]
     tab2=u[which(genehist==1),]
     tab3=u[which(genehist>1),,drop=F]
     if (nrow(tab3)>0){
       tab3u=do.call(rbind,lapply(split(rownames(tab3),rownames(tab3)),function(x){colSums(tab3[x,])}))
       u <- rbind(tab2,tab3u)
     }
     ######
     
     hto <- t(alldat[features$V3=="Custom",])
     rm(alldat)
     gc()
     
     hto_columns <- grep("^HTO",colnames(hto))
     adt_columns=unique(c(grep("^adt_",colnames(hto),val=T),intersect(colnames(hto),adt_list)))
     
   } 
   else{
     
     #####
     if (dir.exists(paste(path,libname,"citeseq",sep="/"))){
       hto=readMM(file =paste(path,libname,"citeseq","matrix.mtx",sep="/"))
       barcodes=read.table(paste(path,libname,"citeseq","barcodes.tsv",sep="/"),stringsAsFactors=F)[,1]
       colnames(hto)=barcodes
       features=read.table(paste(path,libname,"citeseq","features.tsv",sep="/"),stringsAsFactors=F)[,1]
       rownames(hto)=features
       
     }
     else{
       txt_files=grep(list.files(paste(path,"/",libname,sep="")),pattern = "\\.txt$",val=T)
       hto_fn_suf=grep(txt_files,pattern = "adt_hto",val=T)
       if (length(hto_fn_suf)>0){
         hto_fn=paste(path,"/",libname,"/",hto_fn_suf,sep="")
         message("Reading ",hto_fn)
         hto=read.table(hto_fn,header=T,sep="\t",row.names = 1)
         colnames(hto)=gsub("HBC","HTO_",colnames(hto))
         ##########
         # A.L. fix 3/28/19:
         # some early HTO matrices had colnames "HTO" instead of "HTO_"
         if(length(grep("^HTO",colnames(hto)))>0 & length(grep("^HTO_",colnames(hto)))==0){
           colnames(hto)=gsub("HTO","HTO_",colnames(hto))
         }
         ##########
       }
       else{
         hto_fn_suf=grep(txt_files,pattern = "hto",val=T)
         if (length(hto_fn_suf)>0){
           hto_fn=paste(path,"/",libname,"/",hto_fn_suf,sep="")
           adt_fn=paste(path,"/",libname,"/",grep(txt_files,pattern = "adt",val=T),sep="")
           message("Reading ",hto_fn)
           hto=read.table(hto_fn,header=T,sep="\t",row.names = 1)
           colnames(hto)=gsub("HBC","HTO_",colnames(hto))
           ##########
           # A.L. fix 3/28/19:
           # some early HTO matrices had colnames "HTO" instead of "HTO_"
           if(length(grep("^HTO",colnames(hto)))>0 & length(grep("^HTO_",colnames(hto)))==0){
             colnames(hto)=gsub("HTO","HTO_",colnames(hto))
           }
           ##########
           adt=read.table(adt_fn,header=T,sep="\t",row.names = 1)
           adt=adt[intersect(rownames(adt),rownames(hto)),]
           a=merge(adt,hto,by="row.names")
           hto=a[,-1]
           rownames(hto)=a[,1]
           rm(a)
         }
       }  
     }
     u=read_mtx(paste(path,"/",libname,sep=""),type = type,min_umis = min_umi,gene_ID_converter=gene_ID_converter)
     colnames(u)=gsub("-1","",colnames(u))
     if (!is.null(genes)){
       u=u[genes,]
     }
     hto_columns=match(names(hto_to_sample),colnames(hto))
     adt_columns=unique(c(grep("^adt_",colnames(hto),val=T),intersect(colnames(hto),adt_list)))
     
   }
  
  numis_orig=colSums(u)
  unhashed=setdiff(colnames(u),rownames(hto))
  message("Number of barcodes in umitab are not in hto/adt tab : ",length(unhashed))

  # plot HTO vs mRNA UMIs
  numis=colSums(u)
  htomax=apply(hto[,hto_columns,drop=F],1,max,na.rm=T)
  unhashed=c(unhashed,rownames(hto)[htomax<min_umi_maxhto])
  
  htosec=apply(hto[,hto_columns,drop=F],1,function(x){sort(x,decreasing=T)[2]})
  hto_ratio=(htomax)/(htosec)
  barcodes=setdiff(intersect(names(hto_ratio),colnames(u)),unhashed)
  par(mar=c(7,5,3,3))
  plot(log2(1+Matrix::colSums(u)[barcodes]),hto_ratio[barcodes],log="y",pch=".",ylab="#UMIs (max HTO) / #UMIs (2nd HTO) ",xlab="#mRNA UMIs",panel.first=grid(lty=1))
  mtext(text = "1st/2nd HTO",side = 3)
  #  mixmdl = normalmixEM(hto_ratio)
  
  abline(h=hto_UMI_ratio_thresh,col=2,lty=3,lwd=3)
  
  singlets=barcodes[((1+htomax)/(1+htosec))[barcodes]>hto_UMI_ratio_thresh]
  doublets=barcodes[((1+htomax)/(1+htosec))[barcodes]<hto_UMI_ratio_thresh]
  doublets=doublets[!is.na(doublets)]
  singlets=singlets[!is.na(singlets)]
  unhashed=setdiff(colnames(u),c(singlets,doublets))
  
  v=c(mean(names(numis_orig[numis_orig>200])%in%singlets),mean(names(numis_orig[numis_orig>200])%in%doublets))
  v=c(v,1-sum(v))
  
  barplot(v,names.arg = c("Singlets","Multiplet","Not hashed"),las=2,ylim=c(0,1),ylab="Cell-barcode Fraction")
  mtext(text = "Cell-barcode breakdown",side = 3)
  print(round(v,digits=2))
  v=c(singlets=sum(u[,singlets]),multiplets=sum(u[,doublets]),"Not hashed"=sum(u[,unhashed]))
  v=v/sum(v)
  barplot(v,las=2,ylab="UMI fraction",ylim=c(0,1))
  mtext(text = "UMI breakdown",side = 3)
  print(round(v,digits=2))
  bs=colnames(u)[numis>100]
  col=(intersect(bs,rownames(hto))%in%singlets+2*intersect(bs,rownames(hto))%in%doublets+3*intersect(bs,rownames(hto))%in%unhashed)
  plot(log2(1+numis[intersect(bs,rownames(hto))]),log2(1+apply(hto[intersect(bs,rownames(hto)),hto_columns],1,max)),pch=".",xlab="log2 #UMIs (mRNA)",ylab="log2 #UMIs(max HTO)",col=col,panel.first = grid(lty=1))
  mtext(text = "max HTO vs. #UMIs",side = 3)
  plot(log2(1+numis[intersect(bs,rownames(hto))]),log2(1+rowSums(hto[intersect(bs,rownames(hto)),adt_columns])),pch=".",xlab="log2 #UMIs (mRNA)",ylab="log2 #UMIs(adts)",col=col,panel.first=grid(lty=1))
  mtext(text = "ADTs vs mRNAs",side = 3)
  v=log2(1+numis[intersect(bs,rownames(hto))])
  plot(quantile(v[singlets],0:100/100,na.rm=T),0:100/100,type='l',panel.first=grid(lty=1),xlim=range(v),ylim=c(0,1),xlab="log2 #UMIs (mRNA)",ylab="Cumulative Fraction")
  lines(quantile(v[doublets],0:100/100,na.rm=T),0:100/100,type='l',col=2)
  lines(quantile(v[unhashed],0:100/100,na.rm=T),0:100/100,type='l',col=3)
  legend("bottomright",legend=c("Singlets","Multiplet","Not hashed"),col=1:3,lty=1)
  mtext(text = "mRNAs",side = 3)
  
  cell_to_hto=(colnames(hto)[hto_columns])[apply(hto[singlets,hto_columns],1,which.max)]
  names(cell_to_hto)=singlets
  barplot(table(cell_to_hto),las=2,ylab="#Cells")
  mtext(text = "#cells / HTO",side = 3)
  numis=Matrix::colSums(u)
  noise_counts=rowSums(u[,numis>noise_interval[1]&numis<noise_interval[2]])
  
  l=list()
  
  cell_to_sample=hto_to_sample[cell_to_hto]
  names(cell_to_sample)=names(cell_to_hto)
  samples=unique(cell_to_sample)
  for (samp in samples){
    tab=u[,names(cell_to_sample)[cell_to_sample==samp]]
    numis=Matrix::colSums(tab)
    tab=tab[,numis>cell_interval[1]&numis<cell_interval[2]]
    l[[samp]]=list(umitab=tab,noise_counts=noise_counts,doublets_umitab=u[,doublets],unhashed_umitab=u[,unhashed],adttab=t(hto[colnames(tab),adt_columns]),htotab=t(hto[colnames(tab),hto_columns]))
  }
  return(l)
}



#
convert_tab_to_mtx=function(umitab,output_path){
  write.table(rownames(umitab),file=paste(output_path,"/genes.tsv",sep=""),row.names = F,col.names = F,quote = F)
  write.table(colnames(umitab),file=paste(output_path,"/barcodes.tsv",sep=""),row.names = F,col.names = F,quote = F)
  writeMM(Matrix(as.matrix(umitab)),file = paste(output_path,"/matrix.mtx",sep=""))
}



downsample_old=function(u,min_umis){
  non_zero_mask=rowSums(u)!=0
  all_op=1:nrow(u[non_zero_mask,])
  base_tab=rep(0,sum(non_zero_mask))
  names(base_tab)=all_op
  
  downsamp_one=function(x,n){
    tab=base_tab
    tab2=table(sample(rep(all_op,x),size = n,replace=F))
    tab[names(tab2)]=tab2
    return(tab)
  }
  
  cell_mask=colnames(u)[Matrix::colSums(u,na.rm=T)>min_umis]  
  message("Downsampling ", length(cell_mask), " cells to ",min_umis," UMIs")
  
  chunk_size=1000
  breaks=unique(c(seq(from=1,to = length(cell_mask),by = chunk_size),length(cell_mask)+1))
  ds=Matrix(0,nrow =nrow(u),length(cell_mask),dimnames = list(rownames(u),cell_mask))
  for (i in 1:(length(breaks)-1)){
    ds[non_zero_mask,breaks[i]:(breaks[i+1]-1),drop=F]=Matrix(apply(u[non_zero_mask,cell_mask[breaks[i]:(breaks[i+1]-1)],drop=F],2,downsamp_one,min_umis))
  }
  
  return(ds)
}

downsample=function(u,min_umis,parallel=F,chunk_size=100){
  
  ######## 
  # A.L. fix 3/28/19
  # My Windows O.S. does not support the parallel feature
  if(Sys.info()[1]=="Windows"){
    parallel=F
  }
  #########
  
  
  non_zero_mask=Matrix::rowSums(u)!=0
  all_op=1:nrow(u[non_zero_mask,])
  base_tab=rep(0,sum(non_zero_mask))
  names(base_tab)=all_op
  
  downsamp_one=function(x,n){
    tab=base_tab
    tab2=table(sample(rep(all_op,x),size = n,replace=F))
    tab[names(tab2)]=tab2
    return(tab)
  }
  
  cell_mask=colnames(u)[Matrix::colSums(u,na.rm=T)>min_umis]  
  print(paste("Downsampling ", length(cell_mask), " cells to ",min_umis," UMIs",sep=""))
  
  
  breaks=unique(c(seq(from=1,to = length(cell_mask),by = chunk_size),length(cell_mask)+1))
  ds=Matrix(0,nrow =nrow(u),length(cell_mask),dimnames = list(rownames(u),cell_mask))
  
  if (parallel){
    library(parallel)
    numCores=detectCores()
    njobs=min(numCores,length(breaks)-1)
    x=mclapply(1:(length(breaks)-1),function(i){Matrix(apply(u[non_zero_mask,cell_mask[breaks[i]:(breaks[i+1]-1)],drop=F],2,downsamp_one,min_umis))},mc.cores =njobs)
    ds[non_zero_mask,]=do.call(cbind,x)
  }
  else{
    for (i in 1:(length(breaks)-1)){
      ds[non_zero_mask,breaks[i]:(breaks[i+1]-1),drop=F]=Matrix(apply(u[non_zero_mask,cell_mask[breaks[i]:(breaks[i+1]-1)],drop=F],2,downsamp_one,min_umis))
    }
  }
  
  return(ds)
}




split_sparse=function(sparse_umitab,by,rows_or_columns="columns"){
  
  groups=sort(as.character(unique(by)))
  l=list()
  
  if (rows_or_columns=="columns"){
    for (cl in groups){
      l[[cl]]=sparse_umitab[,by==cl,drop=F]
    }
  } 
  else{
    if (rows_or_columns=="rows"){
      for (cl in groups){
        l[[cl]]=sparse_umitab[by==cl,,drop=F]
      }
    }
  }
  return(l)
}



read_and_compile_umitabs=function(sample_IDs,input_paths,output_path,filtered_or_raw="raw",prefix="",noise_interval=c(100,25000),cell_interval=c(200,25000),ds_numis=c(200,500,1000,2000),sample_annots_path=NULL,cite_seq_params=NULL,rename_cell_ids=F,only_load=T,output_qc_path=NULL,gene_ID_converter=NULL,genes=NULL){
  sample_IDs=as.character(sample_IDs)
  umitab=c()
  noise_models=c()
  cell_to_batch=c()
  adt_l=list()
  doublet_umitab=list()
  unhashed_umitab=list()
  hto_l=list()
  if (is.null(cite_seq_params)){
    cite_seq_params=list()
    cite_seq_params$hto_UMI_ratio_thresh=5
    cite_seq_params$min_umi_maxhto=20
  }
  

  path_param_exists=F
  if (!is.null(sample_annots_path)){
   
    
    sample_annots_tab=as.data.frame(read.csv(sample_annots_path,stringsAsFactors = F,row.names = 1))
    annots_tab=sample_annots_tab[sample_IDs,]
    path_param_exists="path"%in%colnames(annots_tab)
    amp_batches=unique(annots_tab[,"amp_batch_ID"])
    
    amp_batche_to_old_libname=annots_tab$old_lib_name[match(amp_batches,annots_tab$amp_batch_ID)]
    names(amp_batche_to_old_libname)=amp_batches
    
    if ("HTO"%in%colnames(annots_tab)){
      cite_seq=(!is.na(annots_tab[,"HTO"]))[match(amp_batches,annots_tab[,"amp_batch_ID"])]
    }
    
    else{
      cite_seq=rep(F,length(amp_batches))
    }
  }
  else{
    cite_seq=rep(F,length(sample_IDs))
    amp_batches=sample_IDs
  }
  
  
  if (length(input_paths)==1){
    if (!path_param_exists){
      message("Assuming all data is in ",input_paths)
      input_paths=rep(input_paths,length(amp_batches))
    }
  }
  else{
    if (length(input_paths)!=length(amp_batches)){
      stop("ERROR! Number of paths should be ",length(input_paths)," or alternatively 'path' column should be added to the sample_annots table.")
    }
  }
  if (only_load){
    message("Reading batches ",paste(amp_batches,collapse=","))
  }
  else{
    message("compiling batches ",paste(amp_batches,collapse=","))
  }
  

  
  for (i in 1:length(amp_batches)){
    
    message("Amp batch ", amp_batches[i])

    suff=strsplit(input_paths[i], split="\\.")[[1]][2]
    if (is.na(suff)){
      suff=""
    }
    if (suff=="h5"){
      require(TENxGenomics)
      group=setdiff(unique(h5ls(input_paths[i])$group),"/")
      tenx <- TENxGenomics(input_paths[i],group)
      v=colnames(tenx)
      cell_to_amp_batch=sapply(strsplit(v,"-"),head,1)
      names(cell_to_amp_batch)=v
      rm(v)
      l=list()
      l[[amp_batches[i]]]=list()
      tab=Matrix(as.matrix(tenx[,cell_to_amp_batch==amp_batches[i]]))
      if (!is.null(gene_ID_converter)){
        converted_genes=gene_ID_converter(rownames(tab))
        tab=tab[!is.na(converted_genes),]
        rownames(tab)=converted_genes[!is.na(converted_genes)]
      }
      tab=pool_rows_by_id(tab)
      if (!is.null(genes)){
        tab=tab[genes,]
      }
      numis=Matrix::colSums(tab)
      noise_counts=Matrix::rowSums(tab[,numis>noise_interval[1]&numis<noise_interval[2]])
      tab=tab[,numis>cell_interval[1]&numis<cell_interval[2]]

      l[[amp_batches[i]]]$umitab=tab
      l[[amp_batches[i]]]$noise_counts=noise_counts
      names(l)=amp_batches[i]
    }
    else if (!cite_seq[i]){
      
      l=read_mtx_and_split_by_10x_batch(path = input_paths[i],libname = amp_batches[i],type = filtered_or_raw,noise_interval =noise_interval,cell_interval = cell_interval,gene_ID_converter=gene_ID_converter,genes=genes) 
      if (exists("annots_tab")){
        sample_ID=rownames(annots_tab)[annots_tab[,"amp_batch_ID"]==amp_batches[i]]
        if (length(l)==1){
          names(l)=sample_ID
        }
        else{
          names(l)=paste(sample_ID,names(l),sep="_")
        }
      }
    }
    else{
      if (!is.null(output_qc_path)){
        pdf(file = paste(output_qc_path,"/",amp_batches[i],".pdf",sep=""),width = 12,height = 8)
        layout(matrix(c(1,1,1,1,2,5,6,7,3,4,8,9),3,4,byrow=T),heights=c(1,5,5))
        par(mar=c(0,0,0,0))
        plot.new()
        if (is.null(annots_tab$old_lib_name)){
          mtext(side=3,text = paste("AB",amp_batches[i]),line=-2)
        }
        else{
          mtext(side=3,text = paste("AB",amp_batches[i],"-",amp_batche_to_old_libname[as.character(amp_batches[i])]),line=-2)
        }
      }
      mask_i=annots_tab[,"amp_batch_ID"]==amp_batches[i]
      hto_list=strsplit(as.character(annots_tab[mask_i,"HTO"]),",")
      hto_to_sample=rep(rownames(annots_tab)[mask_i],sapply(hto_list,length))
      names(hto_to_sample)=paste("HTO",unlist(hto_list),sep="_")
      l=read_mtx_and_split_by_hto(input_paths[i],amp_batches[i],hto_UMI_ratio_thresh=cite_seq_params$hto_UMI_ratio_thresh,noise_interval=noise_interval,cell_interval=cell_interval,min_umi_maxhto=cite_seq_params$min_umi_maxhto,hto_to_sample=hto_to_sample,gene_ID_converter=gene_ID_converter,genes=genes)
      if (!is.null(output_qc_path)){
        dev.off()
      }
      #    l=read_mtx_and_split_by_hto(path = input_paths[i],libname = sample_IDs[i],type = filtered_or_raw,noise_interval =noise_interval,cell_interval = cell_interval,hto_UMI_ratio_thresh = 5,min_umi_maxhto = 20,hto_to_sample=hto_to_sample_list[[i]]) 
    }
    
    for (b in names(l)){
      message(names(b))
      env=new.env()
      env$ds_numis=ds_numis
      env$min_numis=noise_interval[1]
      env$lib_name=b
      env$ds=list()
      env$noise_model=matrix(l[[b]]$noise_counts/sum(l[[b]]$noise_counts),ncol=1)
      rownames(env$noise_model)=rownames(l[[b]]$umitab)
      colnames(env$noise_model)=b
      numis=Matrix::colSums(l[[b]]$umitab)
      numis_mask=numis>=cell_interval[1]&numis<=cell_interval[2]
      
      env$umitab=l[[b]]$umitab[,which(numis_mask)]
      env$adttab=l[[b]]$adttab[,which(numis_mask)]
      env$htotab=l[[b]]$htotab[,which(numis_mask)]
      env$doublets_umitab=l[[b]]$doublets_umitab
      env$unhashed_umitab=l[[b]]$unhashed_umitab
      if (rename_cell_ids){
        colnames(env$umitab)=paste(b,1:ncol(env$umitab),sep="_")
        colnames(env$adttab)=colnames(env$umitab)
        colnames(env$htotab)=colnames(env$umitab)
      }
      
      if (!only_load){
        for (j in 1:length(ds_numis)){
          env$ds[[j]]=downsample(env$umitab,ds_numis[j])
        }
        
        fn=paste(output_path,"/",prefix,"data_",b,".rd",sep="")
        message("writing ",fn)
        save(file=fn,list = names(env),envir = env)
      }
      else{
        cell_to_batch_b=rep(b,ncol(env$umitab))
	colnames(env$umitab)=paste(b,colnames(env$umitab),sep="_")
        names(cell_to_batch_b)=colnames(env$umitab)
        cell_to_batch=c(cell_to_batch,cell_to_batch_b)
        umitab=cBind(umitab,env$umitab)
        noise_models=cbind(noise_models,env$noise_model)
        if (!is.null(env$adttab)){
          adt_l[[b]]=env$adttab
          hto_l[[b]]=env$htotab
        }
        if (!is.null(env$doublets_umitab)){
          doublet_umitab[[b]]=env$doublets_umitab
        }
        if (!is.null(env$unhashed_umitab)){
          unhashed_umitab[[b]]=env$unhashed_umitab
        }
      }
    }
   }
  
  return(list(cell_to_batch=cell_to_batch,umitab=umitab,noise_models=noise_models,adt_l=adt_l,hto_l=hto_l,doublet_umitab=doublet_umitab,unhashed_umitab=unhashed_umitab))
}





# read_multiple_mtx itarates over libnames and calls to read_mtx with dir=path/libnames[i]
# 
# returned list:
# umitabs_by_batch - a list containing the umitab of each library
# cell_to_batch - a vector mapping from cell_names_to_library_name
# umitab_sampled - a umitab containing all the sampled cells
# umitan - a umitab containg all the cells

read_multiple_mtx=function(paths,libnames,sample_annots_path=NULL,noise_interval=c(40,100),cell_interval=c(201,20000),rename_cell_ids=T,  type="raw",cite_seq_params=NULL,gene_ID_converter=NULL){
  if (is.null(cite_seq_params)){
    cite_seq_params=list()
    cite_seq_params$hto_UMI_ratio_thresh=5
    cite_seq_params$min_umi_maxhto=20
  }
  l=list()
  noise_l=list()
  cell_to_batch=c()
  genes_intersection=c()
  
  if (length(paths)==1){
    paths=rep(paths,length(libnames))
  }
  if (!is.null(sample_annots_path)){
    sample_annots_tab=as.data.frame(read.csv(sample_annots_path))
    annots_tab=sample_annots_tab[libnames,]
  }
  for (i in 1:length(libnames)){
    message(libnames[i])
    if (file.exists(paste(paths[i],libnames[i],"hto.txt",sep="/"))){
      l_umitab_i=read_mtx_and_split_by_hto(paths[i],libnames[i],hto_UMI_ratio_thresh=cite_seq_params$hto_UMI_ratio_thresh,noise_interval=noise_interval,cell_interval=cell_interval,min_umi_maxhto=cite_seq_params$min_umi_maxhto,htos_to_load =strsplit(annots_tab[i,]$HTO,",")[[1]],gene_ID_converter=gene_ID_converter)
    }
    else{
      l_umitab_i=read_mtx_and_split_by_10x_batch(paths[i],libnames[i],type,noise_interval=noise_interval,cell_interval=cell_interval,gene_ID_converter=gene_ID_converter)
    }
    l=append(l,l_umitab_i)
  }
  
  for (i in 1:length(l)){
    umitab_i=l[[i]]$umitab
    numis=Matrix::colSums(umitab_i)
    noise_l[[i]]=l[[i]]$noise_counts
    
    
    if (i==1){
      genes_intersection=rownames(umitab_i)
      
    }
    else{
      genes_intersection=intersect(genes_intersection,rownames(umitab_i))
    }
    if (rename_cell_ids){
      colnames(umitab_i)=paste(i,1:ncol(umitab_i),sep="_")
    }
    
    cell_to_batch_i=rep(names(l)[i],ncol(umitab_i))
    print(names(l)[i])
    names(cell_to_batch_i)=colnames(umitab_i)
    cell_to_batch=c(cell_to_batch,cell_to_batch_i)
    l[[i]]=(umitab_i)
    
  }
  umitab=l[[1]][genes_intersection,]
  l[[1]]=l[[1]][genes_intersection,]
  
  
  noise_barcodes=list()
  noise_counts=matrix(noise_l[[1]][genes_intersection],,1)
  
  if (length(l)>1){
    for (i in 2:length(l)){
      
      l[[i]]=l[[i]][genes_intersection,]
      
      umitab=cBind(umitab,l[[i]])
      
      noise_counts=cbind(noise_counts,noise_l[[i]][genes_intersection])
    }
  }
  rownames(noise_counts)=genes_intersection
  colnames(noise_counts)=names(l)
  noise_models=t(t(noise_counts)/colSums(noise_counts))
  return(list(umitabs_by_batch=l,cell_to_batch=cell_to_batch,umitab=umitab[genes_intersection,],noise_models=noise_models,noise_barcodes=noise_barcodes))
}



human_ensmbl_to_geneSymbol_converter=function(ensembl_ids,table_path='/Users/kenige01/Documents/GitHub/scClustering/gene_id_converters/ensembl.txt'){
  if (!exists("human_ensmbl_to_geneSymbol")){
    a=read.delim(table_path,stringsAsFactors = F)
    human_ensmbl_to_geneSymbol<<-a[,1]
    names(human_ensmbl_to_geneSymbol)<<-a[,2]
  }
  return(human_ensmbl_to_geneSymbol[ensembl_ids])
}

