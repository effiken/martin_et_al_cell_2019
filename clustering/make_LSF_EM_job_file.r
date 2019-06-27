make_EM_job_file=function(job_name="scClustering",account= "acc_Chessa01a",queue="premium",ncores=1,memk=8000,wall_time="01:00",headnode="minerva",log_prefix="seed_",path="",data_l_path,model_name,k){
  v=c(
    "#!/bin/bash",
    paste("#BSUB -J", job_name),
    paste("#BSUB -P", account),
    paste("#BSUB -q", queue),
    paste("#BSUB -n", ncores),
    paste("#BSUB -R rusage[mem=",memk,"]",sep=""),
    paste("#BSUB -W", wall_time),
#    paste("#BSUB -m", headnode),
    paste("#BSUB -o log/",log_prefix,"out_%J",sep=""),
    paste("#BSUB -eo log/",log_prefix,"err_%J",sep=""),
    "#BSUB -L /bin/bash",                                                                                                                                                                                                          
    "",
    "module load R",
    paste("Rscript scClustering/run_EM.r",data_l_path,model_name,k)
  )

  writeLines(con=paste(path,"job_",job_name,".sh",sep=""),v)
}
