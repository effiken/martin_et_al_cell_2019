make_seeding_job_file=function(job_name="scClustering",account= "acc_Chessa01a",queue="premium",ncores=1,memk=8000,wall_time="01:00",headnode="minerva",log_prefix="seed_",path="",model_name=""){
  v=c(
    "#!/bin/bash",
    paste("#BSUB -J", job_name),
    paste("#BSUB -P", account),
    paste("#BSUB -q", queue),
    paste("#BSUB -n", ncores),
    paste("#BSUB -R rusage[mem=",memk,"]",sep=""),
    paste("#BSUB -W", wall_time),
#    paste("#BSUB -m", headnode),
    paste("#BSUB -o log/",log_prefix,"out_%J.%I",sep=""),
    paste("#BSUB -eo log/",log_prefix,"err_%J.%I",sep=""),
    "#BSUB -L /bin/bash",                                                                                                                                                                                                          
    "",
    "echo TASK_ID=$LSB_JOBINDEX",
    "module load R",
    paste("Rscript scClustering/run_get_seed.r",model_name,"$LSB_JOBINDEX")
  )

  writeLines(con=paste(path,"job_",job_name,".sh",sep=""),v)
}
