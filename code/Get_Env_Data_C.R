# GET ENV DATA C
#Script does not check for new data
#Script evaluates the most recent version of data available, and then runs ecocast
#Script will not overwrite pre-existing final products (i.e. if final products were created by 3_Get_Env_Data_B.R)

#ONLY RUN SCRIPT ONCE AT END OF DAY. THIS IS THE ONLY SCRIPT THAT WILL CREATE A FINAL PRODUCT DESPITE MISSING DATA

## check to see if a final ecocast product exists, if not, run get_EnvData()
#currently machine accepts "local" or "lasdev"
Get_Env_Data_C=function(machine){
  
  ############ 1. Define directories
  if(machine=="local"){
    path="/home/ehazen/ecocast/EcoCast_Lasdev"
    envdir=paste(path,"/Satellite/",sep="") 
    logodir<-paste0(path,"/logo/",sep="")
    source(paste0(path,"/code/load_libraries.R"),chdir = TRUE)
    source(paste0(path,"/code/predict_CIs.R"),chdir = TRUE)
    source(paste0(path,"/code/plot_EcoCast.R"),chdir = TRUE)
  }
  outdir <- paste(path,"/EcoCastRuns/",sep="")
  logdir=paste(outdir,"logs/",sep="")
  staticdir=paste0(path,"/static_variables/")
  ecocastdir=paste(outdir,"output/",sep="")
  moddir<-paste(path,"/ModRepFiles/",sep="")
  
  ############ 2. Define time and dynamic directories
  get_date=Sys.Date()
  most_recent=as.character(get_date-1) ##change for each user
  get_date_composite=get_date-4
  
  tmpdir=paste(path,"/Real_time_netcdfs_raw/temp_",get_date,sep="")
  finaldir=paste(envdir,get_date,sep="")
  
  logfile = paste(logdir,"log_",get_date,".txt",sep="") 
  sink(logfile, type=c("output","message"),append = TRUE) #set all output to templog
  
  # C. Define species weightings
  namesrisk<-c("Blue shark bycatch","Blue sharks","Sea lions","Leatherbacks","Swordfish")
  ecocastrisk<-c(-0.1,-0.1,-0.05,-0.9,0.9) #upweight swordfish a bit
  bycatchrisk<-c(-0.1,-0.1,-0.1,-0.7,0) #all non-target species
  
  ############ 3. Define functions
  get_EnvData <- function(path,envdir,most_recent,get_date){
    
    ### F. Checks for the most recent layer for a given variable within a three-day window and returns the full path
    move_file=function(final_name){
      if(file.exists(paste(envdir,most_recent,"/",final_name,".grd",sep=""))==TRUE){ # if url isn't successfully found, start checking for older layers, but check mtime of layer to make sure it's within 3 days window
        print(paste(final_name," doesn't exist for ",get_date,", using most recent file instead (",most_recent,").",sep=""))
        file_path=paste(envdir,most_recent,"/",final_name,".grd",sep="")
      }else if (file.exists(paste(envdir,as.Date(most_recent)-1,"/",final_name,".grd",sep=""))==TRUE){
        print(paste(final_name," doesn't exist for ",get_date,", using most recent file instead (",as.Date(most_recent)-1,").",sep=""))
        file_path=paste(envdir,as.Date(most_recent)-1,"/",final_name,".grd",sep="")
      }else if (file.exists(paste(envdir,as.Date(most_recent)-2,"/",final_name,".grd",sep=""))==TRUE){
        print(paste(final_name," doesn't exist for ",get_date,", using most recent file instead (",as.Date(most_recent)-2,").",sep=""))
        file_path=paste(envdir,as.Date(most_recent)-2,"/",final_name,".grd",sep="")
      }else{
        print(paste(final_name," not available within the past three days, EcoCast run for ",get_date," will not include ",final_name,sep=""))
        file_path=NULL
      }
      return(file_path)
    }
    
    ############ 4. Define global objects
    studyarea=readOGR(dsn=staticdir,layer="sa_square_coast3")
    
    print("**************************************************************************************")
    print("********************Starting script Get_Env_Data_C.R********************")
    print("********************Finding most current versions of environmental data for get_date********************")
    
    ############ 11.Get a list of the paths of the env variables for get_date, or the most recent path if missing
    FileList_get_date=list.files(paste(envdir,get_date,sep=""),pattern="*.grd$") # all the files from get_date
    FileList_full=c("analysed_sst.grd","analysed_sst_sd.grd","l.blendChl.grd","l.eke_mean.grd","sla.grd","sla_sd.grd","ywind.grd") # all of the dynamic variables, static ones will always be there
    FileList_missing=setdiff(FileList_full,FileList_get_date) # list of dynamic variables missing from get_date
    FileList_final=list.files(paste(envdir,get_date,sep=""),pattern="*.grd$",full.names = TRUE) # start of final list to pass to preCIs script
    
    for(missing in FileList_missing){ # for each missing dynamic variable
      print(paste(missing," is missing from ",get_date,sep=""))
      final_name=gsub(".grd","",missing) # get rid of .grd to match move_file function (need the grid to differentiate between the means and sds)
      path=move_file(final_name=final_name) # get the pathway of the most recent version of each dynamic variable and..
      FileList_final=unlist(list(FileList_final,path)) # add it to our final path list
    }
    
    a=lapply(FileList_final,function(x)unlist(strsplit(x,"/")))
    available=unlist(lapply(a,function(x)x[length(x)])) # files available for get_date or get_date - 3
    FileList_missing=setdiff(FileList_full,available) # files unavailable for get_date or get_date - 3
    
    return_list=list("FileList_final"=FileList_final,"FileList_missing"=FileList_missing) # return list comprised to two lists: files we have for get date, and files we don't have for get_date
    
    return(return_list) # return the final path list to main script
  }
  
filelist=list.files(paste(path,"/EcoCastRuns/output/mean",sep=""),pattern=as.character(get_date),full.names = TRUE)

if(length(filelist)==0){
  return_list=get_EnvData(path=path,envdir=envdir,most_recent=most_recent,get_date=get_date)
  print(paste("Starting the rest of the EcoCast scripts. Time is ",Sys.time(),sep=""))
  
  print("Running function predict_CIs")
  predCIs_master(get_date=get_date,envdir = envdir,moddir=moddir,outdir = outdir,path = path,final_path_list=return_list)
  print(paste("Finished running function predict_CIs. ","Time is ",Sys.time(),sep=""))
  
  #Now run function to make the plots
  print("Running function plot_EcoCast")
  Run_ecocast(get_date=get_date,moddir=moddir,outdir = outdir,ecocastdir = ecocastdir,namesrisk=namesrisk,ecocastrisk=ecocastrisk,bycatchrisk=bycatchrisk,final_path_list=return_list,logodir=logodir,studyarea=studyarea,staticdir=staticdir)
  print(paste("Finished running function plot_EcoCast. ","Time is ",Sys.time(),sep=""))
  warnings()
  print("**************************************************************************************")
}
}

Get_Env_Data_C(machine = "local")





