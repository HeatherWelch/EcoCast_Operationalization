####### EcoCast script One

# Get data sequence number one : Create final and temporary envdirs, acquire all static variables
# ONLY RUN ONCE AT BEGINNING OF DAY
#currently machine accepts "local" or "lasdev"

Get_Env_Data_A <- function(machine){
  
  ############ 1. Define directories
  if(machine=="local"){
    path="~/EcoCast_Operationalization"
    envdir=paste(path,"/Satellite/",sep="") 
    source(paste0(path,"/code/load_libraries.R"),chdir = TRUE)
  }
  outdir <- paste(path,"/EcoCastRuns/",sep="")
  logdir=paste(outdir,"logs/",sep="")
  staticdir=paste0(path,"/static_variables/")
  ecocastdir=paste(outdir,"output/",sep="")
 
  ############ 2. Define time and dynamic directories
  get_date=Sys.Date()
  most_recent=as.character(Sys.Date()-1) ##change for each user
  
  tmpdir=paste(path,"/Real_time_netcdfs_raw/temp_",get_date,sep="")
  if(!file.exists(tmpdir)){
    dir.create(tmpdir) #create a temp directory to hold unfinished layers, change for each user
  }
  
  finaldir=paste(envdir,get_date,sep="");dir.create(finaldir) #change for each user
  if(!file.exists(finaldir)){
    dir.create(finaldir)
  }
  
  
  # C. Set up logfile
  templogfile = paste(logdir,"log_",get_date,".txt",sep="") 
  logfile <- file(templogfile, open="wt")
  sink(logfile, type=c("output","message")) #set all output to templog
  
  print("**************************************************************************************")
  print(paste("Starting Ecocast run for ",get_date,". Time is ",Sys.time(),sep=""))
  
  ############ 3. Define global objects
  #these are grabbed from sla_mean.grd in /EcoCast_CodeArchive/SpatialPredictions_EnvData/Satellite/2012-08-04
  template=raster() ##create template for resampling
  res(template)=0.2487562
  ncol(template)=201
  nrow(template)=201
  xmin(template)=-149.875
  xmax(template)=-99.875
  ymin(template)=10.125
  ymax(template)=60.125
  projection(template)="+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  
  ################### Acquire static variables
  print("Starting script 3_get_Env_Data_A.R")
  
  ############ 4. bathy and bathy sd
  print(paste("Copying bathymetry files into folder for ",get_date,sep=""))
  files=as.list(c("z.grd","zsd.grd","z_pt25.grd","z.gri","zsd.gri","z_pt25.gri"))
  lapply(files,function(x) file.copy(paste(staticdir,x,sep=""),finaldir))
  
  ############ 5. lunillum
  print(paste("Calculating lunillum for ",get_date,sep=""))
  value <- lunar.illumination(get_date)
  lunar_ras=template
  values(lunar_ras)=value
  writeRaster(lunar_ras,paste(finaldir,"/lunillum",sep=""),overwrite=T)
  
  ############ 6. wipe files in /latest
  mean=list.files(paste0(ecocastdir,"mean/latest/"),full.names = T)
  lapply(mean,function(x)file.remove(x))
  
  se=list.files(paste0(ecocastdir,"se/latest/"),full.names = T)
  lapply(se,function(x)file.remove(x))
  
  #warnings()
  print("**************************************************************************************")
  close(logfile)

}

Get_Env_Data_A(machine="local")



