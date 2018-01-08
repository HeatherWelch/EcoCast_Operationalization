## Script to predict confidence intervals for all species for target data 

predCIs_master <- function(get_date,envdir,moddir,outdir,path,final_path_list){

  ############ 1. Load functions

  predCIs<-function(get_date,spname,modrep,envdir,outdir,baseGrid,modtype="pa",maskCCS=FALSE,varnames=NA){
    
    ## A. PREDICT FROM EACH MODEL OBJECT ON PREDICTION SURFACE RASTER STACK
    ## BRT
    mod_pred10 <- lapply(modrep,FUN=predict.gbm,newdata=predsurfsDF,n.trees=1000,type='response')
    mod_pred10s <- do.call(cbind,lapply(mod_pred10,data.frame,stringsAsFactors=FALSE))
    colnames(mod_pred10s) <- as.character(seq(1,10,by=1))
    ## mean prediction over 10 runs
    meanPred <- apply(mod_pred10s,1,mean)
    ## standard error
    sdPred <- apply(mod_pred10s,1,sd)
    sePred <- sdPred/sqrt(ncol(mod_pred10s))
    ## confidence intervals
    lowPred <- meanPred - 1.96*sePred
    highPred <- meanPred + 1.96*sePred
    ## make rasters 
    meanPredR <- setValues(baseGrid,meanPred)
    lowPredR <-setValues(baseGrid,lowPred)
    highPredR <-setValues(baseGrid,highPred)
    sePredR <- setValues(baseGrid,sePred)
    
    ## PLOTS
    pl <- rbind(c(-116,30), c(-122,30), c(-128,40), c(-128,48), c(-116,48), c(-115,30))
    pl <- SpatialPolygons(list(Polygons(list(Polygon(pl)), 1))); projection(pl) <- projstring
    meanplot <- mask(meanPredR,pl)
    lowplot <- mask(lowPredR,pl)
    highplot <- mask(highPredR,pl)
    seplot <- mask(sePredR,pl)
    
    ## DON'T CUT OUT CORNER OF MAP
    if (!maskCCS) {
      meanplot <- meanPredR
      lowplot <- lowPredR
      highplot <- highPredR
      seplot <- sePredR
    }
    
    # mean
    png(paste(outdir,spname,'_',modtype,'_',get_date,'_mean.png',sep=''),width=960,height=1280,units='px',pointsize=20)
    par(mar=c(3,3,.5,.5),las=1,font=2)
    if (modtype=='pa'){
      image(meanplot,col=cols,xlim=c(-128,-116),ylim=c(30,48),zlim=c(0,1)) 
    } else {
      image(meanplot,col=cols,xlim=c(-128,-116),ylim=c(30,48),zlim=c(0,15)) 
    }
    map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
    box(); dev.off()
    writeRaster(meanplot,filename=paste(outdir,spname,'_',modtype,'_',get_date,'_mean.grd',sep=''),overwrite=TRUE)
    ## low CI
    png(paste(outdir,spname,'_',modtype,'_',get_date,'_lowCI.png',sep=''),width=960,height=1280,units='px',pointsize=20)
    par(mar=c(3,3,.5,.5),las=1,font=2)
    if (modtype=='pa'){
      image(lowplot,col=cols,xlim=c(-128,-116),ylim=c(30,48),zlim=c(0,1)) 
    } else {
      image(lowplot,col=cols,xlim=c(-128,-116),ylim=c(30,48),zlim=c(0,15))  
    }
    map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
    box(); dev.off()
    writeRaster(lowplot,filename=paste(outdir,spname,'_',modtype,'_',get_date,'_lowCI.grd',sep=''),overwrite=TRUE)
    ## high CI
    png(paste(outdir,spname,'_',modtype,'_',get_date,'_highCI.png',sep=''),width=960,height=1280,units='px',pointsize=20)
    par(mar=c(3,3,.5,.5),las=1,font=2)
    if (modtype=='pa'){
      image(highplot,col=cols,xlim=c(-128,-116),ylim=c(30,48),zlim=c(0,1)) 
    } else {
      image(highplot,col=cols,xlim=c(-128,-116),ylim=c(30,48),zlim=c(0,15))
    }
    map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
    box(); dev.off()
    writeRaster(highplot,filename=paste(outdir,spname,'_',modtype,'_',get_date,'_highCI.grd',sep=''),overwrite=TRUE)
    # standard error
    png(paste(outdir,spname,'_',modtype,'_',get_date,'_se.png',sep=''),width=960,height=1280,units='px',pointsize=20)
    par(mar=c(3,3,.5,.5),las=1,font=2)
    if (modtype=='pa'){
      image(seplot,col=cols,xlim=c(-128,-116),ylim=c(30,48),zlim=c(0,0.02)) ## MAX SE PA
      # image(seplot,col=cols,xlim=c(-128,-116),ylim=c(30,48),zlim=c(0,0.01)) ## MAX SE PA
    } else {
      image(seplot,col=cols,xlim=c(-128,-116),ylim=c(30,48),zlim=c(0,0.8)) ## MAX SE DENS
      # image(seplot,col=cols,xlim=c(-128,-116),ylim=c(30,48),zlim=c(0,maxsedens)) ## MAX SE DENS
    }
    map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
    box(); dev.off()
    writeRaster(seplot,filename=paste(outdir,spname,'_',modtype,'_',get_date,'_se.grd',sep=''),overwrite=TRUE)
    print(paste(Sys.time(),get_date, 'plotted'))
    #    }
  }
  
  ############ 2. get projection date
  dates<-as.Date(get_date) ## defined in master script
  
  ############ 3. read in rasters for projection date, create dataframe
  print(paste("Creating raster stack and DF of environmental data for ",get_date,sep=""))
  fl <- final_path_list$FileList_final
  predsurfs <- stack(fl)
  
  nm=lapply(fl,function(x)unlist(strsplit(x,"/")))
  nm1=lapply(nm,function(x)x[length(x)])
  nm2=lapply(nm1,function(x)gsub(".grd","",x[[1]]))
  nm3=unlist(nm2)
  names(predsurfs)=nm3
  
  predsurfsDF <- as.data.frame(predsurfs,stringsAsFactors=FALSE)
  
  # adding in missing variables
  missing=final_path_list$FileList_missing
  missing1=unlist(lapply(missing,function(x)gsub(".grd","",x[[1]])))
  
  predsurfsDF[missing1]=NA
  
   ############ 4. rename colunms to reflect every possible model variable naming convention
    print("Renaming colunms to reflect every possible model variable naming convention")
    predsurfsDF$sla_mean<-predsurfsDF$sla
    predsurfsDF$log_eke<-predsurfsDF$l.eke_mean
    predsurfsDF$logChl<-predsurfsDF$l.blendChl
    predsurfsDF$sst<-predsurfsDF$analysed_sst
    predsurfsDF$sst_sd<-predsurfsDF$analysed_sst_sd
    predsurfsDF$sla_mean<-predsurfsDF$sla
    predsurfsDF$windy_new_mean<-predsurfsDF$ywind
    
    ############ 5. load base grid
    baseGrid <- readRDS(paste(path,'/SpatialPredictions_code/baseGrid_full.rds',sep=""))  # change for each user, or make global
    
    ############ 6. predict the confidence intervals
    CIobjs<-list.files(moddir, glob2rx('*.rds'), full.names=F)

    for(i in 1:length(CIobjs)){
      spname<-unlist(strsplit(CIobjs[[i]],"_"))[2]
      print(paste("Calculating confidence intervals for ",spname,sep=""))
	    CIdir<-paste(outdir,spname,"/predCIs/",sep='') ## change this to date. lets save all CIs in folder with same date
	    modrep<-readRDS(paste(moddir,CIobjs[i],sep=''))
	    cols <- brewer.pal(9,'GnBu')
	    projstring <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
	    modtype='pa'
	    predCIs(get_date=dates,spname=spname,modrep=modrep,envdir=envdir,outdir=CIdir,baseGrid=baseGrid) 
    }
}
    