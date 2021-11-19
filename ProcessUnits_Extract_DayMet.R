# This is sample code used for running stream temperature models for the NorWeST project.
# This is uses the Spokoot processing unit as an example.
# Code originally written March 2013 by SJW.
# Annotated and slightly cleaned October 2013 by SJW.
# Minor corrections made to accommodate changes to the SSN package, April 2017, SJW.

# Processor speed is the limiting factor for fitting models. 
# Memory is the limiting factor for making predictions.
# We run models using a version of R compiled with Intel MKL, which permits 
# some automatic parallel processing of the time-consuming matrix inversions.
# Microsoft OpenR is a pre-compiled version of this, though sometimes wonky. 

# Load the SSN library of functions
library(here)
library(SSN)
library(dplyr)
library(tidyr)
library(daymetr)
library(foreach)
library(doParallel)


# Function to standardize variables
stand <- function(x) { (x-mean(x))/(2*sd(x))}

# Function to standardize a prediction dataset based on the fitting dataset 
stdpreds <- function(newset,originalset) {
  xnames <- colnames(newset)
  sx <- matrix(rep(NA,ncol(newset)*nrow(newset)),nrow=nrow(newset))
  for(i in 1:ncol(newset)) {
    var <- with(originalset,get(xnames[i]))
    sx[,i] <- (newset[,i]-mean(var))/(2*sd(var))
  }
  colnames(sx) <- colnames(newset)
  return(sx)
}

# Function to get fixed effects and SEs from glmssn. Used to make tables of outputs. 
ests <- function(x) {
  means <- round(x$estimates$betahat,3)
  ses <- round(sqrt(diag(x$estimates$covb)),3)
  output <- cbind(means,ses)
  colnames(output) <- c("Estimate","SE")
  return(output)
}

# Set the working directory to the parent folder of the SSN.

# Load the ssn and all sets of prediction points
# Points were divided into three groups to prevent memory errors in ArcGIS.
# This wasnâ€™t necessary for later production units.
AugDayMet <- function(Lat,Lon,ID_1KM){
  # Lat <- LatLonsObs$Lat[x]
  # Lon <- LatLonsObs$Lon[x]
  # locID <- LatLonsObs$locID[x]
  
  df1 <- download_daymet(site = 'DayMet', lat=Lat, lon=Lon, start=1980, end=2020, internal=T, simplify = T) %>%
    filter(measurement %in% c('tmin..deg.c.','tmax..deg.c.')) %>%
    mutate(date=as.Date(paste(year,yday,sep='-'), '%Y-%j')) %>%
    mutate(month=format(date,'%m')) %>% 
    select(year,yday,month,measurement,value) %>% 
    pivot_wider(names_from=measurement, values_from=value) %>% 
    mutate('tmean..deg.c.'= tmin..deg.c. + (tmax..deg.c.-tmin..deg.c.)/2)
  
  AugTemp <- df1 %>% filter(month == '08') %>% group_by(year) %>% summarize(MeanAug=mean(tmean..deg.c.),MeanMaxAug=mean(tmax..deg.c.))
  AugTemp$ID_1KM <-  ID_1KM
    
  return(AugTemp)
}

GetDayMet <<- function(Lat,Lon,ID_1KM){
  
  #Download All Daymet measurements for this site
  df1 <- daymetr::download_daymet(site = 'DayMet', lat=Lat, lon=Lon, start=1980, end=2020, internal=T, simplify = T) %>%
    mutate(date=as.Date(paste(year,yday,sep='-'), '%Y-%j')) %>%
    mutate(month=format(date,'%m')) %>% 
    select(year,yday,month,measurement,value) %>% 
    pivot_wider(names_from=measurement, values_from=value)
  
  #Add id column to link to DayMet
  df1$ID_1KM <- ID_1KM
  return(df1)
}

Units <- c('Clearwater.ssn','Midsnake.ssn','MissouriHW.ssn','Salmon.ssn','SnakeBear.ssn','Spokoot.ssn','UpMissMarias.ssn','UpYellBighorn.ssn')

for(u in c(1)){ #1,2,3,NOT4,5,6,7
  u<-8
  UnitName <- substr(Units[u],1,nchar(Units[u])-4)
  
  if(Units[u]=='Midsnake.ssn'){
    UnitIn <- importSSN(here('Regions',Units[u]),predpts='pred_ne')
    UnitIn <- importPredpts(UnitIn,'pred_se','ssn')
    UnitIn <- importPredpts(UnitIn,'pred_we','ssn')
    
    if(!file.exists(here('Regions',Units[u],'distance'))){
      createDistMat(UnitIn,o.write=T,predpts="pred_ne",amongpreds=T)
      createDistMat(UnitIn,o.write=T,predpts="pred_se",amongpreds=T) 
      createDistMat(UnitIn,o.write=T,predpts="pred_we",amongpreds=T)
    }
    
  }else if(Units[u]=='SnakeBear.ssn'){
    UnitIn <- importSSN(here('Regions',Units[u]),predpts='prednorth')
    UnitIn <- importPredpts(UnitIn,'predsouth','ssn')
    
    if(!file.exists(here('Regions',Units[u],'distance'))){
      createDistMat(UnitIn,o.write=T,predpts="prednorth",amongpreds=T)
      createDistMat(UnitIn,o.write=T,predpts="predsouth",amongpreds=T) 
    }
    
  }else if(Units[u]=='Spokoot.ssn'){
    UnitIn <- importSSN(here('Regions',Units[u]),predpts='prednorth')
    UnitIn <- importPredpts(UnitIn,'predse','ssn')
    UnitIn <- importPredpts(UnitIn,'predsw','ssn')
    
    if(!file.exists(here('Regions',Units[u],'distance'))){
      createDistMat(UnitIn,o.write=T,predpts="prednorth",amongpreds=T)
      createDistMat(UnitIn,o.write=T,predpts="predse",amongpreds=T) 
      createDistMat(UnitIn,o.write=T,predpts="predsw",amongpreds=T) 
    }
    
  }else{
    UnitIn <- importSSN(here('Regions',Units[u]),predpts='preds')
    
    if(!file.exists(here('Regions',Units[u],'distance'))){
      createDistMat(UnitIn,o.write=T,predpts="preds",amongpreds=T)
    }
  
  }
  proj4string<-proj4string(UnitIn)
  #UnitIn@obspoints@SSNPoints[[1]]
  
  #Pull DayMet predictions for all observation locations
  LatLons_Obs <- apply(UnitIn@obspoints@SSNPoints[[1]]@point.coords,1,FUN=function(x){
    p1<-proj4::project(x,proj4string,inverse=TRUE)
  }) %>% t() %>% data.frame()
  colnames(LatLons_Obs) <- c('Lon','Lat')
  LatLons_Obs$ID_1KM <- UnitIn@obspoints@SSNPoints[[1]]@point.data$ID_1KM #as.numeric(levels(UnitIn@obspoints@SSNPoints[[1]]@point.data$ID_1KM)[UnitIn@obspoints@SSNPoints[[1]]@point.data$ID_1KM])
  LatLons_Obs <- LatLons_Obs %>% distinct() %>% select(ID_1KM, Lat, Lon)
  #write.table(LatLons_Obs,file=paste0('Regions/LatLons_Obs_',UnitName,'.csv'),row.names=F,sep=',')
  
  # AugTemp_Obs<-future_lapply(1:nrow(LatLons_Obs),FUN=function(x){
  #   #x <- 1
  #   otp<-AugDayMet(Lon=LatLons_Obs$Lon[x],Lat=LatLons_Obs$Lat[x],ID_1KM = LatLons_Obs$ID_1KM[x])
  #   return(otp)  
  # }) %>% bind_rows()
  
  DayMet_Obs<-lapply(1:nrow(LatLons_Obs),FUN=function(x){
    #x <- 1
    otp<-GetDayMet(Lon=LatLons_Obs$Lon[x],Lat=LatLons_Obs$Lat[x],ID_1KM = LatLons_Obs$ID_1KM[x])
    return(otp)  
  }) %>% bind_rows()
  
  #save(DayMet_Obs,file=paste0('Regions/DayMet/DayMet_Obs',UnitName,'.Rdata'))
  saveRDS(DayMet_Obs,file=here('Regions','DayMet',UnitName,paste0(UnitName,'.Obs','.RDS')))
  
  #Pull DayMet predictions for all prediction sites
  LatLons_Pred <- apply(UnitIn@predpoints@SSNPoints[[1]]@point.coords,1,FUN=function(x){
    p1<-proj4::project(x,proj4string,inverse=TRUE)
  }) %>% t() %>% data.frame()
  colnames(LatLons_Pred) <- c('Lon','Lat')
  LatLons_Pred$ID_1KM <- UnitIn@predpoints@SSNPoints[[1]]@point.data$ID_1KM #as.numeric(levels(UnitIn@predpoints@SSNPoints[[1]]@point.data$ID_1KM)[UnitIn@predpoints@SSNPoints[[1]]@point.data$ID_1KM])
  LatLons_Pred <- LatLons_Pred %>% distinct() %>% select(ID_1KM, Lat, Lon)
  
  #Break Into Groups of 1000
  nGroups <- ceiling(nrow(LatLons_Pred)/1000)
  GroupIndices <- list()
  for(i in 1:nGroups){
    startInd<- 1+1000*(i-1)
    if(i!=nGroups){
      endInd<-i*1000
    }else{
      endInd<-nrow(LatLons_Pred)
    }
    GroupIndices[[i]] <- seq(startInd,endInd)
  }
  
  Start.time<-Sys.time()
  
  #GI<-GroupIndices[1:2]
  #lapply(GroupIndices,FUN=function(x){
  for(i in 94:nGroups){
    x <- GroupIndices[[i]]
    
    # DayMetUnit <- mclapply(xL,FUN=function(y,LatLons=LatLons_Pred){
    #   #y <- xL[[1]]
    #   #LatLons <- LatLons_Pred
    #   otp<-tryCatch({GetDayMet(Lon=LatLons$Lon[y],Lat=LatLons$Lat[y],ID_1KM = LatLons$ID_1KM[y])})
    #   return(otp)  
    # },mc.cores=4) %>% bind_rows()
    
    #Parallel run daymet extractions
    cl <- makeCluster(12)
    registerDoParallel(cl)
    
    DayMetUnit <- foreach(y=x,.packages=c('dplyr','daymetr','tidyr')) %dopar% {
      #y <- xL[[1]]
      otp<-tryCatch({GetDayMet(Lon=LatLons_Pred$Lon[y],Lat=LatLons_Pred$Lat[y],ID_1KM = LatLons_Pred$ID_1KM[y])})
      return(otp)
    } %>% bind_rows()
    
    stopCluster(cl)
    
    saveRDS(DayMetUnit,file=here('Regions','DayMet',UnitName,'PredsBy1000',paste(UnitName,x[1],'RDS',sep='.')))
    print(i)
    print(round(i/nGroups,2))
  }
  
  
  difftime(Sys.time(),Start.time)
  
}
#   
#   UnitIn@predpoints@SSNPoints[[1]]@point.data
#   
#   #Draw Region Specific Torgegram
#   unit_tg <- Torgegram(UnitIn,"STREAM_AUG",nlag=20)
#   jpeg(paste0('Regions/',Units[u],"_rawtorg.jpg"))
#     plot(unit_tg, main = "Raw Data Torgegram") 
#   dev.off()
#   
#   #Get dataframe
#   unit_df <- getSSNdata.frame(UnitIn)
#   
#   unit_df$GLACIER <- as.numeric(unit_df$GLACIER!=-9999)
#   
#   #Standardize covariates
#   continuous <- unit_df %>% select(ELEV,CANOPY,SLOPE,PRECIP,CUMDRAINAG,Y_COORD,NLCD11PC,BFI,Air_Aug,Flow_Aug)
#   cont.s <- apply(continuous,2,stand)
#   colnames(cont.s) <- c("elev","canopy","slope","precip","drainage","lat","water","bfi","airtemp","flow")
#   unit_df.s <- data.frame(unit_df,cont.s)
#   unit_df.s$yearf <- factor(unit_df.s$SAMPLEYEAR)
#   UnitIn_s <- putSSNdata.frame(unit_df.s,UnitIn,"Obs")
#   
#   GlacierBool <- sum(unit_df.s$GLACIER)>0
#   TailwaterBool <- sum(unit_df.s$TAILWATER)>0
#     
#   # A-spatial models with
#   if(GlacierBool & !(TailwaterBool)){
#     unit.aspatial <- lm(STREAM_AUG ~ elev + canopy + slope + precip + GLACIER + drainage + lat + water + bfi + airtemp + flow, data=unit_df.s)
#   }else if(!GlacierBool & (TailwaterBool)){
#     unit.aspatial <- lm(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + drainage + lat + water + bfi + airtemp + flow, data=unit_df.s)
#   }else if(GlacierBool & TailwaterBool){
#     unit.aspatial <- lm(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + GLACIER + drainage + lat + water + bfi + TAILWATER + airtemp + flow, data=unit_df.s)
#   }else{
#     unit.aspatial <- lm(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi  + airtemp + flow, data=unit_df.s)
#   }
#   
#   predictA <- predict(unit.aspatial)
#   sqrt(mean((predictA-unit_df.s$STREAM_AUG)^2))
#   
#   if(GlacierBool & !(TailwaterBool)){ #If glaciers but not tailwaters
#     unit.tu <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.tailup"), addfunccol = "afvArea"))
#     unit.td <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.taildown"), addfunccol = "afvArea"))
#     unit.tu.eu <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.tailup","Exponential.Euclid"), addfunccol = "afvArea"))
#     unit.td.eu <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea"))
#     
#   }else if(!GlacierBool & (TailwaterBool)){ #If tailwaters but not glaciers
#     unit.tu <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.tailup"), addfunccol = "afvArea"))
#     unit.td <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.taildown"), addfunccol = "afvArea"))
#     unit.tu.eu <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.tailup","Exponential.Euclid"), addfunccol = "afvArea"))
#     unit.td.eu <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea"))
#     
#   }else if(GlacierBool & TailwaterBool){ #If tailwaters and glaciers
#     unit.tu <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.tailup"), addfunccol = "afvArea"))
#     unit.td <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.taildown"), addfunccol = "afvArea"))
#     unit.tu.eu <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.tailup","Exponential.Euclid"), addfunccol = "afvArea"))
#     unit.td.eu <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea"))
#     
#   }else{ #If no tailwaters and no glaciers
#     unit.tu <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.tailup"), addfunccol = "afvArea"))
#     unit.td <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.taildown"), addfunccol = "afvArea"))
#     unit.tu.eu <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.tailup","Exponential.Euclid"), addfunccol = "afvArea"))
#     unit.td.eu <- try(glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID","yearf","Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea"))
#   }
#   
#   AICcomp <- c(AIC(unit.tu),
#                AIC(unit.td),
#                AIC(unit.tu.eu),
#                AIC(unit.td.eu))
#   BestMod<-switch(which.min(AICcomp),unit.tu,unit.td,unit.tu.eu,unit.td.eu)
#   BestMod <- update(unit.tu[[1]],EstMeth='ML')
#   
#   assign(paste0('BestMod_',substr(Units[u],1,nchar(Units[u])-4)),BestMod)
#   save(get(paste0('BestMod_',substr(Units[u],1,nchar(Units[u])-4))),paste0('Regions/BestModel_',substr(Units[u],1,nchar(Units[u])-4)))
#   save(ls(),paste0('Regions/Workspace_',substr(Unit[u],1,nchar(Units[u])-4)))
#   
#   BestMod_r <- residuals(BestMod, cross.validation=T)
#   BestMod_rdf <- getSSNdata.frame(BestMod_r)
#   
#   # Torgegram of fitted model residuals
#   BestMod_t <- Torgegram(BestMod_r,"_resid.crossv_",nlag=20)
#   jpeg(paste0('Regions/BestMod_',substr(Units[u],1,nchar(Units[u])-4),"_residtorg.jpg"))
#     plot(BestMod_t, main= "BestModel Residuals Torgegram") 
#   dev.off()
#   
#   #Root mean squared error (RMSE) of cross-validated predictions
#   sqrt(mean((BestMod_rdf[,"_CrossValPred_"]-BestMod_rdf$obsval)^2))
#   #RMSE of fixed effects only
#   sqrt(mean((BestMod_rdf[,"_fit_"]-BestMod_rdf$obsval)^2))
#   # 1.935
#   # Null RMSE
#   sqrt(mean((BestMod_rdf$obsval-mean(BestMod_rdf$obsval))^2))
#   # 2.996
#   
#   #Pseudo-r2 of cross-validated predictions. 
#   cor(BestMod_rdf$obsval,BestMod_rdf[,"_CrossValPred_"])^2
# }
# 
# 
# 
# msh <- importSSN("MissouriHW.ssn",predpts="preds")
# msh@data
# #msh1 <- importPredpts(msh,"predse","ssn")
# #msh <- importPredpts(msh,"predsw","ssn")
# 
# # Create distance matrices
# createDistMat(msh,o.write=T,predpts="preds",amongpreds=T)
# #createDistMat(msh,o.write=T,predpts="predse",amongpreds=T) 
# #createDistMat(msh,o.write=T,predpts="predsw",amongpreds=T)
# 
# 
# # Create raw torgegram
#   msh_tg <- Torgegram(msh,"STREAM_AUG",nlag=20)
#   jpeg("msh_rawtorg.jpg")
#   plot(msh_tg, main = "Raw Data Torgegram") 
# dev.off()
# 
# # Extract dataframe and fit basic aspatial model with un-standardized predictors
# mshdf <- getSSNdata.frame(msh)
# msh.lm <- lm(STREAM_AUG ~ ELEV + CANOPY + SLOPE + PRECIP + CUMDRAINAG + Y_COORD + NLCD11PC + GLACIER + BFI + TAILWATER + Air_Aug + Flow_Aug, data=mshdf)
# summary(msh.lm)
# 
# # Standardize continuous covariates and add factor for year
# continuous <- mshdf[,c(11:21)]
# cont.s <- apply(continuous,2,stand)
# colnames(cont.s) <- c("elev","canopy","slope","precip","drainage","lat","water","glacier","bfi","airtemp","flow")
# mshdf.s <- data.frame(mshdf,cont.s)
# mshdf.s$yearf <- factor(mshdf.s$SAMPLEYEAR)
# mshs <- putSSNdata.frame(mshdf.s,msh,"Obs")
# 
# # Version without glacier + standardized variables
# msh.lm2 <- lm(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + TAILWATER + airtemp + flow, data=mshdf.s)
# 
# # Extract pure aspatial parameter estimates, back-transform, and save both versions of estimates
# # Must be careful about the order of parameters! No attempt is made to match them by name.
# # This WILL need to be adjusted depending if glaciers are or are not used in the model.
# # Ditto for tailwater.
# 
# mshAe <- summary(msh.lm2)$coefficients
# backtrans <- mshAe[-c(1,10),1:2]/(2*sapply(continuous[,-8],sd))
# esttable <- cbind(rbind(mshAe[1,1:2],backtrans[1:8,],mshAe[10,1:2],backtrans[9:10,]),mshAe)
# rownames(esttable) <- rownames(spokAe)
# write.csv(esttable,"aspatialestimates.csv")
# 
# # Aspatial performance
# predictA <- predict(msh.lm2)
# sqrt(mean((predictA-mshdf.s$STREAM_AUG)^2))
# # 1.885
# 
# # Examine correlations
# library(ellipse)
# cor(cont.s)
# jpeg("corrplot.jpg")
# plotcorr(cor(cont.s),type="lower")
# dev.off()
# 
# # Get VIFs from linear model as an indicator of multicollinearity
# library(car)
# vif(msh.lm2)
# #      ELEV     CANOPY      SLOPE     PRECIP CUMDRAINAG    Y_COORD   NLCD11PC 
# #  2.399760   1.584861   1.351117   1.577125   1.909954   2.585873   1.305458 
# #   GLACIER        BFI  TAILWATER    Air_Aug   Flow_Aug 
# #  1.048601   1.255432   1.323322   1.048707   1.051820 
# 
# # Standardize preds based on obs and add factor for year- prednorth
# mshpreddf <- getSSNdata.frame(msh, "preds")
# contpred <- mshpreddf[,c(11:21)]
# contpred.s <- stdpreds(contpred,continuous)
# colnames(contpred.s) <- c("elev","canopy","slope","precip","drainage","lat","water","glacier","bfi","airtemp","flow")
# mshpreddf.s <- data.frame(mshpreddf,contpred.s)
# mshpreddf.s$yearf <- factor(mshpreddf.s$SAMPLEYEAR)
# mshs <- putSSNdata.frame(mshpreddf.s,msh,"preds")
# 
# 
# ###########
# # Model1  #
# ###########
# 
# # We run all models with the same predictor variables, whether significant or not, except:
# # Glacier and Tailwater are omitted as a covariates in units where they are absent 
# # In this unit we use tailwater but not glacier.
# 
# # Note that the correlation model structure was developed after testing of the Salmon and
# # Clearwater production units. Once the model structure was finalized, we decided
# # to run all final models using EstMeth= â€œMLâ€.
# 
# starttime <- Sys.time() # not necessary; just for timing longer runs
# 
# unique(mshs@obspoints@SSNPoints[[1]]@point.data$locID)
# 
# #msh1 <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + TAILWATER + airtemp + flow, mshs, EstMeth= "ML", family="gaussian",CorModels = c("locID","yearf","Exponential.tailup","Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea")
# 
# #get errors for putting tailup and taildown in the same model
# msh.tu <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + TAILWATER + airtemp + flow, mshs, EstMeth= "ML", family="gaussian",CorModels = c("locID","Exponential.tailup"), addfunccol = "afvArea")
# msh.td <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + TAILWATER + airtemp + flow, mshs, EstMeth= "ML", family="gaussian",CorModels = c("locID","yearf","Exponential.taildown"), addfunccol = "afvArea")
# msh.tu.eu <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + TAILWATER + airtemp + flow, mshs, EstMeth= "ML", family="gaussian",CorModels = c("locID","yearf","Exponential.tailup","Exponential.Euclid"), addfunccol = "afvArea")
# #msh.td.eu <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + TAILWATER + airtemp + flow, mshs, EstMeth= "ML", family="gaussian",CorModels = c("Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea") 
# 
# #msh.tu.td.eu <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + TAILWATER + airtemp + flow, mshs, EstMeth= "ML", family="gaussian",CorModels = c("locID","yearf","Exponential.tailup","Exponential.Euclid",'Exponential.taildown'), addfunccol = "afvArea")
# 
# 
# AIC(msh.tu)
# AIC(msh.td)
# AIC(msh.tu.eu)
# 
# msh1 <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + TAILWATER + airtemp + flow, mshs, EstMeth= "ML", family="gaussian",CorModels = c("Exponential.tailup","Exponential.Euclid"), addfunccol = "afvArea")
# 
# elapsed <- Sys.time()-starttime # See above note
# elapsed
# 
# # Extract predictions/ Leave-one-out cross validation predictions
# msh1r <- residuals(msh1, cross.validation=T)
# msh1rrdf <- getSSNdata.frame(msh1r)
# 
# # Torgegram of fitted model residuals
# msh1t <- Torgegram(msh1r,"_resid.crossv_",nlag=20)
# jpeg("msh1residtorg.jpg")
# plot(msh1t, main= "Model1 Residuals Torgegram") 
# dev.off()
# 
# #Root mean squared error (RMSE) of cross-validated predictions
# sqrt(mean((spok1rdf[,"_CrossValPred_"]-spok1rdf$obsval)^2))
# # 0.969
# 
# #RMSE of fixed effects only
# sqrt(mean((spok1rdf[,"_fit_"]-spok1rdf$obsval)^2))
# # 1.935
# 
# # Null RMSE
# sqrt(mean((spok1rdf$obsval-mean(spok1rdf$obsval))^2))
# # 2.996
# 
# #Pseudo-r2 of cross-validated predictions. 
# cor(spok1rdf$obsval,spok1rdf[,"_CrossValPred_"])^2
# # 0.895
# 
# # Get parameter estimates, back-transform, and save both versions of estimates
# # Must be careful about the order of parameters! No attempt is made to match them by name.
# # This WILL need to be adjusted depending if glaciers are or are not used in the model.
# # Ditto for tailwater.
# 
# spok1e <- summary(spok1)$fixed.effects.estimates
# backtrans <- spok1e[-c(1,10),2:3]/(2*sapply(continuous[,c(1:7,9:11)],sd))
# esttable <- cbind(spok1e[,1],rbind(spok1e[1,2:3],backtrans[1:8,],spok1e[10,2:3],backtrans[9:10,]),spok1e[,2:5])
# write.csv(esttable,"spok1estimates.csv")
# 
# # Outlier analysis
# resids <- spok1rdf[,"_resid.student_"]
# outlie <- spok1rdf[abs(resids)>5,c("pid","obsval","_CrossValPred_","_resid.student_")]
# # 18 records have abs studentized residuals >5; 5 records are >7; 1 record >10. 
# 
# 
# ###############
# # Predictions #
# ###############
# 
# 
# # Split prediction datasets into pieces for computers with lowish RAM. 
# # The "chunksof" setting depends on available memory on the computer.
# # 10-12K is about the limit on a machine with 16GB RAM.
# # Used chunks of 7000 for predse and predsw just to make the two batches even, but 10K would have been OK too.
# # Many datasets run without splitting on machines with 32GB RAM; 
# # I think all will run without splitting on machines with 64GB RAM.
# 
# spok1$ssn.object <- splitPredictions(spok1$ssn.object, "prednorth",chunksof=10000)
# spok1$ssn.object <- splitPredictions(spok1$ssn.object, "predse",chunksof=7000)
# spok1$ssn.object <- splitPredictions(spok1$ssn.object, "predsw",chunksof=7000)
# 
# # Make predictions
# spok1p1 <- predict(spok1,"prednorth-1")
# spok1p2 <- predict(spok1,"prednorth-2")
# spok1p3 <- predict(spok1,"prednorth-3")
# spok1p4 <- predict(spok1,"predse-1")
# spok1p5 <- predict(spok1,"predse-2")
# spok1p6 <- predict(spok1,"predsw-1")
# spok1p7 <- predict(spok1,"predsw-2")
# 
# # Extrac prediction data frames
# pred1df <- getSSNdata.frame(spok1p1,"prednorth-1")
# pred2df <- getSSNdata.frame(spok1p2,"prednorth-2")
# pred3df <- getSSNdata.frame(spok1p3,"prednorth-3")
# pred4df <- getSSNdata.frame(spok1p4,"predse-1")
# pred5df <- getSSNdata.frame(spok1p5,"predse-2")
# pred6df <- getSSNdata.frame(spok1p6,"predsw-1")
# pred7df <- getSSNdata.frame(spok1p7,"predsw-2")
# 
# # Reassemble the pieces into one batch.
# allpreds <- rbind(pred1df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")],pred2df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")],pred3df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")],pred4df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")],pred5df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")],pred6df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")],pred7df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")])
# colnames(allpreds) <- c("OBSPREDID","predtemp","predtempse")
# 
# # Export prediction dataset as a csv (good for general use) and dbf (good for GIS)
# write.csv(allpreds,"spok1preds.csv",row.names=F)
# write.dbf(allpreds,"spok1preds.dbf")
# 
# # Compare predicted and observed at fitting sites, and export. 
# predobs <- data.frame(spok1rdf$OBSPREDID,spok1rdf[,"_CrossValPred_"],spok1rdf$obsval)
# colnames(predobs) <- c("obspredid","predicted","observed")
# write.csv(predobs,"predobs.csv",row.names=F)
# 
# 
# ###########
# # Model2  #
# ###########
# 
# # â€œAspatialâ€ model. Same predictors as spatial model. Includes random effects for site and year. 
# 
# starttime <- Sys.time()
# spok2 <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + TAILWATER + airtemp + flow, spoks, EstMeth= "ML", family="gaussian",CorModels = c("locID","yearf"), addfunccol = "afvArea")
# elapsed <- Sys.time()-starttime
# elapsed
# 
# # Extract predictions/ LOO CV predictions
# spok2r <- residuals(spok2,cross.validation=T)
# spok2rdf <- getSSNdata.frame(spok2r)
# 
# # Create Torgegram of residuals
# spok2t <- Torgegram(spok2r,"_resid.crossv_",nlag=20)
# jpeg("spok2residtorg.jpg")
# plot(spok2t, main= "Model2 Residuals Torgegram") 
# dev.off()
# 
# #RMSPE of cross-validated predictions
# sqrt(mean((spok2rdf[,"_CrossValPred_"]-spok2rdf$obsval)^2))
# # 1.176
# 
# #RMSPE of fixed effects only
# sqrt(mean((spok2rdf[,"_fit_"]-spok2rdf$obsval)^2))
# # 1.890
# 
# # Null RMSPE
# sqrt(mean((spok2rdf$obsval-mean(spok2rdf$obsval))^2))
# # 2.996
# 
# #r2 of cross-validated predictions. 
# cor(spok2rdf$obsval,spok2rdf[,"_CrossValPred_"])^2
# # r2 = 0.846
# 
# # Get parameter estimates, back-transform, and save both versions of estimates
# spok2e <- summary(spok2)$fixed.effects.estimates
# backtrans <- spok2e[-c(1,10),2:3]/(2*sapply(continuous[,c(1:7,9:11)],sd))
# esttable <- cbind(spok2e[,1],rbind(spok2e[1,2:3],backtrans[1:8,],spok2e[10,2:3],backtrans[9:10,]),spok2e[,2:5])
# write.csv(esttable,"spok2estimates.csv")
# 
# # Save the full image including all modeling results. Optional.
# save.image(file = â€œspok.RData")
# 
# 
