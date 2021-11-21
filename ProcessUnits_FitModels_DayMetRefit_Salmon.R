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
library(foreign)
library(dplyr)
library(tidyr)
library(doParallel)
library(foreach)
library(lme4)

data_dir <- '../Annual_NorWest_Data/'

#data_dir <- '/Volumes/Cline_USGS/Annual_NorWest_Data/'

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


Units <- c('Clearwater.ssn','Midsnake.ssn','MissouriHW.ssn','Salmon.ssn','SnakeBear.ssn','Spokoot.ssn','UpMissMarias.ssn','UpYellBighorn.ssn')
  u<-4
  UnitName <- substr(Units[u],1,nchar(Units[u])-4)

  UnitIn<-readRDS(paste0(data_dir,paste0('DayMetSSN/',UnitName,'.DayMetSSN.RDS')))
  UnitIn@path <- paste0(data_dir,paste0('Regions/',UnitName,'.ssn'))
  proj4string<-proj4string(UnitIn)
  
  #UnitIn@obspoints@SSNPoints[[1]]
  
  #Get dataframe
  unit_df <- getSSNdata.frame(UnitIn)
  
  if('GLACIER' %in% colnames(unit_df)){
    unit_df$GLACIER[unit_df$GLACIER==-9999] <- 0
  }
  
  #Standardize covariates
  continuous <- unit_df %>% select(ELEV,CANOPY,SLOPE,PRECIP,CUMDRAINAG,Y_COORD,NLCD11PC,BFI,Air_Aug,Flow_Aug)
  cont.s <- apply(continuous,2,stand)
  colnames(cont.s) <- c("elev","canopy","slope","precip","drainage","lat","water","bfi","airtemp","flow")
  unit_df.s <- data.frame(unit_df,cont.s)
  unit_df.s$yearf <- factor(unit_df.s$SAMPLEYEAR)
  
  
  GlacierBool <- ifelse('GLACIER' %in% colnames(unit_df.s) & sum(unit_df.s$GLACIER)>0, 1, 0)
  TailwaterBool <- ifelse('TAILWATER' %in% colnames(unit_df.s) & sum(unit_df.s$TAILWATER)>0,1,0)
  
  #if('TAILWATER' %in% colnames(unit_df.s)){unit_df.s$TAILWATER <- as.factor(unit_df.s$TAILWATER)}
  
  unit_df.s$locID2 <- as.factor(unit_df.s$ID_1KM)
  
  # A-spatial models with
  if(GlacierBool & !(TailwaterBool)){
    unit.aspatial <- lmer(STREAM_AUG ~ elev + canopy + slope + precip + GLACIER + drainage + lat + water + bfi + airtemp + flow + (1|locID2) + (1|yearf), data=unit_df.s)
  }else if(!GlacierBool & (TailwaterBool)){
    unit.aspatial <- lmer(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + drainage + lat + water + bfi + airtemp + flow + (1|locID2) + (1|yearf), data=unit_df.s)
  }else if(GlacierBool & TailwaterBool){
    unit.aspatial <- lmer(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + GLACIER + drainage + lat + water + bfi + airtemp + flow + (1|locID2) + (1|yearf), data=unit_df.s)
  }else{
    unit.aspatial <- lmer(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi  + airtemp + flow + (1|locID2) + (1|yearf), data=unit_df.s)
  }
  
  summary(unit.aspatial)
  predictA <- predict(unit.aspatial)
  sqrt(mean((predictA-unit_df.s$STREAM_AUG)^2))
  
  UnitIn_s <- putSSNdata.frame(unit_df.s,UnitIn,"Obs")
  #options(mc.cores=parallel::detectCores())
  #unit_df.s$date <- paste(unit_df.s$SAMPLEYEAR,'08','15',sep='-')
  #ssnb1<-ssnbayes(STREAM_AUG ~ airtemp,data=unit_df.s,path=here('Regions',Units[u]),time_method = list('ar','date'),space_method = list('use_ssn',c('Exponential.taildown')),iter=1000,warmup=100,chains=3,net=2,addfunccol='afvArea')
  # path <- system.file("extdata/clearwater.ssn", package = "SSNbayes")
  # n <- importSSN(path, predpts = "preds", o.write = TRUE)
  # clear <- readRDS(system.file("extdata/clear_obs.RDS", package = "SSNbayes"))
  #ssn1<-glmssn(STREAM_AUG ~ airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("Exponential.tailup"), addfunccol = "afvArea")
  
  if(GlacierBool & !(TailwaterBool)){ #If glaciers but not tailwaters
    # unit.tu <- tryCatch({glmssn(STREAM_AUG ~ elev + canopy + slope + precip + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.tailup"), addfunccol = "afvArea")},error=function(e){return(NA)})
    # unit.td <- tryCatch({glmssn(STREAM_AUG ~ elev + canopy + slope + precip + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.taildown"), addfunccol = "afvArea")},error=function(e){return(NA)})
    # unit.tu.eu <- tryCatch({glmssn(STREAM_AUG ~ elev + canopy + slope + precip + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.tailup","Exponential.Euclid"), addfunccol = "afvArea")},error=function(e){return(NA)})
    # unit.td.eu <- tryCatch({glmssn(STREAM_AUG ~ elev + canopy + slope + precip + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea")},error=function(e){return(NA)})
    # 
    unit.tu.td.eu <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + GLACIER + drainage + lat + water + bfi + airtemp + flow, 
                            UnitIn_s, EstMeth= "ML", family="gaussian",CorModels = c("locID2","yearf","Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea")
  }else if(!GlacierBool & (TailwaterBool)){ #If tailwaters but not glaciers
    # unit.tu <- tryCatch({glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.tailup"), addfunccol = "afvArea")},error=function(e){message(e);return(NA)})
    # unit.td <- tryCatch({glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.taildown"), addfunccol = "afvArea")},error=function(e){message(e);return(NA)})
    # unit.tu.eu <- tryCatch({glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.tailup","Exponential.Euclid"), addfunccol = "afvArea")},error=function(e){message(e);return(NA)})
    # unit.td.eu <- tryCatch({glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea")},error=function(e){message(e);return(NA)})
    # 
    unit.tu.td.eu <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + drainage + lat + water + bfi + airtemp + flow, 
                            UnitIn_s, EstMeth= "ML", family="gaussian",CorModels = c("locID2","yearf","Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea")
  
  }else if(GlacierBool & TailwaterBool){ #If tailwaters and glaciers
    # unit.tu <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.tailup"), addfunccol = "afvArea")
    # unit.td <- tryCatch({glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.taildown"), addfunccol = "afvArea")},error=function(e){return(NA)})
    # unit.tu.eu <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.tailup",'Exponential.taildown',"Exponential.Euclid"), addfunccol = "afvArea")
    # unit.td.eu <- tryCatch({glmssn(STREAM_AUG ~ elev + canopy + slope + precip + TAILWATER + GLACIER + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea")},error=function(e){return(NA)})
    # 
    unit.tu.td.eu <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + GLACIER + TAILWATER + drainage + lat + water + bfi + airtemp + flow, 
                            UnitIn_s, EstMeth= "ML", family="gaussian",CorModels = c("locID2","yearf","Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea")
  }else{ #If no tailwaters and no glaciers
    # unit.tu <- tryCatch({glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.tailup"), addfunccol = "afvArea")},error=function(e){return(NA)})
    # unit.td <- tryCatch({glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.taildown"), addfunccol = "afvArea")},error=function(e){return(NA)})
    # unit.tu.eu <- tryCatch({glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.tailup","Exponential.Euclid"), addfunccol = "afvArea")},error=function(e){return(NA)})
    # unit.td.eu <- tryCatch({glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + airtemp + flow, UnitIn_s, EstMeth= "REML", family="gaussian",CorModels = c("locID2","yearf","Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea")},error=function(e){return(NA)})
    unit.tu.td.eu <- glmssn(STREAM_AUG ~ elev + canopy + slope + precip + drainage + lat + water + bfi + airtemp + flow, 
                            UnitIn_s, EstMeth= "ML", family="gaussian",CorModels = c("locID2","yearf","Exponential.taildown","Exponential.Euclid"), addfunccol = "afvArea")
  }
  
  # AICcomp <- c(tryCatch({AIC(unit.tu)},error=function(e){return(NA)}),
  #              tryCatch({AIC(unit.td)},error=function(e){return(NA)}),
  #                tryCatch({AIC(unit.tu.eu)},error=function(e){return(NA)}),
  #                  tryCatch({AIC(unit.td.eu)},error=function(e){return(NA)}))
  
  BestMod<-unit.tu.td.eu#switch(which.min(AICcomp),unit.tu,unit.td,unit.tu.eu,unit.td.eu)
  #BestMod <- update(BestMod[[1]],EstMeth='ML')
  
  saveRDS(BestMod,paste(data_dir,'Regions/',paste0('BestModel_',UnitName,'.RDS')))
  save(list=ls(),file=paste(data_dir,'Regions/',paste0('Workspace_',UnitName,'.Rdata')))
  
  BestMod_r <- residuals(BestMod, cross.validation=T)
  BestMod_rdf <- getSSNdata.frame(BestMod_r)
  
  # Torgegram of fitted model residuals
  BestMod_t <- Torgegram(BestMod_r,"_resid.crossv_",nlag=20)
  jpeg(paste0(data_dir,'Regions/BestMod_',substr(Units[u],1,nchar(Units[u])-4),"_residtorg.jpg"))
    plot(BestMod_t, main= "BestModel Residuals Torgegram") 
  dev.off()

#   # #Root mean squared error (RMSE) of cross-validated predictions
#   # sqrt(mean((BestMod_rdf[,"_CrossValPred_"]-BestMod_rdf$obsval)^2))
#   # #RMSE of fixed effects only
#   # sqrt(mean((BestMod_rdf[,"_fit_"]-BestMod_rdf$obsval)^2))
#   # # 1.935
#   # # Null RMSE
#   # sqrt(mean((BestMod_rdf$obsval-mean(BestMod_rdf$obsval))^2))
#   # # 2.996
#   # 
#   # #Pseudo-r2 of cross-validated predictions. 
#   # cor(BestMod_rdf$obsval,BestMod_rdf[,"_CrossValPred_"])^2
#   # 
#   # #Read in all PredPoint DayMet temperatures
#   # DayMet_Preds_AugTemp <- lapply(list.files(here('Regions','DayMet',UnitName,'PredsBy1000')), FUN=function(x){
#   #   #x<-list.files(here('Regions','DayMet',UnitName,'PredsBy1000'))[1]
#   #   t1<-readRDS(here('Regions','DayMet',UnitName,'PredsBy1000',x))
#   #   t1 <- t1 %>% filter(month=='08') %>% 
#   #     select(ID_1KM,year,yday,month,tmax..deg.c.,tmin..deg.c.) %>%
#   #     mutate(tmean..deg.c = tmin..deg.c. + (tmax..deg.c.-tmin..deg.c.)/2) %>%
#   #     group_by(year,ID_1KM) %>% 
#   #     summarise(AugMeanTemp = mean(tmean..deg.c)) %>% 
#   #     rename(SAMPLEYEAR = year)
#   #   return(t1)
#   # }) %>% bind_rows()
#   # 
#   # nrow(DayMet_Preds_AugTemp)
#   # 
#   # Preds_AllYears <- lapply(1980:2020,FUN=function(yy){#foreach(yy = 1980:2020, .packages=c('dplyr','tidyr','SSN')) %do% {
#   #   #Replace air temp with DayMet for this prediction year
#   #   UnitIn_preddf <- getSSNdata.frame(UnitIn, "preds")
#   #   DayMetThisYear <- DayMet_Preds_AugTemp %>% filter(SAMPLEYEAR == yy)
#   #   OriginalAirData <- UnitIn_preddf %>% select(ID_1KM,SAMPLEYEAR,Air_Aug)
#   #   OriginalAirData$SAMPLEYEAR <- yy
#   #   CombiAirTemp <- OriginalAirData %>% left_join(DayMetThisYear,by=c('ID_1KM','SAMPLEYEAR'))
#   #   
#   #   UnitIn_preddf$Air_Aug <- CombiAirTemp$AugMeanTemp
#   #   UnitIn_preddf$GLACIER <- as.numeric(UnitIn_preddf$GLACIER!=-9999)
#   #   
#   #   #Standarized variables and create a new SSN for model predictions
#   #   contpred <- UnitIn_preddf %>% select(ELEV,CANOPY,SLOPE,PRECIP,CUMDRAINAG,Y_COORD,NLCD11PC,BFI,Air_Aug,Flow_Aug)
#   #   contpred.s <- stdpreds(contpred,continuous)
#   #   colnames(contpred.s) <- c("elev","canopy","slope","precip","drainage","lat","water","bfi","airtemp","flow")
#   #   UnitIn_preddf.s <- data.frame(UnitIn_preddf,contpred.s)
#   #   UnitIn_preddf.s$yearf <- factor(UnitIn_preddf.s$SAMPLEYEAR)
#   #   UnitIn_Pred<- putSSNdata.frame(UnitIn_preddf.s,UnitIn,"preds")
#   #   
#   #   
#   #   PredMod <- BestMod
#   #   PredMod$ssn.object@predpoints <- UnitIn_Pred@predpoints
#   #   p1 <- predict(PredMod,'preds')
#   #   pred1<-getSSNdata.frame(p1,'preds') %>% select("OBSPREDID",'STREAM_AUG','STREAM_AUG.predSE')
#   #   allpreds <- pred1
#   #   colnames(allpreds) <- c("OBSPREDID","predtemp","predtempse")
#   #   
#   #   write.csv(x=allpreds,file=here('AnnualPredictions',UnitName,paste0("Pred_AugTemp_",yy,".csv")),row.names=F)
#   #   write.dbf(allpreds,file=here('AnnualPredictions',UnitName,paste0("Pred_AugTemp_",yy,".dbf")))
#   # })
#   # 
# #}
#   
#   # lf<-list.files(here('AnnualPredictions',UnitName),pattern='.csv')
#   # AllTemps<-lapply(1:length(lf),FUN=function(x){
#   #   #x<-1
#   #   t1<-read.csv(here('AnnualPredictions',UnitName,lf[x])) %>% mutate(SAMPLEYEAR = 1979+x)
#   #   return(t1)
#   # }) %>% bind_rows()
#   # head(AllTemps)
#   # xt1<-xtabs(predtemp ~ OBSPREDID + SAMPLEYEAR,AllTemps)
#   # plot(xt1[1,],type='l')
#   # for(ff in 1:100){
#   #   points(xt1[ff,],type='l',col='#00000055')
#   # }
#   
#     # 
#     # 
#     # ?splitPredictions
#     # 
#     # UnitIn_Pred$ssn.object <- splitPredictions(UnitIn_Pred$ssn.object, "prednorth",chunksof=10000)
#     # UnitIn_Pred$ssn.object <- splitPredictions(UnitIn_Pred$ssn.object, "predse",chunksof=7000)
#     # UnitIn_Pred$ssn.object <- splitPredictions(UnitIn_Pred$ssn.object, "predsw",chunksof=7000)
#     # 
#     # # Make predictions
#     # spok1p1 <- predict(spok1,"prednorth-1")
#     # spok1p2 <- predict(spok1,"prednorth-2")
#     # spok1p3 <- predict(spok1,"prednorth-3")
#     # spok1p4 <- predict(spok1,"predse-1")
#     # spok1p5 <- predict(spok1,"predse-2")
#     # spok1p6 <- predict(spok1,"predsw-1")
#     # spok1p7 <- predict(spok1,"predsw-2")
#     # 
#     # # Extrac prediction data frames
#     # pred1df <- getSSNdata.frame(spok1p1,"prednorth-1")
#     # pred2df <- getSSNdata.frame(spok1p2,"prednorth-2")
#     # pred3df <- getSSNdata.frame(spok1p3,"prednorth-3")
#     # pred4df <- getSSNdata.frame(spok1p4,"predse-1")
#     # pred5df <- getSSNdata.frame(spok1p5,"predse-2")
#     # pred6df <- getSSNdata.frame(spok1p6,"predsw-1")
#     # pred7df <- getSSNdata.frame(spok1p7,"predsw-2")
#     # 
#     # # Reassemble the pieces into one batch.
#     # allpreds <- rbind(pred1df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")],pred2df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")],pred3df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")],pred4df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")],pred5df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")],pred6df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")],pred7df[,c("OBSPREDID","STREAM_AUG","STREAM_AUG.predSE")])
#     # colnames(allpreds) <- c("OBSPREDID","predtemp","predtempse")
#     # 
#     # # Export prediction dataset as a csv (good for general use) and dbf (good for GIS)
#     # write.csv(allpreds,"spok1preds.csv",row.names=F)
#     # write.dbf(allpreds,"spok1preds.dbf")
#     # 
#     # # Compare predicted and observed at fitting sites, and export. 
#     # predobs <- data.frame(spok1rdf$OBSPREDID,spok1rdf[,"_CrossValPred_"],spok1rdf$obsval)
#     # colnames(predobs) <- c("obspredid","predicted","observed")
#     # write.csv(predobs,"predobs.csv",row.names=F)
#     # 
#     # 
#   }
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
#   
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
