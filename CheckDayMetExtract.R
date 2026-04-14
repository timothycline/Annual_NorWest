library(here)

Units <- c('Clearwater','Midsnake','MissouriHW','Salmon','SnakeBear','Spokoot','UpMissMarias','UpYellBighorn')

cat("=== DayMet Extraction Check ===\n\n")

for(UnitName in Units){

  cat(sprintf("--- %s ---\n", UnitName))

  # Check obs RDS
  obs_path <- here('Regions','DayMet',UnitName,paste0(UnitName,'.Obs.RDS'))
  if(file.exists(obs_path)){
    obs <- tryCatch(readRDS(obs_path), error=function(e) NULL)
    if(!is.null(obs)){
      cat(sprintf("  Obs RDS:  OK  (%d rows, %d sites)\n", nrow(obs), length(unique(obs$ID_1KM))))
    } else {
      cat("  Obs RDS:  CORRUPT - could not read\n")
    }
  } else {
    cat("  Obs RDS:  MISSING\n")
  }

  # Check pred groups
  preds_dir <- here('Regions','DayMet',UnitName,'PredsBy1000')
  if(dir.exists(preds_dir)){
    files <- list.files(preds_dir, pattern='\\.RDS$')
    if(length(files) > 0){
      # Check first and last file are readable
      first <- tryCatch(readRDS(file.path(preds_dir,files[1])), error=function(e) NULL)
      last  <- tryCatch(readRDS(file.path(preds_dir,files[length(files)])), error=function(e) NULL)
      ok <- !is.null(first) & !is.null(last)
      cat(sprintf("  Pred RDS: %s (%d group files)\n", ifelse(ok,'OK','CORRUPT'), length(files)))
    } else {
      cat("  Pred RDS: MISSING - no files in PredsBy1000/\n")
    }
  } else {
    cat("  Pred RDS: MISSING - PredsBy1000/ directory not found\n")
  }

  cat("\n")
}
