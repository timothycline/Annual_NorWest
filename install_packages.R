user_lib <- path.expand("~/R/library")
dir.create(user_lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(user_lib, .libPaths()))

# Remove any stale lock directories
locks <- list.files(user_lib, pattern="^00LOCK", full.names=TRUE)
if(length(locks) > 0){
  message("Removing stale locks: ", paste(locks, collapse=", "))
  unlink(locks, recursive=TRUE)
}

required_packages <- c('here','sf','terra','spmodel','SSN2',
                       'dplyr','tidyr','daymetr','foreach','doParallel','lme4')

for(pkg in required_packages){
  if(!requireNamespace(pkg, quietly=TRUE)){
    message(Sys.time(), " -- Installing: ", pkg)
    tryCatch(
      install.packages(pkg, repos="https://cloud.r-project.org", lib=user_lib),
      error = function(e) message("ERROR installing ", pkg, ": ", e$message)
    )
    if(requireNamespace(pkg, quietly=TRUE)){
      message(Sys.time(), " -- OK: ", pkg)
    } else {
      message(Sys.time(), " -- FAILED: ", pkg)
    }
  } else {
    message(Sys.time(), " -- Already installed: ", pkg)
  }
}
