user_lib <- path.expand("~/R/library")
dir.create(user_lib, recursive=TRUE, showWarnings=FALSE)
.libPaths(c(user_lib, .libPaths()))

# Remove any stale lock directories
locks <- list.files(user_lib, pattern="^00LOCK", full.names=TRUE)
if(length(locks) > 0){
  message("Removing stale locks: ", paste(locks, collapse=", "))
  unlink(locks, recursive=TRUE)
}

required_packages <- c('here','SSN2','dplyr','tidyr','daymetr','foreach','doParallel','lme4','sf','terra','spmodel')
missing <- required_packages[!sapply(required_packages, requireNamespace, quietly=TRUE)]

if(length(missing) > 0){
  message("Installing: ", paste(missing, collapse=", "))
  install.packages(missing, repos="https://cloud.r-project.org", lib=user_lib)
} else {
  message("All packages already installed.")
}
