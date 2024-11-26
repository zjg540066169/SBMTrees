.onLoad <- function(libname, pkgname) {
  # Path to the R script in inst/scripts
  script_path <- system.file("R_scripts", "Rcpp_function.R", package = pkgname)
  
  # Source the script to load the functions
  source(script_path, local = environment())
}