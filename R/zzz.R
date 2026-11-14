

.onLoad <- function(lib, pkg) {

  library.dynam("SelectionTools", pkg, lib)

  
  st.set.num.threads(2)
  set.info.level(0)

  
  plabsim.init()

  
  
  
}

.onAttach <- function(lib, pkg) {

  ver <- as.character(utils::packageVersion(pkg))
  has_omp <- isTRUE(.Call("st_has_openmp"))

  if (has_omp) {
    packageStartupMessage(
      paste("SelectionTools (OMP)", ver)
    )
  } else {
    packageStartupMessage(
      paste("SelectionTools", ver)
    )
  }
}
