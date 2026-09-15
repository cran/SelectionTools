###

.onLoad <- function(lib, pkg) {

  library.dynam("SelectionTools", pkg, lib)

  ## Starting values
  st.set.num.threads(2)
  set.info.level(0)

  ## Random numbers for simulations
  plabsim.init()

  ## Directory prefixes are not created in .GlobalEnv.
  ## When used, legacy user-defined st.input.dir, st.output.dir and
  ## st.data.dir values are read without modifying the user's workspace.
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
