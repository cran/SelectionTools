ph.linkage.map.create <- function(map, disperse=FALSE, file.disperse=NA,
                                   disperse.factor=100)
    .Call("st_legacy_sim", list(action="map", map=map, disperse=disperse,
          factor=disperse.factor), PACKAGE="SelectionTools")

define.genome <- function(map, length.vec=NULL)
    invisible(.Call("st_legacy_sim", list(action="genome", map=map,
                    lengths=length.vec), PACKAGE="SelectionTools"))

evaluate.allele2 <- function(PopName, Loci, missing=0)
    .Call("st_legacy_sim", list(action="allele", population=PopName,
          loci=Loci, missing=missing,
          background=plabsim$allele.missing.indicator, drop.empty=FALSE,
          frequencies=FALSE), PACKAGE="SelectionTools")

evaluate.allele.freq2 <- function(PopName, Loci, missing=0)
    .Call("st_legacy_sim", list(action="allele", population=PopName,
          loci=Loci, missing=missing,
          background=plabsim$allele.missing.indicator, drop.empty=FALSE,
          frequencies=TRUE), PACKAGE="SelectionTools")

define.effects.df <- function(LocNames) {
  for (locus in LocNames) {
    if (paste0("EFF", locus, 1) %in% list.effects()$effect) {
      warning(paste0("Effects already exist for ", locus, ", using preexisting effects."))
      next
    }
    invisible(.Call("st_legacy_sim", list(action="effect.locus", locus=as.character(locus),
                    name=paste0("EFF", locus, 1)), PACKAGE="SelectionTools"))
  }
}

sel.target.define <- function(criteria.df)
    .Call("st_legacy_sim", list(action="target.classify", criteria=criteria.df),
          PACKAGE="SelectionTools")




select.genotypes <- function(newPop, oldPop, criteria) {
  if (!is.data.frame(criteria) || ncol(criteria) < 3L)
    stop("criteria must contain locus and two allele columns")
  criteria[[1L]] <- as.character(criteria[[1L]])
  invisible(.Call("st_legacy_sim", list(action="select.genotypes",
      source=as.character(oldPop), destination=as.character(newPop),
      criteria=criteria), PACKAGE="SelectionTools"))
}

st.read.marker.data.df <- function(markerdata, format = NULL,
                                   data.set = "default") {
  if (is.null(format)) format <- "list"
  if (!is.character(format) || length(format) != 1L ||
      is.na(format) || !format %in% c("list", "matrix"))
    stop('format must be "list" or "matrix"')
  if (!is.character(data.set) || length(data.set) != 1L ||
      is.na(data.set) || !nzchar(data.set))
    stop("data.set must be one non-empty character string")
  
  if (format == "list") {
    invisible(.Call("gs_read_marker_data_DF_list", markerdata,
                    data.set, PACKAGE = "SelectionTools"))
  } else {
    if (!(is.data.frame(markerdata) || is.matrix(markerdata)))
      stop("matrix format requires a data frame or matrix")
    if (!nrow(markerdata) || !ncol(markerdata))
      stop("markerdata must contain individuals and markers")
    if (is.data.frame(markerdata) && .row_names_info(markerdata, 1L) < 0L)
      stop("supply explicit individual row names")
    markerdata <- as.matrix(markerdata)
    invisible(.Call("gs_read_marker_data_DF_matrix", markerdata,
                    data.set, PACKAGE = "SelectionTools"))
  }
  invisible(NULL)
}

st.read.map.df <- function(map, m.stretch = 1, data.set = "default") {
  if (!(is.data.frame(map) || is.matrix(map)) || ncol(map) != 3L || !nrow(map))
    stop("map must be a non-empty data frame or matrix with three MCP columns")
  if (!is.numeric(m.stretch) || length(m.stretch) != 1L ||
      !is.finite(m.stretch) || m.stretch <= 0)
    stop("m.stretch must be finite and positive")
  if (!is.character(data.set) || length(data.set) != 1L ||
      is.na(data.set) || !nzchar(data.set))
    stop("data.set must be one non-empty character string")
  map <- as.data.frame(map, stringsAsFactors = FALSE)
  
  number <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (!(is.numeric(x) || is.character(x))) stop("invalid numeric map column")
    suppressWarnings(as.double(x))
  }
  if (!(is.character(map[[1L]]) || is.factor(map[[1L]])))
    stop("marker names must be character")
  map <- list(as.character(map[[1L]]), number(map[[2L]]), number(map[[3L]]))
  invisible(.Call("st_read_map_df", map, as.double(m.stretch), data.set,
                  PACKAGE = "SelectionTools"))
}
