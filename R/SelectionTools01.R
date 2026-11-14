



.GS_ALLELE_MAX <- 254L
.GS_STVCF_MAX_ALT <- .GS_ALLELE_MAX - 1L

.st.restore.graphics <- function(settings) {
  graphics::par(settings[!names(settings) %in% c("pin", "plt", "mfg", "fig", "new")])
  graphics::par(mfg=settings$mfg)
  if (!isTRUE(all.equal(graphics::par("fig"), settings$fig, tolerance=1e-12)))
    graphics::par(fig=settings$fig)
  graphics::par(cex=settings$cex, mex=settings$mex)
  graphics::par(new=settings$new)
  invisible(NULL)
}




st.info <-  function( lev, msg )
  {
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="g_r_info",
        arguments=list(as.integer(lev), as.character(msg))), PACKAGE="SelectionTools");
  }

"st.set.num.threads" <- function( active )
{
  if (!is.numeric(active) || length(active) != 1L || is.na(active) ||
      !is.finite(active) || active != floor(active))
    stop("'active' must be one finite integer", call. = FALSE)
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_set_num_threads",
        arguments=list(as.integer(max(1, min(256, active))))), PACKAGE="SelectionTools")
}

"st.set.openblas.threads" <- function( active )
{
  if (!is.numeric(active) || length(active) != 1L || is.na(active) ||
      !is.finite(active) || active != floor(active) || active < 1 || active > 256)
    stop("'active' must be one integer between 1 and 256", call. = FALSE)
  c  <-  .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_set_openblas_threads",
        arguments=list(as.integer(active))), PACKAGE="SelectionTools")
}


"st.get.num.threads" <- function(complete = FALSE)
{
  if (!is.logical(complete) || length(complete) != 1L || is.na(complete))
    stop("'complete' must be TRUE or FALSE", call. = FALSE)
  if (complete) return(.Call("st_thread_control", TRUE))
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_get_num_threads",
        arguments=list(t = 0L),
        result=1L), PACKAGE="SelectionTools")
  c
}





"gs.info" <-  function( lev, msg ){
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_r_info",
        arguments=list(as.integer(lev), as.character(msg))), PACKAGE="SelectionTools");
  }

"gs.set.num.threads" <- function( active )
{
  st.set.num.threads(active)
}

"gs.reset" <-  function()
{
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_reset",
        arguments=list()), PACKAGE="SelectionTools");
  }

"gs.info" <-  function( lev, msg ){
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_r_info",
        arguments=list(as.integer(lev), as.character(msg))), PACKAGE="SelectionTools");
  }

"gs.set.info.level" <- function( level,data.set ="default" )
{
  c  <-  .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_set_info_level_GV",
        arguments=list(as.integer(level), as.character(data.set))), PACKAGE="SelectionTools")
}

"gs.set.all.info.levels" <- function( level )
{
  c  <-  .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_set_all_info_levels",
        arguments=list(as.integer(level))), PACKAGE="SelectionTools")
}


"gs.get.info.level" <- function ( data.set ="default" )
{
  level = 0
  c  <-  .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_get_info_level_GV",
        arguments=list(level = as.integer(level), as.character(data.set)),
        result=1L), PACKAGE="SelectionTools")
  return (c)
}

"gs.test.mp" <- function(n = 0) {

  c  <-  .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_test_mp",
        arguments=list(as.integer(n))), PACKAGE="SelectionTools")
}


"st.set.matrix.ops" <- function(select.s="mt",select.p="st") {
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="st_set_matrix_ops",
        arguments=list(as.character(select.s), as.character(select.p)),
        result=0L), PACKAGE="SelectionTools")
  return(invisible(NULL))
}

"gs.start.timer" <- function (depth=1)
  {
    st.start.timer(depth)
  }


"gs.stop.timer" <- function ( depth=1, info.level=0 )
{
  st.stop.timer(depth,info.level)
 }


"gs.build.Z" <- function(out.filename = "Z.matrix", auxfiles = FALSE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_build_Z_01_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L, as.character(data.set)), retval = 3L, enabled = auxfiles, outputs = list(list(argument = structure(1L, names = ""), 
        header = NA, matrix = 2L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}


"gs.mme.coeff" <- function(out.filename = "mme.coeff", auxfiles = FALSE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_mme_coeff_01_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L, as.character(data.set)), retval = 3L, enabled = auxfiles, outputs = list(list(header = FALSE, 
        argument = structure(1L, names = ""), matrix = 1L, keep.names = FALSE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.build.V" <- function(out.filename = "V.matrix", auxfiles = FALSE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_build_V_01_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L, as.character(data.set)), retval = 3L, enabled = auxfiles, outputs = list(list(header = FALSE, 
        argument = structure(1L, names = ""), matrix = 1L, keep.names = FALSE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.get.Z" <- function ( data.set = "default", names = FALSE)
{
  gs.start.timer()
  Z <- .Call( "st_get_Z_named",
              data.set, names, PACKAGE = "SelectionTools")
  gs.stop.timer (info.level=1)
  return( invisible(Z) )
}


"gs.get.V" <- function ( data.set = "default")
{
  gs.start.timer()
  V <- .Call( "st_get_V",
              data.set)
  gs.stop.timer (info.level=1)
  return( invisible(V) )
}

gs.make.A <- function(Z, coding = "0 1 2")
{
  .Call("gs_make_A", Z, coding, PACKAGE = "SelectionTools")
}

"gs.get.y" <- function ( data.set = "default", names = FALSE)
{
  gs.start.timer()
  y <- .Call( "st_get_y_named",
              data.set, names, PACKAGE = "SelectionTools")
  gs.stop.timer (info.level=1)
  return( invisible(y) )
}

"gs.get.im" <- function ( data.set = "default")
{
    gs.start.timer()
    i <- as.integer(0)
    m <- as.integer(0)
    y <- .Call( "st_get_im",
               i,
               m,
               data.set)
  gs.stop.timer (info.level=1)
  return( invisible( c(i=i,m=m)) )
}


"gs.mme.restcoeff" <- function(out.filename = "mme.coeff", auxfiles = FALSE) {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_mme_restcoeff_01_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L), retval = 3L, enabled = auxfiles, outputs = list(list(header = FALSE, 
        argument = structure(1L, names = ""), matrix = 1L, keep.names = FALSE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.mme.rhs" <- function(out.filename = "mme.rhs", auxfiles = FALSE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_mme_rhs_01_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L, as.character(data.set)), retval = 3L, enabled = auxfiles, outputs = list(list(header = FALSE, 
        argument = structure(1L, names = ""), matrix = 1L, keep.names = FALSE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}


"gs.mme.restrhs" <- function(out.filename = "mme.rhs", auxfiles = FALSE) {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_mme_restrhs_01_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L), retval = 3L, enabled = auxfiles, outputs = list(list(header = FALSE, 
        argument = structure(1L, names = ""), matrix = 1L, keep.names = FALSE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.mme.solve" <- function(out.filename = "mme.solution", auxfiles = FALSE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_mme_solve_01_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L, as.character(data.set)), retval = 3L, enabled = auxfiles, outputs = list(list(header = TRUE, 
        argument = structure(1L, names = ""), matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.mme.invcoeff" <- function(out.filename = "inv.mme", auxfiles = TRUE) {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_mme_invcoeff_01_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L), retval = 3L, enabled = auxfiles, outputs = list(list(header = FALSE, 
        argument = structure(1L, names = ""), matrix = 1L, keep.names = FALSE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}











































































"gs.single.marker.reg" <- function(out.filename = "sm.reg", auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_single_marker_reg_01_GV", 
        arguments = list(NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 3L, enabled = auxfiles, 
        outputs = list(list(header = TRUE, argument = structure(1L, names = ""), matrix = 1L, keep.names = TRUE))), 
        PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.single.marker.aov" <- function(alpha = 1, out.filename = "sm.reg", auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_single_marker_aov_01_GV", 
        arguments = list(as.double(alpha), NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 4L, 
        enabled = auxfiles, outputs = list(list(header = TRUE, argument = structure(2L, names = ""), 
            matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}


"gs.lambda.const" <- function(lambda, out.filename = "lambda", auxfiles = FALSE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_lambda_const_01_GV", arguments = list(as.double(lambda), 
        NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 4L, enabled = auxfiles, outputs = list(list(header = TRUE, 
        argument = structure(2L, names = ""), matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.set.lambda" <- function(lambda, out.filename = "lambda", auxfiles = FALSE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_lambda_var_01_GV", arguments = list(as.double(lambda), 
        as.integer(length(lambda)), NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 5L, 
        enabled = auxfiles, outputs = list(list(header = TRUE, argument = structure(3L, names = ""), 
            matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.lambda.hsq" <- function(hsq, out.filename = "lambda", auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_lambda_const_02_GV", arguments = list(as.double(hsq), 
        NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 4L, enabled = auxfiles, outputs = list(list(header = TRUE, 
        argument = structure(2L, names = ""), matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}



"gs.lambda.reg" <- function(hsq, out.filename = "lambda", auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_lambda_reg_01_GV", arguments = list(as.double(hsq), 
        NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 4L, enabled = auxfiles, outputs = list(list(header = TRUE, 
        argument = structure(2L, names = ""), matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.lambda.aov" <- function(hsq, alpha = 1, out.filename = "lambda", auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_lambda_aov_01_GV", arguments = list(as.double(hsq), 
        as.double(alpha), NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 5L, enabled = auxfiles, 
        outputs = list(list(header = TRUE, argument = structure(3L, names = ""), matrix = 1L, keep.names = TRUE))), 
        PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.lambda.emstep" <- function(constvar = FALSE, out.filename = "lambda", auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_lambda_emstep_01_GV", arguments = list(as.integer(constvar), 
        NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 4L, enabled = auxfiles, outputs = list(list(header = TRUE, 
        argument = structure(2L, names = ""), matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}


"gs.lambda.rmlv" <- function(maxiter = 1000, precision = 0.001, out.filename = "lambda", auxfiles = TRUE, 
    data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_lambda_rmlv_01_GV", arguments = list(as.integer(maxiter), 
        as.double(precision), NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 5L, enabled = auxfiles, 
        outputs = list(list(header = TRUE, argument = structure(3L, names = ""), matrix = 1L, keep.names = TRUE))), 
        PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.lambda.rmlc" <- function(maxiter = 1000, precision = 1e-04, hsq = 0.9, out.filename = "lambda", auxfiles = TRUE, 
    data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_lambda_rmlc_01_GV", arguments = list(as.integer(maxiter), 
        as.double(precision), as.double(hsq), NULL, as.integer(auxfiles), 0L, as.character(data.set)), 
        retval = 6L, enabled = auxfiles, outputs = list(list(header = TRUE, argument = structure(4L, names = ""), 
            matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.lambda.rmlr" <- function(maxiter = 1000, precision = 0.001, hsq = 0.9, out.filename = "lambda", auxfiles = TRUE, 
    data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_lambda_rmlr_01_GV", arguments = list(as.integer(maxiter), 
        as.double(precision), as.double(hsq), NULL, as.integer(auxfiles), 0L, as.character(data.set)), 
        retval = 6L, enabled = auxfiles, outputs = list(list(header = TRUE, argument = structure(4L, names = ""), 
            matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.lambda.rmla" <- function(alpha = 1, maxiter = 1000, precision = 1e-04, hsq = 0.9, out.filename = "lambda", 
    auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_lambda_rmla_01_GV", arguments = list(as.double(alpha), 
        as.integer(maxiter), as.double(precision), as.double(hsq), NULL, as.integer(auxfiles), 0L, as.character(data.set)), 
        retval = 7L, enabled = auxfiles, outputs = list(list(header = TRUE, argument = structure(5L, names = ""), 
            matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}


"gs.set.allele.codes" <- function(    aa = 0,
                                      aA = 1,
                                      AA = 2,
                                      an = 0.5,
                                      nn = 1  ,
                                      nA = 1.5,
                                      data.set ="default" )
{
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_set_allele_codes_01_GV",
        arguments=list(as.double(aa), as.double(aA), as.double(AA), as.double(an), as.double(nn), as.double(nA), as.character(data.set)),
        result=0L), PACKAGE="SelectionTools")
  return( invisible(NULL) )
}


"gs.build.esvs" <- function (es.size         ,
                                base.set       = "default" ,
                                estimation.set = "none"    ,
                                validation.set = "none"    ,
                                es_m = NULL ,
                                es_p = NULL ,
                                vs_m = NULL ,
                                vs_p = NULL ,
                                auxfiles = FALSE                )
{
  st.output.dir <- get0("st.output.dir", envir = .GlobalEnv,
                        inherits = FALSE, ifnotfound = "")

  if (TRUE == auxfiles &&
      (is.null(es_m) || is.null(es_p) || is.null(vs_m) || is.null(vs_p)))
    stop("es_m, es_p, vs_m, and vs_p must be supplied when auxfiles = TRUE")
  if ( "" == st.output.dir ) op <- "" else op <- "/"
    out.es_m    <- if (is.null(es_m)) "" else paste(st.output.dir,op,es_m,sep="")
    out.es_p    <- if (is.null(es_p)) "" else paste(st.output.dir,op,es_p,sep="")
    out.vs_m    <- if (is.null(vs_m)) "" else paste(st.output.dir,op,vs_m,sep="")
    out.vs_p    <- if (is.null(vs_p)) "" else paste(st.output.dir,op,vs_p,sep="")
    gs.start.timer()
    retval <- 0;
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_build_esvs_01_GV",
        arguments=list(as.integer(es.size), as.character(base.set), as.character(estimation.set), as.character(validation.set), as.character(out.es_m), as.character(out.es_p), as.character(out.vs_m), as.character(out.vs_p), as.integer  (auxfiles), retval = as.integer(retval)),
        result=0L), PACKAGE="SelectionTools")
    gs.stop.timer (info.level=1)
  }

"gs.build.tsps" <- function (ts.size         ,
                             base.set       = "default" ,
                             training.set   = "none"    ,
                             prediction.set = "none"    ,
                             es_m = NULL ,
                             es_p = NULL ,
                             vs_m = NULL ,
                             vs_p = NULL ,
                             auxfiles = FALSE                )
{
  st.output.dir <- get0("st.output.dir", envir = .GlobalEnv,
                        inherits = FALSE, ifnotfound = "")

  if (TRUE == auxfiles &&
      (is.null(es_m) || is.null(es_p) || is.null(vs_m) || is.null(vs_p)))
    stop("es_m, es_p, vs_m, and vs_p must be supplied when auxfiles = TRUE")
  if ( "" == st.output.dir ) op <- "" else op <- "/"
    out.es_m    <- if (is.null(es_m)) "" else paste(st.output.dir,op,es_m,sep="")
    out.es_p    <- if (is.null(es_p)) "" else paste(st.output.dir,op,es_p,sep="")
    out.vs_m    <- if (is.null(vs_m)) "" else paste(st.output.dir,op,vs_m,sep="")
    out.vs_p    <- if (is.null(vs_p)) "" else paste(st.output.dir,op,vs_p,sep="")
    gs.start.timer()
    retval <- 0;
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_build_esvs_01_GV",
        arguments=list(as.integer(ts.size), as.character(base.set), as.character(training.set), as.character(prediction.set), as.character(out.es_m), as.character(out.es_p), as.character(out.vs_m), as.character(out.vs_p), as.integer  (auxfiles), retval = as.integer(retval)),
        result=0L), PACKAGE="SelectionTools")
    gs.stop.timer (info.level=1)
  }


"gs.estimate.gv" <- function(estimation.set = "default", prediction.set = "default", out.filename = "breeding.values", 
    auxfiles = TRUE) {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_estimate_bv_01_GV", arguments = list(as.character(estimation.set), 
        as.character(prediction.set), NULL, as.integer(auxfiles), 0L), retval = 5L, enabled = auxfiles, 
        outputs = list(list(header = FALSE, skip = 1, argument = structure(3L, names = ""), matrix = 0L, 
            keep.names = TRUE, names = c("i", "y", "yhat")))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}



"gs.predict.genotypes" <- function(training.set = "default", prediction.set = "default", out.filename = "breeding.values", 
    auxfiles = TRUE) {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_estimate_bv_01_GV", arguments = list(as.character(training.set), 
        as.character(prediction.set), NULL, as.integer(auxfiles), 0L), retval = 5L, enabled = auxfiles, 
        outputs = list(list(header = FALSE, skip = 1, argument = structure(3L, names = ""), matrix = 0L, 
            keep.names = TRUE, names = c("i", "y", "yhat")))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}






"gs.restrict.marker.data.01" <- function ( outfile,
                                          traitfile = NULL,
                                          NoAll.MAX = 999,
                                          MaMis.MAX = 1,
                                          ExHet.MIN = 0,
                                          InMis.MAX = 1,
                                          data.set  = "default" )
{
  st.input.dir <- get0("st.input.dir", envir = .GlobalEnv,
                       inherits = FALSE, ifnotfound = "")
  st.output.dir <- get0("st.output.dir", envir = .GlobalEnv,
                        inherits = FALSE, ifnotfound = "")

    if ( "" == st.input.dir )  ip <- "" else ip <- "/"
    if ( "" == st.output.dir ) op <- "" else op <- "/"
    out.filename    <- paste(st.input.dir,op,outfile,sep="")
    if ( !is.null(traitfile) )
      trait.filename  <- paste(st.input.dir,ip,traitfile,sep="")
    else
      trait.filename = ""
    gs.start.timer()
    retval <- 0;
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_restrict_marker_data_01_GV",
        arguments=list(as.character(out.filename), as.character(trait.filename), as.integer(NoAll.MAX), as.double (MaMis.MAX), as.double (ExHet.MIN), as.double (InMis.MAX), retval = as.integer(retval), as.character(data.set)),
        result=0L), PACKAGE="SelectionTools")
    gs.stop.timer (info.level=1)
  }







"st.reset" <-  function(){
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_reset",
        arguments=list()), PACKAGE="SelectionTools");
  }

"st.data.set.list" <- function(digits=getOption("digits"))
{
  digits <- as.integer(digits)
  if (length(digits) != 1L || is.na(digits) || digits < 1L || digits > 22L)
    stop("digits must be one integer between 1 and 22")
  x <- .Call("st_data_set_list")
  attr(x,"digits") <- digits
  class(x) <- c("st.data.set.list","data.frame")
  x
}

"print.st.data.set.list" <- function(x,digits=attr(x,"digits"),complete=FALSE,...)
{
  if (is.null(digits)) digits <- getOption("digits")
  digits <- as.integer(digits)
  if (length(digits) != 1L || is.na(digits) || digits < 1L || digits > 22L)
    stop("digits must be one integer between 1 and 22")
  if (length(complete) != 1L || is.na(complete))
    stop("complete must be TRUE or FALSE")
  y <- if (isTRUE(complete)) x else
       x[,c("data.set","no.markers","no.individuals","map",
            "BitsPerMarker","A_matrix","Z_matrix","ram.core.MB","ram.work.MB"),drop=FALSE]
  class(y) <- "data.frame"
  print(y,digits=digits,...)
  invisible(x)
}

"st.compact.data.set" <- function(data.set="default")
{
  invisible(.Call("st_compact_data_set",as.character(data.set)))
}

"st.data.set.remove" <- function(data.set)
{
  invisible(.Call("st_data_set_remove",as.character(data.set)))
}

"st.data.set.remove.all" <- function()
{
  invisible(.Call("st_data_set_remove_all"))
}

"st.data.combine" <- function(target.data.set, source.data.sets)
{
   invisible(.Call("st_data_combine", target.data.set, source.data.sets,
                   PACKAGE="SelectionTools"))
}

"st.keep.intersection.markers" <- function(data.set1,data.set2)
{
  invisible(.Call("st_keep_intersection_markers",
                  as.character(data.set1),as.character(data.set2)))
}

"st.info" <-  function( lev, msg ){
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="st_r_info",
        arguments=list(as.integer(lev), as.character(msg))), PACKAGE="SelectionTools");
  }

"st.set.info.level" <- function( level )
{
  c  <-  .Call("st_legacy_marker", list(action="legacy.result",
        routine="st_set_info_level",
        arguments=list(as.integer(level))), PACKAGE="SelectionTools")
}

"st.get.info.level" <- function ()
{
  level = 0
  c  <-  .Call("st_legacy_marker", list(action="legacy.result",
        routine="st_get_info_level",
        arguments=list(level = as.integer(level)),
        result=1L), PACKAGE="SelectionTools")
  return (c)
}

"st.start.timer" <- function(depth=1)
  {
    if (!is.numeric(depth) || length(depth) != 1L || is.na(depth) ||
        !is.finite(depth) || depth != floor(depth) || depth < 1 || depth > 5)
        stop("'depth' must be one integer between 1 and 5", call. = FALSE)
    depth <- as.integer(depth) - 1L

      c  <-  .Call("st_legacy_marker", list(action="legacy.result",
        routine="st_start_timer",
        arguments=list(as.integer(depth), as.double( proc.time()[3] ))), PACKAGE="SelectionTools")
  }

"st.stop.timer" <- function ( depth=1, info.level=0 )
{
    if (!is.numeric(depth) || length(depth) != 1L || is.na(depth) ||
        !is.finite(depth) || depth != floor(depth) || depth < 1 || depth > 5)
        stop("'depth' must be one integer between 1 and 5", call. = FALSE)
    depth <- as.integer(depth) - 1L

    etime = 0
    c  <-  .Call("st_legacy_marker", list(action="legacy.result",
        routine="st_stop_timer",
        arguments=list(as.integer(depth), etime=as.double( etime )),
        result=2L), PACKAGE="SelectionTools")

    time <- proc.time()[3] - c

    if (time > 3600){
      time <- time / 3600
      tunit <- "h"
    } else
    if (time > 60) {
      time <- time / 60
      tunit <- "m"
    } else
    {
      tunit <- "s"
    }
    st.info(info.level,paste("Elapsed time:",paste(round(time,digits=2),tunit,sep="")))
  }



"st.read.marker.data" <- function (filename, format, auxfiles = FALSE, data.set = "default") 
{
    st.start.timer()
    routine <- switch(format, l = "gs_read_marker_data_01_GV", m = "gs_read_marker_data_02_GV", n = "gs_read_marker_data_03_GV", 
        t = "gs_read_marker_data_04_GV", stop("Unsupported marker format"))
    result <- .Call("st_legacy_marker", list(action = "call", routine = routine, arguments = list(as.character(.st.legacy.path(filename)), 
        NULL, NULL, NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 6L, enabled = auxfiles, 
        list = TRUE, outputs = list(list(argument = 2L, header = TRUE, check.names = FALSE), list(argument = 3L, 
            header = TRUE, check.names = FALSE), list(argument = 4L, header = TRUE, check.names = FALSE)), 
        output.names = c("individual.list", "marker.list", "genotypes")), PACKAGE = "SelectionTools")
    st.stop.timer(info.level = 1)
    invisible(result)
}

"st.read.performance.data" <- function (filename, out.filename = "y.vector", auxfiles = TRUE, data.set = "default") 
{
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_read_performance_data_01_GV", 
        arguments = list(as.character(.st.legacy.path(filename)), NULL, as.integer(auxfiles), 0L, as.character(data.set)), 
        retval = 4L, enabled = auxfiles, outputs = list(list(argument = 2L, skip = 1L, names = c("i", 
            "y")))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"st.return.performance.data" <- function (out.filename = "y.vector", auxfiles = TRUE, data.set = "default") 
{
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_return_performance_data_01_GV", 
        arguments = list(NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 3L, enabled = auxfiles, 
        outputs = list(list(argument = 1L, skip = 1L, colClasses = c("character", "numeric"), names = c("i", 
            "y")))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}


"st.read.map" <- function (filename, skip = 0, format = "cpms", m.stretch = 1, auxfiles = TRUE, data.set = "default") 
{
    st.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_read_map_01_GV", arguments = list(as.character(.st.legacy.path(filename)), 
        as.integer(skip), as.character(format), as.double(m.stretch), NULL, as.integer(auxfiles), 0L, 
        as.character(data.set)), retval = 7L, enabled = auxfiles, outputs = list(list(argument = 5L, 
        names = c("Chrom", "Pos", "Name", "Class")))), PACKAGE = "SelectionTools")
    st.stop.timer(info.level = 1)
    invisible(result)
}

.st.valid.identifier <- function(x, max.bytes=255L)
    .Call("st_legacy_marker", list(action="vcf.identifiers", data=x,
          max.bytes=max.bytes), PACKAGE="SelectionTools")

"st.dataframe.to.STvcf" <- function(data)
    .Call("st_legacy_marker", list(action="vcf.from.frame", data=data),
          PACKAGE="SelectionTools")

"st.STvcf.to.dataframe" <- function(STvcf)
    .Call("st_legacy_marker", list(action="vcf.to.frame", data=STvcf),
          PACKAGE="SelectionTools")

"read.vcf" <- function(vcffile, data.set="default", pedigree=TRUE,
                       strictACTG=FALSE)
{
   if (!is.logical(strictACTG) || length(strictACTG) != 1L ||
       is.na(strictACTG))
      stop("strictACTG must be TRUE or FALSE")
   if (length(vcffile) != 1L || is.na(vcffile) ||
       !is.character(vcffile) || !nzchar(vcffile))
      stop("vcffile must be one non-empty character string")
   if (length(data.set) != 1L || is.na(data.set) ||
       !is.character(data.set) || !nzchar(data.set))
      stop("data.set must be one non-empty character string")

   input <- vcffile
   if (grepl("\\.gz$", vcffile, ignore.case=TRUE)) {
      input <- tempfile(fileext=".vcf")
      on.exit(unlink(input), add=TRUE)

      con.in <- gzfile(vcffile, open="rb")
      on.exit(try(close(con.in), silent=TRUE), add=TRUE)
      con.out <- file(input, open="wb")
      on.exit(try(close(con.out), silent=TRUE), add=TRUE)

      repeat {
         z <- readBin(con.in, what="raw", n=1048576L)
         if (!length(z)) break
         writeBin(z, con.out)
      }
      close(con.in)
      close(con.out)
   }

   input <- as.character(input)
   attr(input,"pedigree.enabled") <- isTRUE(pedigree)
   status <- .Call("st_read_vcf", input,
                   as.character(data.set), strictACTG)
   if (!identical(status, 0L))
      stop("read.vcf failed")
   invisible(NULL)
}


"write.vcf" <- function(vcffile, data.set="default", pedigree=TRUE)
{
   if (length(vcffile) != 1L || is.na(vcffile) ||
       !is.character(vcffile) || !nzchar(vcffile))
      stop("vcffile must be one non-empty character string")
   if (length(data.set) != 1L || is.na(data.set) ||
       !is.character(data.set) || !nzchar(data.set))
      stop("data.set must be one non-empty character string")

   output <- vcffile
   compressed <- grepl("\\.gz$", vcffile, ignore.case=TRUE)
   if (compressed) {
      output <- tempfile(fileext=".vcf")
      on.exit(unlink(output), add=TRUE)
   }

   output <- as.character(output)
   attr(output,"pedigree.enabled") <- isTRUE(pedigree)
   status <- .Call("st_write_vcf", output,
                   as.character(data.set))
   if (!identical(status, 0L))
      stop("write.vcf failed")

   if (compressed) {
      con.in <- file(output, open="rb")
      on.exit(try(close(con.in), silent=TRUE), add=TRUE)
      con.out <- gzfile(vcffile, open="wb")
      on.exit(try(close(con.out), silent=TRUE), add=TRUE)

      repeat {
         z <- readBin(con.in, what="raw", n=1048576L)
         if (!length(z)) break
         writeBin(z, con.out)
      }
      close(con.in)
      close(con.out)
   }

   invisible(NULL)
}


"st.load.vcf.data" <- function(STvcf, data.set="default", strictACTG=FALSE)
{
   if (!is.logical(strictACTG) || length(strictACTG) != 1L ||
       is.na(strictACTG))
      stop("strictACTG must be TRUE or FALSE")
   st.start.timer()
   status <- .Call("st_load_vcf_data",
                   STvcf,
                   as.character(data.set), strictACTG)
   st.stop.timer(info.level=1)
   invisible(status)
}


"st.markerdata.to.STvcf" <- function(data.set="default")
{
   st.start.timer()
   x <- .Call("st_markerdata_to_STvcf", as.character(data.set))
   if (!is.null(x))
   {
      map.text <- paste0(seq_along(x$CHROM.levels), "=",
                         x$CHROM.levels, collapse=", ")
      st.info(0, paste("Chromosome coding:", map.text))
   }
   st.stop.timer(info.level=1)
   x
}


"st.load.performance.data" <- function(performance.data,
                                         data.set="default")
{
   if (!(is.data.frame(performance.data) || is.matrix(performance.data)))
      stop("performance.data must be a data frame or matrix")

   d <- as.data.frame(performance.data, stringsAsFactors=FALSE,
                      check.names=FALSE)
   nms <- names(d)
   if (all(c("individual", "performance") %in% nms))
   {
      id <- d[["individual"]]
      y <- d[["performance"]]
   }
   else if (all(c("ii", "yy") %in% nms))
   {
      id <- d[["ii"]]
      y <- d[["yy"]]
   }
   else if (2L == ncol(d))
   {
      id <- d[[1L]]
      y <- d[[2L]]
   }
   else
   {
      stop("performance data require individual and performance columns")
   }

   id <- as.character(id)
   if (any(!.st.valid.identifier(id, max.bytes=Inf)))
      stop("performance individual IDs are invalid for SelectionTools")

   if (is.factor(y))
      y <- as.character(y)

   numeric.y <- suppressWarnings(as.double(y))
   missing.y <- is.na(y) | trimws(as.character(y)) == "NA"
   numeric.y[is.na(numeric.y) & !missing.y] <- NaN

   st.start.timer()
   invisible(.Call("st_load_performance_data",
                   id,
                   numeric.y,
                   as.character(data.set)))
   st.stop.timer(info.level=1)
   invisible(NULL)
}


"st.get.map" <- function(filename = "map", data.set = "default") {
    st.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_get_map_01_GV", arguments = list(NULL, 
        0L, as.character(data.set)), retval = 2L, enabled = TRUE, outputs = list(list(header = FALSE, 
        check.names = FALSE, argument = structure(1L, names = ""), matrix = 0L, keep.names = TRUE, names = c("Chrom", 
            "Pos", "Name", "Class")))), PACKAGE = "SelectionTools")
    st.stop.timer(info.level = 1)
    invisible(result)
}

"st.restrict.marker.data" <- function (ind.list = NULL, ind.file = NULL, mar.list = NULL, mar.file = NULL, 
    NoAll.MAX = 999, MaMis.MAX = 1, ExHet.MIN = 0, InMis.MAX = 1, data.set = "default") 
{
    st.start.timer()
    if (!is.null(ind.file) && !is.null(ind.list)) {
        st.info(-2, "Provide either ind.file or ind.list")
        return(invisible(NULL))
    }
    if (!is.null(mar.file) && !is.null(mar.list)) {
        st.info(-2, "Provide either mar.file or mar.list")
        return(invisible(NULL))
    }
    .Call("st_legacy_marker", list(action = "call", routine = "gs_restrict_marker_data_03_GV", arguments = list(if (is.null(ind.file)) "" else as.character(.st.legacy.path(ind.file)), 
        if (is.null(mar.file)) "" else as.character(.st.legacy.path(mar.file)), as.integer(NoAll.MAX), 
        as.double(MaMis.MAX), as.double(ExHet.MIN), as.double(InMis.MAX), 0L, as.character(data.set)), 
        retval = 7L, inputs = c(if (!is.null(ind.list)) list(list(argument = 1L, data = ind.list)), if (!is.null(mar.list)) list(list(argument = 2L, 
            data = mar.list)))), PACKAGE = "SelectionTools")
    st.stop.timer(info.level = 1)
}

"st.copy.marker.data" <- function (  target.data.set,
                                     source.data.set )
{

  retval <- 0;
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_copy_marker_data_01_GV",
        arguments=list(retval = as.integer(retval), as.character(target.data.set), as.character(source.data.set))), PACKAGE="SelectionTools")
}



"st.write.marker.data" <- function ( format     = "m",
                                     lfilename  = NULL,
                                     mfilename  = NULL,
                                     nfilename  = NULL,
                                     ifilename  = "",
                                     f.ind      = 0,
                                     l.ind      = 0,
                                     add.inum   = ".1",
                                     auxfiles   = FALSE,
                                     data.set   = "default" )
{
  format.code <- as.character(format)[1L]
  if (format.code %in% c("l", "a") && is.null(lfilename))
    stop("lfilename must be supplied for format 'l' or 'a'")
  if (identical(format.code, "m") && is.null(mfilename))
    stop("mfilename must be supplied for format 'm'")
  if (identical(format.code, "n") && is.null(nfilename))
    stop("nfilename must be supplied for format 'n'")
  st.start.timer()
  
  output.filename <- function(filename, extension = "") {
    st.output.dir <- get0("st.output.dir", envir = .GlobalEnv,
                          inherits = FALSE, ifnotfound = "")

    if (is.null(filename)) return("")
    is.absolute <- startsWith(filename, "/") ||
                   startsWith(filename, "\\") ||
                   (nchar(filename) >= 3L &&
                    grepl("^[A-Za-z]$", substr(filename, 1L, 1L)) &&
                    substr(filename, 2L, 2L) == ":" &&
                    substr(filename, 3L, 3L) %in% c("/", "\\"))
    if (is.absolute) {
      paste0(filename, extension)
    } else if ("" == st.output.dir) {
      paste0(filename, extension)
    } else {
      paste(st.output.dir, "/", filename, extension, sep="")
    }
  }
  
  lpo.filename <- output.filename(lfilename, ".lpo")
  
  mmp.filename <- output.filename(mfilename, ".mmp")
  
  mpo.filename <- output.filename(mfilename, ".mpo")
  
  ngp.filename <- output.filename(nfilename, ".ngp")
  
  nmp.filename <- output.filename(nfilename, ".nmp")
  
  npo.filename <- output.filename(nfilename, ".npo")
  
  ind.filename <- ""
  if ( "" != ifilename)
    ind.filename <- output.filename(ifilename)

  retval <- 0;
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_write_marker_data_02_GV",
        arguments=list(as.character(format), as.character(lpo.filename), as.character(mmp.filename), as.character(mpo.filename), as.character(ngp.filename), as.character(nmp.filename), as.character(npo.filename), as.character(ind.filename), as.integer(f.ind), as.integer(l.ind), as.character(add.inum), as.integer(auxfiles), retval = as.integer(retval), as.character(data.set)),
        result=0L), PACKAGE="SelectionTools")

  outp <- NULL
  st.stop.timer (info.level=1)
  return( invisible(outp) )

}

"st.write.map" <- function(mfilename, auxfiles = TRUE, data.set = "default")
{
  st.output.dir <- get0("st.output.dir", envir = .GlobalEnv,
                       inherits = FALSE, ifnotfound = "")
  st.start.timer()
  if ("" == st.output.dir) op <- "" else op <- "/"
  mmp.filename <- paste(st.output.dir,op,mfilename,".mmp",sep="")
  x <- .Call("st_legacy_marker", list(action="map.write", file=as.character(mmp.filename),
      dataset=as.character(data.set), environment=environment(),
      output=list(header=FALSE, names=c("Chrom","Pos","Name","Class"))),
      PACKAGE="SelectionTools")
  st.stop.timer(info.level=1)
  invisible(x)
}



"st.recode.ref" <- function ( reference  = 0,
                              descending = FALSE,
                              data.set   = "default" )
{
  st.start.timer()
  retval <- 0

  if (is.character(reference)) {
    if ((1L != length(reference)) || is.na(reference) || !nzchar(reference))
      stop("reference must be one non-empty individual name or one positive integer position")
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_recode_ref_name_01_GV",
        arguments=list(as.character(reference), as.integer(descending), retval = as.integer(retval), as.character(data.set)),
        result=0L), PACKAGE="SelectionTools")
  } else {
    if ((!is.numeric(reference) && !is.integer(reference)) ||
        (1L != length(reference)) || is.na(reference) ||
        !is.finite(reference) || (reference != floor(reference)) ||
        (reference < 1) || (reference > .Machine$integer.max))
      stop("reference must be one non-empty individual name or one positive integer position")
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_recode_ref_01_GV",
        arguments=list(as.integer(reference), as.integer(descending), retval = as.integer(retval), as.character(data.set)),
        result=0L), PACKAGE="SelectionTools")
  }

  st.stop.timer (info.level=1)
}

"st.def.hblocks" <- function(hap = 5, hap.unit = 1, ld.threshold = 0.8, ld.criterion = "none", tolerance = 0, 
    hap.symbol = "b", out.filename = "blocks", auxfiles = TRUE, data.set = "default") {
    st.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_def_hblocks_01_GV", arguments = list(as.integer(hap), 
        as.integer(hap.unit), as.double(ld.threshold), as.character(ld.criterion), as.integer(tolerance), 
        as.character(hap.symbol), NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 9L, 
        enabled = auxfiles, outputs = list(list(header = FALSE, argument = structure(7L, names = ""), 
            matrix = 0L, keep.names = TRUE, names = c("Chrom", "Pos", "Name", "Class", "Markers")))), 
        PACKAGE = "SelectionTools")
    st.stop.timer(info.level = 1)
    invisible(result)
}

"st.recode.hbc" <- function ( reference  = 1,
                             data.set   = "default" )
{
  st.start.timer()
  retval <- 0;
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_recode_hbc_01_GV",
        arguments=list(as.integer(reference), retval = as.integer(retval), as.character(data.set)),
        result=0L), PACKAGE="SelectionTools")
  st.stop.timer (info.level=1)
}

"st.recode.hil" <- function(out.filename = "blocks", auxfiles = TRUE, data.set = "default") {
    st.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_recode_hil_01_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L, as.character(data.set)), retval = 3L, enabled = auxfiles, outputs = list(list(header = FALSE, 
        argument = structure(1L, names = ""), matrix = 0L, keep.names = TRUE, names = c("Block", "AlleleNr", 
            "AlleleDef")))), PACKAGE = "SelectionTools")
    st.stop.timer(info.level = 1)
    invisible(result)
}

"st.plot.ggt.src" <- function (
                        map   ,
                        pop   ,
                        f.lr  ,
                        f.tb  ,
                        d.h   ,
                        d.v   ,
                        p.t   ,
                        p.s   ,
                        d.t   ,
                        d.map ,
                        c.nme ,
                        i.nme ,
                        color ,
                        z.min ,
                        z.max ,
                        cex.chrom,
                        cex.ind
                       )
{
  m <- map
  g <- pop
  missing.allele <- c("-1", "--", "-", ".", "NA")
  
  
   oldpar <- graphics::par(no.readonly = TRUE)
   on.exit(.st.restore.graphics(oldpar), add = TRUE)
  existing.screens <- split.screen()
  if (!identical(existing.screens, FALSE))
    stop("st.plot.ggt cannot be used while split-screen mode is active")
  par(mfrow=c(1L,1L), mar=c(0,0,0,0))
  
  
  n.ind <- ncol(g)
  
  clr <- color
  
  n.chrom <- max(m[,1])
  
  if ( ( 1 == length(c.nme) ) && ( "" == c.nme[1] ) ) {
    c.nme <- paste("Chr. ",1:n.chrom)
  }
  
  if ( ( 1 == length(i.nme) ) && ( "" == i.nme[1] ) ) {
    i.nme <- names(g)
  }
  
  {
    m.chrom <- 1:n.chrom
    l.chrom <- 1:n.chrom
    b.chrom <- vector("list",n.chrom)
    p.chrom <- vector("list",n.chrom)
    for (c in 1:n.chrom) {
      idx <- m[,1] == c
      m.chrom[c] <- sum(idx)
      m.pos <- m[idx,][,2]
      l.chrom[c] <- max(m.pos) - min(m.pos)
      m.brk <- 0:length(m.pos)
      if (length(m.pos) > 1L)
        for (j in seq.int(2L, length(m.pos))) {
          m.brk[j] <- m.pos[j-1] + (m.pos[j] - m.pos[j-1]) / 2
        }
      m.brk[1+length(m.pos)] <- m.pos[length(m.pos)]
      b.chrom[[c]] <- m.brk
      p.chrom[[c]] <- m.pos
    }
  }
  
  ndc  <- 1 + n.chrom
  layout.length <- l.chrom
  if (any(layout.length <= 0))
    layout.length[layout.length <= 0] <- if (any(layout.length > 0))
      min(layout.length[layout.length > 0]) else 1
  f.disp <- c(f.lr*mean(layout.length),layout.length)
  f.disp <- f.disp / sum (f.disp)
  s.disp <- f.disp
  for (kk in 2:ndc) s.disp[kk] <- s.disp[kk-1] +s.disp[kk]
  s.disp <- c(0,s.disp)
  
  dcol <- matrix(0,ncol=2,nrow=ndc)
  for (kk in 1:ndc)
    {
      dcol[kk,1] = s.disp[kk];
      dcol[kk,2] = s.disp[kk+1];
    }

  
  ndr <- 2 + n.ind
  f.disp <- c ( f.tb, rep(1,n.ind) , f.tb)
  f.disp <- f.disp / sum (f.disp)
  s.disp <- f.disp
  for (kk in 2:ndr) s.disp[kk] <- s.disp[kk-1] +s.disp[kk]
  s.disp <- c(0,s.disp)
  
  drow <- matrix(0,ncol=2,nrow=ndr)
  for (kk in 1:ndr)
    {
      drow[kk,1] = s.disp[kk];
      drow[kk,2] = s.disp[kk+1];
    }

  
  figs <- matrix (ncol=4,nrow=ndr*ndc )
  for (i in 1:ndr) for (c in 1:ndc)
    {
      idx  <- (i-1)*ndc + c
      figs[idx,] <- c(dcol[c,],drow[i,]  )
    }
  
  sve      <- 1 - figs[,3]
  figs[,3] <- 1 - figs[,4]
  figs[,4] <- sve
  
  figs[,1] <- figs[,1] + d.h
  figs[,2] <- figs[,2] - d.h
  figs[,3] <- figs[,3] + d.v
  figs[,4] <- figs[,4] - d.v
  
  figs <- round(figs,digits=5)
  
  neu <- split.screen ( figs  )
  on.exit(close.screen(all.screens=TRUE), add=TRUE, after=FALSE)
  
  
  for (c in 1:n.chrom)
    {
      scr  <- (c+1)
      screen(scr)
      plot(0,0,t="n",axes=FALSE, xlab="", ylab="")
      text (0,0,c.nme[c],cex=cex.chrom)

      scr <- (ndr-1) * ndc + (c+1)
      screen(scr)
      mi <- min(b.chrom[[c]])
      ma <- max(b.chrom[[c]])
     if ( ( p.t == TRUE) || ( p.s == TRUE) )
       {
         plot(0,0, axes=FALSE,t="n", xlab="", ylab="", xlim=c(mi,ma),ylim=c(-d.t,d.t) )
         lines( c(mi,ma),c(0,0))
       }
      if ( p.t == TRUE)
        {
          for (mm in 1:length(p.chrom[[c]]) )
            lines (c(p.chrom[[c]][mm],p.chrom[[c]][mm]),c(0,d.t/2))
        }
      if ( p.s == TRUE)
        {
          t <- seq(0,ma,d.map)
          text( t,rep(-d.t/2,length(t)),t  )
        }

     }
  
  for (i in 1:n.ind)
    {
      
      scr  <- (i)*ndc + (1)
      screen(scr)
      plot(0,0,t="n",axes=FALSE, xlab="", ylab="")
      text (0,0,i.nme[i],cex=cex.ind)
    }
  
  for (i in 1:n.ind) for (c in 1:n.chrom)
    {
      
      scr  <- (i)*ndc + (c+1)
      screen(scr)
     
      p.m <- matrix( 0 ,
                   nrow=m.chrom[c] , ncol=2 ,
                   byrow=FALSE                    )
      
      idx <- m[,1] == c
      m.nme <- as.character ( m[idx,][,3] )
      idx2 <- row.names(g) %in% m.nme
      p.g <- as.character(g[idx2,i])
      
      if (length(p.g) != m.chrom[c])
        stop("Marker map and genotype matrix do not match")
      
      for (p in 1:m.chrom[c]) {
        gt <- p.g[p]
        if (is.na(gt) || gt %in% missing.allele) {
          p.m[p,] <- NA_real_
          next
        }
        splt <- strsplit(gt,"/",fixed=TRUE)[[1]]
        if (2L != length(splt))
          stop(paste("Invalid genotype", gt, "at marker", m.nme[p]))
        if (any(splt %in% missing.allele)) {
          p.m[p,] <- NA_real_
          next
        }
        if (any(!grepl("^[0-9]+$",splt)))
          stop(paste("Invalid genotype", gt, "at marker", m.nme[p]))
        p.m[p,] <- as.numeric(splt)
      }
      image(
            x = b.chrom[[c]]    ,
            y = c(0,1,2)        ,
            z = p.m             ,
            axes = 0            ,
            zlim=c(z.min,z.max),
            col = clr
            )
      box()
      
    }
}


"st.plot.ggt" <- function (data.set = "default", ifilename = "", i.list = "", 
    f.ind = 0, l.ind = 0, f.lr = 2, f.tb = 2, d.h = 0.001, d.v = 0.001, 
    p.t = FALSE, p.s = FALSE, d.t = 50, d.map = 100, c.nme = "", 
    i.nme = "", cex.chrom = 1, cex.ind = 1, color = c("yellow", 
        "blue", "red", "green", "wheat", "skyblue", "tomato", 
        "palegreen", "yellow4", "darkblue", "darkblue", "darkgreen"), 
    z.min = 0, z.max = 12, sort.performance = FALSE, show.performance = FALSE, 
    performance.decreasing = TRUE, performance.format = "%4.2f") 
{
    st.input.dir <- get0("st.input.dir", envir = .GlobalEnv, 
        inherits = FALSE, ifnotfound = "")
    st.start.timer()
    if ("" == st.input.dir) 
        ip <- ""
    else ip <- "/"
    if (!is.logical(sort.performance) || (1L != length(sort.performance)) || 
        is.na(sort.performance)) 
        stop("sort.performance must be TRUE or FALSE")
    if (!is.logical(show.performance) || (1L != length(show.performance)) || 
        is.na(show.performance)) 
        stop("show.performance must be TRUE or FALSE")
    if (!is.logical(performance.decreasing) || (1L != length(performance.decreasing)) || 
        is.na(performance.decreasing)) 
        stop("performance.decreasing must be TRUE or FALSE")
    if (!is.character(performance.format) || (1L != length(performance.format)) || 
        is.na(performance.format) || !nzchar(performance.format)) 
        stop("performance.format must be one non-empty character string")
    performance.data <- NULL
    if (sort.performance || show.performance) {
        performance.data <- st.return.performance.data(data.set = data.set)
        if (is.null(performance.data)) stop("No performance data available")
    }
    tmp.ifilename <- ""
    i.list.supplied <- !((1L == length(i.list)) && ("" == i.list[1]))
    if (i.list.supplied) {
        i.list <- as.character(i.list)
    }
    else if ("" != ifilename) {
        i.list <- read.table(paste(st.input.dir, ip, ifilename, 
            sep = ""), stringsAsFactors = FALSE)[, 1]
        i.list <- as.character(i.list)
        i.list.supplied <- TRUE
    }
    if (sort.performance) {
        if (i.list.supplied) 
            sorted.ind <- i.list
        else sorted.ind <- performance.data$i
        custom.i.nme <- !((1L == length(i.nme)) && ("" == i.nme[1]))
        sorted.i.nme <- NULL
        if (custom.i.nme) {
            if (length(i.nme) != length(sorted.ind)) 
                stop(paste("When sort.performance=TRUE, i.nme must contain one", 
                  "label for each individual before sorting"))
            sorted.i.nme <- i.nme
        }
        idx <- match(sorted.ind, performance.data$i)
        missing.ind <- unique(sorted.ind[is.na(idx)])
        if (0 < length(missing.ind)) 
            stop(paste("Individual(s) not found in data set '", 
                data.set, "': ", paste(missing.ind, collapse = ", "), 
                sep = ""))
        performance <- performance.data$y[idx]
        if (performance.decreasing) 
            ord <- order(-performance, seq_along(performance))
        else ord <- order(performance, seq_along(performance))
        sorted.ind <- sorted.ind[ord]
        if (custom.i.nme) 
            sorted.i.nme <- sorted.i.nme[ord]
        n.ind <- length(sorted.ind)
        f <- as.integer(f.ind)[1L]
        l <- as.integer(l.ind)[1L]
        if (is.na(f) || is.na(l)) 
            stop("f.ind and l.ind must be integers")
        if (n.ind < f) 
            f <- n.ind
        if (n.ind < l) 
            l <- n.ind
        if (0 < f) 
            f <- f - 1L
        else if (0 > f) 
            f <- 0L
        if ((1 > l) || (n.ind < l)) 
            l <- n.ind
        if (f < l) {
            keep <- seq.int(f + 1L, l)
            sorted.ind <- sorted.ind[keep]
            if (custom.i.nme) 
                sorted.i.nme <- sorted.i.nme[keep]
        }
        else {
            sorted.ind <- character(0)
            if (custom.i.nme) 
                sorted.i.nme <- character(0)
        }
        if (0L == length(sorted.ind)) 
            stop("No individuals selected for plotting")
        i.list <- sorted.ind
        if (custom.i.nme) 
            i.nme <- sorted.i.nme
        i.list.supplied <- TRUE
        f.ind <- 0
        l.ind <- 0
    }
    buffers <- .Call("st_legacy_marker", list(action="call", routine="gs_write_marker_data_02_GV",
        arguments=list("m", "", NULL, NULL, "", "", "", "", as.integer(f.ind), as.integer(l.ind), "", 1L, 0L, as.character(data.set)),
        retval=13L, list=TRUE, inputs=if(i.list.supplied) list(list(argument=8L, data=i.list)) else NULL,
        outputs=list(list(argument=3L), list(argument=4L, header=NA, check.names=FALSE))), PACKAGE="SelectionTools")
    if (is.null(buffers)) return(invisible(NULL))
    m <- buffers[[1L]]
    p <- buffers[[2L]]
    plotted.ind <- colnames(p)
    if (i.list.supplied) {
        missing.ind <- unique(i.list[!i.list %in% plotted.ind])
        if (0 < length(missing.ind)) 
            stop(paste("Individual(s) not found in data set '", 
                data.set, "': ", paste(missing.ind, collapse = ", "), 
                sep = ""))
        if (!identical(i.list, plotted.ind)) 
            stop("Requested individual list does not match plotted individuals")
    }
    if (show.performance && ((1L == length(i.nme)) && ("" == 
        i.nme[1]))) {
        idx <- match(plotted.ind, performance.data$i)
        if (anyNA(idx)) 
            stop("Performance data do not match plotted individuals")
        name.width <- max(nchar(plotted.ind, type = "width"))
        i.nme <- paste0(sprintf("%*s", name.width, plotted.ind), 
            ": y=", sprintf(performance.format, performance.data$y[idx]))
    }
    if (!((1 == length(i.nme)) && ("" == i.nme[1])) && (length(i.nme) != 
        ncol(p))) 
        stop("Number of individual names in i.nme does not match plotted individuals")
    n0 <- min(table(m[, 1]))
    if (1 == n0) {
        st.info(-2, "At least two loci are required on each chromosome ")
        return()
    }
    n1 <- max(m[, 1])
    n2 <- length(unique(m[, 1]))
    if (n1 != n2) {
        st.info(-2, "Chromsosomes must be numbered consecutively starting from 1 ")
        return()
    }
    st.plot.ggt.src(map = m, pop = p, f.lr, f.tb, d.h, d.v, p.t, 
        p.s, d.t, d.map, c.nme, i.nme, color, z.min, z.max, cex.chrom, 
        cex.ind)
    st.stop.timer(info.level = 1)
}



"st.marker.data.statistics" <- function (data.set = "default", filename = "marker.stats", ind = TRUE, 
    mar = TRUE, gen = TRUE, auxfiles = TRUE) 
{
    st.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_marker_stats_01_GV", arguments = list(NULL, 
        NULL, NULL, as.integer(ind), as.integer(mar), as.integer(gen), as.integer(auxfiles), 0L, as.character(data.set)), 
        retval = 8L, enabled = auxfiles, list = TRUE, outputs = list(list(argument = 1L, header = TRUE, 
            check.names = FALSE, enabled = ind), list(argument = 2L, header = TRUE, check.names = FALSE, 
            enabled = mar), list(argument = 3L, header = TRUE, check.names = FALSE, enabled = gen)), 
        output.names = c("individual.list", "marker.list", "genotypes")), PACKAGE = "SelectionTools")
    st.stop.timer(info.level = 1)
    invisible(result)
}


"gs.return.effects" <- function(out.filename = "effects", data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_return_effects_01_GV", arguments = list(NULL, 
        0L, as.character(data.set)), retval = 2L, enabled = TRUE, outputs = list(list(header = TRUE, 
        argument = structure(1L, names = ""), matrix = 0L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    result
}

"gs.set.effects" <- function (eff, filename = "seteffects", data.set = "default") 
{
    gs.start.timer()
    .Call("st_legacy_marker", list(action = "call", routine = "gs_set_effects_01_GV", arguments = list(NULL, 
        0L, as.character(data.set)), retval = 2L, inputs = list(list(argument = 1L, data = eff, header = TRUE))), 
        PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
}

"gs.return.pvals" <- function(out.filename = "effects", data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_return_pvals_01_GV", arguments = list(NULL, 
        0L, as.character(data.set)), retval = 2L, enabled = TRUE, outputs = list(list(header = TRUE, 
        argument = structure(1L, names = ""), matrix = 0L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.check.pvals" <- function( data.set ="default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_check_pvals_01_GV",
        arguments=list(retval = as.integer(retval), as.character(data.set)),
        result=1L,
        logical=TRUE), PACKAGE="SelectionTools")
  gs.stop.timer (info.level=1)
  return( c )
}

"gs.plot.effects" <- function (data.set = "default", absolute = TRUE, alpha = 0.05, 
    p.adjust.method = "none", pvals = NULL, xlim = NULL, ylim = NULL, 
    plt = TRUE, col.sp = "blue", pch.sp = 19, col.sn = "blue", 
    pch.sn = 19, col.ns = "blue", pch.ns = 1, ...) 
{
    gs.start.timer()
    {
        pvals.av <- gs.check.pvals(data.set)
        x <- if (pvals.av) gs.return.pvals(data.set=data.set) else gs.return.effects(data.set=data.set)
        if (is.null(x)) return(invisible(NULL))
        effect <- x$effect[-1L]
        if (pvals.av) {
            p <- if ("p" %in% names(x)) x$p else x[["pvalue"]]
            p <- p[-1L]
        }
        if (!is.null(pvals)) { p <- pvals; pvals.av <- TRUE }
        mp <- st.get.map(data.set=data.set)
        if (is.null(mp)) return(invisible(NULL))
        names(mp) <- c("chrom", "pos", "name", "class")
        mp$pos <- mp$pos / 100
        cl <- as.numeric(tapply(mp$pos, mp$chrom, function(x) as.numeric(sprintf("%.4f", .01 + x[length(x)]))))
        if (TRUE == absolute) 
            mp$effect <- abs(effect)
        else mp$effect <- effect
        if (TRUE == pvals.av) {
            mp$p <- p
            mp$p.adjusted <- p.adjust(mp$p, p.adjust.method)
        }
        c.pos <- cl
        if (length(c.pos) > 1L)
            for (i in seq.int(2L, length(c.pos))) c.pos[i] <- c.pos[i] + c.pos[i - 1L]
        for (i in seq_len(nrow(mp))) if (mp$chrom[i] > 1) 
            mp$pos[i] <- mp$pos[i] + c.pos[mp$chrom[i] - 1]
    }
    if (TRUE == plt) {
        oldpar <- graphics::par(no.readonly = TRUE)
        on.exit(.st.restore.graphics(oldpar), add = TRUE)
        par(mar = c(4, 4, 1, 1))
        if (is.null(ylim)) 
            ylim = c(min(c(min(mp$effect), 0)), max(mp$effect))
        if (is.null(xlim)) 
            xlim = c(min(c(min(mp$pos), 0)), max(mp$pos))
        plot(mp$pos, mp$effect, type = "n", xaxt = "n", xlab = "Chromosome", 
            ylab = "Effect size", ylim = ylim, xlim = xlim, ...)
        c.pos <- c(0, c.pos)
        axis(1, at = c.pos, labels = c(1:length(cl), " "), tck = -0.01)
        grid()
        for (i in 1:(length(c.pos))) {
            x <- c.pos[i]
            lines(c(x, x), ylim)
        }
        if (FALSE == absolute) 
            lines(xlim, c(0, 0))
        tol = 1e-08
        if (TRUE == pvals.av) {
            mp1 <- mp[mp$p.adjusted < alpha, ]
            mp1ap <- mp1[mp1$effect > tol, ]
            mp1an <- mp1[mp1$effect < -tol, ]
            mp1b <- mp1[abs(mp1$effect) <= tol, ]
            mp2 <- mp[mp$p.adjusted >= alpha, ]
            points(mp1b$pos, mp1b$effect, type = "p", lwd = 2, 
                col = col.ns, pch = pch.ns)
            points(mp2$pos, mp2$effect, type = "p", lwd = 2, 
                col = col.ns, pch = pch.ns)
            points(mp1ap$pos, mp1ap$effect, type = "p", lwd = 2, 
                col = col.sp, pch = pch.sp)
            points(mp1an$pos, mp1an$effect, type = "p", lwd = 2, 
                col = col.sn, pch = pch.sn)
        }
        else {
            points(mp$pos, mp$effect, type = "p", lwd = 2, col = "blue")
        }
    }
    ret.val <- data.frame(chrom = mp$chrom, pos = mp$pos, name = mp$name, 
        cpos = mp$pos, effect = mp$effect)
    if (TRUE == pvals.av) 
        ret.val$p <- mp$p
    gs.stop.timer(info.level = 1)
    return(invisible(ret.val))
}


"gs.compare.effects" <- function (data.set.a, data.set.b, label.a = data.set.a, label.b = data.set.b) 
{
    gs.start.timer()
    x <- gs.return.effects(data.set = data.set.a)
    if (is.null(x)) return(invisible(NULL))
    effects.a <- x$effect[-1L]
    x <- gs.return.effects(data.set = data.set.b)
    if (is.null(x)) return(invisible(NULL))
    effects.b <- x$effect[-1L]
    max <- max(c(abs(effects.a), abs(effects.b)))
    min <- -max
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(.st.restore.graphics(oldpar), add = TRUE)
    par(mar = c(4, 4, 1, 1))
    plot(effects.b ~ effects.a, type = "p", lwd = 2, col = "blue", 
        ylim = c(min, max), xlim = c(min, max), xlab = label.a, 
        ylab = label.b)
    lines(c(min, max), c(min, max))
    gs.stop.timer(info.level = 1)
    return(invisible(NULL))
}


gs.plot.model.fit <- function ( training.set ="default",
                                title        = ""         )
{
  c1 <- gs.predict.genotypes ( training.set   = training.set,
                               prediction.set = training.set )
  if ( is.null(c1) ) return (invisible(NULL) )
  
  c1 <- c1[,2:3]
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(.st.restore.graphics(oldpar), add = TRUE)
  par(mar=c(4,4,1,1) )
  
  l <- c(min(c1),max(c1))
  p1 <- l[1] + 0.65*(l[2]-l[1] ) ; p2 <- l[1] + 0.10*(l[2]-l[1] )
  p3 <- l[1] + 0.05*(l[2]-l[1] ) ; p4 <- l[1] + 0.90*(l[2]-l[1] )
  
  plot(c1,ylim=l,xlim=l,col="blue")
  text(p1,p2,sprintf("r = %4.2f\nrho = %4.2f",cor(c1)[2],
                   cor(c1,method="spearman")[2]),pos=4)
  text(p3,p4,sprintf("%s\nModel fit",title),pos=4)
  lines(l,l)
  lines ( c(mean(c1[,1]),mean(c1[,1])), c(min(c1[,2]),max(c1[,2])) )
  lines ( c(min(c1[,1]),max(c1[,1])), c(mean(c1[,2]),mean(c1[,2])) )
}

gs.plot.validation <- function ( training.set ,
                                 prediction.set ,
                                 title=""          )
{
  
  c1 <- gs.predict.genotypes ( training.set   = training.set,
                               prediction.set = training.set )
  c2 <- gs.predict.genotypes ( training.set   = training.set,
                               prediction.set = prediction.set )
  
  if ( 3 == ncol(c1) ) c1 <- c1[,2:3] else return (invisible(NULL) )
  if ( 3 == ncol(c2) ) c2 <- c2[,2:3] else
   {
     st.info (-2,"No phenotypic data in the prediction set")
     return (invisible(NULL) )
   }
  
  if ( is.null(c1) || is.null(c1) ) return (invisible(NULL) )
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(.st.restore.graphics(oldpar), add = TRUE)
  par( mfcol=c(1,2) , mar=c(4,4,1,1) )
  
  l <- c( min(c1,c2),max(c1,c2) )
  
  p1 <- l[1] + 0.65*(l[2]-l[1] ) ; p2 <- l[1] + 0.10*(l[2]-l[1] )
  p3 <- l[1] + 0.05*(l[2]-l[1] ) ; p4 <- l[1] + 0.90*(l[2]-l[1] )
  
  plot(c1,ylim=l,xlim=l,col="blue")
  text(p1,p2,sprintf("r = %4.2f\nrho = %4.2f",cor(c1)[2],
                     cor(c1,method="spearman")[2]),pos=4)
  text(p3,p4,sprintf("%s\nModel fit",title),pos=4)
  lines(l,l)
  lines ( c(mean(c1[,1]),mean(c1[,1])), c(min(c1[,2]),max(c1[,2])) )
  lines ( c(min(c1[,1]),max(c1[,1])), c(mean(c1[,2]),mean(c1[,2])) )
  
  plot(c2,ylim=l,xlim=l,col="blue")
  text(p1,p2,sprintf("r = %4.2f\nrho = %4.2f",cor(c2)[2],
                     cor(c2,method="spearman")[2]),pos=4)
  text(p3,p4,sprintf("%s\nValidation",title),pos=4)
  lines(l,l)
  lines ( c(mean(c2[,1]),mean(c2[,1])), c(min(c2[,2]),max(c2[,2])) )
  lines ( c(min(c2[,1]),max(c2[,1])), c(mean(c2[,2]),mean(c2[,2])) )
}

st.plot.corr <- function ( x , y,
                          title        = "" ,
                          pch          = 16 ,
                          cex          = 1  ,
                          xlab         = "x",
                          ylab         = "y"  )
{
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(.st.restore.graphics(oldpar), add = TRUE)
  par ( mar=c(4,4,1,1) )
  
  lx <- c(min(x),max(x))
  ly <- c(min(y),max(y))
  p1 <- lx[1] + 0.65*(lx[2]-lx[1] ) ; p2 <- ly[1] + 0.10*(ly[2]-ly[1] )
  p3 <- lx[1] + 0.05*(lx[2]-lx[1] ) ; p4 <- ly[1] + 0.90*(ly[2]-ly[1] )
  
  plot(x,y,xlab=xlab,ylab=ylab,type="n")
  points(x,y,pch=pch,cex=cex);
  text(p1,p2,sprintf("r = %4.2f\nrho = %4.2f",cor(x,y),
                   cor(x,y,method="spearman")),pos=4)

  text(p3,p4,sprintf("%s",title),pos=4)
  lines ( c(lx[1],lx[2]) , c(ly[1],ly[2]))
  lines ( c(mean(x),mean(x)), c(min(y),max(y)) )
  lines ( c(min(x),max(x)), c(mean(y),mean(y)) )
}

st.plot.corr.l <- function ( x , y,
                          title        = "" ,
                          pch          = 16 ,
                          cex          = 1  ,
                          xlab         = "x",
                          ylab         = "y",
                          ylim,
                          xlim )
{
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(.st.restore.graphics(oldpar), add = TRUE)
  par ( mar=c(4,4,1,1) )
  
  lx <- c(min(x),max(x))
  ly <- c(min(y),max(y))
  p1 <- lx[1] + 0.65*(lx[2]-lx[1] ) ; p2 <- ly[1] + 0.10*(ly[2]-ly[1] )
  p3 <- lx[1] + 0.05*(lx[2]-lx[1] ) ; p4 <- ly[1] + 0.90*(ly[2]-ly[1] )
  
  plot(x,y,xlab=xlab,ylab=ylab,type="n",xlim=xlim,ylim=ylim)
  points(x,y,pch=pch,cex=cex);
  text(p1,p2,sprintf("r = %4.2f\nrho = %4.2f",cor(x,y),
                   cor(x,y,method="spearman")),pos=4)

  text(p3,p4,sprintf("%s",title),pos=4)
  lines ( c(lx[1],lx[2]) , c(ly[1],ly[2]))
  lines ( c(mean(x),mean(x)), c(min(y),max(y)) )
  lines ( c(min(x),max(x)), c(mean(y),mean(y)) )
}


st.genetic.distances.02 <- function ( measure = "mrd", 
                                      format  = "l",   
                                      filename  = NULL,
                                      auxfiles  = FALSE,
                                      data.set  = "default" )
{
  st.output.dir <- get0("st.output.dir", envir = .GlobalEnv,
                        inherits = FALSE, ifnotfound = "")

   if (TRUE == auxfiles && is.null(filename))
     stop("filename must be supplied when auxfiles = TRUE")
   out.filename <- NULL
   if (TRUE == auxfiles) {
     if ( "" == st.output.dir ) op <- "" else op <- "/"
     out.filename <- paste(st.output.dir,op,filename,".gdi",sep="")
   }
  st.start.timer()
  distance <- .Call("st_legacy_sim", list(action="marker.distance",
      data.set=as.character(data.set), measure=measure, format=format,
      auxiliary=auxfiles), PACKAGE="SelectionTools")
  st.stop.timer(info.level=1)
  if  (auxfiles)
    {
      write.table ( distance$table,
                    file=out.filename,
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE        )
      distance <- distance$result
    }

  st.stop.timer (info.level=1)
  invisible(distance)
}

st.genetic.distances <- function ( measure = "mrd", 
                                   format  = "l",   
                                   filename  = NULL,
                                   auxfiles  = FALSE,
                                   data.set  = "default" )
{
  st.output.dir <- get0("st.output.dir", envir = .GlobalEnv,
                        inherits = FALSE, ifnotfound = "")

   if (TRUE == auxfiles && is.null(filename))
     stop("filename must be supplied when auxfiles = TRUE")
   out.filename <- NULL
   if (TRUE == auxfiles) {
     if ( "" == st.output.dir ) op <- "" else op <- "/"
     out.filename <- paste(st.output.dir,op,filename,".gdi",sep="")
   }

   gs.cross.eval.gd ( dist=measure , data.set=data.set )
   c <- gs.cross.info ( data.set=data.set )
   distance <- data.frame (OTU1=c$P1Name,OTU2=c$P2Name,Measure=c$gd)


  if  (auxfiles)
    {
      write.table ( distance,
                    file=out.filename,
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE        )
    }

  st.stop.timer (info.level=1)
   if ("m" == format)
     {
       c <- st.marker.data.statistics ( ind=TRUE,
                                        mar=FALSE,
                                        gen=FALSE,
                                        auxfiles=TRUE,
                                        data.set=data.set)
       n <- c$individual.list$Name
       aa <- distance[,3]
       bb <- length(n)
       cc <- n
       dd <- structure(aa, Size = bb, Labels = cc, Diag = FALSE,
                       Upper = TRUE, method = "euclidean", class = "dist")
       return ( invisible(dd) )
     }
   else
     {
       return ( invisible(distance) )
     }
}


"st.indir" <- function ( filename )
{
  st.input.dir <- get0("st.input.dir", envir = .GlobalEnv,
                       inherits = FALSE, ifnotfound = "")

  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  paste(st.input.dir,ip,filename,sep="")
}

"st.outdir" <- function ( filename )
{
  st.output.dir <- get0("st.output.dir", envir = .GlobalEnv,
                        inherits = FALSE, ifnotfound = "")

  if ( "" == st.output.dir ) op <- "" else op <- "/"
  paste(st.output.dir,op,filename,sep="")
}

"st.datadir" <- function ( filename )
{
  st.data.dir <- get0("st.data.dir", envir = .GlobalEnv,
                      inherits = FALSE, ifnotfound = "")

   if ( "" == st.data.dir ) dp <- "" else dp <- "/"
   paste(st.data.dir,dp,filename,sep="")
}







"st.id" <- function ( filename )
{
  st.input.dir <- get0("st.input.dir", envir = .GlobalEnv,
                       inherits = FALSE, ifnotfound = "")

  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  paste(st.input.dir,ip,filename,sep="")
}

"st.od" <- function ( filename )
{
  st.output.dir <- get0("st.output.dir", envir = .GlobalEnv,
                        inherits = FALSE, ifnotfound = "")

  if ( "" == st.output.dir ) op <- "" else op <- "/"
  paste(st.output.dir,op,filename,sep="")
}

"st.dd" <- function ( filename )
{
  st.data.dir <- get0("st.data.dir", envir = .GlobalEnv,
                      inherits = FALSE, ifnotfound = "")

   if ( "" == st.data.dir ) dp <- "" else dp <- "/"
   paste(st.data.dir,dp,filename,sep="")
}








"gs.cross.eval.gd" <- function ( dist     = "rd",
                                 data.set = "default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_cross_eval_gd_01_GV",
        arguments=list(as.character(dist), retval = as.integer(retval), as.character(data.set)),
        result=0L), PACKAGE="SelectionTools")
  gs.stop.timer (info.level=1)
}


"gs.cross.info.gd" <- function(out.filename = "cross.info.gd", data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_cross_info_gd_01_GV", arguments = list(NULL, 
        0L, as.character(data.set)), retval = 2L, enabled = TRUE, outputs = list(list(header = TRUE, 
        argument = structure(1L, names = ""), matrix = 0L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}


"gs.cross.eval.mu" <- function ( data.set = "default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_cross_eval_mu_01_GV",
        arguments=list(retval = as.integer(retval), as.character(data.set)),
        result=0L), PACKAGE="SelectionTools")
  gs.stop.timer (info.level=1)
}

"gs.cross.eval.mi" <- function ( data.set = "default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_cross_eval_mi_01_GV",
        arguments=list(retval = as.integer(retval), as.character(data.set)),
        result=0L), PACKAGE="SelectionTools")
  gs.stop.timer (info.level=1)
}

"gs.cross.eval.ma" <- function ( data.set = "default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_cross_eval_ma_01_GV",
        arguments=list(retval = as.integer(retval), as.character(data.set)),
        result=0L), PACKAGE="SelectionTools")
  gs.stop.timer (info.level=1)
}


"gs.cross.eval.va" <- function ( data.set     = "default"  ,
                                 pop.type     = "unlinked" ,
                                 t            = 0          ,
                                 map.function = "Haldane"    )
{
  gs.start.timer()
  .Call("st_legacy_marker", list(action="cross.va", dataset=as.character(data.set),
        pop.type=as.character(pop.type), t=as.double(t), map.function=as.character(map.function)),
        PACKAGE="SelectionTools")
  gs.stop.timer (info.level=1)
}

"gs.cross.eval.es" <- function (data.set = "default", alpha = 0.1, N = 0, G = 0) 
{
    gs.start.timer()
    .Call("st_legacy_marker", list(action = "cross.es", dataset = as.character(data.set), alpha = as.double(alpha), 
        N = as.double(N), G = as.double(G)), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
}


st.cross.allocate <- function(plan, scores, total.progeny, n.selected = 10,
                              min.progeny = n.selected, max.progeny = total.progeny,
                              accuracy = 1, control = list(),
                              selection.method = c("burrows", "order.statistics"),
                              constraints = NULL, genetic.model = NULL) {
    selection.method <- match.arg(selection.method)
    .Call("st_legacy_marker", list(action = "cross.allocate", plan = plan,
        scores = scores, total.progeny = total.progeny, n.selected = n.selected,
        min.progeny = min.progeny, max.progeny = max.progeny,
        accuracy = accuracy, control = control, selection.method = selection.method,
        constraints = constraints, genetic.model = genetic.model),
        PACKAGE = "SelectionTools")
}

st.cross.evaluate <- function(plan, scores = NULL, n.selected = NULL,
                               accuracy = NULL, K = NULL, selection.method = NULL,
                               constraints = NULL, genetic.model = NULL) {
    .Call("st_legacy_marker", list(action = "cross.evaluate", plan = plan,
        scores = scores, n.selected = n.selected, accuracy = accuracy, K = K,
        selection.method = selection.method, constraints = constraints,
        genetic.model = genetic.model),
        PACKAGE = "SelectionTools")
}

print.st.cross.population <- function(x, digits = NULL, top = NULL, bottom = NULL,
                                      char.per.line = NULL, ...) {
    .Call("st_legacy_marker", list(action = "cross.evaluate.print", object = x,
        options = st.options(), digits = digits, top = top, bottom = bottom,
        char.per.line = char.per.line), PACKAGE = "SelectionTools")
    invisible(x)
}

print.st.cross.allocation <- function(x, digits = NULL, top = NULL, bottom = NULL,
                                      char.per.line = NULL, ...) {
    .Call("st_legacy_marker", list(action = "cross.allocate.print", object = x,
        options = st.options(), digits = digits, top = top, bottom = bottom,
        char.per.line = char.per.line), PACKAGE = "SelectionTools")
    invisible(x)
}

"gs.cross.info" <- function(bestn = 0, sortby = "index", out.filename = "cross.info", data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "cross.info", dataset = as.character(data.set),
        bestn = as.double(bestn), sortby = as.character(sortby)), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"st.chrom.stats" <- function(filename = "chrom.stats", auxfiles = TRUE, data.set = "default") {
    st.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_chrom_stats_01_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L, as.character(data.set)), retval = 3L, enabled = auxfiles, outputs = list(list(header = FALSE, 
        argument = structure(1L, names = ""), matrix = 0L, keep.names = TRUE, names = c("Chrom", "NLoci", 
            "Length")))), PACKAGE = "SelectionTools")
    st.stop.timer(info.level = 1)
    invisible(result)
}

"st.calc.rf" <- function(map.function = "Haldane", auxfiles = TRUE, data.set = "default") {
    st.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_calc_rf_01_GV", arguments = list(as.character(map.function), 
        NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 4L, enabled = auxfiles, outputs = list(list(header = FALSE, 
        argument = structure(2L, names = ""), matrix = 0L, keep.names = TRUE, names = c("Chrom", "Locus1", 
            "Locus2", "RecFreq")))), PACKAGE = "SelectionTools")
    st.stop.timer(info.level = 1)
    invisible(result)
}

"st.calc.q" <- function(pop.type = "DH", t = 0, auxfiles = TRUE, data.set = "default") {
    st.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_calc_q_01_GV", arguments = list(as.character(pop.type), 
        as.integer(t), NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 5L, enabled = auxfiles, 
        outputs = list(list(header = FALSE, argument = structure(3L, names = ""), matrix = 0L, keep.names = TRUE, 
            names = c("Chrom", "Locus1", "Locus2", "q")))), PACKAGE = "SelectionTools")
    st.stop.timer(info.level = 1)
    invisible(result)
}

"st.get.simpop" <- function ( pop.name ,
                              data.set ="default" )
  {
   if ((1L != length(pop.name)) || is.na(pop.name) ||
       nchar(pop.name, type="bytes") < 2L ||
       nchar(pop.name, type="bytes") >= 256L ||
       !grepl("^[A-Za-z0-9]+$", pop.name))
      stop("pop.name must contain 2 to 255 ASCII letters or digits")
   bg <- plabsim$allele.missing.indicator
   retval <- 0
   c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="st_get_simpop_GV",
        arguments=list(as.character        (pop.name), as.integer          (bg), retval = as.integer (retval), as.character        (data.set))), PACKAGE="SelectionTools")
}



"st.get.simpop2" <- function ( pop.name ,
                              data.set ="default" )
  {
   if ((1L != length(pop.name)) || is.na(pop.name) ||
       nchar(pop.name, type="bytes") < 2L ||
       nchar(pop.name, type="bytes") >= 256L ||
       !grepl("^[A-Za-z0-9]+$", pop.name))
      stop("pop.name must contain 2 to 255 ASCII letters or digits")
   bg <- plabsim$allele.missing.indicator
   retval <- 0
   c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="st_get_simpop2_GV",
        arguments=list(as.character        (pop.name), as.integer          (bg), retval = as.integer (retval), as.character        (data.set))), PACKAGE="SelectionTools")
}





il.eval.lines <- function(allele=1, chrom=0, lower=0, upper=0, exclude=0,
                          data.set="default") {
  if (length(allele) != 1L || is.na(allele) || allele != as.integer(allele) ||
      allele < 1L || allele > .GS_ALLELE_MAX)
    stop("allele must be an integer from 1 to 254")
  .Call("st_legacy_marker", list(action="il.lines", source=data.set,
        allele=as.integer(allele), chromosome=as.integer(chrom),
        lower=as.double(lower), upper=as.double(upper), exclude=as.integer(exclude)),
        PACKAGE="SelectionTools")
}


il.ideal.library <- function(n.c=5, l.c=100, s.l=20, data.set="default")
    invisible(.Call("st_legacy_marker", list(action="il.create", source=data.set,
          chromosomes=n.c, length=l.c, segment=s.l, overlap=FALSE),
          PACKAGE="SelectionTools"))

il.overlapping.library <- function(n.c=5, l.c=100, s.l=20, data.set="default")
    invisible(.Call("st_legacy_marker", list(action="il.create", source=data.set,
          chromosomes=n.c, length=l.c, segment=s.l, overlap=TRUE),
          PACKAGE="SelectionTools"))


"il.eval.library" <- function ( rp.allele = 1,
                                dp.allele = 2,
                                chrom     = 0,
                                lower     = 0,
                                upper     = 0,
                                data.set  = "default"  )
{
  if ((length(rp.allele) != 1L) || is.na(rp.allele) ||
      (rp.allele != as.integer(rp.allele)) ||
      (rp.allele < 1L) || (rp.allele > .GS_ALLELE_MAX) ||
      (length(dp.allele) != 1L) || is.na(dp.allele) ||
      (dp.allele != as.integer(dp.allele)) ||
      (dp.allele < 1L) || (dp.allele > .GS_ALLELE_MAX))
    stop("rp.allele and dp.allele must be integers from 1 to 254")

  .Call("st_legacy_marker", list(action="il.summary", source=data.set,
        donor=as.integer(dp.allele), recipient=as.integer(rp.allele),
        chromosome=as.integer(chrom), lower=as.double(lower),
        upper=as.double(upper)), PACKAGE="SelectionTools")

}

  st.select.phen <- function(p=NULL, n=0, t=0, decreasing=TRUE,
                           nsmall=1, data.set=NULL) {
  if (is.null(p) && is.null(data.set))
    stop("either p or data.set must be supplied")
  p <- .Call("st_legacy_marker", list(action="select.phen", data=p, source=data.set,
             n=n, threshold=t, decreasing=decreasing), PACKAGE="SelectionTools")
  i.col <- if ("i" %in% names(p)) "i" else "ii"
  y.col <- if ("y" %in% names(p)) "y" else "yy"
  p$descr <- paste(p[[i.col]], format(p[[y.col]], nsmall=nsmall))
  p
}

st.mark.alleles <- function(markers=NULL, data.set="default", effect.set=NULL,
                            target.set, hap.lst=NULL, alpha=0.05,
                            p.adjust.method="none", replace="3")
    invisible(.Call("st_legacy_marker", list(action="mark.alleles",
                    markers=markers, source=data.set, effect=effect.set,
                    target=target.set, haplotypes=hap.lst, alpha=alpha,
                    adjust=p.adjust.method, replace=as.character(replace)),
                    PACKAGE="SelectionTools"))


gs.write.pseff <- function(e, beta0=NA, file) {
  e <- as.data.frame(e)
  output <- .Call("st_legacy_marker", list(action="pseff", data=e, intercept=beta0),
                  PACKAGE="SelectionTools")
  write(output[[1L]], file)
  write.table(output[[2L]], file, quote=FALSE, row.names=FALSE,
              col.names=FALSE, append=TRUE)
}

"gs.neg.effall" <- function( data.set ="default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_neg_effall_01_GV",
        arguments=list(retval = as.integer(retval), as.character(data.set)),
        result=0L), PACKAGE="SelectionTools")
  gs.stop.timer (info.level=1)
}

"gs.pos.effall" <- function( data.set ="default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_pos_effall_01_GV",
        arguments=list(retval = as.integer(retval), as.character(data.set)),
        result=0L), PACKAGE="SelectionTools")
  gs.stop.timer (info.level=1)
}





"gs.mmet.coeff" <- function(out.filename = "mme.coeff", auxfiles = FALSE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_mmet_coeff_01_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L, as.character(data.set)), retval = 3L, enabled = auxfiles, outputs = list(list(header = FALSE, 
        argument = structure(1L, names = ""), matrix = 1L, keep.names = FALSE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.mmet.coeff.3" <- function(out.filename = "mme.coeff", auxfiles = FALSE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_mmet_coeff_03_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L, as.character(data.set)), retval = 3L, enabled = auxfiles, outputs = list(list(header = FALSE, 
        argument = structure(1L, names = ""), matrix = 1L, keep.names = FALSE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}


"gs.mmet.rhs" <- function(out.filename = "mme.rhs", auxfiles = FALSE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_mmet_rhs_01_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L, as.character(data.set)), retval = 3L, enabled = auxfiles, outputs = list(list(header = FALSE, 
        argument = structure(1L, names = ""), matrix = 1L, keep.names = FALSE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.mmet.solve" <- function(out.filename = "mme.solution", auxfiles = FALSE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_mmet_solve_01_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L, as.character(data.set)), retval = 3L, enabled = auxfiles, outputs = list(list(header = TRUE, 
        argument = structure(1L, names = ""), matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}


"gs.lambda.rmla.02" <- function(alpha = 1, maxiter = 1000, precision = 0.001, hsq = 0.9, out.filename = "lambda", 
    auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_lambda_rmla_02_GV", arguments = list(as.double(alpha), 
        as.integer(maxiter), as.double(precision), as.double(hsq), NULL, as.integer(auxfiles), 0L, as.character(data.set)), 
        retval = 7L, enabled = auxfiles, outputs = list(list(header = TRUE, argument = structure(5L, names = ""), 
            matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.w.aov" <- function(alpha = 1, out.filename = "W", auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_w_aov_01_GV", arguments = list(as.double(alpha), 
        NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 4L, enabled = auxfiles, outputs = list(list(header = TRUE, 
        argument = structure(2L, names = ""), matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}


"gs.esteff.rr" <- function(method = "BLUP", maxiter = 1000, precision = 1e-04, hsq = 0.8, alpha = 1, 
    out.filename = "eff.rr", auxfiles = FALSE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_esteff_rr_01_GV", arguments = list(as.character(method), 
        as.integer(maxiter), as.double(precision), as.double(hsq), as.double(alpha), NULL, as.integer(auxfiles), 
        0L, as.character(data.set)), retval = 8L, enabled = auxfiles, outputs = list(list(header = TRUE, 
        argument = structure(6L, names = ""), matrix = 1L, keep.names = TRUE)), warn.nonconvergence = TRUE), 
        PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.cross.validation" <- function (estimation.method, n.ts, n.runs, hsq = 0.8, alpha = 1, maxiter = 100, 
    precision = 1e-04, out.filename.r = "crosval", out.filename.u = "effects", auxfiles = TRUE, data.set = "default") 
{
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_cross_validation_02_GV", 
        arguments = list(as.character(estimation.method), as.integer(n.ts), as.integer(n.runs), as.double(hsq), 
            as.double(alpha), as.integer(maxiter), as.double(precision), NULL, NULL, as.integer(auxfiles), 
            0L, as.character(data.set)), retval = 11L, enabled = auxfiles, list = TRUE, warn.nonconvergence = TRUE, 
        outputs = list(list(argument = 8L, header = TRUE, matrix = 1L), list(argument = 9L, header = TRUE, 
            matrix = 1L)), output.names = c("cor", "eff")), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}









"st.calc.ld" <- function(ld.measure = "r2", auxfiles = TRUE, data.set = "default") {
    st.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_calc_ld_01_GV", arguments = list(as.character(ld.measure), 
        NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 4L, enabled = auxfiles, outputs = list(list(header = FALSE, 
        argument = structure(2L, names = ""), matrix = 0L, keep.names = TRUE, names = c("Chrom", "Locus1", 
            "Locus2", "Name1", "Name2", "LD")))), PACKAGE = "SelectionTools")
    st.stop.timer(info.level = 1)
    invisible(result)
}


"st.calc.ld.2" <- function (ld.measure = "r2", auxfiles = TRUE, keep.ldfile = FALSE, data.set = "default", 
    ld.filename = NULL) 
{
    if (!keep.ldfile) 
        return(invisible(st.calc.ld(ld.measure = ld.measure, auxfiles = auxfiles, data.set = data.set)))
    if (auxfiles && is.null(ld.filename)) 
        stop("ld.filename must be supplied when keep.ldfile = TRUE and auxfiles = TRUE")
    st.start.timer()
    .Call("st_legacy_marker", list(action = "call", routine = "gs_calc_ld_01_GV", arguments = list(as.character(ld.measure), 
        if (auxfiles) as.character(.st.legacy.path(ld.filename, TRUE)) else "", as.integer(auxfiles), 
        0L, as.character(data.set)), retval = 4L), PACKAGE = "SelectionTools")
    st.stop.timer(info.level = 1)
    invisible(NULL)
}

st.LDplot.ld <- function (ld,chrom)
  {
    l <- ld[ld$Chrom==chrom,]
    n <- nrow(l)
    d <- (1+sqrt(1+8*n))/2
    LD <- matrix(nrow=d,ncol=d)
    for (i in seq_len(n)) LD[ l[i,2],l[i,3] ] <- l[i,6]
    return(LD)
  }

st.LDplot.map <- function(chrom, data.set)
  {
    mp <- st.get.map(data.set=data.set)
    m <- mp[mp$Chrom==chrom,]$Pos
    return(m)
  }

"st.set.hblocks" <- function (haplotype.list, hap.symbol = "s", out.filename = "blocks", auxfiles = TRUE, 
    data.set = "default") 
{
    st.start.timer()
    result <- .Call("st_legacy_marker", list(action = "hblocks.set", markers = as.character(haplotype.list$Markers), 
        symbol = as.character(hap.symbol), auxfiles = as.logical(auxfiles), dataset = as.character(data.set), 
        output = list(names = c("Chrom", "Pos", "Name", "Class", "Markers"))), PACKAGE = "SelectionTools")
    st.stop.timer(info.level = 1)
    invisible(result)
}











st.genetic.distances.fct <- function ( measure = "mrd", 
                                       format  = "l",   
                                       split   = 1,     
                                       filename  = NULL,
                                       auxfiles  = FALSE,
                                       data.set  = "default" )
{
  st.output.dir <- get0("st.output.dir", envir = .GlobalEnv,
                        inherits = FALSE, ifnotfound = "")

   if (TRUE == auxfiles && is.null(filename))
     stop("filename must be supplied when auxfiles = TRUE")
   out.filename <- NULL
   if (TRUE == auxfiles) {
     if ( "" == st.output.dir ) op <- "" else op <- "/"
     out.filename <- paste(st.output.dir,op,filename,".gdi",sep="")
   }

   gs.cross.eval.gd.fct ( dist=measure , split=split, data.set=data.set )
   c <- gs.cross.info ( data.set=data.set )
   distance <- data.frame (OTU1=c$P1Name,OTU2=c$P2Name,Measure=c$gd)


  if  (auxfiles)
    {
      write.table ( distance,
                    file=out.filename,
                    quote = FALSE,
                    row.names = FALSE,
                    col.names = FALSE        )
    }

  st.stop.timer (info.level=1)
   if ("m" == format)
     {
       c <- st.marker.data.statistics ( ind=TRUE,
                                        mar=FALSE,
                                        gen=FALSE,
                                        auxfiles=TRUE,
                                        data.set=data.set)
       n <- c$individual.list$Name
       aa <- distance[,3]
       bb <- length(n)
       cc <- n
       dd <- structure(aa, Size = bb, Labels = cc, Diag = FALSE,
                       Upper = TRUE, method = "euclidean", class = "dist")
       return ( invisible(dd) )
     }
   else
     {
       return ( invisible(distance) )
     }
}

"gs.cross.eval.gd.fct" <- function ( dist     = "rd",
                                     split    = 1   ,
                                     data.set = "default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_cross_eval_gd_fct_01_GV",
        arguments=list(as.character(dist), as.integer(split), retval = as.integer(retval), as.character(data.set)),
        result=0L), PACKAGE="SelectionTools")
  gs.stop.timer (info.level=1)
}









"st.recode.ref.2" <- function ( r1 = 0,
                                r2 = 1,
                                descending = FALSE,
                                data.set   = "default" )
{
  st.start.timer()
  retval <- 0;
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_recode_ref_02_GV",
        arguments=list(as.integer(r1), as.integer(r2), as.integer(descending), retval = as.integer(retval), as.character(data.set)),
        result=0L), PACKAGE="SelectionTools")
  st.stop.timer (info.level=1)

}










"gs.lambda.rrblup" <- function(maxiter = 1000, precision = 1e-04, hsq = 0.8, out.filename = "lambda", 
    auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_lambda_rrblup_01_GV", arguments = list(as.integer(maxiter), 
        as.double(precision), as.double(hsq), NULL, as.integer(auxfiles), 0L, as.character(data.set)), 
        retval = 6L, enabled = auxfiles, outputs = list(list(header = TRUE, argument = structure(4L, names = ""), 
            matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.vc.rrblup" <- function(maxiter = 1000, precision = 1e-04, hsq = 0.8, out.filename = "lambda", auxfiles = TRUE, 
    data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_vc_rrblup_01_GV", arguments = list(as.integer(maxiter), 
        as.double(precision), as.double(hsq), NULL, as.integer(auxfiles), 0L, as.character(data.set)), 
        retval = 6L, enabled = auxfiles, outputs = list(list(header = TRUE, argument = structure(4L, names = ""), 
            matrix = 1L, keep.names = TRUE)), warn.nonconvergence = TRUE), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}












"gs.esteff.lsq" <- function(out.filename = "eff.lsq", auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_esteff_lsq_01_GV", arguments = list(NULL, 
        as.integer(auxfiles), 0L, as.character(data.set)), retval = 3L, enabled = auxfiles, outputs = list(list(header = TRUE, 
        argument = structure(1L, names = ""), matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}


"gs.testeff.lsq" <- function(nperm = 0, out.filename = "test.lsq", auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_testeff_lsq_01_GV", arguments = list(as.integer(nperm), 
        NULL, as.integer(auxfiles), 0L, as.character(data.set)), retval = 4L, enabled = auxfiles, outputs = list(list(header = TRUE, 
        argument = structure(2L, names = ""), matrix = 0L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.esteff.rmlc" <- function(maxiter = 1000, precision = 1e-04, hsq = 0.9, out.filename = "eff.rmlc", 
    auxfiles = FALSE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_esteff_rmlc_01_GV", arguments = list(as.integer(maxiter), 
        as.double(precision), as.double(hsq), NULL, as.integer(auxfiles), 0L, as.character(data.set)), 
        retval = 6L, enabled = auxfiles, outputs = list(list(header = TRUE, argument = structure(4L, names = ""), 
            matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.testeff.rmlc" <- function(nperm = 1000, maxiter = 1000, precision = 0.001, hsq = 0.9, out.filename = "test.rmlc", 
    auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_testeff_rmlc_01_GV", arguments = list(as.integer(nperm), 
        as.integer(maxiter), as.double(precision), as.double(hsq), NULL, as.integer(auxfiles), 0L, as.character(data.set)), 
        retval = 7L, enabled = auxfiles, outputs = list(list(header = TRUE, argument = structure(5L, names = ""), 
            matrix = 0L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}


"gs.esteff.rmlv" <- function(maxiter = 1000, precision = 0.001, hsq = 0.9, out.filename = "eff.rmlv", 
    auxfiles = FALSE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_esteff_rmlv_01_GV", arguments = list(as.integer(maxiter), 
        as.double(precision), as.double(hsq), NULL, as.integer(auxfiles), 0L, as.character(data.set)), 
        retval = 6L, enabled = auxfiles, outputs = list(list(header = TRUE, argument = structure(4L, names = ""), 
            matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.testeff.rmlv" <- function(nperm = 1000, maxiter = 1000, precision = 0.001, hsq = 0.9, out.filename = "test.rmlv", 
    auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_testeff_rmlv_01_GV", arguments = list(as.integer(nperm), 
        as.integer(maxiter), as.double(precision), as.double(hsq), NULL, as.integer(auxfiles), 0L, as.character(data.set)), 
        retval = 7L, enabled = auxfiles, outputs = list(list(header = TRUE, argument = structure(5L, names = ""), 
            matrix = 0L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}


"gs.esteff.rmlr" <- function(maxiter = 1000, precision = 0.001, hsq = 0.9, out.filename = "eff.rmlr", 
    auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_esteff_rmlr_01_GV", arguments = list(as.integer(maxiter), 
        as.double(precision), as.double(hsq), NULL, as.integer(auxfiles), 0L, as.character(data.set)), 
        retval = 6L, enabled = auxfiles, outputs = list(list(header = TRUE, argument = structure(4L, names = ""), 
            matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}

"gs.esteff.rmla" <- function(alpha = 1, maxiter = 1000, precision = 1e-04, hsq = 0.9, out.filename = "eff.rmla", 
    auxfiles = TRUE, data.set = "default") {
    gs.start.timer()
    result <- .Call("st_legacy_marker", list(action = "call", routine = "gs_esteff_rmla_01_GV", arguments = list(as.double(alpha), 
        as.integer(maxiter), as.double(precision), as.double(hsq), NULL, as.integer(auxfiles), 0L, as.character(data.set)), 
        retval = 7L, enabled = auxfiles, outputs = list(list(header = TRUE, argument = structure(5L, names = ""), 
            matrix = 1L, keep.names = TRUE))), PACKAGE = "SelectionTools")
    gs.stop.timer(info.level = 1)
    invisible(result)
}









"st.set.sim.gp" <- function (data.set="default")
    {
        r <- 0
        c <-  .Call("st_legacy_marker", list(action="legacy.result",
        routine="st_set_sim_gp",
        arguments=list(r = as.integer(r), as.character(data.set))), PACKAGE="SelectionTools")
    }

"st.set.sim.mp" <- function (data.set="default")
    {
        r <- 0
        c <-  .Call("st_legacy_marker", list(action="legacy.result",
        routine="st_set_sim_mp",
        arguments=list(r = as.integer(r), as.character(data.set))), PACKAGE="SelectionTools")
    }

"st.set.sim.ef" <- function (name="effect", data.set="default")
    {
        r <- 0
        c <-  .Call("st_legacy_marker", list(action="legacy.result",
        routine="st_set_sim_ef",
        arguments=list(r = as.integer(r), as.character(name), as.character(data.set))), PACKAGE="SelectionTools")
    }

"st.set.sim.pp" <- function (pop.name="default", data.set="default")
    {
        if ((1L != length(pop.name)) || is.na(pop.name) ||
            nchar(pop.name, type="bytes") < 2L ||
            nchar(pop.name, type="bytes") >= 256L ||
            !grepl("^[A-Za-z0-9]+$", pop.name))
           stop("pop.name must contain 2 to 255 ASCII letters or digits")
        r <- 0
        c <-  .Call("st_legacy_marker", list(action="legacy.result",
        routine="st_set_sim_pp",
        arguments=list(r = as.integer(r), as.character(pop.name), as.character(data.set))), PACKAGE="SelectionTools")
    }

"st.set.simpop" <- function ( pop.name ,
                              data.set ="default",
                              effect.name="effect" )
{
  if ((1L != length(pop.name)) || is.na(pop.name) ||
      nchar(pop.name, type="bytes") < 2L ||
      nchar(pop.name, type="bytes") >= 256L ||
      !grepl("^[A-Za-z0-9]+$", pop.name))
     stop("pop.name must contain 2 to 255 ASCII letters or digits")

  reset.all()

  c <- st.set.sim.gp(data.set=data.set)

  if (is.null(c$r) ||
      0L != c$r)
    return(invisible(c))

  c <- st.set.sim.mp(data.set=data.set)

  if (is.null(c$r) ||
      0L != c$r)
    return(invisible(c))

  c <- st.set.sim.ef(name=effect.name, data.set=data.set)

  if (is.null(c$r) ||
      0L != c$r)
    return(invisible(c))

  c <- st.set.sim.pp(pop.name, data.set=data.set)

  invisible(c)
}





st.plot.gene.diversity <- function (
   data.set ="default",
   plt.pdf      = FALSE,
   plt.fname    = NULL,
   plt.width    = 10,
   plt.height   = 6,
   plt.ptsize   = 14
)
{
  st.output.dir <- get0("st.output.dir", envir = .GlobalEnv,
                        inherits = FALSE, ifnotfound = "")

    if (TRUE == plt.pdf && is.null(plt.fname))
        stop("plt.fname must be supplied when plt.pdf = TRUE")

    ms <- st.marker.data.statistics(ind=FALSE,mar=TRUE, gen=TRUE,
                                    data.set=data.set)$marker.list
    mp <- st.get.map(data.set=data.set)
    mp <- merge (mp,ms,sort=FALSE)
    mp <- data.frame(chrom=mp$Chrom,pos=mp$Pos,name=mp$Name,gene.div=mp$ExHet)
    c.pos <- st.chrom.stats(data.set=data.set)$Length
    if (length(c.pos) > 1L)
        for (i in seq.int(2L, length(c.pos))) c.pos[i] <- c.pos[i] + c.pos[i-1L]
    mp$cpos <- mp$pos
    for (i in seq_len(nrow(mp)))
        if (mp$chrom[i] > 1)
            mp$cpos[i] <- mp$cpos[i] + c.pos[ mp$chrom[i]-1]
    cl <- st.chrom.stats(data.set=data.set)$Length
      
    if ("" == st.output.dir) op <- "" else op <- "/"
    if (TRUE == plt.pdf)
        out.filename <- paste(st.output.dir, op, plt.fname, ".pdf", sep="")
    if (TRUE == plt.pdf)
        {
            pdf(
                file  = out.filename,
                width = plt.width,
                height=plt.height,
                pointsize=plt.ptsize)
            plot.device <- grDevices::dev.cur()
            oldpar <- graphics::par(no.readonly = TRUE)
            on.exit({
                dl <- grDevices::dev.list()
                if (!is.null(dl) && plot.device %in% dl) {
                    if (grDevices::dev.cur() != plot.device)
                        grDevices::dev.set(plot.device)
                    try(.st.restore.graphics(oldpar), silent = TRUE)
                    grDevices::dev.off(plot.device)
                }
            }, add = TRUE)
        }
    else
        {
            oldpar <- graphics::par(no.readonly = TRUE)
            on.exit(.st.restore.graphics(oldpar), add = TRUE)
       
        }

    par(mar=c(4,4,1,1))
    plot(
        mp$cpos,
        mp$gene.div,
        type="l",
        lwd=2,
        col="blue",
        ylim=c(0,1),
        xaxt="n",
        xlab="Chromosome",
        ylab="Gene diversity"
        )

    c.pos <- c(0,c.pos)
    axis( 1,
         at=c.pos,
         labels=c(1:length(cl)," "),
         tck=-.01)



    for (i in 1:(length(c.pos))){
        x <- c.pos[i]
        y0 <- 0
        y1  <- 1
        lines(c(x,x),c(y0,y1))
    }



    gdiv.chr <- tapply(mp$gene.div,mp$chrom,mean)
    for (i in 1: (length(c.pos)-1) ){
        text((c.pos[i]+c.pos[i+1])/2,0.95, sprintf("%4.2f",gdiv.chr[i]))
    }

    ret.val <- data.frame( chrom = mp$chrom,
                          pos   = mp$pos,
                          name  = mp$name,
                          cpos  = mp$cpos,
                          gene.div=mp$gene.div);

    return(invisible(ret.val))

}





"gs.set.performance.data" <- function( y,
                                       data.set = "default" )
{
  y <- as.double(y)
  if (any(is.infinite(y)))
    stop("y must not contain infinite values")
  retval <- 0;
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="gs_set_performance_data_01_GV",
        arguments=list(as.integer(length(y)), as.double(y), retval = as.integer(retval), as.character(data.set)),
        result=3L,
        naok=TRUE), PACKAGE="SelectionTools")

  return(invisible(c))

}

st.get.simpop.perfdata <- function(pop.name,data.set,effect)
{
    st.get.simpop (pop.name, data.set)
    genotype.population (pop.name)
    evaluate.population (pop.name,effect)
    invisible(.Call("st_legacy_sim", list(action="performance.copy",
        population=as.character(pop.name), data.set=as.character(data.set)),
        PACKAGE="SelectionTools"))
}










st.LDheatmap <- function(LD, map = NULL,
                         distances = "genetical",
                         LDmeasure = "r",
                         title = "",
                         col = grDevices::heat.colors(20),
                         zlim = c(0, 1),
                         add.map = TRUE,
                         add.key = TRUE,
                         depth = Inf,
                         max.dist = Inf,
                         write.ld = NULL,
                         write.snp.id = FALSE,
                         text = FALSE,
                         digits = 2,
                         geneMapLocation = 0.15,
                         geneMapLabelX = NULL,
                         geneMapLabelY = NULL,
                         cex.title = 0.9,
                         cex.map = 0.7,
                         cex.key = 0.6,
                         cex.ld = NULL,
                         cex.snp = NULL)
{
    is.flag <- function(x)
        is.logical(x) && length(x) == 1L && !is.na(x)

    check.cex <- function(x, name)
    {
        if (!is.numeric(x) || length(x) != 1L ||
            !is.finite(x) || x <= 0)
            stop(sprintf("'%s' must be a positive finite number", name),
                 call. = FALSE)
    }

    if (!is.flag(add.map))
        stop("'add.map' must be TRUE or FALSE", call. = FALSE)
    if (!is.flag(add.key))
        stop("'add.key' must be TRUE or FALSE", call. = FALSE)
    if (!is.flag(write.snp.id))
        stop("'write.snp.id' must be TRUE or FALSE", call. = FALSE)
    if (!is.flag(text))
        stop("'text' must be TRUE or FALSE", call. = FALSE)

    if (!is.character(title) || length(title) != 1L || is.na(title))
        stop("'title' must be a single character string", call. = FALSE)

    if (!is.character(LDmeasure) || length(LDmeasure) != 1L ||
        is.na(LDmeasure))
        stop("'LDmeasure' must be a single character string", call. = FALSE)

    distances <- match.arg(distances,
                           c("genetical", "genetic", "physical"))
    if (distances == "genetic")
        distances <- "genetical"

    if (!is.numeric(geneMapLocation) || length(geneMapLocation) != 1L ||
        !is.finite(geneMapLocation) || geneMapLocation < 0)
        stop("'geneMapLocation' must be a non-negative finite number",
             call. = FALSE)

    if (!is.null(geneMapLabelX) &&
        (!is.numeric(geneMapLabelX) || length(geneMapLabelX) != 1L ||
         !is.finite(geneMapLabelX)))
        stop("'geneMapLabelX' must be NULL or a finite number",
             call. = FALSE)

    if (!is.null(geneMapLabelY) &&
        (!is.numeric(geneMapLabelY) || length(geneMapLabelY) != 1L ||
         !is.finite(geneMapLabelY)))
        stop("'geneMapLabelY' must be NULL or a finite number",
             call. = FALSE)

    check.cex(cex.title, "cex.title")
    check.cex(cex.map, "cex.map")
    check.cex(cex.key, "cex.key")
    if (!is.null(cex.ld))
        check.cex(cex.ld, "cex.ld")
    if (!is.null(cex.snp))
        check.cex(cex.snp, "cex.snp")

    if (!is.numeric(depth) || length(depth) != 1L || is.na(depth) ||
        depth < 0 || (!is.infinite(depth) && depth != floor(depth)))
        stop("'depth' must be a non-negative integer or Inf",
             call. = FALSE)

    if (!is.numeric(max.dist) || length(max.dist) != 1L ||
        is.na(max.dist) || max.dist < 0)
        stop("'max.dist' must be a non-negative number or Inf",
             call. = FALSE)

    if (!is.numeric(zlim) || length(zlim) != 2L ||
        any(!is.finite(zlim)) || zlim[1L] >= zlim[2L])
        stop("'zlim' must contain two increasing finite numbers",
             call. = FALSE)

    if (length(col) < 2L || anyNA(col))
        stop("'col' must contain at least two valid colors",
             call. = FALSE)
    tryCatch(grDevices::col2rgb(col),
             error = function(e)
                 stop("'col' contains an invalid color", call. = FALSE))

    if (!is.null(write.ld) && !is.function(write.ld))
        stop("'write.ld' must be NULL or a function", call. = FALSE)

    if (!is.numeric(digits) || length(digits) != 1L || is.na(digits) ||
        digits < 0 || digits != floor(digits))
        stop("'digits' must be a non-negative integer", call. = FALSE)

    LD <- as.matrix(LD)

    if (!is.numeric(LD))
        stop("'LD' must be a numeric matrix", call. = FALSE)

    if (length(dim(LD)) != 2L || nrow(LD) != ncol(LD))
        stop("'LD' must be a square matrix", call. = FALSE)

    n <- nrow(LD)
    if (n < 2L)
        stop("at least two markers are required", call. = FALSE)

    if (any(!is.finite(LD) & !is.na(LD)))
        stop("'LD' contains non-finite values", call. = FALSE)

    if (add.map || is.finite(max.dist)) {
        if (is.null(map))
            stop("'map' is required when 'add.map = TRUE' or 'max.dist' is finite",
                 call. = FALSE)
    }

    if (!is.null(map)) {
        if (!is.numeric(map) || length(map) != n || any(!is.finite(map)))
            stop("'map' must be a finite numeric vector with one position per marker",
                 call. = FALSE)
        if (any(diff(map) < 0))
            stop("'map' positions must be in non-decreasing order",
                 call. = FALSE)
    }

    plot.LD <- LD

    for (j in seq.int(2L, n)) {
        i <- seq_len(j - 1L)
        missing.upper <- is.na(plot.LD[i, j]) & !is.na(plot.LD[j, i])
        if (any(missing.upper))
            plot.LD[i[missing.upper], j] <- plot.LD[j, i[missing.upper]]
    }

    tol <- sqrt(.Machine$double.eps)

    for (j in seq.int(2L, n)) {
        i <- seq_len(j - 1L)
        values <- plot.LD[i, j]
        present <- !is.na(values)

        if (any(present)) {
            x <- values[present]
            if (any(x < zlim[1L] - tol) || any(x > zlim[2L] + tol))
                stop("LD values outside the range specified by 'zlim'",
                     call. = FALSE)
            values[present] <- pmax(zlim[1L], pmin(zlim[2L], x))
            plot.LD[i, j] <- values
        }
    }

    for (j in seq_len(n))
        plot.LD[j:n, j] <- NA_real_

    if (is.finite(depth) || is.finite(max.dist)) {
        for (j in seq.int(2L, n)) {
            i <- seq_len(j - 1L)
            omit <- rep(FALSE, length(i))

            if (is.finite(depth))
                omit <- omit | ((j - i) > depth)

            if (is.finite(max.dist))
                omit <- omit | ((map[j] - map[i]) > max.dist)

            if (any(omit))
                plot.LD[i[omit], j] <- NA_real_
        }
    }

    if (is.null(write.ld) && text) {
        format.ld <- paste0("%.", digits, "f")
        write.ld <- function(x) sprintf(format.ld, x)
    }

    snp.id <- rownames(LD)
    if (is.null(snp.id))
        snp.id <- colnames(LD)
    if (is.null(snp.id))
        snp.id <- as.character(seq_len(n))
    snp.id <- as.character(snp.id)
    bad.id <- is.na(snp.id) | !nzchar(snp.id)
    if (any(bad.id))
        snp.id[bad.id] <- as.character(which(bad.id))

    snp <- ((seq_len(n) - 1) + 0.5) / n
    plot.col <- rev(col)
    breaks <- seq(zlim[1L], zlim[2L], length.out = length(plot.col) + 1L)

    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(.st.restore.graphics(oldpar), add = TRUE)
    graphics::par(mar = c(0.5, 1.0, 1.5, 1.0))

    graphics::plot(NA_real_, NA_real_,
                   type = "n",
                   axes = FALSE,
                   ann = FALSE,
                   xlim = c(0, 1),
                   ylim = c(0, 1),
                   xaxs = "i",
                   yaxs = "i",
                   asp = 1)

    graphics::image(snp, snp, plot.LD,
                    zlim = zlim,
                    col = plot.col,
                    breaks = breaks,
                    axes = FALSE,
                    xlab = "",
                    ylab = "",
                    add = TRUE)

    if (!is.null(write.ld)) {
        ind <- which(!is.na(plot.LD), arr.ind = TRUE)

        if (nrow(ind) > 0L) {
            label <- vapply(plot.LD[ind],
                            function(x) {
                                ans <- write.ld(x)
                                if (length(ans) != 1L)
                                    stop("'write.ld' must return one value for each LD value",
                                         call. = FALSE)
                                if (is.na(ans)) "" else as.character(ans)
                            },
                            character(1L))

            show <- nzchar(label)

            if (any(show)) {
                if (is.null(cex.ld)) {
                    lab <- label[show]
                    cell.size <- 1 / n
                    w <- max(graphics::strwidth(lab, cex = 1))
                    h <- max(graphics::strheight(lab, cex = 1))
                    w <- max(w, .Machine$double.eps)
                    h <- max(h, .Machine$double.eps)
                    cex.ld.use <- min(0.8,
                                      0.75 * cell.size / w,
                                      0.75 * cell.size / h)
                } else {
                    cex.ld.use <- cex.ld
                }

                graphics::text(snp[ind[show, 1L]],
                               snp[ind[show, 2L]],
                               labels = label[show],
                               cex = cex.ld.use)
            }
        }
    }

    regionx <- NULL
    regiony <- NULL

    if (add.map) {
        min.dist <- min(map)
        max.dist.map <- max(map)
        total.dist <- max.dist.map - min.dist

        map.span <- snp[n] - snp[1L]
        if (geneMapLocation >= map.span)
            stop("'geneMapLocation' is too large for the number of markers",
                 call. = FALSE)

        seq.x <- c(snp[1L] + geneMapLocation, snp[n])
        seq.y <- c(snp[1L], snp[n] - geneMapLocation)

        if (total.dist > 0)
            p <- (map - min.dist) / total.dist
        else
            p <- rep(0.5, n)

        regionx <- seq.x[1L] + p * (seq.x[2L] - seq.x[1L])
        regiony <- seq.y[1L] + p * (seq.y[2L] - seq.y[1L])

        graphics::segments(snp, snp, regionx, regiony, lwd = 0.5)
        graphics::segments(seq.x[1L], seq.y[1L], seq.x[2L], seq.y[2L])

        if (distances == "physical") {
            map.label <- paste0("Physical Length:",
                                format(round(total.dist / 1000, 1), trim = TRUE),
                                "kb")
        } else {
            map.label <- paste0("Genetic Map Length:",
                                format(round(total.dist, 1), trim = TRUE),
                                "cM")
        }

        if (is.null(geneMapLabelX))
            geneMapLabelX <- 0.52
        if (is.null(geneMapLabelY))
            geneMapLabelY <- 0.24

        graphics::text(geneMapLabelX, geneMapLabelY,
                       labels = map.label,
                       adj = c(0, 0.5),
                       cex = cex.map)
    }

    if (write.snp.id) {
        if (is.null(cex.snp)) {
            label.length <- max(nchar(snp.id), 1L)
            cex.snp.use <- min(0.7, 12 / (n * sqrt(label.length)))
        } else {
            cex.snp.use <- cex.snp
        }

        if (add.map) {
            graphics::text(regionx, regiony,
                           labels = paste0(" ", snp.id),
                           srt = -45,
                           adj = c(0, 0.5),
                           cex = cex.snp.use,
                           xpd = TRUE)
        } else {
            offset <- 0.01
            graphics::text(snp + offset, snp - offset,
                           labels = paste0(" ", snp.id),
                           srt = -45,
                           adj = c(0, 0.5),
                           cex = cex.snp.use,
                           xpd = TRUE)
        }
    }

    if (add.key) {
        key.x1 <- 0.58
        key.x2 <- 0.95
        key.y1 <- 0.050
        key.y2 <- 0.085

        xx <- seq(key.x1, key.x2, length.out = length(plot.col) + 1L)

        graphics::rect(xx[-length(xx)], key.y1,
                       xx[-1L], key.y2,
                       col = plot.col,
                       border = NA)
        graphics::rect(key.x1, key.y1, key.x2, key.y2,
                       border = "black")

        ticks <- seq(zlim[1L], zlim[2L], length.out = 6L)
        tick.x <- key.x1 +
            (ticks - zlim[1L]) / diff(zlim) * (key.x2 - key.x1)

        graphics::segments(tick.x, key.y1,
                           tick.x, key.y1 - 0.010)
        graphics::text(tick.x, key.y1 - 0.012,
                       labels = format(signif(ticks, 3), trim = TRUE),
                       adj = c(0.5, 1),
                       cex = cex.key)

        measure <- tolower(LDmeasure)
        if (measure %in% c("r", "r2", "r^2"))
            key.label <- expression(R^2 ~ "Color Key")
        else if (measure %in% c("dp", "dprime", "d'", "d"))
            key.label <- "D' Color Key"
        else
            key.label <- paste(LDmeasure, "Color Key")

        graphics::text((key.x1 + key.x2) / 2,
                       key.y2 + 0.018,
                       labels = key.label,
                       cex = cex.key)
    }

    if (nzchar(title))
        graphics::mtext(title,
                        side = 3,
                        line = 0,
                        adj = 0.5,
                        cex = cex.title)

    invisible(list(LDmatrix = plot.LD,
                   map = map,
                   snp.position = snp,
                   map.x = regionx,
                   map.y = regiony,
                   distances = distances,
                   depth = depth,
                   max.dist = max.dist,
                   zlim = zlim,
                   col = col))
}











st.calculate.enrichment <- function(base.pool, enrichment.pool,
                                    haplotype.options=NULL, blocks=NULL,
                                    exclude.compatible=TRUE) {
    .Call("st_enrichment", base.pool, enrichment.pool, haplotype.options,
          blocks, exclude.compatible, NULL, "best.subset", PACKAGE="SelectionTools")
}

st.select.enrichment <- function(base.pool, enrichment.pool, n,
                                 haplotype.options=NULL, blocks=NULL,
                                 exclude.compatible=TRUE,
                                 method="best.subset") {
    if (is.null(n)) stop("n must be a positive integer", call.=FALSE)
    .Call("st_enrichment", base.pool, enrichment.pool, haplotype.options,
          blocks, exclude.compatible, n, method, PACKAGE="SelectionTools")
}


st.detect.core <- function(pool, n.max,
                           haplotype.options=NULL, blocks=NULL,
                           exclude.compatible=TRUE, method="forward") {
    .Call("st_detect_core", pool, n.max, haplotype.options, blocks,
          exclude.compatible, method, PACKAGE="SelectionTools")
}
.st.legacy.path <- function (filename, output = FALSE) 
{
    prefix <- get0(if (output) 
        "st.output.dir"
    else "st.input.dir", envir = .GlobalEnv, inherits = FALSE, 
        ifnotfound = "")
    mab.df(prefix, filename)
}
