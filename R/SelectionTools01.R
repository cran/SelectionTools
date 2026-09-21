###################################################################################
# SelectionTools01.R
###################################################################################

.GS_ALLELE_MAX <- 254L
.GS_STVCF_MAX_ALT <- .GS_ALLELE_MAX - 1L


## Startup message

st.info <-  function( lev, msg )
  {
    c <- .C( "g_r_info",
             as.integer(lev),
             as.character(msg) );
  }

"st.set.num.threads" <- function( active )
{
  c  <-  .C( "gs_set_num_threads", as.integer(active) )
}

"st.set.openblas.threads" <- function( active )
{
  c  <-  .C( "gs_set_openblas_threads", as.integer(active) )
}


"st.get.num.threads" <- function()
{
  active <- 0;
  c  <-  .C( "gs_get_num_threads", t=as.integer(active) )
  return( c$t )
}

###################################################################################
# Genomic selection functions: gs.
###################################################################################

"gs.info" <-  function( lev, msg ){
    c <- .C("gs_r_info",
            as.integer(lev),
            as.character(msg)
            );
  }

"gs.set.num.threads" <- function( active )
{
  c  <-  .C( "gs_set_num_threads", as.integer(active) )
}

"gs.reset" <-  function()
{
    c <- .C("gs_reset");
  }

"gs.info" <-  function( lev, msg ){
    c <- .C("gs_r_info",
            as.integer(lev),
            as.character(msg)
            );
  }

"gs.set.info.level" <- function( level,data.set ="default" )
{
  c  <-  .C( "gs_set_info_level_GV",
            as.integer(level),
            as.character(data.set)   )
}

"gs.set.all.info.levels" <- function( level )
{
  c  <-  .C( "gs_set_all_info_levels",
            as.integer(level)   )
}


"gs.get.info.level" <- function ( data.set ="default" )
{
  level = 0
  c  <-  .C(  "gs_get_info_level_GV",
              level = as.integer(level),
              as.character(data.set)      )
  return (c$level)
}

"gs.test.mp" <- function(n = 0) {
  
  c  <-  .C(
    "gs_test_mp",
    as.integer(n)
  )
}


"st.set.matrix.ops" <- function(select.s="mt",select.p="st") {
  c <- .C( "st_set_matrix_ops",
           as.character(select.s),
           as.character(select.p)   )
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


"gs.build.Z" <- function (
    out.filename = "Z.matrix" ,
    auxfiles     = FALSE,
    data.set     = "default"
    )
{
  gs.start.timer()
  X.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(X.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_build_Z_01_GV"          ,
          as.character(X.filename)   ,
          as.integer(auxfiles)       ,
          retval = as.integer(retval),
          as.character(data.set)      
          )
  Z <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    X <- read.table(file=X.filename)
    Z <- as.matrix(X)
    unlink(X.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(Z) )
}


"gs.mme.coeff" <- function(
                              out.filename = "mme.coeff" ,
                              auxfiles     = FALSE,
                              data.set     = "default"
                              )
{
  gs.start.timer()
  A.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(A.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_mme_coeff_01_GV",
          as.character(A.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)      
          )
  A <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    A <- data.matrix( read.table(file=A.filename,header=FALSE) )
    unlink(A.filename)
    colnames(A) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(A) )
}

"gs.build.V" <- function(
                              out.filename = "V.matrix" ,
                              auxfiles     = FALSE,
                              data.set     = "default"
                              )
{
  gs.start.timer()
  A.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(A.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_build_V_01_GV",
          as.character(A.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)      
          )
  A <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    A <- data.matrix( read.table(file=A.filename,header=FALSE) )
    unlink(A.filename)
    colnames(A) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(A) )
}

"gs.get.Z" <- function ( data.set = "default")
{
  gs.start.timer()
  Z <- .Call( "st_get_Z",
              data.set)
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

"gs.get.y" <- function ( data.set = "default")
{
  gs.start.timer()
  y <- .Call( "st_get_y",
              data.set)
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


"gs.mme.restcoeff" <- function( out.filename = "mme.coeff" ,
                                  auxfiles     = FALSE       )
{
  gs.start.timer()
  A.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(A.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_mme_restcoeff_01_GV",
          as.character(A.filename),
          as.integer(auxfiles),
          retval = as.integer(retval)
          )
  A <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    A <- data.matrix( read.table(file=A.filename,header=FALSE) )
    unlink(A.filename)
    colnames(A) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(A) )
}

"gs.mme.rhs" <- function ( out.filename = "mme.rhs" ,
                              auxfiles     = FALSE      ,
                              data.set     = "default" )
{
  gs.start.timer()
  b.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(b.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_mme_rhs_01_GV",
          as.character(b.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)        )
  b <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    b <- data.matrix( read.table(file=b.filename,header=FALSE) )
    unlink(b.filename)
    colnames(b) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(b) )
}


"gs.mme.restrhs" <- function ( out.filename = "mme.rhs" ,
                                 auxfiles = FALSE              )
{
  gs.start.timer()
  b.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(b.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_mme_restrhs_01_GV",
          as.character(b.filename),
          as.integer(auxfiles),
          retval = as.integer(retval)  )
  b <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    b <- data.matrix( read.table(file=b.filename,header=FALSE) )
    unlink(b.filename)
    colnames(b) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(b) )
}

"gs.mme.solve" <- function( out.filename = "mme.solution" ,
                               auxfiles     = FALSE           ,
                               data.set     = "default"      )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_mme_solve_01_GV",
          as.character(x.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)        )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.mme.invcoeff" <- function( out.filename = "inv.mme" ,
                                  auxfiles = TRUE              )
{
  gs.start.timer()
  o.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(o.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_mme_invcoeff_01_GV",
          as.character(o.filename),
          as.integer(auxfiles),
          retval = as.integer(retval)  )
  A <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    A <- data.matrix( read.table(file=o.filename,header=FALSE) )
    unlink(o.filename)
    colnames(A) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(A) )
}

## "gs.check.model.fit" <- function ( cor.filename = "breeding.values"  ,
##                                       auxfiles = TRUE                   )
## {
##   gs.start.timer()
##   o.filename <- paste(st.output.dir,op,cor.filename,sep="")
##   retval <- 0;
##   c <- .C ( "gs_check_model_fit_01_GV"  ,
##             as.character(o.filename)    ,
##             as.integer(auxfiles)        ,
##             retval = as.integer(retval) )
##   x <- NULL
##   if  ( (auxfiles) && (-2 != c$retval ) ) {
##     x <- data.matrix( read.table(file=o.filename,header=TRUE) )
##     unlink(o.filename)
##   }
##   gs.stop.timer (info.level=1)
##   return( invisible(x) )
##}

## "gs.estimate.breeding.values" <- function (
##           mar.filename ,
##           cor.filename = "breeding.values",
##           auxfiles = TRUE                          )
## {
##   gs.start.timer()
##   if ( "" == st.input.dir )  ip <- "" else ip <- "/"
##   if ( "" == st.output.dir ) op <- "" else op <- "/"
##   m.filename <- paste(st.input.dir,ip,mar.filename,sep="")
##   o.filename <- paste(st.output.dir,op,cor.filename,sep="")
##   retval <- 0;
##   c <- .C ( "gs_estimate_breeding_values_01_GV",
##             as.character(m.filename),
##             as.character(o.filename),
##             as.integer(auxfiles),
##             retval = as.integer(retval))
##   x <- NULL
##   if  ( (auxfiles) && (-2 != c$retval ) ) {
##     x <- read.table(file=o.filename,header=FALSE,skip=1,stringsAsFactors=FALSE) 
##     unlink(o.filename)
##     if (2 == ncol(x)) names(x) <- c("i","yhat")
##     if (3 == ncol(x)) names(x) <- c("i","y","yhat")
##   }
##   gs.stop.timer (info.level=1)
##   return( invisible(x) )
## }


## "gs.validation" <- function ( mar.filename ,
##                                 per.filename ,
##                                 cor.filename = "predicted",
##                                 auxfiles = TRUE )
## {
##   gs.start.timer()
##   if ( "" == st.input.dir )  ip <- "" else ip <- "/"
##   if ( "" == st.output.dir ) op <- "" else op <- "/"
##   m.filename <- paste(st.input.dir,ip,mar.filename,sep="")
##   p.filename <- paste(st.input.dir,ip,per.filename,sep="")
##   o.filename <- paste(st.output.dir,op,cor.filename,sep="")
##   retval <- 0;
##   c <- .C ( "gs_validation_01_GV",
##             as.character(m.filename),
##             as.character(p.filename),
##             as.character(o.filename),
##             as.integer(auxfiles),
##             retval = as.integer(retval))
##   x <- NULL
##   if  ( (auxfiles) && (-2 != c$retval ) ) {
##     x <- data.matrix( read.table(file=o.filename,header=TRUE) )
##     unlink(o.filename)
##   }
##   gs.stop.timer (info.level=1)
##   return( invisible(x) )
## }

"gs.single.marker.reg" <- function(  out.filename = "sm.reg" ,
                                        auxfiles     = TRUE      ,
                                        data.set     = "default"  )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_single_marker_reg_01_GV",
           as.character(x.filename)    ,
           as.integer(auxfiles)        ,
           retval = as.integer(retval) ,
           as.character(data.set)        )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.single.marker.aov" <- function(  alpha = 1 ,
                                     out.filename = "sm.reg" ,
                                     auxfiles     = TRUE      ,
                                     data.set     = "default"  )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_single_marker_aov_01_GV",
            as.double(alpha)           ,
            as.character(x.filename)    ,
            as.integer(auxfiles)        ,
            retval = as.integer(retval) ,
            as.character(data.set)        )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}


"gs.lambda.const" <- function( lambda                  ,
                              out.filename = "lambda"  ,
                              auxfiles = FALSE         ,
                              data.set     = "default" )
{
  st.input.dir <- get0("st.input.dir", envir = .GlobalEnv,
                       inherits = FALSE, ifnotfound = "")

  gs.start.timer()
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
 c <- .C( "gs_lambda_const_01_GV",
           as.double(lambda),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval) ,
           as.character(data.set)         )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.set.lambda" <- function( lambda                  ,
                            out.filename = "lambda"  ,
                            auxfiles = FALSE         ,
                            data.set     = "default" )
{
  st.input.dir <- get0("st.input.dir", envir = .GlobalEnv,
                       inherits = FALSE, ifnotfound = "")

  gs.start.timer()
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
 c <- .C( "gs_lambda_var_01_GV",
         as.double(lambda),
         as.integer(length(lambda)),
         as.character(x.filename),
         as.integer(auxfiles),
         retval = as.integer(retval) ,
         as.character(data.set)         )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.lambda.hsq" <- function( hsq,
                             out.filename = "lambda"  ,
                             auxfiles     = TRUE      ,
                             data.set     = "default" )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_lambda_const_02_GV",
           as.double(hsq),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval) ,
           as.character(data.set)         )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}



"gs.lambda.reg" <- function( hsq                      ,
                                out.filename = "lambda"  ,
                                auxfiles = TRUE          ,
                                data.set     = "default"  )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_lambda_reg_01_GV"        ,
           as.double(hsq)               ,
           as.character(x.filename)     ,
           as.integer(auxfiles)         ,
           retval = as.integer(retval)  ,
           as.character(data.set)         )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.lambda.aov" <- function( hsq                      ,
                             alpha = 1                ,
                             out.filename = "lambda"  ,
                             auxfiles = TRUE          ,
                             data.set     = "default"  )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_lambda_aov_01_GV"        ,
           as.double(hsq)               ,
           as.double(alpha)             ,
           as.character(x.filename)     ,
           as.integer(auxfiles)         ,
           retval = as.integer(retval)  ,
           as.character(data.set)         )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.lambda.emstep" <- function (constvar     = FALSE        ,
                                   out.filename = "lambda" ,
                                   auxfiles     = TRUE,
                                   data.set     = "default"     )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_lambda_emstep_01_GV"    ,
           as.integer(constvar)        ,
           as.character(x.filename)    ,
           as.integer(auxfiles)        ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}


"gs.lambda.rmlv" <- function( maxiter     = 1000       ,
                                 precision   = 0.001      ,
                                 out.filename = "lambda"  ,
                                 auxfiles = TRUE          ,
                                 data.set     = "default" )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_lambda_rmlv_01_GV"     ,
           as.integer(maxiter)        ,
           as.double(precision)       ,
           as.character(x.filename)   ,
           as.integer(auxfiles)       ,
           retval = as.integer(retval),
           as.character(data.set)     )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.lambda.rmlc" <- function( maxiter     = 1000       ,
                              precision   = 0.0001      ,
                              hsq         = 0.9        ,
                              out.filename = "lambda"  ,
                              auxfiles = TRUE          ,
                              data.set     = "default" )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_lambda_rmlc_01_GV"     ,
           as.integer(maxiter)        ,
           as.double(precision)       ,
           as.double(hsq)             ,
           as.character(x.filename)   ,
           as.integer(auxfiles)       ,
           retval = as.integer(retval),
           as.character(data.set)     )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
} 

"gs.lambda.rmlr" <- function( maxiter     = 1000       ,
                              precision   = 0.001      ,
                              hsq         = 0.9        ,
                              out.filename = "lambda"  ,
                              auxfiles = TRUE          ,
                              data.set     = "default" )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_lambda_rmlr_01_GV"     ,
           as.integer(maxiter)        ,
           as.double(precision)       ,
           as.double(hsq)             ,
           as.character(x.filename)   ,
           as.integer(auxfiles)       ,
           retval = as.integer(retval),
           as.character(data.set)     )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
} 

"gs.lambda.rmla" <- function( alpha        = 1        ,
                              maxiter      = 1000     ,
                              precision    = 0.0001   ,
                              hsq          = 0.9      ,
                              out.filename = "lambda" ,
                              auxfiles     = TRUE     ,
                              data.set     = "default" )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_lambda_rmla_01_GV"     ,
           as.double(alpha)           ,
           as.integer(maxiter)        ,
           as.double(precision)       ,
           as.double(hsq)             ,
           as.character(x.filename)   ,
           as.integer(auxfiles)       ,
           retval = as.integer(retval),
           as.character(data.set)     )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
} 


"gs.set.allele.codes" <- function(    aa = 0,
                                      aA = 1,
                                      AA = 2,
                                      an = 0.5,
                                      nn = 1  ,
                                      nA = 1.5,
                                      data.set ="default" )
{
  c <- .C( "gs_set_allele_codes_01_GV",
           as.double(aa),
           as.double(aA),
           as.double(AA),
           as.double(an),
           as.double(nn),
           as.double(nA),
           as.character(data.set)   )
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
    c <- .C("gs_build_esvs_01_GV"           ,
            as.integer(es.size)          ,
            as.character(base.set)       ,            
            as.character(estimation.set) ,            
            as.character(validation.set) ,            
            as.character(out.es_m)       ,
            as.character(out.es_p)       ,
            as.character(out.vs_m)       ,
            as.character(out.vs_p)       ,
            as.integer  (auxfiles)       ,
            retval = as.integer(retval)   )
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
    c <- .C("gs_build_esvs_01_GV"           ,
            as.integer(ts.size)          ,
            as.character(base.set)       ,            
            as.character(training.set) ,            
            as.character(prediction.set) ,            
            as.character(out.es_m)       ,
            as.character(out.es_p)       ,
            as.character(out.vs_m)       ,
            as.character(out.vs_p)       ,
            as.integer  (auxfiles)       ,
            retval = as.integer(retval)   )
    gs.stop.timer (info.level=1)
  }


"gs.estimate.gv" <- function ( estimation.set ="default",
                               prediction.set ="default",
                               out.filename = "breeding.values",
                               auxfiles = TRUE                    )
{
  gs.start.timer()
  o.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(o.filename), add = TRUE)
  retval <- 0;
  c <- .C ( "gs_estimate_bv_01_GV",
            as.character(estimation.set)  ,
            as.character(prediction.set)  ,
            as.character(o.filename),
            as.integer(auxfiles)    ,
            retval = as.integer(retval)     )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- read.table(file=o.filename,header=FALSE,skip=1,stringsAsFactors=FALSE)
    unlink(o.filename)
    if (3 == ncol(x))  names(x) <- c("i","y","yhat")
    if (2 == ncol(x))  names(x) <- c("i","yhat")
  }
    gs.stop.timer (info.level=1)
  return( invisible(x) )
}
                                        


"gs.predict.genotypes" <- function ( training.set ="default",
                                     prediction.set ="default",
                                     out.filename = "breeding.values",
                                     auxfiles = TRUE                    )
{
  gs.start.timer()
  o.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(o.filename), add = TRUE)
  retval <- 0;
  c <- .C ( "gs_estimate_bv_01_GV",
            as.character(training.set)  ,
            as.character(prediction.set)  ,
            as.character(o.filename),
            as.integer(auxfiles)    ,
            retval = as.integer(retval)     )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- read.table(file=o.filename,header=FALSE,skip=1,stringsAsFactors=FALSE)
    unlink(o.filename)
    if (3 == ncol(x))  names(x) <- c("i","y","yhat")
    if (2 == ncol(x))  names(x) <- c("i","yhat")
  }
    gs.stop.timer (info.level=1)
  return( invisible(x) )
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
    c <- .C("gs_restrict_marker_data_01_GV"    ,
            as.character(out.filename)   ,
            as.character(trait.filename) ,
            as.integer(NoAll.MAX),
            as.double (MaMis.MAX),
            as.double (ExHet.MIN),
            as.double (InMis.MAX),
            retval = as.integer(retval),
            as.character(data.set)       )
    gs.stop.timer (info.level=1)
  }


###################################################################################
# Selection tools functions: st.
###################################################################################


"st.reset" <-  function(){
    c <- .C("gs_reset");
  }

"st.info" <-  function( lev, msg ){
    c <- .C(  "st_r_info",
              as.integer(lev),
              as.character(msg)    );
  }

"st.set.info.level" <- function( level )
{
  c  <-  .C( "st_set_info_level",
            as.integer(level)  )
}

"st.get.info.level" <- function ()
{
  level = 0
  c  <-  .C(  "st_get_info_level",
              level = as.integer(level)   )
  return (c$level)
}

"st.start.timer" <- function(depth=1)
  {
    if ( 0 > depth ) depth = 0 
    if ( depth > 5 ) depth = 4 else depth = depth-1
    
      c  <-  .C( "st_start_timer",
                 as.integer(depth),
                 as.double( proc.time()[3] ) )
  }

"st.stop.timer" <- function ( depth=1, info.level=0 )
{
    if ( 0 > depth ) depth = 0 
    if ( depth > 5 ) depth = 4 else depth = depth-1

    etime = 0
    c  <-  .C( "st_stop_timer",
               as.integer(depth),
               etime=as.double( etime ) )

    time <- proc.time()[3] - c$etime

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

#######

"st.read.marker.data" <- function ( filename,
                                    format,
                                    auxfiles  = FALSE,
                                    data.set  = "default" )
{
  st.input.dir <- get0("st.input.dir", envir = .GlobalEnv,
                       inherits = FALSE, ifnotfound = "")

  st.start.timer()
  is.absolute <- startsWith(filename, "/") ||
                 startsWith(filename, "\\") ||
                 (nchar(filename) >= 3L &&
                  grepl("^[A-Za-z]$", substr(filename, 1L, 1L)) &&
                  substr(filename, 2L, 2L) == ":" &&
                  substr(filename, 3L, 3L) %in% c("/", "\\"))
  if (is.absolute) {
    dta.filename <- filename
  } else {
    if ( "" == st.input.dir ) ip <- "" else ip <- "/"
    dta.filename <- paste(st.input.dir, ip, filename, sep="")
  }
  tmp.base <- tempfile(
    pattern = paste0("SelectionTools-", basename(filename), "-"),
    tmpdir = tempdir()
  )
  ind.filename <- paste0(tmp.base, ".ind")
  mar.filename <- paste0(tmp.base, ".mar")
  gen.filename <- paste0(tmp.base, ".gen")
  on.exit(unlink(c(ind.filename, mar.filename, gen.filename)), add = TRUE)
  retval <- 0;

  if ("l" == format)
    {
      c <- .C("gs_read_marker_data_01_GV",
              as.character(dta.filename),
              as.character(ind.filename),
              as.character(mar.filename),
              as.character(gen.filename),
              as.integer(auxfiles),
              retval = as.integer(retval),
              as.character(data.set)
              )
    }

    if ("m" == format)
    {
      c <- .C("gs_read_marker_data_02_GV",
              as.character(dta.filename),
              as.character(ind.filename),
              as.character(mar.filename),
              as.character(gen.filename),
              as.integer(auxfiles),
              retval = as.integer(retval),
              as.character(data.set)
              )
    }

    if ("n" == format)
    {
      c <- .C("gs_read_marker_data_03_GV",
              as.character(dta.filename),
              as.character(ind.filename),
              as.character(mar.filename),
              as.character(gen.filename),
              as.integer(auxfiles),
              retval = as.integer(retval),
              as.character(data.set)
              )
    }

    if ("t" == format)
    {
      c <- .C("gs_read_marker_data_04_GV",
              as.character(dta.filename),
              as.character(ind.filename),
              as.character(mar.filename),
              as.character(gen.filename),
              as.integer(auxfiles),
              retval = as.integer(retval),
              as.character(data.set)
              )
    }
  
  outp <- NULL

  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      individual.list <- read.table(file=ind.filename,
                                    check.names = FALSE,
                                    header=TRUE,
                                    stringsAsFactors=FALSE)
      outp$individual.list <- individual.list
      marker.list      <- read.table ( file=mar.filename,
                                       check.names = FALSE,
                                       header=TRUE,
                                       stringsAsFactors=FALSE  )
      outp$marker.list <- marker.list
      genotypes        <- read.table ( file=gen.filename,
                                      check.names = FALSE,
                                      header=TRUE,
                                      stringsAsFactors=FALSE)
      outp$genotypes   <- genotypes
      unlink(ind.filename)
      unlink(mar.filename)
      unlink(gen.filename)
    }
  
  st.stop.timer (info.level=1)
  return(invisible(outp))

}

"st.read.performance.data" <- function ( filename,
                                         out.filename = "y.vector",
                                         auxfiles     = TRUE,
                                         data.set     ="default"     )
{
  st.input.dir <- get0("st.input.dir", envir = .GlobalEnv,
                       inherits = FALSE, ifnotfound = "")

  gs.start.timer()
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  dta.filename <- paste ( st.input.dir,  ip, filename,sep="" )
  per.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(per.filename), add = TRUE)
  unlink(per.filename)
  retval <- 0;
  c <- .C("gs_read_performance_data_01_GV",
          as.character(dta.filename),
          as.character(per.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )
  y <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      y <- read.table(file=per.filename,header=FALSE,skip=1)
      unlink(per.filename)
      names(y) <- c("i","y")
    }
  gs.stop.timer (info.level=1)
  return( invisible(y) )
}

"st.return.performance.data" <- function ( out.filename = "y.vector",
                                           auxfiles     = TRUE,
                                           data.set     ="default"     )
{
  st.input.dir <- get0("st.input.dir", envir = .GlobalEnv,
                       inherits = FALSE, ifnotfound = "")

  gs.start.timer()
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  per.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(per.filename), add = TRUE)
  unlink(per.filename)
  retval <- 0;
  c <- .C("gs_return_performance_data_01_GV",
          as.character(per.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )
  y <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      y <- read.table(file=per.filename, header=FALSE, skip=1,
                      colClasses=c("character","numeric"),
                      stringsAsFactors=FALSE)
      unlink(per.filename)
      names(y) <- c("i","y")
    }
  gs.stop.timer (info.level=1)
  return( invisible(y) )
}


"st.read.map" <- function ( filename,
                            skip = 0,
                            format     = "cpms", # or "mcp" 
                            m.stretch  = 1,
                            auxfiles   = TRUE,
                            data.set   = "default" )
{
  st.input.dir <- get0("st.input.dir", envir = .GlobalEnv,
                       inherits = FALSE, ifnotfound = "")

  st.start.timer()
  is.absolute <- startsWith(filename, "/") ||
                 startsWith(filename, "\\") ||
                 (nchar(filename) >= 3L &&
                  grepl("^[A-Za-z]$", substr(filename, 1L, 1L)) &&
                  substr(filename, 2L, 2L) == ":" &&
                  substr(filename, 3L, 3L) %in% c("/", "\\"))
  if (is.absolute) {
    map.filename <- filename
  } else {
    if ( "" == st.input.dir ) ip <- "" else ip <- "/"
    map.filename <- paste(st.input.dir, ip, filename, sep="")
  }
  out.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(filename), "-"),
    tmpdir = tempdir(),
    fileext = ".srt"
  )
  on.exit(unlink(out.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_read_map_01_GV",
          as.character(map.filename),
          as.integer(skip),
          as.character(format),
          as.double (m.stretch),
          as.character(out.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=out.filename,header=FALSE)
      unlink(out.filename)
      names(x) <- c("Chrom","Pos","Name","Class")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}

"st.dataframe.to.STvcf" <- function(data)
{
   if (!is.data.frame(data))
      stop("data must be a data frame")

   nms <- names(data)
   if (anyNA(nms) || any(!nzchar(nms)))
      stop("data frame column names must be non-missing and non-empty")
   if (ncol(data) < 6L)
      stop("data must contain fixed fields and at least one sample column")

   if (!(nms[1L] %in% c("CHROM", "#CHROM")) ||
       !identical(nms[2:5], c("POS", "ID", "REF", "ALT")))
      stop("first five columns must be CHROM/#CHROM, POS, ID, REF, ALT")

   nr <- nrow(data)
   if (0L == nr)
      stop("data contains no markers")

   chrom.text <- as.character(data[[1L]])
   if (anyNA(chrom.text) || any(!nzchar(chrom.text)))
      stop("CHROM contains missing or empty values")

   chrom.levels <- unique(chrom.text)
   CHROM <- match(chrom.text, chrom.levels)
   map.text <- paste0(seq_along(chrom.levels), "=", chrom.levels,
                      collapse=", ")
   st.info(0, paste("Chromosome coding:", map.text))

   format.candidates <- which(nms == "FORMAT" & seq_along(nms) > 5L)
   valid.format <- logical(length(format.candidates))
   for (k in seq_along(format.candidates))
   {
      z <- as.character(data[[format.candidates[k]]])
      valid.format[k] <- all(!is.na(z) & nzchar(z) & z != "." &
         !endsWith(z, ":") &
         vapply(strsplit(z, ":", fixed=TRUE),
                function(x) all(nzchar(x)) && sum(x == "GT") == 1L,
                logical(1L)))
   }
   format.candidates <- format.candidates[valid.format]
   if (length(format.candidates) > 1L)
      stop("more than one column can be interpreted as FORMAT")
   format.pos <- if (length(format.candidates)) {
      format.candidates[1L]
   } else {
      NA_integer_
   }

   cm.pos <- NA_integer_
   if (!is.na(format.pos))
   {
      meta.idx <- if (format.pos > 6L) 6L:(format.pos - 1L) else integer(0)
      if (length(meta.idx))
      {
         if (any(!(nms[meta.idx] %in% c("QUAL", "FILTER", "INFO", "CM"))))
            stop("unrecognized column before FORMAT")
         if (any(table(nms[meta.idx]) > 1L))
            stop("VCF metadata columns before FORMAT must be unique")
         cm.idx <- meta.idx[nms[meta.idx] == "CM"]
         if (length(cm.idx))
            cm.pos <- cm.idx[1L]
      }
      if (format.pos == ncol(data))
         stop("FORMAT is present but no sample columns follow it")
      sample.idx <- (format.pos + 1L):ncol(data)
   }
   else
   {
      rest <- 6L:ncol(data)
      cm.idx <- rest[nms[rest] == "CM"]
      if (length(cm.idx))
      {
         cm.valid <- logical(length(cm.idx))
         for (k in seq_along(cm.idx))
         {
            z <- data[[cm.idx[k]]]
            if (is.factor(z))
               z <- as.character(z)
            z <- suppressWarnings(as.numeric(z))
            cm.valid[k] <- all(!is.na(z) & is.finite(z) & z >= 0)
         }
         cm.idx <- cm.idx[cm.valid]
         if (length(cm.idx) > 1L)
            stop("more than one column can be interpreted as CM")
         if (length(cm.idx))
            cm.pos <- cm.idx[1L]
      }
      sample.idx <- setdiff(rest, cm.pos)
   }

   if (0L == length(sample.idx))
      stop("data contains no sample genotype columns")
   sample.names <- nms[sample.idx]
   if (anyDuplicated(sample.names))
      stop("sample names must be unique")
   if (any(!grepl("^[A-Za-z0-9]+$", sample.names)))
      stop("sample names may contain only A-Z, a-z, and 0-9")
   if (any(nchar(sample.names, type="bytes") >= 256L))
      stop("sample names must contain fewer than 256 bytes")

   pos.input <- if (!is.na(cm.pos)) data[[cm.pos]] else data[[2L]]
   pos.name <- if (!is.na(cm.pos)) "CM" else "POS"
   if (is.factor(pos.input))
      pos.input <- as.character(pos.input)
   POS <- suppressWarnings(as.numeric(pos.input))
   if (anyNA(POS) || any(!is.finite(POS)))
      stop(paste(pos.name, "contains invalid positions"))
   if (any(POS < 0))
      stop(paste(pos.name, "contains negative genetic positions"))

   ID <- as.character(data[[3L]])
   if (anyNA(ID) || any(!nzchar(ID)) || anyDuplicated(ID))
      stop("marker IDs are required and must be unique")
   if (any(!grepl("^[A-Za-z0-9]+$", ID)))
      stop("marker IDs may contain only A-Z, a-z, and 0-9")
   if (any(nchar(ID, type="bytes") >= 256L))
      stop("marker IDs must contain fewer than 256 bytes")

   REF <- as.character(data[[4L]])
   if (anyNA(REF) || any(!nzchar(REF)) || any(REF == ".") ||
       any(grepl(",", REF, fixed=TRUE)))
      stop("REF must contain one non-missing allele label per marker")

   ALT <- as.character(data[[5L]])
   if (anyNA(ALT) || any(!nzchar(ALT)))
      stop("ALT may be '.' but may not be missing or empty")

   alt.list <- strsplit(ALT, ",", fixed=TRUE)
   n.alt <- integer(nr)
   for (i in seq_len(nr))
   {
      if ("." == ALT[i])
      {
         n.alt[i] <- 0L
      }
      else
      {
         aa <- alt.list[[i]]
         if (endsWith(ALT[i], ",") ||
             any(!nzchar(aa)) || any(aa == ".") ||
             anyDuplicated(aa) || any(aa == REF[i]))
            stop(paste("invalid or duplicated ALT field at marker", ID[i]))
         n.alt[i] <- length(aa)
      }
   }

   if (any(n.alt > .GS_STVCF_MAX_ALT))
      stop("STvcf supports at most 253 ALT alleles per marker")

   ns <- length(sample.idx)
   allele1 <- matrix(as.raw(255L), nrow=nr, ncol=ns)
   allele2 <- matrix(as.raw(255L), nrow=nr, ncol=ns)

   for (i in seq_len(nr))
   {
      gt <- vapply(data[i, sample.idx, drop=FALSE],
                   function(z) as.character(z[[1L]]),
                   character(1L), USE.NAMES=FALSE)

      if (!is.na(format.pos))
      {
         fmt <- as.character(data[[format.pos]][i])
         fields <- strsplit(fmt, ":", fixed=TRUE)[[1L]]
         gt.pos <- match("GT", fields)
         if (is.na(gt.pos))
            stop(paste("FORMAT does not contain GT at marker", ID[i]))

         for (j in seq_len(ns))
         {
            if (is.na(gt[j]) || !nzchar(gt[j]) || "." == gt[j])
            {
               gt[j] <- NA_character_
            }
            else
            {
               x <- strsplit(gt[j], ":", fixed=TRUE)[[1L]]
               if ((length(x) < gt.pos) || !nzchar(x[gt.pos]))
                  stop(paste("sample FORMAT field does not contain GT at",
                             "marker", ID[i], "sample", sample.names[j]))
               gt[j] <- x[gt.pos]
            }
         }
      }

      a1 <- rep(NA_integer_, ns)
      a2 <- rep(NA_integer_, ns)

      for (j in seq_len(ns))
      {
         x <- gt[j]
         if (is.na(x) || !nzchar(x) || "." == x || "-1" == x ||
             "./." == x || ".|." == x || "-1/-1" == x ||
             "-1|-1" == x || "NA" == x)
            next

         m <- regexec("^([^/|]*)[/|]([^/|]*)$", x)
         p <- regmatches(x, m)[[1L]]
         if (0L == length(p))
            stop(paste("non-diploid or invalid GT at marker", ID[i],
                       "sample", sample.names[j]))

         p <- p[2:3]
         if (any(!nzchar(p)))
            stop(paste("invalid GT at marker", ID[i],
                       "sample", sample.names[j]))

         present <- !(p %in% c(".", "-1"))
         if (any(present) && !all(grepl("^[0-9]+$", p[present])))
            stop(paste("invalid allele index at marker", ID[i],
                       "sample", sample.names[j]))

         a <- rep(NA_integer_, 2L)
         if (any(present))
            a[present] <- suppressWarnings(as.integer(p[present]))
         if (anyNA(a[present]) || any(a[present] < 0L) ||
             any(a[present] > n.alt[i]))
            stop(paste("allele index exceeds REF/ALT definition at marker",
                       ID[i], "sample", sample.names[j]))

         if (!all(present))
            next

         a1[j] <- a[1L]
         a2[j] <- a[2L]
      }

      observed <- unique(c(a1[!is.na(a1)], a2[!is.na(a2)]))
      n.states <- length(observed) + as.integer(any(is.na(a1) | is.na(a2)))
      if (n.states > (.GS_ALLELE_MAX + 1L))
         stop(paste("too many observed allele states at marker", ID[i]))

      x1 <- ifelse(is.na(a1), 255L, a1)
      x2 <- ifelse(is.na(a2), 255L, a2)
      allele1[i, ] <- as.raw(x1)
      allele2[i, ] <- as.raw(x2)
   }

   dimnames(allele1) <- NULL
   dimnames(allele2) <- NULL

   structure(
      list(
         CHROM=as.integer(CHROM),
         POS=as.double(POS),
         ID=ID,
         REF=REF,
         ALT=ALT,
         individual=sample.names,
         allele1=allele1,
         allele2=allele2,
         CHROM.levels=chrom.levels
      ),
      class="STvcf"
   )
}


"st.STvcf.to.dataframe" <- function(STvcf)
{
   required <- c("CHROM", "POS", "ID", "REF", "ALT", "individual",
                 "allele1", "allele2", "CHROM.levels")
   nms <- names(STvcf)
   if (!is.list(STvcf) || (length(STvcf) != length(required)) ||
       is.null(nms) || anyNA(nms) || anyDuplicated(nms) ||
       !setequal(nms, required))
      stop("STvcf must contain exactly the nine required components")

   nr <- length(STvcf$ID)
   ns <- length(STvcf$individual)
   if ((0L == nr) || (0L == ns))
      stop("STvcf contains no markers or no individuals")

   if (!is.integer(STvcf$CHROM) || length(STvcf$CHROM) != nr ||
       !is.double(STvcf$POS) || length(STvcf$POS) != nr ||
       !is.character(STvcf$ID) || length(STvcf$ID) != nr ||
       !is.character(STvcf$REF) || length(STvcf$REF) != nr ||
       !is.character(STvcf$ALT) || length(STvcf$ALT) != nr ||
       !is.character(STvcf$individual) ||
       length(STvcf$individual) != ns ||
       !is.character(STvcf$CHROM.levels))
      stop("invalid STvcf component types or lengths")

   if (anyNA(STvcf$POS) || any(!is.finite(STvcf$POS)) ||
       any(STvcf$POS < 0))
      stop("invalid genetic position in STvcf")

   if (anyNA(STvcf$ID) || any(!nzchar(STvcf$ID)) ||
       anyDuplicated(STvcf$ID) ||
       any(!grepl("^[A-Za-z0-9]+$", STvcf$ID)) ||
       any(nchar(STvcf$ID, type="bytes") >= 256L))
      stop("invalid or duplicated marker ID in STvcf")

   if (anyNA(STvcf$individual) || any(!nzchar(STvcf$individual)) ||
       anyDuplicated(STvcf$individual) ||
       any(!grepl("^[A-Za-z0-9]+$", STvcf$individual)) ||
       any(nchar(STvcf$individual, type="bytes") >= 256L))
      stop("invalid or duplicated individual ID in STvcf")

   if (anyNA(STvcf$REF) || any(!nzchar(STvcf$REF)) ||
       any(STvcf$REF == ".") ||
       any(grepl(",", STvcf$REF, fixed=TRUE)) ||
       anyNA(STvcf$ALT) || any(!nzchar(STvcf$ALT)))
      stop("invalid REF or ALT field in STvcf")

   d1 <- dim(STvcf$allele1)
   d2 <- dim(STvcf$allele2)
   if (!identical(d1, c(nr, ns)) || !identical(d2, c(nr, ns)))
      stop("STvcf allele matrices have invalid dimensions")
   if (!((is.raw(STvcf$allele1) && is.raw(STvcf$allele2)) ||
         (is.integer(STvcf$allele1) && is.integer(STvcf$allele2))))
      stop(paste0("STvcf allele matrices must both use raw or both use ",
                  "integer storage"))

   nchrom <- length(STvcf$CHROM.levels)
   if ((0L == nchrom) || anyNA(STvcf$CHROM) ||
       any(STvcf$CHROM < 1L) || any(STvcf$CHROM > nchrom) ||
       anyNA(STvcf$CHROM.levels) ||
       any(!nzchar(STvcf$CHROM.levels)) ||
       anyDuplicated(STvcf$CHROM.levels))
      stop("invalid chromosome coding in STvcf")
   if (!identical(sort(unique(STvcf$CHROM)), seq_len(nchrom)))
      stop("STvcf chromosome codes must be consecutive and used")

   alt.list <- strsplit(STvcf$ALT, ",", fixed=TRUE)
   n.alt <- integer(nr)
   for (i in seq_len(nr))
   {
      if (STvcf$ALT[i] == ".")
      {
         n.alt[i] <- 0L
      }
      else
      {
         aa <- alt.list[[i]]
         if (endsWith(STvcf$ALT[i], ",") ||
             any(!nzchar(aa)) || any(aa == ".") ||
             anyDuplicated(aa) || any(aa == STvcf$REF[i]))
            stop(paste("invalid or duplicated ALT field in STvcf marker",
                       STvcf$ID[i]))
         n.alt[i] <- length(aa)
      }
   }
   if (any(n.alt > .GS_STVCF_MAX_ALT))
      stop("STvcf supports at most 253 ALT alleles per marker")

   gt <- matrix("./.", nrow=nr, ncol=ns)
   for (i in seq_len(nr))
   {
      if (is.raw(STvcf$allele1))
      {
         a1 <- as.integer(STvcf$allele1[i, ])
         miss1 <- (255L == a1)
      }
      else
      {
         a1 <- STvcf$allele1[i, ]
         miss1 <- is.na(a1)
      }
      if (is.raw(STvcf$allele2))
      {
         a2 <- as.integer(STvcf$allele2[i, ])
         miss2 <- (255L == a2)
      }
      else
      {
         a2 <- STvcf$allele2[i, ]
         miss2 <- is.na(a2)
      }
      bad1 <- !miss1 & ((a1 < 0L) | (a1 > n.alt[i]))
      bad2 <- !miss2 & ((a2 < 0L) | (a2 > n.alt[i]))
      if (any(bad1) || any(bad2))
         stop(paste("invalid allele index in STvcf marker", STvcf$ID[i]))

      miss <- miss1 | miss2
      observed <- unique(c(a1[!miss], a2[!miss]))
      n.states <- length(observed) + as.integer(any(miss))
      if (n.states > (.GS_ALLELE_MAX + 1L))
         stop(paste("too many observed allele states at marker",
                    STvcf$ID[i]))

      ok <- !miss
      gt[i, ok] <- paste0(a1[ok], "/", a2[ok])
   }

   out <- data.frame(
      CHROM=STvcf$CHROM.levels[STvcf$CHROM],
      POS=STvcf$POS,
      ID=STvcf$ID,
      REF=STvcf$REF,
      ALT=STvcf$ALT,
      FORMAT=rep("GT", nr),
      stringsAsFactors=FALSE,
      check.names=FALSE
   )
   gtdf <- as.data.frame(gt, stringsAsFactors=FALSE,
                         check.names=FALSE)
   names(gtdf) <- STvcf$individual
   ans <- cbind(out, gtdf)
   names(ans) <- c(names(out), STvcf$individual)
   ans
}


"read.vcf" <- function(vcffile, data.set="default")
{
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

   status <- .Call("st_read_vcf", as.character(input),
                   as.character(data.set))
   if (!identical(status, 0L))
      stop("read.vcf failed")
   invisible(NULL)
}


"write.vcf" <- function(vcffile, data.set="default")
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

   status <- .Call("st_write_vcf", as.character(output),
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


"st.load.vcf.data" <- function(STvcf, data.set="default")
{
   st.start.timer()
   invisible(.Call("st_load_vcf_data",
                   STvcf,
                   as.character(data.set)))
   st.stop.timer(info.level=1)
   invisible(NULL)
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
   if (anyNA(id) || any(!nzchar(id)) ||
       any(!grepl("^[A-Za-z0-9]+$", id)) ||
       any(nchar(id, type="bytes") >= 256L))
      stop("performance individual IDs are invalid for SelectionTools")

   if (is.factor(y))
      y <- as.character(y)

   st.start.timer()
   invisible(.Call("st_load_performance_data",
                   id,
                   as.double(y),
                   as.character(data.set)))
   st.stop.timer(info.level=1)
   invisible(NULL)
}


"st.get.map" <- function (  filename   = "map",
                            data.set   = "default" )
{
  st.start.timer()
  out.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(filename), "-"),
    tmpdir = tempdir(),
    fileext = ".srt"
  )
  on.exit(unlink(out.filename), add = TRUE)
  retval <- 0
  c <- .C("gs_get_map_01_GV",
          as.character(out.filename),
          retval = as.integer(retval),
          as.character(data.set)
          )
  x <- NULL
  if    (-2 != c$retval ) 
    {
      x <-  read.table(file=out.filename,header=FALSE,check.names=FALSE)
      unlink(out.filename)
      names(x) <- c("Chrom","Pos","Name","Class")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}

"st.restrict.marker.data" <- function ( ind.list  = NULL,
                                        ind.file  = NULL,
                                        mar.list  = NULL,
                                        mar.file  = NULL,
                                        NoAll.MAX = 999,
                                        MaMis.MAX = 1,
                                        ExHet.MIN = 0,
                                        InMis.MAX = 1,
                                        data.set  = "default" )
{
  st.input.dir <- get0("st.input.dir", envir = .GlobalEnv,
                       inherits = FALSE, ifnotfound = "")

  st.start.timer()

  if ( "" == st.input.dir ) ip <- "" else ip <- "/"
  
  if (  ! is.null (ind.file) && ! is.null (ind.list) ) {
      st.info(-2,"Provide either ind.file or ind.list")
      return (invisible(NULL));
    }

  if (  ! is.null (mar.file) && ! is.null (mar.list) ) {
      st.info(-2,"Provide either mar.file or mar.list")
      return (invisible(NULL));
    }
  
  ind.filename <- ""

  if ( ! is.null (ind.file) ) {
    ind.filename <- paste(st.input.dir,ip,ind.file,sep="")
  }

  if (  ! is.null( ind.list) ) {
      ind.filename <- tempfile(
        pattern = "SelectionTools-individuals-",
        tmpdir = tempdir()
      )
      on.exit(unlink(ind.filename), add = TRUE)
      write.table (ind.list,file=ind.filename,quote=FALSE,row.names=FALSE,
                   col.names=FALSE)
    }
  
  mar.filename <- ""
  
  if ( ! is.null( mar.file) ) {
    mar.filename <- paste(st.input.dir,ip,mar.file,sep="")
  }

  if (  ! is.null( mar.list) ) {
    mar.filename <- tempfile(
      pattern = "SelectionTools-markers-",
      tmpdir = tempdir()
    )
    on.exit(unlink(mar.filename), add = TRUE)
    write.table (mar.list,file=mar.filename,quote=FALSE,row.names=FALSE,
                 col.names=FALSE)
  }
  
  retval <- 0;
  c <- .C("gs_restrict_marker_data_03_GV"    ,
          as.character(ind.filename),
          as.character(mar.filename),
          as.integer(NoAll.MAX),
          as.double (MaMis.MAX),
          as.double (ExHet.MIN),
          as.double (InMis.MAX),
          retval = as.integer(retval),
          as.character(data.set)       )
  
  if (  ! is.null( ind.list) ) unlink(ind.filename)
  if (  ! is.null( mar.list) ) unlink(mar.filename)

  st.stop.timer (info.level=1)

}

"st.copy.marker.data" <- function (  target.data.set,
                                     source.data.set )
{

  retval <- 0;
  c <- .C("gs_copy_marker_data_01_GV"   ,
          retval = as.integer(retval)   ,
          as.character(target.data.set) ,
          as.character(source.data.set)    )
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
  #
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
  # Data in list format #
  lpo.filename <- output.filename(lfilename, ".lpo")
  # Map for graphical genotypes #
  mmp.filename <- output.filename(mfilename, ".mmp")
  # Data matrix for graphical genotypes #
  mpo.filename <- output.filename(mfilename, ".mpo")
  # Genome paramteters for Plabsim #
  ngp.filename <- output.filename(nfilename, ".ngp")
  # Map for Plabsim #
  nmp.filename <- output.filename(nfilename, ".nmp")
  # Data matrix for Plabsim #
  npo.filename <- output.filename(nfilename, ".npo")
  # Individual list for sorting and selecting #
  ind.filename <- ""
  if ( "" != ifilename)
    ind.filename <- output.filename(ifilename)
  
  retval <- 0;
  c <- .C("gs_write_marker_data_02_GV",
          as.character(format),
          as.character(lpo.filename),
          as.character(mmp.filename),
          as.character(mpo.filename),
          as.character(ngp.filename),
          as.character(nmp.filename),
          as.character(npo.filename),
          as.character(ind.filename),
          as.integer(f.ind),
          as.integer(l.ind),
          as.character(add.inum),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )

  outp <- NULL
  st.stop.timer (info.level=1)
  return( invisible(outp) )

}

"st.write.map" <- function (         mfilename,
                                     auxfiles   = TRUE,
                                     data.set   = "default" )
{
  st.output.dir <- get0("st.output.dir", envir = .GlobalEnv,
                        inherits = FALSE, ifnotfound = "")

  st.start.timer()
  if ( "" == st.output.dir ) op <- "" else op <- "/"
  mmp.filename <- paste(st.output.dir,op,mfilename,".mmp",sep="")

  retval <- 0;
  c <- .C( "gs_write_map_01_GV",
           as.character(mmp.filename),
           retval = as.integer(retval),
           as.character(data.set)          )

  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=mmp.filename,header=FALSE)
      names(x) <- c("Chrom","Pos","Name","Class")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
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
    c <- .C("gs_recode_ref_name_01_GV",
            as.character(reference),
            as.integer(descending),
            retval = as.integer(retval),
            as.character(data.set)
            )
  } else {
    if ((!is.numeric(reference) && !is.integer(reference)) ||
        (1L != length(reference)) || is.na(reference) ||
        !is.finite(reference) || (reference != floor(reference)) ||
        (reference < 1) || (reference > .Machine$integer.max))
      stop("reference must be one non-empty individual name or one positive integer position")
    c <- .C("gs_recode_ref_01_GV",
            as.integer(reference),
            as.integer(descending),
            retval = as.integer(retval),
            as.character(data.set)
            )
  }

  st.stop.timer (info.level=1)
}

"st.def.hblocks" <- function ( hap          =  5,
                               hap.unit     =  1,
                               ld.threshold = 0.8,
                               ld.criterion = "none",
                               tolerance    = 0,
                               hap.symbol   = "b",
                               out.filename = "blocks"  ,
                               auxfiles     = TRUE,
                               data.set     = "default" )
{
  st.start.timer()
  out.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir(),
    fileext = ".hmap"
  )
  on.exit(unlink(out.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_def_hblocks_01_GV",
          as.integer(hap),
          as.integer(hap.unit),
          as.double(ld.threshold),
          as.character(ld.criterion),
          as.integer(tolerance),
          as.character(hap.symbol),
          as.character(out.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=out.filename,header=FALSE)
      unlink(out.filename)
      names(x) <- c("Chrom","Pos","Name","Class","Markers")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}

"st.recode.hbc" <- function ( reference  = 1,
                             data.set   = "default" )
{
  st.start.timer()
  retval <- 0;
  c <- .C("gs_recode_hbc_01_GV",
          as.integer(reference),
          retval = as.integer(retval),
          as.character(data.set)
          )
  st.stop.timer (info.level=1)
}

"st.recode.hil" <- function ( out.filename = "blocks"  ,
                              auxfiles     = TRUE,
                              data.set   = "default" )
{
  st.start.timer()
  out.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir(),
    fileext = ".hall"
  )
  on.exit(unlink(out.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_recode_hil_01_GV",
          as.character(out.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=out.filename,header=FALSE)
      unlink(out.filename)
      names(x) <- c("Block","AlleleNr","AlleleDef")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
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
  #
  #
   oldmar <- graphics::par("mar")
   on.exit(graphics::par(mar=oldmar), add = TRUE)
  existing.screens <- split.screen()
  if (!identical(existing.screens, FALSE))
    stop("st.plot.ggt cannot be used while split-screen mode is active")
  par (mar=c(0,0,0,0) )
  #
  # Number of individuals
  n.ind <- ncol(g)
  #
  clr <- color
  # Number of chromosomes
  n.chrom <- max(m[,1])
  # Chromosome names
  if ( ( 1 == length(c.nme) ) && ( "" == c.nme[1] ) ) {
    c.nme <- paste("Chr. ",1:n.chrom)
  }
  # Individual names
  if ( ( 1 == length(i.nme) ) && ( "" == i.nme[1] ) ) {
    i.nme <- names(g)
  }
  # Length of chromosomes
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
      for ( j in 2:length(m.pos) ) {
        m.brk[j] <- m.pos[j-1] + (m.pos[j] - m.pos[j-1]) / 2
      }
      m.brk[1+length(m.pos)] <- m.pos[length(m.pos)]
      b.chrom[[c]] <- m.brk
      p.chrom[[c]] <- m.pos
    }
  }
  # Width of displays
  ndc  <- 1 + n.chrom
  f.disp <- c(f.lr*mean(l.chrom),l.chrom)
  f.disp <- f.disp / sum (f.disp)
  s.disp <- f.disp
  for (kk in 2:ndc) s.disp[kk] <- s.disp[kk-1] +s.disp[kk]
  s.disp <- c(0,s.disp)
  #
  dcol <- matrix(0,ncol=2,nrow=ndc)
  for (kk in 1:ndc)
    {
      dcol[kk,1] = s.disp[kk];
      dcol[kk,2] = s.disp[kk+1];
    }

  # Height of displays
  ndr <- 2 + n.ind
  f.disp <- c ( f.tb, rep(1,n.ind) , f.tb)
  f.disp <- f.disp / sum (f.disp)
  s.disp <- f.disp
  for (kk in 2:ndr) s.disp[kk] <- s.disp[kk-1] +s.disp[kk]
  s.disp <- c(0,s.disp)
  #
  drow <- matrix(0,ncol=2,nrow=ndr)
  for (kk in 1:ndr)
    {
      drow[kk,1] = s.disp[kk];
      drow[kk,2] = s.disp[kk+1];
    }

  # Define the figures
  figs <- matrix (ncol=4,nrow=ndr*ndc )
  for (i in 1:ndr) for (c in 1:ndc)
    {
      idx  <- (i-1)*ndc + c
      figs[idx,] <- c(dcol[c,],drow[i,]  )
    }
  #
  sve      <- 1 - figs[,3]
  figs[,3] <- 1 - figs[,4]
  figs[,4] <- sve
  # 
  figs[,1] <- figs[,1] + d.h
  figs[,2] <- figs[,2] - d.h
  figs[,3] <- figs[,3] + d.v
  figs[,4] <- figs[,4] - d.v
  #
  figs <- round(figs,digits=5)
  #
  neu <- split.screen ( figs  )
  on.exit(close.screen(all.screens=TRUE), add=TRUE, after=FALSE)
  #
  # Plot chromosome.names and legend
  for (c in 1:n.chrom)
    {
      scr  <- (c+1)
      screen(scr)
      plot(0,0,t="n",axes=FALSE, xlab="", ylab="")
      text (0,0,c.nme[c],cex=cex.chrom)
#
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
  # Plot individual labels
  for (i in 1:n.ind)
    {
      # Activate screen
      scr  <- (i)*ndc + (1)
      screen(scr)
      plot(0,0,t="n",axes=FALSE, xlab="", ylab="")
      text (0,0,i.nme[i],cex=cex.ind)
    }
  # Plot all chromosomes
  for (i in 1:n.ind) for (c in 1:n.chrom)
    {
      # Activate screen
      scr  <- (i)*ndc + (c+1)
      screen(scr)
     #    
      p.m <- matrix( 0 ,
                   nrow=m.chrom[c] , ncol=2 , 
                   byrow=FALSE                    )
      #
      idx <- m[,1] == c
      m.nme <- as.character ( m[idx,][,3] )
      idx2 <- row.names(g) %in% m.nme
      p.g <- as.character(g[idx2,i])
      #
      if (length(p.g) != m.chrom[c])
        stop("Marker map and genotype matrix do not match")
      #
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
      #
    }
}


"st.plot.ggt" <- function (
                        data.set   = "default",
                        ifilename  = "",
                        i.list     = "",
                        f.ind      = 0,
                        l.ind      = 0,
                        f.lr  = 2,      # Width of left description
                        f.tb  = 2,      # Hight of top/bottom description
                        d.h   = 0.001,  # Distance between chromosomes
                        d.v   = 0.001,  # Distance between individuals
                        p.t   = FALSE    ,  # Plot marker tickmarks
                        p.s   = FALSE    ,  # Plot cm scale
                        d.t   = 50   ,  # Lengh of marker tickmarks
                        d.map = 100  ,  # Distance between marker labels
                        c.nme = ""   ,  # Chromosome names
                        i.nme = ""   ,  # Individual names
                        cex.chrom = 1,  # Scaling factor chromosome names
                        cex.ind   = 1,  # Scaling factor individual names
                        color = c ("yellow","blue","red","green",
                              "wheat", "skyblue","tomato","palegreen",
                              "yellow4","darkblue","darkblue","darkgreen" ),
                        z.min        = 0,
                        z.max        = 12,
                        sort.performance       = FALSE,
                        show.performance       = FALSE,
                        performance.decreasing = TRUE,
                        performance.format     = "%4.2f"
                       )
{
  st.input.dir <- get0("st.input.dir", envir = .GlobalEnv,
                       inherits = FALSE, ifnotfound = "")

  st.start.timer()
#
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"

  mfilename <- tempfile(
    pattern = "SelectionTools-ggt-",
    tmpdir = tempdir()
  )
#
   mmp.filename <- paste0(mfilename, ".mmp")
   mpo.filename <- paste0(mfilename, ".mpo")
   on.exit(unlink(c(mmp.filename, mpo.filename)), add = TRUE)

  if (!is.logical(sort.performance) || (1L != length(sort.performance)) ||
      is.na(sort.performance))
    stop("sort.performance must be TRUE or FALSE")
  if (!is.logical(show.performance) || (1L != length(show.performance)) ||
      is.na(show.performance))
    stop("show.performance must be TRUE or FALSE")
  if (!is.logical(performance.decreasing) ||
      (1L != length(performance.decreasing)) ||
      is.na(performance.decreasing))
    stop("performance.decreasing must be TRUE or FALSE")
  if (!is.character(performance.format) ||
      (1L != length(performance.format)) ||
      is.na(performance.format) || !nzchar(performance.format))
    stop("performance.format must be one non-empty character string")

  performance.data <- NULL
  if (sort.performance || show.performance) {
    per.filename <- tempfile(
      pattern = "SelectionTools-ggt-performance-",
      tmpdir = tempdir()
    )
    on.exit(unlink(per.filename), add = TRUE)
    retval <- 0
    c.per <- .C("gs_return_performance_data_01_GV",
                as.character(per.filename),
                as.integer(TRUE),
                retval = as.integer(retval),
                as.character(data.set))
    if ((-2 == c.per$retval) || !file.exists(per.filename))
      stop(paste("No performance data available for data set '",
                 data.set,"'",sep=""))
    performance.data <- read.table(
      file=per.filename,
      header=FALSE,
      skip=1,
      colClasses=c("character","numeric"),
      stringsAsFactors=FALSE
    )
    unlink(per.filename)
    names(performance.data) <- c("i","y")
  }

  tmp.ifilename <- ""
  i.list.supplied <- !((1L == length(i.list)) && ("" == i.list[1]))

  if (i.list.supplied) {
    i.list <- as.character(i.list)
  } else if ("" != ifilename) {
    i.list <- read.table( paste(st.input.dir,ip,ifilename,sep=""),
                          stringsAsFactors=FALSE)[,1]
    i.list <- as.character(i.list)
    i.list.supplied <- TRUE
  }

  if (sort.performance) {
    if (i.list.supplied)
      sorted.ind <- i.list
    else
      sorted.ind <- performance.data$i

    custom.i.nme <- !((1L == length(i.nme)) && ("" == i.nme[1]))
    sorted.i.nme <- NULL
    if (custom.i.nme) {
      if (length(i.nme) != length(sorted.ind))
        stop(paste("When sort.performance=TRUE, i.nme must contain one",
                   "label for each individual before sorting"))
      sorted.i.nme <- i.nme
    }

    idx <- match(sorted.ind,performance.data$i)
    missing.ind <- unique(sorted.ind[is.na(idx)])
    if (0 < length(missing.ind))
      stop(paste("Individual(s) not found in data set '",data.set,"': ",
                 paste(missing.ind,collapse=", "),sep=""))

    performance <- performance.data$y[idx]
    if (performance.decreasing)
      ord <- order(-performance,seq_along(performance))
    else
      ord <- order(performance,seq_along(performance))
    sorted.ind <- sorted.ind[ord]
    if (custom.i.nme) sorted.i.nme <- sorted.i.nme[ord]

    n.ind <- length(sorted.ind)
    f <- as.integer(f.ind)[1L]
    l <- as.integer(l.ind)[1L]
    if (is.na(f) || is.na(l))
      stop("f.ind and l.ind must be integers")
    if (n.ind < f) f <- n.ind
    if (n.ind < l) l <- n.ind
    if (0 < f) f <- f - 1L else if (0 > f) f <- 0L
    if ((1 > l) || (n.ind < l)) l <- n.ind

    if (f < l) {
      keep <- seq.int(f + 1L,l)
      sorted.ind <- sorted.ind[keep]
      if (custom.i.nme) sorted.i.nme <- sorted.i.nme[keep]
    } else {
      sorted.ind <- character(0)
      if (custom.i.nme) sorted.i.nme <- character(0)
    }

    if (0L == length(sorted.ind))
      stop("No individuals selected for plotting")

    i.list <- sorted.ind
    if (custom.i.nme) i.nme <- sorted.i.nme
    i.list.supplied <- TRUE
    f.ind <- 0
    l.ind <- 0
  }

  if (i.list.supplied) {
    ind.filename <- tempfile(
      pattern = "SelectionTools-ggt-individuals-",
      tmpdir = tempdir()
    )
    on.exit(unlink(ind.filename), add = TRUE)
    write.table ( i.list,
                  file=ind.filename,
                  quote=FALSE,
                  row.names = FALSE, col.names = FALSE)
    tmp.ifilename <- ind.filename
  }
  
  st.write.marker.data ( format     = "m",
                         mfilename  = mfilename,
                         ifilename  = tmp.ifilename,
                         f.ind      = f.ind,
                         l.ind      = l.ind,
                         data.set   = data.set  )
                                        #
  if ( !file.exists(mmp.filename))
    {
      st.info(-2,"No linkage map loaded");
      return ();
    }

  plotted.ind <- scan(mpo.filename,what=character(),nlines=1,quiet=TRUE)
  if (i.list.supplied) {
    missing.ind <- unique(i.list[!i.list %in% plotted.ind])
    if (0 < length(missing.ind))
      stop(paste("Individual(s) not found in data set '",data.set,"': ",
                 paste(missing.ind,collapse=", "),sep=""))
    if (!identical(i.list,plotted.ind))
      stop("Requested individual list does not match plotted individuals")
  }

  retval <- 0;

  m <- read.table(mmp.filename,stringsAsFactors=FALSE, check.names = FALSE,)
  p <- read.table(mpo.filename, check.names = FALSE)

  if (show.performance &&
      ((1L == length(i.nme)) && ("" == i.nme[1]))) {
    idx <- match(plotted.ind,performance.data$i)
    if (anyNA(idx))
      stop("Performance data do not match plotted individuals")
    name.width <- max(nchar(plotted.ind, type="width"))
    i.nme <- paste0(sprintf("%*s",name.width,plotted.ind),": y=",
                    sprintf(performance.format,performance.data$y[idx]))
  }

  if ( !((1 == length(i.nme)) && ("" == i.nme[1])) &&
       (length(i.nme) != ncol(p)) )
    stop("Number of individual names in i.nme does not match plotted individuals")

  unlink(mmp.filename)
  unlink(mpo.filename)

  n0 <-  min( table (m[,1]) ) 
  if ( 1 == n0 )
    {
      st.info(-2,"At least two loci are required on each chromosome ");
      return () ;
    }

  n1 <- max(m[,1])             # largest chrom number
  n2 <- length(unique(m[,1]))  # number of chromosomes
  if ( n1 != n2 )
    {
      st.info(-2,"Chromsosomes must be numbered consecutively starting from 1 ");
      return () ;
    }
  
  st.plot.ggt.src  (
                        map = m , 
                        pop = p , 
                        f.lr,  
                        f.tb,  
                        d.h,   
                        d.v,   
                        p.t,   
                        p.s,   
                        d.t,   
                        d.map, 
                        c.nme, 
                        i.nme, 
                        color,
                        z.min,
                        z.max,
                        cex.chrom,
                        cex.ind
               )
  st.stop.timer (info.level=1)
}


"st.marker.data.statistics" <- function ( data.set ="default",
                                          filename ="marker.stats",
                                          ind = TRUE,
                                          mar = TRUE,
                                          gen = TRUE,
                                          auxfiles = TRUE )
{
  st.start.timer()
  tmp.base <- tempfile(
    pattern = paste0("SelectionTools-", basename(filename), "-"),
    tmpdir = tempdir()
  )
  ind.filename <- paste0(tmp.base, ".ind")
  mar.filename <- paste0(tmp.base, ".mar")
  gen.filename <- paste0(tmp.base, ".gen")
  on.exit(unlink(c(ind.filename, mar.filename, gen.filename)), add = TRUE)
  retval <- 0;
  c <- .C("gs_marker_stats_01_GV",
          as.character(ind.filename),
          as.character(mar.filename),
          as.character(gen.filename),
          as.integer(ind),
          as.integer(mar),
          as.integer(gen),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )

  outp <- NULL

  if  ( (auxfiles) && (-2 != c$retval ) )
    {
       if (ind) {
        individual.list <- read.table ( file=ind.filename,
                                        header=TRUE,
                                        check.names=FALSE,                                        
                                        stringsAsFactors=FALSE)
        outp$individual.list <- individual.list
      }
      if (mar) {
        marker.list      <- read.table ( file=mar.filename,
                                         check.names=FALSE,                                        
                                         header=TRUE,
                                         stringsAsFactors=FALSE  )
      outp$marker.list <- marker.list
      }
      if (gen) {
        genotypes        <- read.table ( file=gen.filename,
                                         check.names=FALSE,                                        
                                         header=TRUE,
                                         stringsAsFactors=FALSE)
        outp$genotypes   <- genotypes
      }
      unlink(ind.filename)
      unlink(mar.filename)
      unlink(gen.filename)
     }
  
  st.stop.timer (info.level=1)
  return(invisible(outp))

}


"gs.return.effects" <- function( out.filename ="effects"    ,
                                 data.set     ="default" )
{
  st.input.dir <- get0("st.input.dir", envir = .GlobalEnv,
                       inherits = FALSE, ifnotfound = "")

  gs.start.timer()
  if ( "" == st.input.dir )  ip <- "" else ip <- "/"
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_return_effects_01_GV"      ,
           as.character(x.filename)    ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  gs.stop.timer (info.level=1)
  if  (-2 == c$retval ) return(invisible(NULL))
  x <- read.table(file=x.filename,header=TRUE,stringsAsFactors=FALSE)
  unlink(x.filename)
  return (x)
}

"gs.set.effects" <- function( eff,
                              filename ="seteffects" ,
                              data.set ="default"  )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  write.table(eff,file=x.filename,quote=FALSE,row.names=FALSE)
  retval <- 0;
  c <- .C( "gs_set_effects_01_GV"      ,
           as.character(x.filename)    ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  gs.stop.timer (info.level=1)
}

"gs.return.pvals" <- function( out.filename ="effects"    ,
                               data.set     ="default" )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_return_pvals_01_GV"      ,
           as.character(x.filename)    ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  gs.stop.timer (info.level=1)
  if  (-2 == c$retval ) return(invisible(NULL))
  x <- read.table(file=x.filename,header=TRUE,stringsAsFactors=FALSE)
  unlink(x.filename)
  return(x)
}

"gs.check.pvals" <- function( data.set ="default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C( "gs_check_pvals_01_GV"      ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  gs.stop.timer (info.level=1)
  return( as.logical(c$retval) )
}

gs.plot.effects <- function (data.set = "default",
                             absolute = TRUE,
                             alpha    = 0.05,
                             p.adjust.method = "none",
                             pvals = NULL,
                             xlim  = NULL,
                             ylim  = NULL,
                             plt   = TRUE,
                             col.sp = "blue", pch.sp = 19 ,
                             col.sn = "blue", pch.sn = 19  , 
                             col.ns = "blue", pch.ns = 1 , 
                             ... ) 
{
  gs.start.timer()
  {
    auxfiles <- TRUE
    x.filename <- tempfile(
      pattern = "SelectionTools-effects-",
      tmpdir = tempdir()
    )
    on.exit(unlink(x.filename), add = TRUE)
    retval <- 0
    pvals.av <- gs.check.pvals(data.set)
    if ( TRUE == pvals.av )
      {
        c <- .C("gs_return_pvals_01_GV", as.character(x.filename), 
                retval = as.integer(retval), as.character(data.set))
      }
    else
      {
        c <- .C("gs_return_effects_01_GV", as.character(x.filename), 
                retval = as.integer(retval), as.character(data.set))
      }
    if (-2 == c$retval)  return(invisible(NULL))
    x <- read.table(file = x.filename, header = TRUE, 
                                check.names = FALSE)
    if (TRUE == pvals.av)
      { x <- data.matrix(x)
        effect <- x[2:nrow(x),1]
        pvalue <- x[2:nrow(x),2]
      } else {
        effect <- x$effect[2:nrow(x)]
      }
    if (!is.null(pvals) )
      { 
        pvalue <- pvals
        pvals.av <- TRUE
      }
    unlink(x.filename)
    {
      nfilename <- tempfile(
        pattern = "SelectionTools-effects-marker-data-",
        tmpdir = tempdir()
      )
      smp.file <- paste0(nfilename, ".nmp")
      spo.file <- paste0(nfilename, ".npo")
      sgp.file <- paste0(nfilename, ".ngp")
      on.exit(unlink(c(smp.file, spo.file, sgp.file)), add = TRUE)
      st.write.marker.data(format = "n", nfilename = nfilename, 
                           data.set = data.set)
      if (!file.exists(sgp.file)) {
        st.info(-2, "No linkage map loaded")
        return (invisible(NULL))
      }
      reset.all()
      cl <- scan(sgp.file, quiet = TRUE)
      set.genome.par(no.chrom = length(cl), no.hom = 2, chrom.len = cl)
      linkage.map.load(smp.file)
      unlink(smp.file)
      unlink(spo.file)
      unlink(sgp.file)
      mp <- get.map()
    }
    if (TRUE==absolute) mp$effect <- abs(effect) else  mp$effect <- effect
    if (TRUE == pvals.av) {
      mp$pvalue   <- pvalue
      mp$pvalue.a <- p.adjust(mp$pvalue,p.adjust.method) 
    }
    c.pos <- cl
    for (i in 2:length(c.pos)) c.pos[i] <- c.pos[i] + c.pos[i-1]
    for (i in 1:nrow(mp)) if (mp$chrom[i] > 1) 
      mp$pos[i] <- mp$pos[i] + c.pos[mp$chrom[i] - 1]
  }
  if (TRUE == plt)
      {
        oldpar <- graphics::par(no.readonly = TRUE)
        on.exit(graphics::par(oldpar), add = TRUE)
        par(mar = c(4, 4, 1, 1))
        if ( is.null(ylim) )
          ylim = c( min( c(min(mp$effect),0) ), max(mp$effect) )
        if ( is.null(xlim) )
          xlim = c( min( c(min(mp$pos),0) ), max(mp$pos) )
        plot(mp$pos, mp$effect, type = "n", 
             xaxt = "n", xlab = "Chromosome", ylab = "Effect size",
             ylim=ylim, xlim=xlim,...)
        c.pos <- c(0, c.pos)
        axis(1, at = c.pos, labels = c(1:length(cl), " "), tck = -0.01)
        grid()
        for (i in 1:(length(c.pos))) {
          x <- c.pos[i]
          lines( c(x,x), ylim )
        }
        if (FALSE == absolute) lines ( xlim,c(0,0))
        tol = 1e-8
        if ( TRUE == pvals.av )
          {
            mp1   <- mp[mp$pvalue.a < alpha, ]
            mp1ap <- mp1[mp1$effect > tol, ]
            mp1an <- mp1[mp1$effect < -tol, ]
            mp1b  <- mp1[abs(mp1$effect) <= tol, ]
            mp2   <- mp[mp$pvalue.a >= alpha, ]
            points(mp1b$pos,  mp1b$effect, type = "p", lwd = 2, col=col.ns, pch=pch.ns)
            points(mp2$pos,   mp2$effect,  type = "p", lwd = 2, col=col.ns, pch=pch.ns)
            points(mp1ap$pos, mp1ap$effect, type = "p", lwd = 2, col=col.sp, pch=pch.sp)
            points(mp1an$pos, mp1an$effect, type = "p", lwd = 2, col=col.sn, pch=pch.sn)
          }
        else
          {
            points(mp$pos, mp$effect, type = "p", lwd = 2, col = "blue")
          }
      }
    ret.val <- data.frame(chrom = mp$chrom, pos = mp$pos, name = mp$name, 
        cpos = mp$pos, effect = mp$effect)
    if ( TRUE == pvals.av ) ret.val$pvalue <- mp$pvalue
    gs.stop.timer(info.level = 1)
    return(invisible(ret.val))
}

gs.compare.effects <- function ( data.set.a,
                                 data.set.b,
                                 label.a      = data.set.a,
                                 label.b      = data.set.b   )
{
  gs.start.timer()
  auxfiles <- TRUE
  x.filename <- tempfile(
    pattern = "SelectionTools-effects-",
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)

  retval <- 0;
  c <- .C( "gs_return_effects_01_GV"      ,
           as.character(x.filename)    ,
           retval = as.integer(retval) ,
           as.character(data.set.a)      )

  if (-2 == c$retval) return (invisible(NULL))
  
  x <- read.table(file=x.filename,header=TRUE) 
  effects.a <- x$effect[2:nrow(x)]

  retval <- 0;
  c <- .C( "gs_return_effects_01_GV"      ,
           as.character(x.filename)    ,
           retval = as.integer(retval) ,
           as.character(data.set.b)      )

  if (-2 == c$retval) return (invisible(NULL))
  
  x <- read.table(file=x.filename,header=TRUE) 
  effects.b <- x$effect[2:nrow(x)]
  
  unlink(x.filename)

  max <- max ( c( abs(effects.a),abs(effects.b) ) ) 
  min <- -max;
  
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  par(mar=c(4,4,1,1))
  plot( effects.b ~effects.a ,
        type="p",
        lwd=2,
        col="blue",
        ylim = c(min,max),
        xlim = c(min,max),
        xlab=label.a,
        ylab=label.b
  )
  lines (c(min,max),c(min,max))
#
  gs.stop.timer (info.level=1)
  return (invisible(NULL))
  
 }

gs.plot.model.fit <- function ( training.set ="default",
                                title        = ""         )
{
  c1 <- gs.predict.genotypes ( training.set   = training.set,
                               prediction.set = training.set )
  if ( is.null(c1) ) return (invisible(NULL) )
  #
  c1 <- c1[,2:3]
  #
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  par(mar=c(4,4,1,1) )
  #
  l <- c(min(c1),max(c1))
  p1 <- l[1] + 0.65*(l[2]-l[1] ) ; p2 <- l[1] + 0.10*(l[2]-l[1] )
  p3 <- l[1] + 0.05*(l[2]-l[1] ) ; p4 <- l[1] + 0.90*(l[2]-l[1] )
  #
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
  #
  c1 <- gs.predict.genotypes ( training.set   = training.set,
                               prediction.set = training.set )
  c2 <- gs.predict.genotypes ( training.set   = training.set,
                               prediction.set = prediction.set )
  #
  if ( 3 == ncol(c1) ) c1 <- c1[,2:3] else return (invisible(NULL) )
  if ( 3 == ncol(c2) ) c2 <- c2[,2:3] else
   {
     st.info (-2,"No phenotypic data in the prediction set")
     return (invisible(NULL) )
   }
  #
  if ( is.null(c1) || is.null(c1) ) return (invisible(NULL) )
  #
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)
  par( mfcol=c(1,2) , mar=c(4,4,1,1) )
  #
  l <- c( min(c1,c2),max(c1,c2) )
  #
  p1 <- l[1] + 0.65*(l[2]-l[1] ) ; p2 <- l[1] + 0.10*(l[2]-l[1] )
  p3 <- l[1] + 0.05*(l[2]-l[1] ) ; p4 <- l[1] + 0.90*(l[2]-l[1] )
  #
  plot(c1,ylim=l,xlim=l,col="blue")
  text(p1,p2,sprintf("r = %4.2f\nrho = %4.2f",cor(c1)[2],
                     cor(c1,method="spearman")[2]),pos=4)
  text(p3,p4,sprintf("%s\nModel fit",title),pos=4)
  lines(l,l)
  lines ( c(mean(c1[,1]),mean(c1[,1])), c(min(c1[,2]),max(c1[,2])) )
  lines ( c(min(c1[,1]),max(c1[,1])), c(mean(c1[,2]),mean(c1[,2])) )
  #
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
  on.exit(graphics::par(oldpar), add = TRUE)
  par ( mar=c(4,4,1,1) )
  #
  lx <- c(min(x),max(x))
  ly <- c(min(y),max(y))
  p1 <- lx[1] + 0.65*(lx[2]-lx[1] ) ; p2 <- ly[1] + 0.10*(ly[2]-ly[1] )
  p3 <- lx[1] + 0.05*(lx[2]-lx[1] ) ; p4 <- ly[1] + 0.90*(ly[2]-ly[1] )
  #
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
  on.exit(graphics::par(oldpar), add = TRUE)
  par ( mar=c(4,4,1,1) )
  #
  lx <- c(min(x),max(x))
  ly <- c(min(y),max(y))
  p1 <- lx[1] + 0.65*(lx[2]-lx[1] ) ; p2 <- ly[1] + 0.10*(ly[2]-ly[1] )
  p3 <- lx[1] + 0.05*(lx[2]-lx[1] ) ; p4 <- ly[1] + 0.90*(ly[2]-ly[1] )
  #
  plot(x,y,xlab=xlab,ylab=ylab,type="n",xlim=xlim,ylim=ylim)
  points(x,y,pch=pch,cex=cex);
  text(p1,p2,sprintf("r = %4.2f\nrho = %4.2f",cor(x,y),
                   cor(x,y,method="spearman")),pos=4)
  
  text(p3,p4,sprintf("%s",title),pos=4)
  lines ( c(lx[1],lx[2]) , c(ly[1],ly[2]))
  lines ( c(mean(x),mean(x)), c(min(y),max(y)) )
  lines ( c(min(x),max(x)), c(mean(y),mean(y)) )
}


st.genetic.distances.02 <- function ( measure = "mrd", # "mrd" "rd" "euc"
                                      format  = "l",   # "l" "m"
                                      filename  = NULL,
                                      auxfiles  = FALSE,
                                      data.set  = "default" )
{
  st.output.dir <- get0("st.output.dir", envir = .GlobalEnv,
                        inherits = FALSE, ifnotfound = "")

   if (TRUE == auxfiles && is.null(filename))
     stop("filename must be supplied when auxfiles = TRUE")
   nfilename <- tempfile(
     pattern = "SelectionTools-genetic-distances-",
     tmpdir = tempdir()
   )
#
   smp.file <- paste0(nfilename, ".nmp")
   spo.file <- paste0(nfilename, ".npo")
   sgp.file <- paste0(nfilename, ".ngp")
   on.exit(unlink(c(smp.file, spo.file, sgp.file)), add = TRUE)
   out.filename <- NULL
   if (TRUE == auxfiles) {
     if ( "" == st.output.dir ) op <- "" else op <- "/"
     out.filename <- paste(st.output.dir,op,filename,".gdi",sep="")
   }
#
  st.write.marker.data ( format     = "n",
                         nfilename  = nfilename,
                         data.set   = data.set  )
  dta  <- read.table(spo.file)
#
  unlink(smp.file)
  unlink(spo.file)
  unlink(sgp.file)
#   
  distance  <- gd.genetic.distance (dta,measure=measure) 
#
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
       p <- gd.data.parameters(dta)
       aa <- distance[,3]
       bb <- p$no.pop
       cc <- p$pop.list
       dd <- structure(aa, Size = bb, Labels = cc, Diag = FALSE,
                       Upper = TRUE, method = "euclidean", class = "dist")
       return ( invisible(dd) )
     }
   else
     {
       return ( invisible(distance) )
     }
}

st.genetic.distances <- function ( measure = "mrd", # "mrd" "rd" "euc"
                                   format  = "l",   # "l" "m"
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
#   
   gs.cross.eval.gd ( dist=measure , data.set=data.set )
   c <- gs.cross.info ( data.set=data.set )
   distance <- data.frame (OTU1=c$P1Name,OTU2=c$P2Name,Measure=c$gd)

#
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

## "st.scriptdir" <- function ( filename )
## {
##    if ( "" == st.script.dir ) sp <- "" else sp <- "/"
##    paste(st.script.dir,sp,filename,sep="")
## }

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

## "st.sd" <- function ( filename )
## {
##    if ( "" == st.script.dir ) sp <- "" else sp <- "/"
##    paste(st.script.dir,sp,filename,sep="")
## }


"gs.cross.eval.gd" <- function ( dist     = "rd",
                                 data.set = "default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C("gs_cross_eval_gd_01_GV"       ,
          as.character(dist)          ,   
          retval = as.integer(retval) ,
          as.character(data.set)       )
  gs.stop.timer (info.level=1)
}


"gs.cross.info.gd" <- function( out.filename ="cross.info.gd",
                                data.set     ="default"        )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_cross_info_gd_01_GV"    ,
           as.character(x.filename)    ,
           retval = as.integer(retval) ,
           as.character(data.set)        )
  x <- NULL
  if  (-2 != c$retval ) {
    x <-  read.table(file=x.filename,header=TRUE,stringsAsFactors=FALSE)
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}


"gs.cross.eval.mu" <- function ( data.set = "default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C("gs_cross_eval_mu_01_GV"    ,
          retval = as.integer(retval) ,
          as.character(data.set)        )
  gs.stop.timer (info.level=1)
}

"gs.cross.eval.mi" <- function ( data.set = "default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C("gs_cross_eval_mi_01_GV"    ,
          retval = as.integer(retval) ,
          as.character(data.set)        )
  gs.stop.timer (info.level=1)
}

"gs.cross.eval.ma" <- function ( data.set = "default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C("gs_cross_eval_ma_01_GV"    ,
          retval = as.integer(retval) ,
          as.character(data.set)        )
  gs.stop.timer (info.level=1)
}


"gs.cross.eval.va" <- function ( data.set     = "default"  ,
                                 pop.type     = "unlinked" ,
                                 t            = 0          ,
                                 map.function = "Haldane"    )
{
  gs.start.timer()
  retval <- 0;
  c <- .C("gs_cross_eval_va_01_GV"    ,
          as.character(pop.type)      ,
          as.integer  (t)             ,
          as.character(map.function)  ,
          retval = as.integer(retval) ,
          as.character(data.set)        )
  gs.stop.timer (info.level=1)
}

"gs.cross.eval.es" <- function ( data.set = "default",
                                 alpha    = 0.1      ,
                                 N        = 0        ,
                                 G        = 0          )
{
  gs.start.timer()

  if ( (0==N) || (0==G) ) {
    i <- dnorm(qnorm(1-alpha,0,1),0,1)/alpha
  } else {
    alpha   <- N/G
    i.alpha <- dnorm(qnorm(1-alpha,0,1),0,1)/alpha
    i       <- i.alpha - ((G-N)/(2*N*(G+1)*i.alpha))
  }

  retval <- 0;
  c <- .C("gs_cross_eval_es_01_GV"       ,
          as.double (i)                  ,
          retval = as.integer(retval)    ,
          as.character(data.set)          )
  gs.stop.timer (info.level=1)
}


"gs.cross.info" <- function( bestn  = 0,
                             sortby = "index",
                             out.filename ="cross.info" ,
                             data.set     ="default" )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_cross_info_01_GV"       ,
           as.integer  (bestn)         ,
           as.character(sortby)        ,
           as.character(x.filename)    ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  x <- NULL
  if  (-2 != c$retval ) {
    x <- read.table(file=x.filename, header=TRUE,
                    colClasses=c("integer","integer","character","character",
                                 rep("numeric",6L)),
                    stringsAsFactors=FALSE)
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"st.chrom.stats" <- function ( filename = "chrom.stats",
                               auxfiles = TRUE,
                               data.set ="default" )
{
  st.start.timer()
  out.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(filename), "-"),
    tmpdir = tempdir(),
    fileext = ".chr"
  )
  on.exit(unlink(out.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_chrom_stats_01_GV",
          as.character(out.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=out.filename,header=FALSE)
      unlink(out.filename)
      names(x) <- c("Chrom","NLoci","Length")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}

"st.calc.rf" <- function ( map.function = "Haldane",
                           auxfiles     = TRUE     ,
                           data.set     = "default" )
{
  st.start.timer()
  out.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(map.function), "-"),
    tmpdir = tempdir(),
    fileext = ".rf"
  )
  on.exit(unlink(out.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_calc_rf_01_GV"         ,
          as.character(map.function) ,
          as.character(out.filename) ,
          as.integer(auxfiles)       ,
          retval = as.integer(retval),
          as.character(data.set)       )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=out.filename,header=FALSE)
      unlink(out.filename)
      names(x) <- c("Chrom","Locus1","Locus2","RecFreq")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}

"st.calc.q" <- function ( pop.type = "DH",
                          t        = 0   , 
                          auxfiles = TRUE,
                          data.set ="default" )
{
  st.start.timer()
  out.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(pop.type), "-"),
    tmpdir = tempdir(),
    fileext = ".q"
  )
  on.exit(unlink(out.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_calc_q_01_GV",
          as.character(pop.type),
          as.integer(t),
          as.character(out.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=out.filename,header=FALSE)
      unlink(out.filename)
      names(x) <- c("Chrom","Locus1","Locus2","q")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
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
   c <- .C("st_get_simpop_GV",
            as.character        (pop.name),
            as.integer          (bg),
            retval = as.integer (retval)  ,
            as.character        (data.set)   )
}

#//+DEBUG++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
   c <- .C("st_get_simpop2_GV",
            as.character        (pop.name),
            as.integer          (bg),
            retval = as.integer (retval)  ,
            as.character        (data.set)   )
}

#//+DEBUG++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



"il.eval.lines" <- function ( allele    = 1,
                              chrom     = 0,
                              lower     = 0,
                              upper     = 0,
                              exclude   = 0,
                              data.set  = "default"  )

{
  if ((length(allele) != 1L) || is.na(allele) ||
      (allele != as.integer(allele)) ||
      (allele < 1L) || (allele > .GS_ALLELE_MAX))
    stop("allele must be an integer from 1 to 254")
 
  #
  mar.stats <- st.marker.data.statistics(data.set=data.set,ind=TRUE,mar=FALSE,gen=FALSE)
  #
  NoLines <- 0
  c  <- .C( "il_get_population_size_GV",
             NoLines = as.integer (NoLines) , 
             as.character         (data.set)   )
  #
  retval <- 0
  abs  <- numeric(c$NoLines)
  tot  <- numeric(c$NoLines)
  nseg <- numeric(c$NoLines)
  lseg <- numeric(c$NoLines)

  if ( !( ( length(chrom) == length(lower) ) &&
          ( length(lower) == length(upper) )    )  )
      {
        st.info(-2,"Incorrect target region borders")
        return ( invisible(NULL) )
      }

  if ( ( length(chrom) < c$NoLines ) ||
       ( length(lower) < c$NoLines ) ||
       ( length(upper) < c$NoLines )    )
      {
        chrom <- rep (chrom[1],c$NoLines)
        lower <- rep (lower[1],c$NoLines)
        upper <- rep (upper[1],c$NoLines)
      }
  
  c <- .C ( "il_eval_lines_GV",
           as.integer (allele),
           as.integer (chrom) ,
           as.numeric (lower) ,
           as.numeric (upper) ,
           as.integer(exclude),
           abs = as.numeric (abs),
           tot = as.numeric (tot),
           nseg = as.numeric (nseg),
           lseg = as.numeric (lseg),
           retval = as.integer (retval), 
           as.character        (data.set) )
  
  if ( -2 == c$retval ) { return ( invisible(NULL) )}

  r <- data.frame ( name=mar.stats$individual.list$Name, nseg=c$nseg, lseg=c$lseg,
                   abs=c$abs,tot=c$tot,rel=(c$abs/c$tot) )
    
    return ( r )

}


il.ideal.library <- function ( n.c =  5  ,            # number of chromosomes 
                               l.c = 100 ,            # length of chromosomes
                               s.l =  20 ,            # segment length   
                               data.set = "default" )
{
  m <- NULL
  for (c in 1:n.c) 
  m <- rbind ( m,data.frame( chrom = rep(c,l.c),
                             pos   = 1:l.c,
                             name  = paste ("C",c,"MA",1:l.c,sep=""),
                             class = rep(paste ("C",c,sep=""),l.c) ))
  m.name <- tempfile(
    pattern = "SelectionTools-ideal-library-map-",
    tmpdir = tempdir()
  )
  p.name <- tempfile(
    pattern = "SelectionTools-ideal-library-population-",
    tmpdir = tempdir()
  )
  on.exit(unlink(c(m.name, p.name)), add = TRUE)
  write.table (m,file=m.name,quote=FALSE,row.names=FALSE,col.names=FALSE,sep =" ")
  x <- matrix("AA",nrow=n.c*l.c,ncol=n.c*l.c/s.l)
  colnames(x) <- paste ( "IL", 1:(n.c*l.c/s.l),  sep="" )
  rownames(x) <- m$name
  i <- 1
  for (p in 1:(n.c*l.c)) {
    x[p,i] <- "CC"
    if (p/s.l==round(p/s.l)) i <- i+1
  }
  write.table (x,file=p.name,quote=FALSE, sep =" ")
  st.read.marker.data (p.name, format="m",auxfiles = FALSE,data.set=data.set)
  st.read.map (m.name, data.set=data.set, auxfiles=FALSE)
  unlink(m.name)
  unlink(p.name)
}

il.overlapping.library <- function ( n.c =  5  ,            # number of chromosomes 
                                     l.c = 100 ,            # length of chromosomes
                                     s.l =  20 ,            # segment length   
                                     data.set = "default" )
{
  m <- NULL
  for (c in 1:n.c ) 
  m <- rbind ( m,data.frame( chrom = rep(c,l.c),
                             pos   = 1:l.c,
                             name  = paste ("C",c,"MA",1:l.c,sep=""),
                             class = rep(paste ("C",c,sep=""),l.c) ))
  m.name <- tempfile(
    pattern = "SelectionTools-overlapping-library-map-",
    tmpdir = tempdir()
  )
  p.name <- tempfile(
    pattern = "SelectionTools-overlapping-library-population-",
    tmpdir = tempdir()
  )
  on.exit(unlink(c(m.name, p.name)), add = TRUE)
  write.table (m,file=m.name,quote=FALSE,row.names=FALSE,col.names=FALSE,sep =" ")
  x <- matrix("AA",nrow=n.c*l.c,ncol=n.c*(1+l.c/s.l))
  colnames(x) <- paste ( "IL", 1:(n.c*(1+l.c/s.l)),  sep="" )
  rownames(x) <- m$name
  i <- 1
  for (p in 1:(n.c*l.c)) {
    x[p,i] <- "CC"
    if (p/s.l==round(p/s.l)) i <- i+1
    if (p/l.c==round(p/l.c)) i <- i+1
  }
  i <- 2
  for (p in 1:(n.c*l.c)) {
    x[p,i] <- "CC"
    if (p/s.l==round(p/s.l)) i <- i+1
    if (p/l.c==round(p/l.c)) i <- i+1
  }

  write.table (x,file=p.name,quote=FALSE, sep =" ")
  st.read.marker.data (p.name, format="m",auxfiles = FALSE,data.set=data.set)
  st.read.map (m.name, data.set=data.set, auxfiles=FALSE)
  unlink(m.name)
  unlink(p.name)
}


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

  retval <- 0
  cov    <- 0
  dep    <- 0
  seglib <- 0
  c <- .C ( "il_eval_library_GV",
            cov    = as.numeric (cov   ),
            dep    = as.numeric (dep   ),
            seglib = as.numeric (seglib),
                     as.integer (dp.allele),
            retval = as.integer (retval), 
            as.character        (data.set) )
  if ( -2 == c$retval ) { return ( invisible(NULL) )}
  r <- c ( cov=c$cov , dep=c$dep , seglib=c$seglib )
  #
  NoLines <- 0
  c  <- .C( "il_get_population_size_GV",
           NoLines = as.integer (NoLines) , 
           as.character         (data.set)   )
  abs  <- numeric(c$NoLines)
  tot  <- numeric(c$NoLines)
  nseg <- numeric(c$NoLines)
  lseg <- numeric(c$NoLines)
  retval <- 1
  #
  c0 <- .C ( "il_eval_lines_GV",
           as.integer (dp.allele),            # Donor allele
           as.integer ( rep(0,c$NoLines) ) ,  # genome wide
           as.numeric ( rep(0,c$NoLines) ) ,
           as.numeric ( rep(0,c$NoLines) ) ,
           as.integer ( 0 ) , 
           abs = as.numeric (abs),
           tot = as.numeric (tot),
           nseg = as.numeric (nseg),
           lseg = as.numeric (lseg),
           retval = as.integer (retval), 
           as.character        (data.set) )
  if ( -2 == c0$retval ) { return ( invisible(NULL) )}
  
  r <- c(r,
         dpg  = mean (c0$abs/c0$tot),
         ndps = mean (c0$nseg)      ,
         ldps = mean (c0$lseg)          )
  if ( 0 == chrom[1] ) {  return ( r ) }

  # Analysis if target regions are defined
  if ( ( length(chrom) != c$NoLines ) ||
       ( length(lower) != c$NoLines ) ||
       ( length(upper) != c$NoLines )    )
    {
      st.info(-2,"Incorrect length of target region borders")
      return ( invisible(NULL) )
    }
  #
  c1 <- .C ( "il_eval_lines_GV",
           as.integer (dp.allele),    # Donor allele
           as.integer (chrom) ,
           as.numeric (lower) ,
           as.numeric (upper) ,
           as.integer ( 0 ),          # in the target regions
           abs = as.numeric (abs),
           tot = as.numeric (tot),
           nseg = as.numeric (nseg),
           lseg = as.numeric (lseg),
           retval = as.integer (retval), 
           as.character        (data.set) )
  if ( -2 == c1$retval ) { return ( invisible(NULL) )}
  #
  c2 <- .C ( "il_eval_lines_GV",
           as.integer (rp.allele),    # Recipient allele
           as.integer (chrom) ,
           as.numeric (lower) ,
           as.numeric (upper) ,
           as.integer ( 1 ),          # outside the target regions
           abs = as.numeric (abs),
           tot = as.numeric (tot),
           nseg = as.numeric (nseg),
           lseg = as.numeric (lseg),
           retval = as.integer (retval), 
           as.character        (data.set) )
  if ( -2 == c2$retval ) { return ( invisible(NULL) )}
  #
  c3 <- .C ( "il_eval_lines_GV",
           as.integer (dp.allele),    # Donor allele
           as.integer (chrom) ,
           as.numeric (lower) ,
           as.numeric (upper) ,
           as.integer ( 1 ),          # outside the target regions
           abs = as.numeric (abs),
           tot = as.numeric (tot),
           nseg = as.numeric (nseg),
           lseg = as.numeric (lseg),
           retval = as.integer (retval), 
           as.character        (data.set) )
  if ( -2 == c3$retval ) { return ( invisible(NULL) )}
  #
  r <- c(r,
         tr.dpg  = mean (c1$abs/c1$tot),
         nt.rpg  = mean (c2$abs/c2$tot),
         nt.ndps = mean (c3$nseg)      ,
         nt.ldps = mean (c3$lseg)          )
  #
  return ( r )

}

  st.select.phen <- function (p=NULL, n=0, t=0, decreasing=TRUE,
                                  nsmall=1, data.set=NULL) {
    if (is.null(p)) {
      if (is.null(data.set))
        stop("either p or data.set must be supplied")
      p <- st.return.performance.data(data.set=data.set)
      if (is.null(p))
        stop(paste("No performance data available for data set '",
                   data.set,"'",sep=""))
    }

    if (!is.data.frame(p))
      stop("p must be a data frame")

    nms <- names(p)
    if ("i" %in% nms) {
      i.col <- "i"
    } else if ("ii" %in% nms) {
      i.col <- "ii"
    } else {
      stop("p must contain an individual column named 'i' or 'ii'")
    }

    if ("y" %in% nms) {
      y.col <- "y"
    } else if ("yy" %in% nms) {
      y.col <- "yy"
    } else {
      stop("p must contain a performance column named 'y' or 'yy'")
    }

    if (!is.numeric(p[[y.col]]))
      stop("performance values must be numeric")

    p <- p[order(p[[y.col]], decreasing=decreasing), , drop=FALSE]
    if (t>0) {
      if (decreasing) {
        p <- p[p[[y.col]]>t, , drop=FALSE]
      } else {
        p <- p[p[[y.col]]<t, , drop=FALSE]
      }
    }
    if ((n>0)&&(n<nrow(p)))
      p <- p[seq_len(n), , drop=FALSE]

    p$descr <- paste(p[[i.col]], format(p[[y.col]], nsmall=nsmall))
    return(p)
  }

st.mark.alleles <- function ( markers = NULL,      # Provide markers directly
                              data.set="default", 
                              effect.set=NULL,     # Significant in another data set
                              target.set,          
                              hap.lst = NULL,      # Second data set uses haplotypes  
                              alpha=0.05,           
                              p.adjust.method="none",
                              replace = "3"
                             )
{
  replace.value <- suppressWarnings(as.integer(replace))
  if ((length(replace) != 1L) || is.na(replace.value) ||
      !grepl("^[0-9]+$", as.character(replace)) ||
      (replace.value < 1L) || (replace.value > .GS_ALLELE_MAX))
    stop("replace must be an allele code from 1 to 254")
  replace <- as.character(replace.value)

  gs.start.timer()
 # Determine marker if not provided 
  if ( is.null(markers) ) 
    {
    if ( is.null(effect.set) ) effect.set <- data.set
    #
    auxfiles <- TRUE
    x.filename <- tempfile(
      pattern = "SelectionTools-pvalues-",
      tmpdir = tempdir()
    )
    on.exit(unlink(x.filename), add = TRUE)
    retval <- 0
    pvals.av <- gs.check.pvals(effect.set)
    if ( FALSE == pvals.av ){
      st.info(-2,"No p-values estimated")
      return(invisible(NULL))
    }
    c <- .C("gs_return_pvals_01_GV", as.character(x.filename), 
            retval = as.integer(retval), as.character(effect.set))
    x <- read.table(file = x.filename, header = TRUE, 
        check.names=FALSE, stringsAsFactors=FALSE)
    unlink(x.filename)
    x <- x[-1,]
    tol = 1e-8
    x <- x[abs(x$effect) > tol, ]
    x$pvalue <- p.adjust(x$pvalue,p.adjust.method)
    markers <- substring(rownames(x[x$pvalue < alpha, ]), 1, 7)
    #
    mark.trans <- NULL
    if ( !is.null(hap.lst) ) {
      for (i in 1:length(markers)) {
        m = as.character(hap.lst$Markers[ markers[i] == hap.lst$Name ])
        mark.trans <- c(mark.trans,unlist(strsplit(m,";")))
      }
      markers <- mark.trans
    }
  }
# Replace marker alleles 
    m.file <- tempfile(
      pattern = "SelectionTools-mark-alleles-",
      tmpdir = tempdir()
    )
    mmp.file <- paste0(m.file, ".mmp")
    fpop <- paste0(m.file, ".mpo")
    on.exit(unlink(c(mmp.file, fpop)), add = TRUE)
    #
    st.write.marker.data ( format = "m", mfilename=m.file,data.set=data.set )
    la <- read.table (fpop, stringsAsFactors=FALSE )
    #
    hom.rep = paste(replace,"/",replace,sep="")
    het.rep = paste("1","/",replace,sep="")
    #
    idx <- row.names(la) %in% markers
    for (i in 1:length(idx)) 
      if (TRUE == idx[i]){
        i2 <-  la[i,] == "2/2";        la[i,i2] <- hom.rep
        i2 <-  la[i,] == "1/2";        la[i,i2] <- het.rep
        }
    #     
    write.table(la,file=fpop,quote=FALSE)
# Load data set
    st.read.marker.data ( fpop,
                          format="m",data.set=target.set,
                          auxfiles=FALSE)
    st.read.map ( mmp.file,
                  data.set=target.set,auxfiles=FALSE)
# Remove temp files
    unlink(fpop)
    unlink(mmp.file)
}


gs.write.pseff <- function (e,beta0=NA,file) {
  marker.effects <- e[2:nrow(e),4]
  if ( is.na(beta0) )
    beta0 <- e[1,4]
  write ( beta0 + sum(marker.effects),file )
  ps.effects <- rbind ( cbind (e[2:nrow(e),c(1,2)],
                               effect=e[2:nrow(e),4] / 2 ) ,
                        cbind ( marker= e[2:nrow(e),1],
                                allele= e[2:nrow(e),3],
                                effect=-e[2:nrow(e),4]/2 ) )
  write.table (ps.effects,file,quote=FALSE,
               row.names=FALSE, col.names=FALSE , append=TRUE )
}

"gs.neg.effall" <- function( data.set ="default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C( "gs_neg_effall_01_GV"      ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  gs.stop.timer (info.level=1)
}

"gs.pos.effall" <- function( data.set ="default" )
{
  gs.start.timer()
  retval <- 0;
  c <- .C( "gs_pos_effall_01_GV"      ,
           retval = as.integer(retval) ,
           as.character(data.set)      )
  gs.stop.timer (info.level=1)
}

################################################################################
# BEGIN: Ridge regression after transformation
################################################################################

"gs.mmet.coeff" <- function(
                              out.filename = "mme.coeff" ,
                              auxfiles     = FALSE,
                              data.set     = "default"
                              )
{
  gs.start.timer()
  A.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(A.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_mmet_coeff_01_GV",
          as.character(A.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)      
          )
  A <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    A <- data.matrix( read.table(file=A.filename,header=FALSE) )
    unlink(A.filename)
    colnames(A) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(A) )
}

"gs.mmet.coeff.3" <- function(
                              out.filename = "mme.coeff" ,
                              auxfiles     = FALSE,
                              data.set     = "default"
                              )
{
  gs.start.timer()
  A.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(A.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_mmet_coeff_03_GV",
          as.character(A.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)      
          )
  A <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    A <- data.matrix( read.table(file=A.filename,header=FALSE) )
    unlink(A.filename)
    colnames(A) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(A) )
}


"gs.mmet.rhs" <- function ( out.filename = "mme.rhs" ,
                              auxfiles     = FALSE      ,
                              data.set     = "default" )
{
  gs.start.timer()
  b.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(b.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_mmet_rhs_01_GV",
          as.character(b.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)        )
  b <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    b <- data.matrix( read.table(file=b.filename,header=FALSE) )
    unlink(b.filename)
    colnames(b) <- NULL
  }
  gs.stop.timer (info.level=1)
  return( invisible(b) )
}

"gs.mmet.solve" <- function( out.filename = "mme.solution" ,
                               auxfiles     = FALSE        ,
                               data.set     = "default"      )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_mmet_solve_01_GV",
          as.character(x.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)        )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}


"gs.lambda.rmla.02" <- function( alpha        = 1        ,
                              maxiter      = 1000     ,
                              precision    = 0.001    ,
                              hsq          = 0.9       ,
                              out.filename = "lambda" ,
                              auxfiles     = TRUE     ,
                              data.set     = "default" )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_lambda_rmla_02_GV"     ,
           as.double(alpha)            ,
           as.integer(maxiter)         ,
           as.double(precision)        ,
           as.double(hsq)              ,
           as.character(x.filename)    ,
           as.integer(auxfiles)        ,
           retval = as.integer(retval) ,
           as.character(data.set)       )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
} 

"gs.w.aov" <- function( alpha        = 1        ,
                        out.filename = "W" ,
                        auxfiles     = TRUE     ,
                        data.set     = "default" )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_w_aov_01_GV"     ,
           as.double(alpha)            ,
           as.character(x.filename)    ,
           as.integer(auxfiles)        ,
           retval = as.integer(retval) ,
           as.character(data.set)       )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
} 


"gs.esteff.rr" <- function( method       = "BLUP",  
                            maxiter      = 1000        ,
                            precision    = 0.0001      ,
                            hsq          = 0.80        ,
                            alpha        = 1           ,
                            out.filename = "eff.rr"    ,
                            auxfiles     = FALSE       ,
                            data.set     = "default"     )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_esteff_rr_01_GV"        ,
           as.character(method)        ,
           as.integer(maxiter)         ,
           as.double(precision)        ,
           as.double(hsq)              ,
           as.double(alpha)            ,
           as.character(x.filename)    ,
           as.integer(auxfiles)        ,
           retval = as.integer(retval) ,
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.cross.validation" <- function ( estimation.method          ,
                                    n.ts                       ,
                                    n.runs                     ,
                                    hsq           = 0.8        ,
                                    alpha         = 1          , 
                                    maxiter       = 100        ,
                                    precision     = 0.0001     ,
                                    out.filename.r  = "crosval"  ,
                                    out.filename.u  = "effects"  ,
                                    auxfiles      = TRUE       ,
                                    data.set      = "default"  )
{
  gs.start.timer()
  r.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename.r), "-"),
    tmpdir = tempdir()
  )
  u.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename.u), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(c(r.filename, u.filename)), add = TRUE)
  retval <- 0;
  c <- .C( "gs_cross_validation_02_GV" ,
           as.character(estimation.method)        ,
           as.integer  (n.ts)  ,
           as.integer  (n.runs)        ,
           as.double   (hsq)           ,
           as.double   (alpha)         ,
           as.integer  (maxiter)       ,
           as.double   (precision)     ,
           as.character(r.filename)    ,
           as.character(u.filename)    ,
           as.integer  (auxfiles)      ,
           retval = as.integer(retval) ,
           as.character(data.set)       )
  outp <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <- data.matrix( read.table(file=r.filename,header=TRUE) )
      unlink(r.filename)
      outp$cor <- x
      x <- data.matrix( read.table(file=u.filename,header=TRUE) )
      unlink(u.filename)
      outp$eff <- x
  }
  
  gs.stop.timer (info.level=1)
  return( invisible(outp) )
}

################################################################################
# END: Ridge regression after transformation
################################################################################

##############################################################################
# BEGIN: LD based haplotypes                                           2015-04
##############################################################################
                                        
"st.calc.ld" <- function ( ld.measure   = "r2" ,
                           auxfiles     = TRUE ,
                           data.set     = "default" )
{
  st.start.timer()
  out.filename <- tempfile(
    pattern = "SelectionTools-ld-",
    tmpdir = tempdir()
  )
  on.exit(unlink(out.filename), add = TRUE)
  retval <- 0;
  c <- .C("gs_calc_ld_01_GV"         ,
          as.character(ld.measure) ,
          as.character(out.filename) ,
          as.integer(auxfiles)       ,
          retval = as.integer(retval),
          as.character(data.set)       )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=out.filename,header=FALSE)
      unlink(out.filename)
      names(x) <- c("Chrom","Locus1","Locus2","Name1","Name2","LD")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}

# For problems with large LD files
"st.calc.ld.2" <- function ( ld.measure   = "r2" ,
                             auxfiles     = TRUE ,
                             keep.ldfile  = FALSE,
                             data.set     = "default",
                             ld.filename  = NULL )
{
  st.output.dir <- get0("st.output.dir", envir = .GlobalEnv,
                        inherits = FALSE, ifnotfound = "")

  if (TRUE == keep.ldfile && TRUE == auxfiles && is.null(ld.filename))
    stop("ld.filename must be supplied when keep.ldfile = TRUE and auxfiles = TRUE")
  st.start.timer()
  if (TRUE == keep.ldfile && TRUE == auxfiles) {
    if ( "" == st.output.dir ) op <- "" else op <- "/"
    out.filename <- paste(st.output.dir,op,ld.filename,sep="")
  } else if (TRUE == keep.ldfile) {
    out.filename <- ""
  } else {
    out.filename <- tempfile(
      pattern = "SelectionTools-ld-",
      tmpdir = tempdir()
    )
    on.exit(unlink(out.filename), add = TRUE)
  }
  retval <- 0;
  c <- .C("gs_calc_ld_01_GV"         ,
          as.character(ld.measure) ,
          as.character(out.filename) ,
          as.integer(auxfiles)       ,
          retval = as.integer(retval),
          as.character(data.set)       )
  x <- NULL
  if (TRUE==keep.ldfile)  {return(invisible(NULL))}
  if  ( (auxfiles) && (-2 != c$retval ) )
      {
          y <- scan(file=out.filename,
                    what=list(integer(),integer(),integer(),character(),character(),double()),
                    quiet = TRUE)
          x <- data.frame(y)
          unlink(out.filename)
          names(x) <- c("Chrom","Locus1","Locus2","Name1","Name2","LD")
      }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}

st.LDplot.ld <- function (ld,chrom)
  {
    l <- ld[ld$Chrom==chrom,]
    n <- nrow(l)
    d <- (1+sqrt(1+8*n))/2
    LD <- matrix(nrow=d,ncol=d)
    for (i in 1:n) LD[ l[i,2],l[i,3] ] <- l[i,6]
    return(LD)
  }

st.LDplot.map <- function(chrom, data.set)
  {
    mp <- st.get.map(data.set=data.set)
    m <- mp[mp$Chrom==chrom,]$Pos
    return(m)
  }

"st.set.hblocks" <- function ( haplotype.list ,
                               hap.symbol   = "s",
                               out.filename = "blocks"  ,
                               auxfiles     = TRUE,
                               data.set     = "default" )
{
  st.start.timer()
  out.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir(),
    fileext = ".hmap"
  )
  on.exit(unlink(out.filename), add = TRUE)
  retval <- 0;
  #
    {
      s <- strsplit(as.character(haplotype.list$Markers),";")
      for (i in 1:length(s)) s[[i]] <-  rep(i,length(s[[i]]))
      hl <- unlist(s)
    } 
  #
  c <- .C("gs_set_hblocks_01_GV",
          as.integer(length(hl)),
          as.integer(hl),
          as.character(hap.symbol),
          as.character(out.filename),
          as.integer(auxfiles),
          retval = as.integer(retval),
          as.character(data.set)
          )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) )
    {
      x <-  read.table(file=out.filename,header=FALSE)
      unlink(out.filename)
      names(x) <- c("Chrom","Pos","Name","Class","Markers")
    }
  st.stop.timer (info.level=1)
  return( invisible(x) )
}


##############################################################################
# END:   LD based haplotypes                                           2015-04
##############################################################################

##############################################################################
# Begin:   Genetid distances for factoials                             2015-11
##############################################################################


st.genetic.distances.fct <- function ( measure = "mrd", # "mrd" "rd" "euc"
                                       format  = "l",   # "l" "m"
                                       split   = 1,     # Last individual of the first set
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
#   
   gs.cross.eval.gd.fct ( dist=measure , split=split, data.set=data.set )
   c <- gs.cross.info ( data.set=data.set )
   distance <- data.frame (OTU1=c$P1Name,OTU2=c$P2Name,Measure=c$gd)

#
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
  c <- .C("gs_cross_eval_gd_fct_01_GV"    ,
          as.character(dist)          ,
          as.integer(split)           ,
          retval = as.integer(retval) ,
          as.character(data.set)       )
  gs.stop.timer (info.level=1)
}

##############################################################################
# End:  Genetic distances for factoials                                2015-11
##############################################################################

##############################################################################
# Begin: Recode reference for two parents                              2017-09
##############################################################################

"st.recode.ref.2" <- function ( r1 = 0,
                                r2 = 1,
                                descending = FALSE,
                                data.set   = "default" )
{
  st.start.timer()
  retval <- 0;
  c <- .C("gs_recode_ref_02_GV",
          as.integer(r1),
          as.integer(r2),
          as.integer(descending),
          retval = as.integer(retval),
          as.character(data.set)
          )
  st.stop.timer (info.level=1)

}

##############################################################################
# End: Recode reference for two parents                              2017-09
##############################################################################

##############################################################################
# Begin: New ridge regression procedures                             2017-11
##############################################################################


st.mxd <- function(
    y,      
    X,      
    Z = NULL,      
    V = NULL,      
    estimates = FALSE, 
    maxcyc = 100,
    precision = 1e-4
    )
{
    if ( nrow(X) != length(y) ) {
        st.info (-2,"Incorrect dimensons of model matrix")
        return();   }
    
    if ( is.null(Z) && is.null(V) ) {             
        st.info (-2,"Provide either Z or V")
        return();   }

    if ( !is.null(Z) && !is.null(V) ) {             
        st.info (-2,"Provide either Z or V")
        return ();   }

    
    if ( !is.null(Z) && !is.list(Z) ) Z <- list(Z) 
    if ( !is.null(V) && !is.list(V) ) V <- list(V)

    n.vc  <- 1 + ( if ( is.null(Z) ) { length(V) } else {length(Z)} )
    sigma = rep ( 0.0, n.vc )

    b = 0.0
    if (estimates) { b = rep(1.0,ncol(X)) } 
    
    u = 0.0
    if ( (! is.null(Z) )  && estimates ) { 
        nu <- 0
        for (i in 1:length(Z)) nu = nu + ncol(Z[[i]])
        u <- rep(1.0,nu)}
    
    maxcyc <- as.integer(maxcyc)
    
    r.lst = list ( sigma=sigma, b=b, u=u)
    retval = as.integer(0);
    
    .Call("st_mxd_interface",
          y,      
	  X,      
	  Z,      
	  V,      
	  maxcyc,
	  precision,
          r.lst,
          retval)

    if (-2 != retval) return(r.lst) else return(NULL)
    
}


"gs.lambda.rrblup" <- function( maxiter     = 1000       ,
                                precision   = 1e-4       ,
                                hsq         = 0.8        ,
                                out.filename = "lambda"  ,
                                auxfiles = TRUE          ,
                                data.set     = "default" )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_lambda_rrblup_01_GV"   ,
           as.integer(maxiter)        ,
           as.double(precision)       ,
           as.double(hsq)             ,
           as.character(x.filename)   ,
           as.integer(auxfiles)       ,
           retval = as.integer(retval),
           as.character(data.set)     )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
} 

"gs.vc.rrblup" <- function( maxiter     = 1000       ,
                                precision   = 1e-4       ,
                                hsq         = 0.8        ,
                                out.filename = "lambda"  ,
                                auxfiles = TRUE          ,
                                data.set     = "default" )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_vc_rrblup_01_GV"   ,
           as.integer(maxiter)        ,
           as.double(precision)       ,
           as.double(hsq)             ,
           as.character(x.filename)   ,
           as.integer(auxfiles)       ,
           retval = as.integer(retval),
           as.character(data.set)     )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
} 



##############################################################################
# Begin: New ridge regression procedures                             2017-11
##############################################################################


##############################################################################
# Begin: EM based genomic prediction                                      2012
##############################################################################

"gs.esteff.lsq" <- function ( out.filename = "eff.lsq" ,
                              auxfiles = TRUE           ,
                              data.set ="default"       )
{
  gs.start.timer()
  #
  Z <- gs.build.Z (data.set=data.set,auxfiles=TRUE)
  r <-qr(Z)$rank
  n <- ncol(Z)-1
  if (r<n) { 
    m <- sprintf("Design matrix is not of full rank (%i < %i)",r,n)
    st.info(-2,m)
    gs.stop.timer (info.level=1)
    return (invisible(NULL))
  }
  #
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_esteff_lsq_01_GV",
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}


"gs.testeff.lsq" <- function ( nperm = 0,
                               out.filename = "test.lsq" ,
                               auxfiles = TRUE           ,
                               data.set ="default"       )
{
  gs.start.timer()
  #
  Z <- gs.build.Z (data.set=data.set,auxfiles=TRUE)
  r <-qr(Z)$rank
  n <- ncol(Z)-1
  if (r<n) { 
    m <- sprintf("Design matrix is not of full rank (%i < %i)",r,n)
    st.info(-2,m)
    gs.stop.timer (info.level=1)
    return (invisible(NULL))
  }
  #
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_testeff_lsq_01_GV",
           as.integer(nperm),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <-  read.table(file=x.filename,header=TRUE)
#   x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.esteff.rmlc" <- function( maxiter     = 1000        ,
                              precision   = 0.0001      ,
                              hsq         = 0.9         ,
                              out.filename = "eff.rmlc" ,
                              auxfiles = FALSE          ,
                              data.set ="default"       )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_esteff_rmlc_01_GV",
           as.integer(maxiter),
           as.double(precision),
           as.double(hsq),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}

"gs.testeff.rmlc" <- function( nperm       = 1000         ,
                               maxiter     = 1000         ,
                               precision   = 0.001        ,
                               hsq         = 0.9          ,
                               out.filename = "test.rmlc" ,
                               auxfiles = TRUE            ,
                               data.set ="default"         )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_testeff_rmlc_01_GV",
           as.integer(nperm),
           as.integer(maxiter),
           as.double(precision),
           as.double(hsq),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <-  read.table(file=x.filename,header=TRUE)
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}


"gs.esteff.rmlv" <- function(  maxiter     = 1000        ,
                               precision   = 0.001       ,
                               hsq         = 0.9         ,
                               out.filename = "eff.rmlv" ,
                               auxfiles = FALSE          ,
                               data.set ="default"       )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_esteff_rmlv_01_GV",
           as.integer(maxiter),
           as.double(precision),
           as.double(hsq),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ){
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
  
}

"gs.testeff.rmlv" <- function( nperm       = 1000        ,
                               maxiter     = 1000        ,
                               precision   = 0.001       ,
                               hsq         = 0.9         ,
                               out.filename = "test.rmlv" ,
                               auxfiles = TRUE           ,
                               data.set ="default"       )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_testeff_rmlv_01_GV",
           as.integer(nperm),
           as.integer(maxiter),
           as.double(precision),
           as.double(hsq),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <-  read.table(file=x.filename,header=TRUE)
#   x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
}


"gs.esteff.rmlr" <- function( maxiter     = 1000        ,
                              precision   = 0.001       ,
                              hsq         = 0.9         ,       
                              out.filename = "eff.rmlr" ,
                              auxfiles = TRUE           ,
                              data.set ="default"       )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_esteff_rmlr_01_GV",
           as.integer(maxiter),
           as.double(precision),
           as.double(hsq),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
  
}

"gs.esteff.rmla" <- function( alpha       = 1        ,
                              maxiter     = 1000        ,
                              precision   = 0.0001       ,
                              hsq         = 0.9         ,
                              out.filename = "eff.rmla" ,
                              auxfiles = TRUE           ,
                              data.set ="default"       )
{
  gs.start.timer()
  x.filename <- tempfile(
    pattern = paste0("SelectionTools-", basename(out.filename), "-"),
    tmpdir = tempdir()
  )
  on.exit(unlink(x.filename), add = TRUE)
  retval <- 0;
  c <- .C( "gs_esteff_rmla_01_GV",
           as.double(alpha),
           as.integer(maxiter),
           as.double(precision),
           as.double(hsq),
           as.character(x.filename),
           as.integer(auxfiles),
           retval = as.integer(retval),
           as.character(data.set)   )
  x <- NULL
  if  ( (auxfiles) && (-2 != c$retval ) ) {
    x <- data.matrix( read.table(file=x.filename,header=TRUE) )
    unlink(x.filename)
  }
  gs.stop.timer (info.level=1)
  return( invisible(x) )
  
}

##############################################################################
# End:  EM based genomic prediction                                       2012
##############################################################################

##############################################################################
# New interface routines to simulations                       Jun 2019
##############################################################################

"st.set.sim.gp" <- function (data.set="default")
    {
        r <- 0
        c <-  .C( "st_set_sim_gp",
                  r = as.integer(r),  
                  as.character(data.set) )
    }

"st.set.sim.mp" <- function (data.set="default")
    {
        r <- 0
        c <-  .C( "st_set_sim_mp",
                  r = as.integer(r),  
                  as.character(data.set) )
    }

"st.set.sim.ef" <- function (name="effect", data.set="default")
    {
        r <- 0
        c <-  .C( "st_set_sim_ef",
                  r = as.integer(r),
                  as.character(name),
                  as.character(data.set) )
    }

"st.set.sim.pp" <- function (pop.name="default", data.set="default")
    {
        if ((1L != length(pop.name)) || is.na(pop.name) ||
            nchar(pop.name, type="bytes") < 2L ||
            nchar(pop.name, type="bytes") >= 256L ||
            !grepl("^[A-Za-z0-9]+$", pop.name))
           stop("pop.name must contain 2 to 255 ASCII letters or digits")
        r <- 0
        c <-  .C( "st_set_sim_pp",
                  r = as.integer(r),
                  as.character(pop.name),
                  as.character(data.set) )
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

##############################################################################
# New implementation of st.plot.genediversity                         Jun 2019
##############################################################################

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
    for ( i in 2:length(c.pos) ) c.pos[i] <- c.pos[i] + c.pos[i-1]
    mp$cpos <- mp$pos
    for ( i in 1:nrow(mp) )
        if (mp$chrom[i] > 1)
            mp$cpos[i] <- mp$cpos[i] + c.pos[ mp$chrom[i]-1]
    cl <- st.chrom.stats(data.set=data.set)$Length
      #
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
                    try(graphics::par(oldpar), silent = TRUE)
                    grDevices::dev.off(plot.device)
                }
            }, add = TRUE)
        }
    else
        {
            oldpar <- graphics::par(no.readonly = TRUE)
            on.exit(graphics::par(oldpar), add = TRUE)
       # dev.new()
        }
# Plot
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
#
    c.pos <- c(0,c.pos)
    axis( 1,
         at=c.pos,
         labels=c(1:length(cl)," "),
         tck=-.01)
                                        
#
# Chromosome borders
    for (i in 1:(length(c.pos))){
        x <- c.pos[i]
        y0 <- 0
        y1  <- 1
        lines(c(x,x),c(y0,y1))
    }
                                        
#
# Gene diversity per chromosome
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

##############################################################################
# Begin: Write performance data to a marker data set,                     2022
##############################################################################

"gs.set.performance.data" <- function( y,
                                       data.set = "default" )
{
  retval <- 0;
  c <- .C( "gs_set_performance_data_01_GV",
          as.integer(length(y)),
          as.double(y),
          retval = as.integer(retval) ,
          as.character(data.set)         )

  return(invisible(c$retval))

}

st.get.simpop.perfdata <- function(pop.name,data.set,effect)
{
    st.get.simpop (pop.name, data.set)
    genotype.population (pop.name)
    evaluate.population (pop.name,effect)
    gv <- get.population.gvalue (pop.name)$gvalue
    gs.set.performance.data (gv,data.set)
}


##############################################################################
# END: Write performance data to a marker data set,                     2022
##############################################################################

##############################################################################
# LD Heatmap replacement                                                 2026
##############################################################################

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
    on.exit(graphics::par(oldpar), add = TRUE)
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


##############################################################################
# END LDheatmap replacement                                              2026
##############################################################################



##############################################################################
# Begin: Genomic prediction                                              2026
##############################################################################
st_mixed <- function(y,
                     id = NULL,
                     Z,
                     Z_ids = NULL,
                     # hybrid-only (optional)
                     pool1 = NULL,
                     pool2 = NULL,
                     hybrid_id = NULL,
                     include_sca = FALSE,
                     ridge = 1e-8,
                     # shared
                     method = c("std", "vanraden", "sommer"),
                     coding = c("detect", "0 0.5 1", "-1 0 1", "0 1 2"),
                     varcmp = NULL,
                     maxcyc = 100L,
                     tol = 1e-10,
                     calcVC = 1,
                     model_type = c("kernel", "marker"),
                     bootstrap = list(nObs=NULL, 
                                      nMarkers=NULL, 
                                      nRep=1, 
                                      seed=NULL))
{
  method <- match.arg(method)
  coding <- match.arg(coding)
  model_type <- match.arg(model_type)

  # normalize bootstrap
  if (is.null(bootstrap) || !is.list(bootstrap)) bootstrap <- list()
  if (is.null(bootstrap$nObs)) bootstrap$nObs <- NULL
  if (is.null(bootstrap$nMarkers)) bootstrap$nMarkers <- NULL
  if (is.null(bootstrap$nRep)) bootstrap$nRep <- 1L
  if (is.null(bootstrap$seed)) bootstrap$seed <- NULL
  
  bootstrap$nRep <- as.integer(bootstrap$nRep)[1]
  if (is.na(bootstrap$nRep) || bootstrap$nRep < 1L) bootstrap$nRep <- 1L
  
  calcVC <- as.logical(calcVC)[1]
  if (is.na(calcVC)) calcVC <- TRUE
  
  # y: accept vector or n x 1 matrix with rownames
  if (is.matrix(y)) {
    if (ncol(y) != 1) stop("y matrix must be n x 1")
    if (is.null(id)) id <- rownames(y)
    y <- y[, 1]
  }
  if (!is.numeric(y)) stop("y must be numeric")
  n <- length(y)
  
  is_hybrid <- !is.null(pool1) || !is.null(pool2) || is.list(Z)
  
  if (!is_hybrid) {
    # ---- single path (unchanged behavior; stubs passed through) ----
    if (is.null(id)) {
      id <- names(y)
      if (is.null(id)) stop("id is NULL but y has no names; provide id or rownames(y)")
    }
    if (!is.character(id)) id <- as.character(id)
    if (length(id) != n) stop("length(id) must equal length(y)")
    
    if (!is.matrix(Z)) stop("Z must be a matrix")
    if (!is.numeric(Z)) Z <- matrix(as.numeric(Z), nrow(Z), ncol(Z), dimnames = dimnames(Z))
    
    if (is.null(Z_ids)) {
      Z_ids <- rownames(Z)
      if (is.null(Z_ids)) stop("Z_ids is NULL and Z has no rownames; provide Z_ids")
    }
    if (!is.character(Z_ids)) Z_ids <- as.character(Z_ids)
    if (length(Z_ids) != nrow(Z)) stop("length(Z_ids) must equal nrow(Z)")
    
    args <- list(
      mode       = "single",
      y          = as.numeric(y),
      id         = id,
      Z          = Z,
      Z_ids      = Z_ids,
      method     = method,
      coding     = coding,
      varcmp     = varcmp,
      maxcyc     = as.integer(maxcyc),
      tol        = as.numeric(tol),
      calcVC     = as.logical(calcVC),
      model_type = model_type,   # stub for later
      bootstrap  = bootstrap     # stub for later
    )

    return(.Call("st_mixed", args))
  }
  
  # ---- hybrid path ----
  if (is.null(pool1) || is.null(pool2)) stop("hybrid mode requires pool1 and pool2")
  if (length(pool1) != n || length(pool2) != n) stop("pool1/pool2 must be same length as y")
  pool1 <- as.character(pool1)
  pool2 <- as.character(pool2)
  
  if (is.null(hybrid_id)) hybrid_id <- paste(pool1, pool2, sep=":")
  if (length(hybrid_id) != n) stop("hybrid_id must be length(y)")
  hybrid_id <- as.character(hybrid_id)
  
  if (anyDuplicated(hybrid_id)) {
    dups <- unique(hybrid_id[duplicated(hybrid_id)])
    stop(sprintf("Duplicated hybrids are not allowed. First few: %s",
                 paste(head(dups, 10), collapse = ", ")))
  }
  
  if (!is.list(Z) || length(Z) != 2) stop("Z must be list(Z_pool1, Z_pool2)")
  Z1 <- Z[[1]]; Z2 <- Z[[2]]
  if (!is.matrix(Z1) || !is.matrix(Z2)) stop("Z[[1]] and Z[[2]] must be matrices")
  if (!is.numeric(Z1)) Z1 <- matrix(as.numeric(Z1), nrow(Z1), ncol(Z1), dimnames=dimnames(Z1))
  if (!is.numeric(Z2)) Z2 <- matrix(as.numeric(Z2), nrow(Z2), ncol(Z2), dimnames=dimnames(Z2))
  
  if (is.null(Z_ids)) {
    Z1_ids <- rownames(Z1); Z2_ids <- rownames(Z2)
    if (is.null(Z1_ids) || is.null(Z2_ids)) stop("Z_ids is NULL but Z matrices lack rownames")
  } else {
    if (!is.list(Z_ids) || length(Z_ids) != 2) stop("Z_ids must be list(ids1, ids2)")
    Z1_ids <- as.character(Z_ids[[1]])
    Z2_ids <- as.character(Z_ids[[2]])
  }
  
  if (!(is.numeric(ridge) && length(ridge) == 1 && ridge >= 0)) stop("ridge must be a single numeric >= 0")
  
  # For deterministic bootstrap in C, set R RNG seed here
  if (identical(model_type, "marker") && isTRUE(include_sca) && !is.null(bootstrap$seed)) {
    set.seed(as.integer(bootstrap$seed)[1])
  }
  
  args <- list(
    mode        = "hybrid",
    y           = as.numeric(y),
    hybrid_id   = hybrid_id,
    pool1       = pool1,
    pool2       = pool2,
    Z1          = Z1,
    Z1_ids      = as.character(Z1_ids),
    Z2          = Z2,
    Z2_ids      = as.character(Z2_ids),
    include_sca = as.logical(include_sca),
    method      = method,
    coding      = coding,
    varcmp      = varcmp,
    maxcyc      = as.integer(maxcyc),
    tol         = as.numeric(tol),
    ridge       = as.numeric(ridge),
    calcVC      = as.logical(calcVC),
    model_type  = model_type,
    bootstrap   = bootstrap
  )
  
  .Call("st_mixed", args)
}

##############################################################################
# END: Genomic prediction                                                 2026
##############################################################################

