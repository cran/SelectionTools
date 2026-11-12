
## Functions

"sm.start.timer" <- function(depth=1)
  {
    st.start.timer(depth)
  }

"sm.stop.timer" <- function(depth=1,info.level=0)
{
  st.stop.timer(depth,info.level)
}

mab.df <- function(d, f){
  if (!is.character(f) || length(f) != 1L || is.na(f) || !nzchar(f))
    stop("f must be a non-empty file name")
  is.absolute <- startsWith(f, "/") ||
                 startsWith(f, "\\") ||
                 (nchar(f) >= 3L &&
                  grepl("^[A-Za-z]$", substr(f, 1L, 1L)) &&
                  substr(f, 2L, 2L) == ":" &&
                  substr(f, 3L, 3L) %in% c("/", "\\"))
  if (is.absolute || is.null(d) || !nzchar(d)) f else file.path(d, f)
}

mab.simulate <- function (
   simulation.name  ,
   simulation.run   ,
   repetitions      = 1000,      
   recurrent.parent = "H",
   donor.parent     = "H",
   linkage.map      ,
   target.loci      = NULL,
   flanking.loci    = NULL,
   recipient.loci   = NULL,
   gen.type         = NULL,
   population.size  = NULL,
   sel.strategy     = NULL,
   no.selected      = NULL,
   no.preselected   = NULL,
   reg.chr          = NULL,                            
   reg.begin        = NULL,                            
   reg.end          = NULL,                            
   missing.allele   = "9",
   result.file      ,
   success.factor   = 1   ,
   recode.infiles   = FALSE,
   recode.file      = NULL
){
  st.input.dir <- get0("st.input.dir", envir=.GlobalEnv, inherits=FALSE, ifnotfound="")
  st.output.dir <- get0("st.output.dir", envir=.GlobalEnv, inherits=FALSE, ifnotfound="")
  st.data.dir <- get0("st.data.dir", envir=.GlobalEnv, inherits=FALSE, ifnotfound="")
  if (!is.character(result.file) || length(result.file) != 1L ||
      is.na(result.file) || !nzchar(result.file))
    stop("result.file must be supplied as a non-empty file name")
  if (recode.infiles && (!is.character(recode.file) || length(recode.file) != 1L ||
                        is.na(recode.file) || !nzchar(recode.file)))
    stop("recode.file must be supplied when recode.infiles = TRUE")
  map <- read.table(mab.df(st.input.dir, linkage.map), stringsAsFactors=FALSE)
  parent1 <- parent2 <- parent.labels <- NULL
  if (recurrent.parent != "H" && donor.parent != "H") {
    parent1 <- read.table(mab.df(st.input.dir, donor.parent), skip=1,
                          colClasses="character", stringsAsFactors=FALSE)
    parent2 <- read.table(mab.df(st.input.dir, recurrent.parent), skip=1,
                          colClasses="character", stringsAsFactors=FALSE)
    parent.labels <- c(scan(mab.df(st.input.dir, donor.parent), what="character", n=1, quiet=TRUE),
                       scan(mab.df(st.input.dir, recurrent.parent), what="character", n=1, quiet=TRUE))
  }
  sm.start.timer()
  outp <- .Call("st_legacy_sim", list(action="mab", simulation.name=simulation.name,
      simulation.run=simulation.run, repetitions=repetitions,
      recurrent.parent=recurrent.parent, donor.parent=donor.parent,
      linkage.map=linkage.map, target.loci=target.loci, flanking.loci=flanking.loci,
      recipient.loci=recipient.loci, gen.type=gen.type, population.size=population.size,
      sel.strategy=as.character(sel.strategy), no.selected=no.selected,
      no.preselected=no.preselected, reg.chr=reg.chr, reg.begin=reg.begin, reg.end=reg.end,
      missing.allele=as.character(missing.allele), success.factor=success.factor,
      map=map, parent1=parent1, parent2=parent2, parent.labels=parent.labels,
      recode.path=if (recode.infiles) mab.df(st.output.dir, recode.file) else NULL),
      PACKAGE="SelectionTools")
  if (!is.null(outp)) save(outp, file=mab.df(st.data.dir, result.file))
  sm.stop.timer()
  invisible(outp)
}


mab.load.data <- function(
   simulation.name  ,
   simulation.run   ,
   result.file
){
  st.data.dir <- get0("st.data.dir", envir = .GlobalEnv,
                      inherits = FALSE, ifnotfound = "")

   if (!is.character(result.file) || length(result.file) != 1L ||
       is.na(result.file) || !nzchar(result.file))
     stop("result.file must be supplied as a non-empty file name")
   outp <- NULL
   load(file=mab.df(st.data.dir,result.file))
   outp
}  



mab.tabulate<- function(
               simulation.names.runs,
               par.summary       = TRUE,
               ext.summary       = FALSE,
               RPG.gw.mean       = FALSE,
               TargetAlleles     = FALSE,
               RPG.gw.Q10        = FALSE,
               RPG.reg.Q10       = FALSE,
               RPG.reg.mean      = FALSE,
               DonorSegments.gw  = FALSE,
               DonorSegments.reg = FALSE,
               LinkageDrag       = FALSE,
               MarkerAnalyses    = FALSE,
               digits            = 1,
               percent           = TRUE,
               result.files
                            )
{
  if (percent) percent <- TRUE else percent <- FALSE
  snr <- matrix(simulation.names.runs, ncol=2, byrow=TRUE)
  n.sc <- nrow(snr)
  if (!is.character(result.files) || length(result.files) != n.sc ||
      anyNA(result.files) || any(!nzchar(result.files)))
    stop("result.files must contain one non-empty file name per scenario")
  outp <- vector("list", n.sc)
  for (i in 1:n.sc)
    outp[[i]] <- mab.load.data(simulation.name=snr[i,1],
                              simulation.run=snr[i,2],
                              result.file=result.files[i])
  if (n.sc > 1)
    for (i in 2:n.sc)
      if (!identical(outp[[i]]$target.loci, outp[[i-1]]$target.loci)) {
        info(-2, "Target loci are not identical in all scenarios")
        return(invisible(NULL))
      }
  if (RPG.reg.mean || RPG.reg.Q10)
    if (n.sc > 1)
      for (i in 2:n.sc)
        if (!(identical(outp[[i]]$reg.chr, outp[[i-1]]$reg.chr) &&
              identical(outp[[i]]$reg.begin, outp[[i-1]]$reg.begin) &&
              identical(outp[[i]]$reg.begin, outp[[i-1]]$reg.begin))) {
          info(-2, "Chromosome regions are not identical in all scenarios")
          return(invisible(NULL))
        }
  .Call("st_legacy_sim", list(action="mab.tabulate", data=outp, scenarios=snr,
      par.summary=par.summary, ext.summary=ext.summary,
      RPG.gw.mean=RPG.gw.mean, TargetAlleles=TargetAlleles, RPG.gw.Q10=RPG.gw.Q10,
      RPG.reg.Q10=RPG.reg.Q10, RPG.reg.mean=RPG.reg.mean,
      DonorSegments.gw=DonorSegments.gw, DonorSegments.reg=DonorSegments.reg,
      LinkageDrag=LinkageDrag, MarkerAnalyses=MarkerAnalyses,
      digits=digits, percent=percent), PACKAGE="SelectionTools")
}


mab.save.inds <- function(pop, pname, file.prefix){
  st.output.dir <- get0("st.output.dir", envir=.GlobalEnv,
                        inherits=FALSE, ifnotfound="")
  if (!is.character(file.prefix) || length(file.prefix) != 1L ||
      is.na(file.prefix) || !nzchar(file.prefix))
    stop("file.prefix must be supplied as a non-empty file name prefix")
  if (length(pop) > 1L) pop <- paste(pop, collapse=" ")
  invisible(.Call("st_legacy_sim", list(action="mab.individuals",
        population=as.character(pop), background=as.integer(plabsim$allele.missing.indicator),
        name=pname, prefix=file.prefix, directory=st.output.dir,
        path.function=mab.df), PACKAGE="SelectionTools"))
}


mab.compare <- function (
   a.simulation.name  ,
   a.simulation.run   ,
   a.repetitions      = 1000,      
   a.recurrent.parent = "H",
   a.donor.parent     = "H",
   a.linkage.map      ,
   a.target.loci      = NULL,
   a.flanking.loci    = NULL,
   a.recipient.loci   = NULL,
   a.gen.type         = NULL,
   a.population.size  = NULL,
   a.sel.strategy     = NULL,
   a.no.selected      = NULL,
   a.no.preselected   = NULL,
   a.reg.chr          = NULL,                            
   a.reg.begin        = NULL,                            
   a.reg.end          = NULL,                            
   a.missing.allele   = "9" ,
   a.success.factor   = 1   ,                      
#
   b.simulation.name =NULL, c.simulation.name =NULL, d.simulation.name =NULL,
   b.simulation.run  =NULL, c.simulation.run  =NULL, d.simulation.run  =NULL,
   b.repetitions     =NULL, c.repetitions     =NULL, d.repetitions     =NULL,      
   b.recurrent.parent=NULL, c.recurrent.parent=NULL, d.recurrent.parent=NULL,
   b.donor.parent    =NULL, c.donor.parent    =NULL, d.donor.parent    =NULL,
   b.linkage.map     =NULL, c.linkage.map     =NULL, d.linkage.map     =NULL,
   b.target.loci     =NULL, c.target.loci     =NULL, d.target.loci     =NULL,
   b.flanking.loci   =NULL, c.flanking.loci   =NULL, d.flanking.loci   =NULL,
   b.recipient.loci  =NULL, c.recipient.loci  =NULL, d.recipient.loci  =NULL,
   b.gen.type        =NULL, c.gen.type        =NULL, d.gen.type        =NULL,
   b.population.size =NULL, c.population.size =NULL, d.population.size =NULL,
   b.sel.strategy    =NULL, c.sel.strategy    =NULL, d.sel.strategy    =NULL,
   b.no.selected     =NULL, c.no.selected     =NULL, d.no.selected     =NULL,
   b.no.preselected  =NULL, c.no.preselected  =NULL, d.no.preselected  =NULL,
   b.reg.chr         =NULL, c.reg.chr         =NULL, d.reg.chr         =NULL,
   b.reg.begin       =NULL, c.reg.begin       =NULL, d.reg.begin       =NULL,
   b.reg.end         =NULL, c.reg.end         =NULL, d.reg.end         =NULL,
   b.missing.allele  =NULL, c.missing.allele  =NULL, d.missing.allele  =NULL,
   b.success.factor  =NULL, c.success.factor  =NULL, d.success.factor  =NULL,
#
   e.simulation.name =NULL, f.simulation.name =NULL, g.simulation.name =NULL,
   e.simulation.run  =NULL, f.simulation.run  =NULL, g.simulation.run  =NULL,
   e.repetitions     =NULL, f.repetitions     =NULL, g.repetitions     =NULL,      
   e.recurrent.parent=NULL, f.recurrent.parent=NULL, g.recurrent.parent=NULL,
   e.donor.parent    =NULL, f.donor.parent    =NULL, g.donor.parent    =NULL,
   e.linkage.map     =NULL, f.linkage.map     =NULL, g.linkage.map     =NULL,
   e.target.loci     =NULL, f.target.loci     =NULL, g.target.loci     =NULL,
   e.flanking.loci   =NULL, f.flanking.loci   =NULL, g.flanking.loci   =NULL,
   e.recipient.loci  =NULL, f.recipient.loci  =NULL, g.recipient.loci  =NULL,
   e.gen.type        =NULL, f.gen.type        =NULL, g.gen.type        =NULL,
   e.population.size =NULL, f.population.size =NULL, g.population.size =NULL,
   e.sel.strategy    =NULL, f.sel.strategy    =NULL, g.sel.strategy    =NULL,
   e.no.selected     =NULL, f.no.selected     =NULL, g.no.selected     =NULL,
   e.no.preselected  =NULL, f.no.preselected  =NULL, g.no.preselected  =NULL,
   e.reg.chr         =NULL, f.reg.chr         =NULL, g.reg.chr         =NULL,
   e.reg.begin       =NULL, f.reg.begin       =NULL, g.reg.begin       =NULL,
   e.reg.end         =NULL, f.reg.end         =NULL, g.reg.end         =NULL,
   e.missing.allele  =NULL, f.missing.allele  =NULL, g.missing.allele  =NULL,
   e.success.factor  =NULL, f.success.factor  =NULL, g.success.factor  =NULL,
#
   h.simulation.name =NULL, i.simulation.name =NULL,
   h.simulation.run  =NULL, i.simulation.run  =NULL,
   h.repetitions     =NULL, i.repetitions     =NULL,     
   h.recurrent.parent=NULL, i.recurrent.parent=NULL,
   h.donor.parent    =NULL, i.donor.parent    =NULL,
   h.linkage.map     =NULL, i.linkage.map     =NULL,
   h.target.loci     =NULL, i.target.loci     =NULL,
   h.flanking.loci   =NULL, i.flanking.loci   =NULL,
   h.recipient.loci  =NULL, i.recipient.loci  =NULL,
   h.gen.type        =NULL, i.gen.type        =NULL,
   h.population.size =NULL, i.population.size =NULL,
   h.sel.strategy    =NULL, i.sel.strategy    =NULL,
   h.no.selected     =NULL, i.no.selected     =NULL,
   h.no.preselected  =NULL, i.no.preselected  =NULL, 
   h.reg.chr         =NULL, i.reg.chr         =NULL,
   h.reg.begin       =NULL, i.reg.begin       =NULL,
   h.reg.end         =NULL, i.reg.end         =NULL,
   h.missing.allele  =NULL, i.missing.allele  =NULL,
   h.success.factor  =NULL, i.success.factor  =NULL,
   result.files
)
{
  n.scn <- 0

  prefix <- c("a.","b.","c.","d.","e.","f.","g.","h.","i.")
  suffix <- c("simulation.name", 
              "simulation.run",  
              "repetitions",     
              "recurrent.parent",
              "donor.parent",    
              "linkage.map",     
              "target.loci",     
              "flanking.loci",
              "recipient.loci",
              "gen.type",        
              "population.size", 
              "sel.strategy",    
              "no.selected", 
              "no.preselected", 
              "reg.chr",     
              "reg.begin",  
              "reg.end",    
              "missing.allele",
              "success.factor" )

  an <- c( "sn",
           "sr",
           "re",
           "rp",
           "dp",
           "li",
           "tl",
           "fl",
           "rl",
           "gt",
           "ps",
           "ss",
           "ns",
           "np",
           "gh",
           "gb",
           "ge",
           "am",
           "sf"
          )

 ### FIX: make variables visible to R CMD check (no behavior change)
  sn <- sr <- re <- rp <- dp <- li <- tl <- fl <- rl <- gt <- ps <- ss <- ns <- np <- gh <-
      gb <- ge <- am <- sf <- vector("list", 9)


  
 for (i in 1:length(an)){
   cmd <- paste(an[i], " <- vector('list',9)",sep="")
   eval(parse(text=cmd))
 }


  for (i in 1:9)
    {
      if (
          eval(parse(text=paste(
           " !is.null(",
           prefix[i],"simulation.name",
          ")" ,
           sep="" )))
          )
        {
  
          n.scn   <- n.scn+1

          sn[[i]] <-
            eval(parse(text=paste(
                         prefix[i],"simulation.name",
                         sep="" )))
          
          for (k in 2:length(an)){
            cmd <- paste( "if( !is.null( a.",suffix[k],") ) ",
                         an[k],"[[",i,"]]","<- a.", suffix[k],
                         sep="" )
            eval(parse(text=cmd))
            cmd <- paste( "if( !is.null(",prefix[i],suffix[k],") ) ", 
                         an[k],"[[",i,"]]", "<-", prefix[i],suffix[k],
                         sep="" )
            eval(parse(text=cmd))
          }
          
        }
    }
 
  
  if (!is.character(result.files) || length(result.files) != n.scn ||
      anyNA(result.files) || any(!nzchar(result.files)))
    stop("result.files must contain one non-empty file name per active scenario")
  if (anyDuplicated(result.files))
    stop("result.files must contain unique file names")

  for (i in 1:n.scn){
    mab.simulate(
                   simulation.name   =  sn[[i]] , 
                   simulation.run    =  sr[[i]] ,
                   repetitions       =  re[[i]] ,
                   recurrent.parent  =  rp[[i]] ,
                   donor.parent      =  dp[[i]] ,
                   linkage.map       =  li[[i]] ,
                   target.loci       =  tl[[i]] ,
                   flanking.loci     =  fl[[i]] ,
                   recipient.loci    =  rl[[i]] ,
                   gen.type          =  gt[[i]] ,
                   population.size   =  ps[[i]] ,
                   sel.strategy      =  ss[[i]] ,
                   no.selected       =  ns[[i]] ,
                   no.preselected    =  np[[i]] ,
                   reg.chr           =  gh[[i]] ,                            
                   reg.begin         =  gb[[i]] ,                            
                   reg.end           =  ge[[i]] ,                            
                   missing.allele    =  am[[i]] ,
                   result.file       =  result.files[i],
                   success.factor    =  sf[[i]] 
                   )
  }
}

