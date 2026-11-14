


plabsim.R.version <- "19.3"
plabsim.R.date    <- "2019/09/11"

plabsim <- list()
plabsim$deprecated.info <- FALSE
plabsim$allele.missing.indicator <- 65535


"info" <-  function(lev,msg){
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="r_info",
        arguments=list(as.integer(lev), as.character(msg))), PACKAGE="SelectionTools");
  }

"info.cat" <- function(lev, msg)
  invisible(.Call("st_legacy_sim", list(action="info.cat", environment=environment()),
                  PACKAGE="SelectionTools"))


"write.version.2" <-
  function(){
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="print_version_2",
        arguments=list(msg = as.character(plabsim.R.date))), PACKAGE="SelectionTools");
  }

"plabsim.init" <- function() {
  rng.init()
}


"population.append" <- function(
  NameP1,
  NameP2
) {
  c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="append_population",
        arguments=list(as.character(NameP1), as.character(NameP2))), PACKAGE="SelectionTools")
}

"population.concat" <- function(
  NameP1,
  NameP2
) {
  c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="concat_population",
        arguments=list(as.character(NameP1), as.character(NameP2))), PACKAGE="SelectionTools")
}

"population.copy" <- function(
  NameP1,
  NameP2,
  start=1,
  n=-1
) {
  c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="copy_population",
        arguments=list(as.character(NameP1), as.character(NameP2), format(start,scientific=FALSE), format(n,scientific=FALSE))), PACKAGE="SelectionTools")
}

"population.divide" <- function(
  NameP1,
  NameP2,
  NoI=-1
) {
  c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="divide_population",
        arguments=list(as.character(NameP1), as.character(NameP2), format(NoI,scientific=FALSE))), PACKAGE="SelectionTools")
}

"population.exist" <- function(PopName) {
   invisible(.Call("st_legacy_sim", list(action="population.exists",
       name=PopName), PACKAGE="SelectionTools"))
}

"population.info.get" <- function(PopName, ind) {
  if ((ind > 0) && (ind <= population.size.get(PopName)$NoInds))
    .Call("st_legacy_sim", list(action="population.info", population=as.character(PopName),
          individual=format(ind,scientific=FALSE)), PACKAGE="SelectionTools")
}


"population.name.swap" <- function(
  NameP1,
  NameP2
) {
  c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="swap_population_name",
        arguments=list(as.character(NameP1), as.character(NameP2))), PACKAGE="SelectionTools")
}

"population.optimize" <-
  function(PopNames){
    c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="optimize_population",
        arguments=list(as.character(PopNames))), PACKAGE="SelectionTools")
  }

"population.remove" <-  function(
  PopNames
) {
  if (length(PopNames) > 1) {
    PopNames <- paste(PopNames, collapse=" ")
  }
  c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="remove_population",
        arguments=list(as.character(PopNames))), PACKAGE="SelectionTools")
}

"population.remove.all" <- function() {
  invisible(.Call("st_legacy_sim", list(action="remove.all"), PACKAGE="SelectionTools"))
}


"population.rename" <- function(
  OldName,
  NewName
) {
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="rename_population",
        arguments=list(as.character(OldName), as.character(NewName))), PACKAGE="SelectionTools")
}

"population.resize" <- function(
  PopName,
  newSize
) {
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="resize_population",
        arguments=list(as.character(PopName), format(newSize,scientific=FALSE))), PACKAGE="SelectionTools")
}

"population.sample" <-
  function(NameP1, NameP2, size=-1, replace=FALSE) {
    if (replace) {
      rep <- 1
    } else {
      rep <- 0
    }
    c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="sample_population",
        arguments=list(as.character(NameP1), as.character(NameP2), format(size,scientific=FALSE), as.integer(rep))), PACKAGE="SelectionTools")
  }

"population.size.get" <- function (PopName) 
.Call("st_legacy_sim", list(action = "output", operation = "size", population = as.character(PopName)), 
    PACKAGE = "SelectionTools")

population.transfer <- function(
  NameP1,
  NameP2,
  population.NameP2.delete=TRUE
) {
  if (population.NameP2.delete) {
    population.concat(NameP1, NameP2)
  } else {
    population.append(NameP1, NameP2)
  }
}

"be.quiet" <-
  function(){
    set.info.level(0);
  }

"cross" <- function(
  NamePg,
  NameP1,
  NameP2,
  NoPg = -1,
  self = 0,
  mode = 1,
  crossing.scheme = 0
) {
  significant <- -1
  c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="cross",
        arguments=list(as.character(NamePg), as.character(NameP1), as.character(NameP2), as.double(self), format(NoPg,scientific=FALSE), as.double(significant), as.integer(mode), as.integer(crossing.scheme))), PACKAGE="SelectionTools")
}

st.cross.plan.execute <- function(plan, population, name, scheme = "DH") invisible(.Call("st_legacy_sim",
    list(action = "cross.plan.execute", plan = plan, population = population, name = name,
        scheme = scheme), PACKAGE = "SelectionTools"))

"define.effects" <- function (name, description) 
invisible(.Call("st_legacy_sim", list(action = "define.effects", name = as.character(name), description = as.character(description)), 
    PACKAGE = "SelectionTools"))

"define.map" <- function (description) 
invisible(.Call("st_legacy_sim", list(action = "define.map", description = as.character(description)), 
    PACKAGE = "SelectionTools"))



"dh" <- function (
  NamePg,
  NameP1,
  NoPg
) {
  significant <- -1;
  c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="dh",
        arguments=list(as.character(NamePg), as.character(NameP1), format(NoPg,scientific=FALSE), as.double(significant))), PACKAGE="SelectionTools")
}

"dh" <- function (
  NamePg,
  NameP1,
  NoPg=-1
) {
  significant <- -1;
  c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="dh",
        arguments=list(as.character(NamePg), as.character(NameP1), format(NoPg,scientific=FALSE), as.double(significant))), PACKAGE="SelectionTools")
}

"evaluate.allele" <- function(PopName, Loci, missing=0)
    .Call("st_legacy_sim", list(action="allele", population=PopName,
          loci=unlist(strsplit(Loci, " ")), missing=missing,
          background=plabsim$allele.missing.indicator, drop.empty=TRUE,
          frequencies=FALSE), PACKAGE="SelectionTools")

"evaluate.allele.freq" <- function(PopName, Loci, missing=0)
{
    if (length(Loci) > 1L) {
        info(-2, "Only available for single loci")
        return(invisible(NULL))
    }
    .Call("st_legacy_sim", list(action="allele", population=PopName,
          loci=unlist(strsplit(Loci, " ")), missing=missing,
          background=plabsim$allele.missing.indicator, drop.empty=TRUE,
          frequencies=TRUE), PACKAGE="SelectionTools")
}

evaluate.genome <- function(PopName, allele, chromosome=NULL, begin=NULL,
                            end=NULL, verbose=FALSE)
    .Call("st_legacy_sim", list(action="genome.evaluate", population=PopName,
          allele=allele, chromosome=chromosome, begin=begin, end=end, verbose=verbose),
          PACKAGE="SelectionTools")

"evaluate.genotype" <- function(PopName, Loci, alleles=0, mode=2, bg=0)
    .Call("st_legacy_sim", list(action="genotype", population=PopName,
          loci=unlist(strsplit(Loci, " ")), alleles=alleles, mode=mode,
          omit=bg, background=plabsim$allele.missing.indicator),
          PACKAGE="SelectionTools")

"evaluate.haplotype" <- function (PopName, Loci, missing = 0, bg = plabsim$allele.missing.indicator) 
{
    genotype.population(PopName)
    y <- .Call("evaluate_haplotype", as.character(PopName), as.character(Loci), as.integer(missing), 
        as.integer(bg), PACKAGE = "SelectionTools")
    if (is.null(y)) 
        stop("incorrect number of dimensions")
    y
}


"evaluate.hap.init" <-
  function(PopName, doubleHetero=0, missing=0){




    bg <- plabsim$allele.missing.indicator
    genotype.population(PopName);
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="haplotype_init",
        arguments=list(as.character(PopName), as.integer(missing), as.integer(doubleHetero), as.integer(bg))), PACKAGE="SelectionTools");
  }


"evaluate.hap.calc" <- function(locusA, locusB)
  .Call("st_legacy_sim", list(action="hap.next", locusA=as.integer(locusA),
        locusB=as.integer(locusB)), PACKAGE="SelectionTools")

"evaluate.hap.calc.d" <- function(locusA, locusB)
  invisible(.Call("st_legacy_sim", list(action="hap.d"), PACKAGE="SelectionTools"))


"evaluate.hap.return.genotype" <- function () 
.Call("haplotype_return_genotype", PACKAGE = "SelectionTools")

"evaluate.hap.return.table" <- function () 
.Call("haplotype_return_table", PACKAGE = "SelectionTools")

"evaluate.hap.free" <-
  function(){
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="haplotype_free",
        arguments=list()), PACKAGE="SelectionTools");
  }

"evaluate.locus" <- function (PopName, eLoci) 
.Call("st_legacy_sim", list(action = "output", operation = "locus", population = as.character(PopName), 
    loci = as.character(eLoci), background = as.integer(plabsim$allele.missing.indicator)), PACKAGE = "SelectionTools")

"evaluate.all.loci" <- function(
  PopName,
  loci = paste(as.vector(unlist(linkage.map.get()[,3])), collapse=" ")
) {
  evaluate.locus(
    PopName,
    loci
  )
}

evaluate.mdp <- function(NamePg, NameP1, NameP2, effectfile=NULL)
    invisible(.Call("st_legacy_sim", list(action="mdp", progeny=NamePg,
              parent1=NameP1, parent2=NameP2, effect=effectfile),
              PACKAGE="SelectionTools"))

"evaluate.population" <- function(
  PopNames,
  effName=NULL
) {
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="evaluate_population",
        arguments=list(as.character(PopNames), as.character(effName))), PACKAGE="SelectionTools")
}

"phenotype.population" <- function(
  PopName,
  hsq
) {
  effects <- list.effects()
  n.effects <- nrow(effects)
  hsq <- as.double(hsq)

  if (0 == n.effects) {
    info(-2, "No effects defined")
    return(invisible(NULL))
  }

  if (length(hsq) != n.effects) {
    info(-2, "One heritability must be supplied for each effect")
    return(invisible(NULL))
  }

  if (any(!is.finite(hsq)) ||
      any(hsq <= 0) ||
      any(hsq > 1)) {
    info(-2, "Heritabilities must be in (0,1]")
    return(invisible(NULL))
  }

  genotype.population(PopName)
  evaluate.population(PopName)

  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="phenotype_population",
        arguments=list(as.character(PopName), as.double(hsq), as.integer(length(hsq)), retval=as.integer(-2))), PACKAGE="SelectionTools")
}


"generate.effect.file" <- function(
  fName,
  description
) {
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="generate_effect_file",
        arguments=list(as.character(fName), as.character(description))), PACKAGE="SelectionTools")
}


"generate.map.file" <- function(
  fName,
  description
) {
  if (length(fName) != 1 ||
      is.na(fName)) {
    info(-2, "Invalid map file name")
    return(invisible(list(retval=-2L)))
  }

  if (length(description) < 1 ||
      anyNA(description)) {
    info(-2, "Invalid map description")
    return(invisible(list(retval=-2L)))
  }

  if (length(description) > 1) {
    description <- paste(
      description,
      collapse=", "
    )
  }

  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="generate_map_file",
        arguments=list(as.character(fName), as.character(description), retval=as.integer(-2))), PACKAGE="SelectionTools")

  invisible(c)
}



generate.population <- function(dta, backcross=FALSE)
    invisible(.Call("st_legacy_sim", list(action="data", operation="generate",
          data=dta, backcross=backcross, background=plabsim$allele.missing.indicator),
          PACKAGE="SelectionTools"))

"genome.contribution" <- function(pops, allele, chromosome=NULL,
                                      begin=NULL, end=NULL)
    .Call("st_legacy_sim", list(action="summary", kind="ratio",
          pops=unlist(strsplit(pops, " ")), allele=allele,
          chromosome=chromosome, begin=begin, end=end), PACKAGE="SelectionTools")

"genotype.population" <-
  function(PopNames){
    if (length(PopNames) > 1) {
      PopNames <- paste(PopNames, collapse=" ")
    }
    c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="genotype_population",
        arguments=list(as.character(PopNames))), PACKAGE="SelectionTools")
  }


"genome.parameter.get" <- function () 
.Call("st_legacy_sim", list(action = "output", operation = "genome"), PACKAGE = "SelectionTools")

"linkage.map.get" <- function () 
.Call("st_legacy_sim", list(action = "output", operation = "map"), PACKAGE = "SelectionTools")

"get.mdp" <- function () 
.Call("st_legacy_sim", list(action = "output", operation = "mdp"), PACKAGE = "SelectionTools")

"set.mdp" <- function(
  MDP=0
) {
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="set_mdp",
        arguments=list(format(MDP,scientific=FALSE))), PACKAGE="SelectionTools")
}

"get.population" <- function (PopName) 
.Call("st_legacy_sim", list(action = "output", operation = "population", population = as.character(PopName)), 
    PACKAGE = "SelectionTools")

"get.population.gvalue" <- function (name, EffName = NULL) 
.Call("st_legacy_sim", list(action = "output", operation = "gvalue", population = as.character(name), 
    effect = EffName), PACKAGE = "SelectionTools")


"get.population.pvalue" <- function (name, EffName = NULL) 
.Call("st_legacy_sim", list(action = "output", operation = "pvalue", population = as.character(name), 
    effect = EffName), PACKAGE = "SelectionTools")


"get.score" <-
  function(pop,effectfile=NULL){
    genotype.population(pop);
    evaluate.population(pop);
    g <- get.population.gvalue(pop,effectfile)
    g
  }


homozygote <- function(allele, NoInd=1)
    .Call("st_legacy_sim", list(action="plant", genotype=allele, n=NoInd,
          constant=TRUE), PACKAGE="SelectionTools")

plant <- function(genotype, NoInd=1)
    .Call("st_legacy_sim", list(action="plant", genotype=genotype, n=NoInd,
          constant=FALSE), PACKAGE="SelectionTools")

"init.population" <- function(name, data, delete=TRUE) {
  invisible(.Call("st_legacy_sim", list(action="population.init",
      population=as.character(name), data=data, delete=delete),
      PACKAGE="SelectionTools"))
}

"list.effects" <- function () 
.Call("st_legacy_sim", list(action = "output", operation = "effects"), PACKAGE = "SelectionTools")

"population.list" <- function () 
.Call("st_legacy_sim", list(action = "output", operation = "populations"), PACKAGE = "SelectionTools")

"load.effmap" <- function(
  name,
  file = NA
) {
  if (is.na(file)) {
    file <- name
  }
  allPops <- paste(as.vector(population.list()$PopName),collapse = " ")
  remove.evaluate.population(allPops)
  c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="load_effmap",
        arguments=list(as.character(file), as.character(name))), PACKAGE="SelectionTools")
}

"load.internal.effmap" <-
  function(name, spec, weight){
    c <- population.list();
    y <- as.vector((c$PopName));
    allPops <- paste(y,collapse = " ");
    remove.evaluate.population(allPops)
    b  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="load_effmap_mod",
        arguments=list(as.character(name), as.character(spec), as.double(weight))), PACKAGE="SelectionTools");
  }

"linkage.map.load" <-
  function(
    file,
    disperse=TRUE,
    file.disperse=NA,
    disperse.factor=100
  ) {
    if (length(file) != 1 ||
        is.na(file)) {
      info(-2, "Invalid linkage map file")
      return(invisible(list(retval=-2L)))
    }

    if (length(disperse) != 1 ||
        is.na(disperse)) {
      info(-2, "Invalid disperse argument")
      return(invisible(list(retval=-2L)))
    }

    if (length(file.disperse) != 1) {
      info(-2, "Invalid dispersed map file name")
      return(invisible(list(retval=-2L)))
    }

    if (length(disperse.factor) != 1) {
      info(-2, "Invalid disperse.factor")
      return(invisible(list(retval=-2L)))
    }

    if (is.na(file.disperse)) {
      file.disperse <- ""
    } else {
      file.disperse <- as.character(file.disperse)
    }

    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="linkage_map_load",
        arguments=list(as.character(file), as.integer(as.logical(disperse)), as.character(file.disperse), as.double(disperse.factor), retval=as.integer(-2))), PACKAGE="SelectionTools")

    invisible(c)
  }

"linkage.map.save" <- function(file, dta.map=linkage.map.get()) {
  utils::write.table(
    dta.map,
    file,
    row.names = FALSE,
    col.names = FALSE,
    quote     = FALSE
  )
}

"population.individual.remove" <- function(PopName, IndPos) {
  invisible(.Call("st_legacy_sim", list(action="population.individual.remove",
      population=as.character(PopName), position=IndPos), PACKAGE="SelectionTools"))
}

"effmap.remove" <-
  function(Name){
    c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="effmap_remove",
        arguments=list(as.character(Name))), PACKAGE="SelectionTools")
  }

"effmap.remove.all" <-
  function(){
    c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="effmap_remove_all",
        arguments=list()), PACKAGE="SelectionTools");
  }

"remove.effmaps" <-
  function(){
    effmap.remove.all()
  }

"remove.evaluate.population" <-
  function(PopNames){
    c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="remove_evaluate_population",
        arguments=list(as.character(PopNames))), PACKAGE="SelectionTools")
  }

"remove.genotype.population" <-
  function(PopNames){
    c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="remove_genotype_population",
        arguments=list(as.character(PopNames))), PACKAGE="SelectionTools")
  }


"remove.map" <-
  function(){
    c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="remove_map",
        arguments=list()), PACKAGE="SelectionTools")
  }

"reset.all" <-
  function(){
    population.remove.all();
    remove.effmaps();
    remove.map();
    reset.mdp();
  }

"reset.mdp" <-
  function(){
    c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="reset_mdp",
        arguments=list()), PACKAGE="SelectionTools");
  }

"resources" <-
  function(expr) {
    expr <- substitute(expr)
    stime <- proc.time()
    w <- eval(expr)
    etime <- proc.time()
    on.exit()
    time <- etime - stime
    time[3] <- max(time[3], time[1] + time[2])
    info(0, c(CPU = time[1] + time[2],
            Elapsed = time[3],
            "% CPU" = round((100 * (time[1] + time[2]))/time[3], 1),
           Child = time[4] + time[5]))
    invisible(w)
  }


"return.population" <- function(PopNames, missing.rm = TRUE) {
    if (1 < length(PopNames)) {
        PopNames <- paste(PopNames, collapse = " ")
    }
    if (missing.rm) {
        bg <- plabsim$allele.missing.indicator
    }
    else {
        bg <- 3334
    }
    genotype.population(PopNames)
    y <- .Call("return_population", as.character(PopNames), as.integer(bg))
    if (is.null(y)) 
        return(invisible(NULL))
    y
}

"select.all.best" <- function(
  newPop,
  oldPop,
  effect = NULL,
  x = 1,
  decreasing = TRUE,
  selection.criterion = "G"
) {
  selection.criterion <- match.arg(
    selection.criterion,
    c("G", "P")
  )

  popsize <- population.size.get(oldPop)$NoInds
  if (is.null(popsize) || (0 == popsize)) {
    info(-2,paste("Population",oldPop,"is of size 0"))
    return(invisible(NULL))
  }

  if ("G" == selection.criterion) {
    genotype.population(oldPop)
    evaluate.population(oldPop,effect)

    c <- population.sort(
      oldPop,
      decreasing=decreasing,
      selection.criterion="G"
    )
    if (is.null(c$retval) || 0L != c$retval)
      return(invisible(c))

    s <- select.all.best.intern(
      oldPop,
      x,
      selection.criterion="G"
    )
  } else {
    c <- population.sort(
      oldPop,
      decreasing=decreasing,
      selection.criterion="P",
      effect=effect
    )
    if (is.null(c$retval) || 0L != c$retval)
      return(invisible(c))

    s <- select.all.best.intern(
      oldPop,
      x,
      selection.criterion="P",
      effect=effect
    )
  }

  if (is.null(s))
    return(invisible(NULL))

  population.divide(newPop,oldPop,s$n)
  invisible(s)
}


"select.all.best.intern" <- function (NamePop, n = 1, selection.criterion = "G", effect = NULL) 
.Call("st_legacy_sim", list(action = "output", operation = "select.best", population = as.character(NamePop), 
    n = as.integer(n), criterion = match.arg(selection.criterion, c("G", "P")), effect = effect), PACKAGE = "SelectionTools")


"select.n.best" <- function(
  newPop,
  oldPop,
  effect=NULL,
  n=1,
  decreasing=TRUE,
  selection.criterion="G"
) {
  selection.criterion <- match.arg(selection.criterion,c("G","P"))
  invisible(.Call("st_population_select_n_best",enc2native(as.character(newPop)),
      enc2native(as.character(oldPop)),as.character(effect),n,format(n,scientific=FALSE),
      as.integer(decreasing),selection.criterion,PACKAGE="SelectionTools"))
}

"set.co.freq" <-
  function(cofreq=1){
    c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="set_co_freq",
        arguments=list(as.double(cofreq))), PACKAGE="SelectionTools")
  }

"set.NoLociInit" <-
  function(NoLociInit){
    c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="set_nolociinit",
        arguments=list(as.integer(NoLociInit))), PACKAGE="SelectionTools")
  }


"effect.weight.set.all" <- function(
  weight
) {
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="effect_weight_set_all",
        arguments=list(as.double(weight))), PACKAGE="SelectionTools")
}


"set.eff.weight" <- function(
  fname,
  weight
) {
  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="set_eff_weight",
        arguments=list(as.character(fname), as.double(weight))), PACKAGE="SelectionTools")
}

"genome.parameter.set" <- function(
  no.chrom,
  no.hom=2,
  chrom.len
) {
  if (length(chrom.len) != no.chrom) {
    info(-2, "More chromosomes than chromosome lenghts specified")
  } else {
    population.remove.all();
    remove.map();
    c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="set_genome_par",
        arguments=list(as.integer(no.chrom), as.integer(no.hom), as.double(chrom.len))), PACKAGE="SelectionTools")
  }
}

"set.info.level" <- function(
  level
) {
  c  <-  .Call("st_legacy_marker", list(action="legacy.result",
        routine="set_info_level",
        arguments=list(as.integer(level))), PACKAGE="SelectionTools")
}

"set.population.gvalue" <- function(
  name,
  gvalue
) {
  popsize <- population.size.get(name)
  gvalue <- as.double(gvalue)

  if (is.null(popsize)) {
    info(-2,paste("Population",name,"is not defined"))
    return(invisible(NULL))
  }

  if (length(gvalue) != popsize$NoInds) {
    info(-2,"Length of gvalue does not match population size")
    return(invisible(NULL))
  }

  if (any(!is.finite(gvalue))) {
    info(-2,"Genotypic values must be finite")
    return(invisible(NULL))
  }

  if (0 < nrow(list.effects())) {
    genotype.population(name)
    evaluate.population(name)
  }

  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="set_population_gvalue",
        arguments=list(as.character(name), as.double(gvalue), as.integer(length(gvalue)), retval=as.integer(-2))), PACKAGE="SelectionTools")
}


"population.info.set" <-
  function(PopName,ind,info){
    if ((ind > 0) && (ind <= population.size.get(PopName)$NoInds)) {
      c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="set_population_info",
        arguments=list(as.character(PopName), format(ind,scientific=FALSE), as.character(info))), PACKAGE="SelectionTools")
    }
  }

single.cross <- function(PgName, PopName, classSize=1)
    invisible(.Call("st_legacy_sim", list(action="single.cross", target=PgName,
              source=PopName, size=classSize), PACKAGE="SelectionTools"))

"population.sort" <- function(
  PopName,
  decreasing=TRUE,
  selection.criterion="G",
  effect=NULL
) {
  selection.criterion <- match.arg(
    selection.criterion,
    c("G", "P")
  )
  effect.name <- if (is.null(effect)) "" else as.character(effect)

  c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="sort_population",
        arguments=list(as.character(paste(PopName,sep=" ")), as.integer(decreasing), as.character(selection.criterion), effect.name, retval=as.integer(-2))), PACKAGE="SelectionTools")
}


"ssd.mating" <- function(
  NamePg,
  NameP,
  NoPg=1,
  maxcycles=1
) {
  significant <- -1;
  c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="ssd_mating",
        arguments=list(as.character(NamePg), as.character(NameP), format(NoPg,scientific=FALSE), as.integer(maxcycles), as.double(significant))), PACKAGE="SelectionTools")
}

"summarize.gvalue" <- function(pops, effectfile=NULL)
    .Call("st_legacy_sim", list(action="summary", kind="gvalue",
          pops=unlist(strsplit(pops, " ")), effect=effectfile),
          PACKAGE="SelectionTools")

"talk.to.me" <-
  function(){
    set.info.level(1);
  }

data.params <- function(dta) gd.data.parameters(dta)


"splitdt" <-
  function(m.a) {
    .Call("splitdt", as.character(m.a) )
  }

"plabsim.version" <- function() {
  write.version.2()
}


"population.matrix.load" <- function(
  file,
  backcross=FALSE,
  keep.IndID=TRUE
) {
  dta <- utils::read.table(file,header=TRUE)
  generate.population(dta,backcross=backcross)

  if (keep.IndID) {
    dta.param <- gd.data.parameters(dta)
    line <- 1
    for (p in 1:dta.param$no.pop) {
      for (i in 1:dta.param$no.ind[1,p]) {
        set.population.info(
          dta.param$pop.list[p],
          i,
          dta.param$colhead[line,2]
        )
        line <- line + 1
      } 
    } 
  }
}

"population.matrix.save" <- function (file,
                                      PopNames,
                                      missing.rm=TRUE,
                                      keep.IndID=TRUE) {
  dta <- return.population(PopNames, missing.rm)

  if (keep.IndID) {
    dta.param <- gd.data.parameters(dta)
    col.names <- NULL
    for (p in 1:dta.param$no.pop) {
      for (i in 1:dta.param$no.ind[1, p]) {
        indid <- as.integer(population.info.get(dta.param$pop.list[p], i))
        if (is.na(indid)) {
          info(-2,"IndID not numeric")
        } else {
        col.names <- c(col.names,
                       paste(dta.param$pop.list[p],
                             indid,
                             sep=".")
                       )
      }
      } 
    } 
    colnames(dta) <- col.names
  }

  utils::write.table(dta,
              file = file,
              quote = FALSE)
}

"population.plabsim.load" <- function(file, PopName=NA) {
  dta.description <- utils::read.table(
    file,
    comment.char = "",
    nrows        = 1
  )
  PopName <- ifelse(
    is.na(PopName),
    as.character(dta.description[1,2]),
    PopName
  )
  dta.info <- utils::read.table(
    file,
    comment.char = "",
    skip         = 1,
    nrows        = dta.description[1,3]
  )
  dta.info <- as.character(dta.info[,3])
  col.classes <- c("integer",      
                   "integer",      
                   "integer",      
                   "numeric",      
                   "integer")      

  dta <- utils::read.table(
    file,
    colClasses=col.classes,
    skip = 1+dta.description[1,3]
  )
  colnames(dta) <- c("ind","chrom","hom","pos","all")
  init.population(PopName,dta)
  for (ind in 1:dta.description[1,3]) {
    set.population.info(PopName,ind,dta.info[ind])
  } 
}


"population.plabsim.save" <- function(file, PopName) {
  dta <- get.population(PopName)
  write(
    file = file,
    paste("#",PopName,population.size.get(PopName)$NoInds,sep=" ")
  )
  for (i in 1:population.size.get(PopName)$NoInds) {
    write(
      file   = file,
      paste("# ",i," \"",population.info.get(PopName,i),"\"",sep=""),
      append = TRUE
    )
  } 
  utils::write.table(
    file     = file,
    dta,
    quote    = FALSE,
    col.name = FALSE,
    row.name = FALSE,
    append   = TRUE
  )
}







"rng.choose" <-
  function(rngName="") {
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="rng_choose",
        arguments=list(as.character(rngName))), PACKAGE="SelectionTools");
  }

"rng.info" <-
  function() {
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="rng_info",
        arguments=list()), PACKAGE="SelectionTools")
  }

"rng.list" <-
  function() {
    c <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="rng_list",
        arguments=list()), PACKAGE="SelectionTools")
  }

"rng.init" <-
  function(seed = NULL) {
    if (is.null(seed))
      return(invisible(.Call("st_legacy_marker", list(action="legacy.result",
        routine="rng_init",
        arguments=list()), PACKAGE="SelectionTools")))
    if (!is.numeric(seed) || is.complex(seed) || length(seed) != 1L || is.na(seed) ||
        !is.finite(seed) || seed != trunc(seed) ||
        seed < -.Machine$integer.max || seed > .Machine$integer.max)
      stop("seed must be NULL or a single finite integer in R's seed range")
    set.seed(as.integer(seed))
    result <- .Call("st_legacy_sim", list(action="rng.init"),
                    PACKAGE="SelectionTools")
    if (result$retval != 0L)
      stop("Could not initialize internal random-number generators")
    invisible(result)
  }







"deprecated" <-  function(NameOld,NameNew) {}

"append.population" <-
  function(NameP1, NameP2){
    deprecated("append.population()","population.append()")
    population.append(NameP1, NameP2)
  }

"concat.population" <-
  function(NameP1, NameP2){
    deprecated("concat.population()","population.concat()")
    population.concat(NameP1, NameP2)
  }

"copy.population" <-
  function(NameP1, NameP2, start=1, n=-1) {
    deprecated("copy.population()","population.copy()")
    population.copy(NameP1, NameP2, start, n)
  }

"divide.population" <-
  function(NameP1, NameP2, NoI=-1) {
    deprecated("divide.population()","population.divide()")
    population.divide(NameP1, NameP2, NoI)
  }

"get.genome.par" <-
  function() {
    deprecated("get.genome.par()","genome.parameter.get()")
    genome.parameter.get()
  }

"get.map" <-
  function(){
    deprecated("get.map()","linkage.map.get()")
    linkage.map.get()
  }

"get.population.info" <-
  function(PopName,ind){
    deprecated("get.population.info()","population.info.get()")
    population.info.get(PopName,ind)
  }

"get.population.size" <-
  function(PopName){
    deprecated("get.population.size()","population.size.get()")
    population.size.get(PopName)
  }

"load.linkage.map" <-
  function(file, disperse=TRUE, file.disperse=NA, disperse.factor=100){
    deprecated("load.linkage.map()","linkage.map.load()")
    linkage.map.load(file, disperse, file.disperse, disperse.factor)
  }

"list.populations" <-
  function(){
    deprecated("list.populations()","population.list()")
    population.list()
  }

"optimize.population" <-
  function(PopNames){
    deprecated("optimize.population()","population.optimize()")
    population.optimize(PopNames)

  }

"population.swap.name" <-
  function(NameP1,NameP2){
    deprecated("population.swap.name()","population.name.swap()")
    population.name.swap(NameP1,NameP2)
  }

"remove.all.populations" <-
  function(){
    deprecated("remove.all.populations()","population.remove.all()")
    population.remove.all()
  }

"remove.population" <-
  function(PopNames) {
    deprecated("remove.population()","population.remove()")
    population.remove(PopNames)
  }

"rename.population" <-
  function(OldName, NewName){
    deprecated("rename.population()","population.rename()")
    population.rename(OldName, NewName)
  }

"resize.population" <-
  function(PopName, newSize){
    deprecated("rename.population()","population.rename()")
    population.resize(PopName, newSize)
  }

"sample.population" <-
  function(NameP1, NameP2, NoI=-1, rep=1){
    deprecated("sample.population()","population.sample()")
    population.sample(NameP1, NameP2, NoI, rep)
  }

"save.linkage.map" <- function(file, dta.map=get.map()) {
    deprecated("save.linkage.map()","linkage.map.save()")
    linkage.map.save(file, dta.map)
  }


"set.genome.par" <-
  function(no.chrom,no.hom,chrom.len) {
    deprecated("set.genome.par()","genome.parameter.set()")
    genome.parameter.set(no.chrom,no.hom,chrom.len)
  }

"set.population.info" <-
  function(PopName,ind,info){
    deprecated("set.population.info()","population.info.set()")
    population.info.set(PopName,ind,info)
  }









"swap.population.name" <-
  function(NameP1,NameP2) {
    deprecated("swap.population.name()","population.swap.name()")
    population.swap.name(NameP1,NameP2)
  }






genome.segments <- function(pops, allele, chromosome=NULL,
                            begin=NULL, end=NULL)
    .Call("st_legacy_sim", list(action="summary", kind="segments",
          pops=unlist(strsplit(pops, " ")), allele=allele,
          chromosome=chromosome, begin=begin, end=end), PACKAGE="SelectionTools")

select.n.best.segments <- function(newPop, oldPop, n=1, allele=1, decreasing=FALSE)
    invisible(.Call("st_legacy_sim", list(action="select.segments", target=newPop,
              source=oldPop, n=n, allele=allele, decreasing=decreasing),
              PACKAGE="SelectionTools"))

linkage.drag <- function(pops, allele, chrom, pos)
    .Call("st_legacy_sim", list(action="summary", kind="drag",
          pops=unlist(strsplit(pops, " ")), allele=allele,
          chromosome=chrom, pos=pos), PACKAGE="SelectionTools")

"evaluate.ld" <- function (PopName, allele, chrom, pos) 
.Call("st_legacy_sim", list(action = "output", operation = "ld", population = as.character(PopName), 
    allele = as.integer(allele), chromosome = as.integer(chrom), position = as.double(pos)), PACKAGE = "SelectionTools")









"gd.allow.zero.frequencies" <-
  function(allow=0){
    c  <- .Call("st_legacy_marker", list(action="legacy.result",
        routine="set_allow_zeros",
        arguments=list(as.integer(allow))), PACKAGE="SelectionTools")
  }

"gd.status.zero.frequencies" <- function()
  .Call("st_legacy_sim", list(action="allow.zeros"), PACKAGE="SelectionTools")

gd.allele.frequencies <- function(dta)
    .Call("st_legacy_sim", list(action="data", operation="frequencies", data=dta),
          PACKAGE="SelectionTools")

gd.correct.missing <- function(dta)
    .Call("st_legacy_sim", list(action="data", operation="correct", data=dta),
          PACKAGE="SelectionTools")

gd.list.missing <- function(dta)
    .Call("st_legacy_sim", list(action="data", operation="missing", data=dta),
          PACKAGE="SelectionTools")

gd.list.irregular <- function(dta1)
    .Call("st_legacy_sim", list(action="data", operation="irregular", data=dta1),
          PACKAGE="SelectionTools")

gd.distance.similarity <- function(dta, measure, par, boot=0,
        no.pop.u=0, pop.u=0, no.ind.u=0, ind.u=0, no.mar.u=0, mar.u=0)
    .Call("st_legacy_sim", list(action="data", operation="distance", data=dta,
          measure=measure, par=par, boot=boot, no.pop.u=no.pop.u, pop.u=pop.u,
          no.ind.u=no.ind.u, ind.u=ind.u, no.mar.u=no.mar.u, mar.u=mar.u),
          PACKAGE="SelectionTools")

"gd.similarity.coefficient" <-
  function(dta,measure="jac",par="dist",boot=0,
           no.pop.u=0,pop.u=0,no.ind.u=0,ind.u=0,no.mar.u=0,mar.u=0,
           resample.primers=FALSE){

    if ( (measure != "dic") && (measure != "jac") && (measure != "sma") ){
      info(-2,"Distance measure not defined")
      return (invisible(NULL))
    }
    p <- gd.data.parameters(dta)
    if ( sum(p$no.all)==p$no.mar){
      if ( gd.status.zero.frequencies()==0 ){
        info (-1,"Dominant marker data?")
        info (-1,"Consider 'gd.allow.zero.frequencies(1)'")
      }
    }
    if ( sum(p$no.all)>p$no.mar){
      if ( gd.status.zero.frequencies()==1 ){
        info (-1,"Co-dominant marker data?")
        info (-1,"Consider 'gd.allow.zero.frequencies(0)'")
      }
    }
    if ( sum(p$no.ind)>p$no.pop){
      info (-1,"Each population should consist of only one individual")
      info (-1,"when using similarity coefficients.")
    }
    if (resample.primers)
      .Call("st_legacy_sim", list(action="data", operation="distance", data=dta,
            measure=measure, par=par, boot=boot, no.pop.u=no.pop.u, pop.u=pop.u,
            no.ind.u=no.ind.u, ind.u=ind.u, no.mar.u=no.mar.u, mar.u=mar.u,
            primers=TRUE), PACKAGE="SelectionTools")
    else
      gd.distance.similarity(dta,measure,par,boot,
             no.pop.u,pop.u,no.ind.u,ind.u,no.mar.u,mar.u)
  }

"gd.genetic.distance" <-
  function(dta,measure="euc",par="dist",boot=0,
           no.pop.u=0,pop.u=0,no.ind.u=0,ind.u=0,no.mar.u=0,mar.u=0){

    if ( (measure != "euc") && (measure != "rd") && (measure != "mrd") ){
      info(-2,"Distance measure not defined")
      return (invisible(NULL))
    }

    p <- gd.data.parameters(dta)

    if ( sum(p$no.all)==p$no.mar){
      if ( gd.status.zero.frequencies()==0 ){
        info (-1,"Dominant marker data?")
        info (-1,"Consider 'gd.allow.zero.frequencies(1)'")
      }
      if (measure!="euc") {
        info (-1,"Consider the Euclidean distance for dominant data")
      }
      if ( sum(p$no.ind)>p$no.pop){
        info (-1,"analysis on population level for dominant data")
      }
    }

    if ( sum(p$no.all)>p$no.mar){
      if ( gd.status.zero.frequencies()==1 ){
        info (-1,"Co-dominant marker data?")
        info (-1,"Consider 'gd.allow.zero.frequencies(0)'")
      }
    }
    gd.distance.similarity(dta,measure,par,boot,
           no.pop.u,pop.u,no.ind.u,ind.u,no.mar.u,mar.u)
  }

gd.data.parameters <- function(dta)
    .Call("st_legacy_sim", list(action="data", operation="parameters", data=dta),
          PACKAGE="SelectionTools")

"gd.splitdot" <-
  function(m.a) {
    .Call("splitdot", as.character(m.a) )
  }

gd.mk.matr <- function(dta, distance)
    .Call("st_legacy_sim", list(action="data", operation="matrix", data=dta,
          distance=distance[,3]), PACKAGE="SelectionTools")

"gd.pcoa" <- function(dta, k=3)
    .Call("st_legacy_sim", list(action="pcoa", data=dta, k=k),
          PACKAGE="SelectionTools")





sdev <- function(x, na.rm=FALSE)
    .Call("st_legacy_sim", list(action="sdev", data=x, na.rm=na.rm),
          PACKAGE="SelectionTools")


population.pedigree <- function(name, source=c("auto", "simulation", "marker"))
    .Call("st_population_pedigree", as.character(name), match.arg(source), PACKAGE="SelectionTools")

population.cross <- function(name, pool1, pool2, n)
    .Call("st_population_cross", as.character(name), as.character(pool1), as.character(pool2),
          as.double(n), PACKAGE="SelectionTools")

population.marker.map <- function(n.markers)
    invisible(.Call("st_population_marker_map", n.markers, PACKAGE="SelectionTools"))


extend.marker.map <- function(n.markers)
    invisible(.Call("st_extend_marker_map", n.markers, PACKAGE="SelectionTools"))
