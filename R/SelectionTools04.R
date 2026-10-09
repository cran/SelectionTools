ph.linkage.map.create <-
  function(map, disperse=FALSE, file.disperse=NA, disperse.factor=100){
    allPops <- paste(as.vector((population.list()$PopName)),collapse = " ");
    remove.evaluate.population(allPops);  
    remove.effmaps();
    remove.genotype.population(allPops);
    remove.map();
    if (disperse) {
      left    <- 1
      right   <- 1;
      counter <- 0;
      if (dim(map)[1]>1) { # more than one marker locus
        for (line in 2:dim(map)[1]) {
          
          # Wir zaehlen bis entweder Marker oder Chromosom unterschiedlich sind
          # Unterschiedlich fuehrt zu TRUE und wir beginnen mit disperse fuer alle gezaehlten Marker
          if ((map[line-1,2] != map[line,2])||((map[line-1,1] != map[line,1]))) {
            
            # Haben wir mehr als einen Marker? Wenn ja, disperse
            if (left<right) {
              counter <- counter + right - left # Anzahl der zu aendernden Marker
              
              # Wenn Chromosom gleich ist, muss Position anders sein:
              if ((map[line-1,1] == map[line,1])) { #Marker Loci befinden sich nicht am Ende des Chromosoms
                
                # Position des ersten anderen Markers - Position Beginn des identischen Blocks
                delta <- (map[line,2]-map[left,2])/(disperse.factor*(right-left+1))
                # Erster Marker bleibt danach sequenz 1:x *delta fuer x identische Marker
                map[left:right,2] <- map[left:right,2] + (0:(right-left))*delta
              } else {
                # Ist Position des Markers gleich der maximalen Laenge des Chromosoms
                if (map[left,2] == get.genome.par()$chrom.len[map[line-1,1]]) { #vorher | marker loci sitzen am Ende des Chromosomes | hinter ihnen ist kein Platz
                  
                  # Es ist nicht der erste Marker und es existiert noch ein 
                  # Marker auf dem Chromosom vor Marker left
                  if ((left>1)&&(map[left,1]==map[left-1,1])) { 
                    # Wir duerfen unser Delta/Disperse nach links ausrichten
                    delta <- (map[left,2]-map[left-1,2])/(disperse.factor*(right-left+1))
                    map[left:right,2] <- map[left:right,2] - ((right-left):0)*delta
                  } else { # Marker ist erster Marker auf dem Chromosom
                    if (0 < map[left,2]) { #zwischen Marker left und chromosomanfang ist noch "platz"
                      delta <- (map[left,2])/(disperse.factor*(right-left+1))
                      map[left:right,2] <- map[left:right,2] - ((right-left):0)*delta
                    } else {
                      info(-1,paste("Dispersion not possible for marker loci:",paste(as.vector(map[left:right,3]), collapse =" "),"!"))
                    }
                  }
                } else {# marker left right sind am ende des Chromosoms, doch dahinter ist noch platz
                  delta <- (get.genome.par()$chrom.len[map[line-1,1]]-map[left,2])/(disperse.factor*(right-left+1))
                  map[left:right,2] <- map[left:right,2] + (0:(right-left))*delta
                }
              }
            }
            left <- line; # Wenn wir einmal disperse durchgefuehrt haben, wird aktuelle line left
          }
          right <- line; # Solange die Marker gleich sind, wird right nach rechts verschoben
        } #line
      }
      
      # Falls Loop auf identischem Block endet noch einmal disperse ausfuehren:
      if (left<right) {
        counter <- counter + right - left
        if (map[left,2] == get.genome.par()$chrom.len[map[line-1,1]]) {
          if ((left>1)&&(map[left,1]==map[left-1,1])) {
            delta <- (map[left,2]-map[left-1,2])/(disperse.factor*(right-left+1))
            map[left:right,2] <- map[left:right,2] - ((right-left):0)*delta
          } else {
            if (0 < map[left,2]) {
              delta <- (map[left,2])/(disperse.factor*(right-left+1))
              map[left:right,2] <- map[left:right,2] - ((right-left):0)*delta
            } else {
              info(-1,paste("Dispersion not possible for marker loci:",paste(as.vector(map[left:right,3]), collapse =" "),"!"))
            }
          }
        } else {
          delta <- (get.genome.par()$chrom.len[map[line-1,1]]-map[left,2])/(disperse.factor*(right-left+1))
          map[left:right,2] <- map[left:right,2] + (0:(right-left))*delta
        }
      }
      info(0,paste(counter,"loci positions modified"))
      # file.disperse <- ifelse(is.na(file.disperse),paste(file,".disperse",sep=""),file.disperse)
      # info(0,paste("Dispersed linkage map saved to file:",file.disperse))
      # write.table(map,file=file.disperse,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
      # file <- file.disperse
    }
    # Das c kann wahrscheinlich raus
    c  <- .Call("define_map_df", map)
    
    # Ob das geht?
    if(disperse) return(map)
  }

##### Define genome + map
define.genome <- function(map, length.vec=NULL){
  # Check format 
  if (!is.data.frame(map)) {info(-2,"Input is not a data.frame");return(invisible(NULL))}
  if(nrow(map)==0){
    info(-2,"Dataframe must contain at least one row")
    return(invisible(NULL))
  }
  
  # No NAs allowed
  if (any(is.na(map))) {info(-2,"Dataframe contains NA");return(invisible(NULL))}
  
  # Data type checks  
  if(!is.numeric(map[,1]))    {info(-2,"Chromosome must be numeric"); return(invisible(NULL))}
  if(any(map[,1] %% 1 != 0))  {info(-2,"Chromosome must be whole number"); return(invisible(NULL))} else {
    map[,1] <- as.integer(map[,1])
  }
  if(!is.numeric(map[,2]))    {info(-2,"Position must be numeric"); return(invisible(NULL))}
  if(any(map[,2] < 0))        {info(-2,"Position cannot be less than 0"); return(invisible(NULL))}
  if(!is.character(map[,3]))  {info(-2,"Name column (3) must be character"); return(invisible(NULL))}
  
  # IF class is given: IS class character?
  if(ncol(map)==4){
    if(!is.character(map[,4]))  {info(-2,"Class column (4) must be character"); return(invisible(NULL))}
  } else {
    map[,4] <- paste0("Class",map[,3])
  }
  
  
  # Genomeparameter vorbereiten
  if(!is.null(length.vec)){
    if(!is.numeric(length.vec))    {info(-2,"Vector provided for chromsome length must be numeric"); return(invisible(NULL))}
    no.chrom <- length(length.vec)
  } else {
    no.chrom <- length(unique(map[,1]))
  }
  no.hom <- 2
  
  ## Checks aus C nach R
  # Falls Leute nicht sequentielle indices eingeben: neu faktorisieren
  # Bsp: 1,4,5 wird zu 1,2,3
  if( !all((1:no.chrom) == sort(unique(map[,1])))){
    map[,1] <- match(map[,1], sort(unique(map[,1])))
    info(-1, "Refactored chromosome indices")
  }
  map <- map[order(map[,1], map[,2]),] #  order by Chromosome, Pos
  
  
  
  if (!is.null(length.vec)) {
    if (length(length.vec) != length(unique(map[, 1]))) {
      info(-2,
           "Vector provided for chromsome length must match number of chromosomes")
      return(invisible(NULL))
    }
    chrom.len <- length.vec
  } else {
    # what if Chom.length = 0
    chrom.len <- as.vector(tapply(map[, 2], map[, 1], max))
  }
  if (any(chrom.len == 0)) {
    chrom.len[chrom.len == 0] <- 0.01
    warning("Chromosome length cannot be 0; increased to 0.01", call. = FALSE)
  }
  
  genome.parameter.set(no.chrom, no.hom, chrom.len)
  
  
  ph.linkage.map.create(map)
}

evaluate.allele2 <- function(
    PopName,
    Loci,
    missing=0
) {
  # Das if statement funktioniert noch nicht wegen viel zu viel blabla
  # da muss man nochmal ran
  # if (1==missing) {
  #   evallo <- evaluate.locus(PopName,Loci)
  #   summe <- evallo[,1] + evallo[,2]
  #   locivec <- unlist(strsplit(Loci," "))
  #   Loci <- paste(locivec[which(summe>0)],collapse=" ")
  # }
  if (0 < length(Loci)) {
    bg <- plabsim$allele.missing.indicator
    result <- .Call("allele_evaluate",
                    as.character(PopName),
                    as.character(Loci),
                    as.integer(missing),
                    as.integer(bg));
    if ( is.null(result) ) return (invisible( NULL) )
    result <- as.data.frame(result)
    result <- result[,-(dim(result)[2])]
  } else {
    return(invisible( NULL))
  }
  result
}

evaluate.allele.freq2 <- function( PopName,
                                   Loci,
                                   missing=0 )
{
  # Quatsch, oder?
  # if ( 1 <length(Loci) ){
  #   info(-2,"Only available for single loci")
  #   return(invisible(NULL))
  # }
  e.a <- evaluate.allele2(PopName,Loci,missing)
  if (is.null(e.a)) return(invisible(NULL))
  e.a$freq <-e.a$count/sum(e.a$count)
  return(e.a)
}


define.effects.df <- function(LocNames) {
  for (lnames in LocNames) {
    eff.def <- data.frame(as.character(lnames), as.integer(1), as.numeric(1))
    
    
    n.eff <- sum(grepl(lnames, list.effects()))
    if (n.eff == 0) {
      eff.name <- paste0("EFF", lnames, 1)
    } else {
      warning(paste0("Effects already exist for ",lnames,", using preexisting effects."))
      return(invisible(NULL))
    }
    
    result <- .Call("define_effects_df", eff.def, eff.name)
  }
}

sel.target.define <- function(criteria.df){
  # H1 = Homozygot 1
  # H2 = Homozygot 2 
  if(nrow(criteria.df)==1){
    row.sum <- criteria.df[1,2] + criteria.df[1,3] 
    if( row.sum == 2){
      sel.crit <- "H1"
    } else if (row.sum == 3){
      sel.crit <- "MIX"
    } else if (row.sum == 4){
      sel.crit <- "H2"
    } else {
      info(-2,paste("Faulty selection criterion for ", criteria.df[1,1]))
      return (NULL)
    }
  }
  
  if(nrow(criteria.df)==2){
    row.sum1 <- criteria.df[1,2] + criteria.df[1,3] 
    row.sum2 <- criteria.df[2,2] + criteria.df[2,3] 
    tot.sum <- row.sum1 + row.sum2
    
    if(row.sum1 == row.sum2){
      sel.crit <- sel.target.define(criteria.df[1,])
    } else  if( tot.sum == 5){
      sel.crit <- c("H1","MIX")
    } else if (tot.sum == 6){
      sel.crit <- c("H2","H1")
    } else if (tot.sum == 7){
      sel.crit <- c("MIX","H1")
    } else {
      info(-2,paste("Faulty selection criterion for ", criteria.df[1,1]))
      return (NULL)
    }
  }
  
  if(nrow(criteria.df)==3){
    sel.crit <- "All"
  }
  if(nrow(criteria.df)>3){
    info(-2,paste("More than three possible selection criteria for ", criteria.df[1,1]))
    return (NULL)
  }
  
  return(sel.crit)
  
}



# Das neue select all best
select.genotypes <- function(newPop, oldPop, criteria) {
  # Criteria is a dataframe with: Locusname | Allel1 | Allel2
  # Wir brauchen: Effekte fuer alle unique(locusname)
  
  
  # We now use the old select best function
  # But with some hacks to get it right
  popsize <- population.size.get(oldPop)$NoInds
  if (is.null(popsize) || (0 == popsize)) {
    info(-2, paste("Population", oldPop, "is of size 0"))
    return (invisible(NULL))
  }
  
  names <- unique(criteria[, 1])
  for (n in names) {
    tmpPopName <- paste0("tmp",n)    
    define.effects.df(n)
    
    # For first trait
    # By definition, allele 1 gets value 1,
    # allele 2 gets value 0
    criterion.target <- sel.target.define(criteria[criteria[,1] == n, ])
    
    eff.name <- paste0("EFF", n, 1)
    genotype.population(oldPop)
    evaluate.population(oldPop, eff.name)
    
    # Alle sollen ausgewaehlt werden
    if (criterion.target[1] == "All") {
      population.sort(oldPop, decreasing = TRUE)
      c <- select.all.best.intern(oldPop, 3)
      population.divide(tmpPopName, oldPop, c$n)
      # Zwei Auswahlkriterien
    } else if (length(criterion.target) == 2) {
      # Kritierien enthalten heterozygoten G
      if ("MIX" %in% criterion.target) {
        crit2 <-  setdiff(criterion.target, "MIX")
        
        switch(crit2, H1 = {
          population.sort(oldPop, decreasing = TRUE)
          c <- select.all.best.intern(oldPop, 2)
          population.divide(tmpPopName, oldPop, c$n)
        }, H2 = {
          population.sort(oldPop, decreasing = FALSE)
          c <- select.all.best.intern(oldPop, 2)
          population.divide(tmpPopName, oldPop, c$n)
        })
        # Kriterien enthalten nur die homozygoten
      } else if (all(c("H1", "H2") %in% criterion.target)) {
        population.sort(oldPop, decreasing = TRUE)
        c <- select.all.best.intern(oldPop, 1)
        population.divide(tmpPopName, oldPop, c$n)
        
        population.sort(oldPop, decreasing = FALSE)
        c <- select.all.best.intern(oldPop, 1)
        population.divide("tmp", oldPop, c$n)
        population.concat(tmpPopName, "tmp")
      }
      # Wenn nur ein Kriterium gilt
    } else if (length(criterion.target) == 1) {
      switch(criterion.target[1],
             
             H1 = {
               population.sort(oldPop, decreasing = TRUE)
               c <- select.all.best.intern(oldPop, 1)
               population.divide(tmpPopName, oldPop, c$n)
             },
             MIX = {
               population.sort(oldPop, decreasing = TRUE)
               c <- select.all.best.intern(oldPop, 2)
               population.divide("tempPop", oldPop, c$n)
               
               population.sort("tempPop", decreasing = FALSE)
               c <- select.all.best.intern(oldPop, 1)
               population.divide(tmpPopName, "tempPop", c$n)
               population.concat(oldPop, "tempPop")
             },
             H2 = {
               population.sort(oldPop, decreasing = FALSE)
               c <- select.all.best.intern(oldPop, 1)
               population.divide(tmpPopName, oldPop, c$n)
             })
      
    }
    population.copy(newPop, tmpPopName)
    oldPop <- tmpPopName
  }
  population.remove(paste0("tmp",names))
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
  # Matrices containing marker names have character storage throughout.
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
