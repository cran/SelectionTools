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



evaluate.genotype2 <- function(PopName, Loci, alleles=0, mode=2, bg=0) {
  result <- .Call("genotype_evaluate",
                  as.character(PopName),
                  Loci,
                  as.integer(mode));
  if ( is.null(result) ) return( invisible(NULL) )
  if (bg == 1) {
    col.name <- colnames(result)
    res <- matrix(result[,1:(dim(result)[2]-2)],ncol=dim(result)[2]-2)
    result <- matrix(result[(0==rowSums(plabsim$allele.missing.indicator == res)),],ncol=dim(result)[2])
    result[,dim(result)[2]] <- result[,dim(result)[2]-1]/sum(result[,dim(result)[2]-1])
    colnames(result) <- col.name
  } else {
    result[,dim(result)[2]] <- result[,dim(result)[2]-1]/result[,dim(result)[2]]
  }
  if (alleles[1] != 0) {
    ng <- nrow(result)
    na <- length(alleles)
    passt <- 0; i <- 0
    while ( (i < ng) && (passt==0) ){
      passt <- 1 ; i <- i+1
      for (j in 1:na){
        if (result[i,j] != alleles[j]) {passt <- 0}
      }
      if (passt==1) {select <- i} 
    }
    if (passt==1) {result <- result[select,]}
  }  
  as.data.frame(result)
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
st.read.marker.data.df <- function(
    markerdata,
    format = NULL,
    data.set = "default"
) {
  
  # Hier chekcs einfuegen um format selbst zu bestimmen
  # zB IF df has three column FORMAT is LIST
  
  
  # Oder einfach alles in matrix verwandeln?
  
  # Dann C Funktion die Daten einliest und data.set erstellt
  C <- .Call("gs_read_marker_data_01_DF",
             markerdata,
             as.character(data.set))
  return(invisible(NULL))
}
