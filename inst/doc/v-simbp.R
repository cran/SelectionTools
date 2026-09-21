### R code from vignette source 'v-simbp.Rnw'

###################################################
### code chunk number 1: v-simbp.Rnw:207-212
###################################################
library ("SelectionTools")

vignette.data <- new.env(parent = emptyenv())
data("v-tropmaize-vcf", package="SelectionTools", envir=vignette.data)
data("v-tropmaize-phe", package="SelectionTools", envir=vignette.data)


###################################################
### code chunk number 2: v-simbp.Rnw:218-219
###################################################
st.STvcf.to.dataframe ( vignette.data$v.tropmaize.vcf ) [ 1:10, 1:12 ]


###################################################
### code chunk number 3: v-simbp.Rnw:222-223
###################################################
vignette.data$v.tropmaize.phe [1:10,]


###################################################
### code chunk number 4: v-simbp.Rnw:226-228
###################################################
st.load.vcf.data(vignette.data$v.tropmaize.vcf)
st.load.performance.data(vignette.data$v.tropmaize.phe)


###################################################
### code chunk number 5: v-simbp.Rnw:231-232
###################################################
x <- st.marker.data.statistics ()      # Overview of marker data


###################################################
### code chunk number 6: v-simbp.Rnw:234-235
###################################################
x$genotypes[1:10,1:5]


###################################################
### code chunk number 7: v-simbp.Rnw:238-239
###################################################
x$individual.list[1:5,]


###################################################
### code chunk number 8: v-simbp.Rnw:242-243
###################################################
x$marker.list[1:5,]


###################################################
### code chunk number 9: v-simbp.Rnw:267-271
###################################################
st.restrict.marker.data (NoAll.MAX=2  ) # Max. no. of alleles
st.restrict.marker.data (MaMis.MAX=0.1) # Max. missing at a marker
st.restrict.marker.data (ExHet.MIN=0.1) # Minimum gene diversity
st.restrict.marker.data (InMis.MAX=0.1) # Max. missing per individual


###################################################
### code chunk number 10: v-simbp.Rnw:275-276
###################################################
x <- st.marker.data.statistics ()      


###################################################
### code chunk number 11: v-simbp.Rnw:278-279
###################################################
x$marker.list[1:5,]


###################################################
### code chunk number 12: v-simbp.Rnw:287-288
###################################################
st.copy.marker.data ( "C0", source.data.set="default" )


###################################################
### code chunk number 13: v-simbp.Rnw:297-298
###################################################
gs.esteff.rr ( method="BLUP", data.set="C0") 


###################################################
### code chunk number 14: v-simbp.Rnw:303-304
###################################################
gs.esteff.rr ( method="BLUP", hsq=0.85, maxiter=0, data.set="C0") 


###################################################
### code chunk number 15: v-simbp.Rnw:307-308
###################################################
gs.plot.model.fit ( training.set="C0" )


###################################################
### code chunk number 16: v-simbp.Rnw:332-334 (eval = FALSE)
###################################################
## reset.all()               
## st.set.info.level(0)      


###################################################
### code chunk number 17: v-simbp.Rnw:337-338 (eval = FALSE)
###################################################
## st.set.simpop ( pop.name="C0", data.set="C0" )


###################################################
### code chunk number 18: v-simbp.Rnw:346-347 (eval = FALSE)
###################################################
## list.populations()                  # List all defined populations


###################################################
### code chunk number 19: v-simbp.Rnw:355-356 (eval = FALSE)
###################################################
## list.effects()                      # List all defined effects


###################################################
### code chunk number 20: v-simbp.Rnw:386-389 (eval = FALSE)
###################################################
## genotype.population("C0")         # Construct marker data matrix 
## evaluate.population("C0")         # Calculate genotypic values 
## g <- get.population.gvalue("C0")$gvalue          # Save in a variable 


###################################################
### code chunk number 21: v-simbp.Rnw:393-410
###################################################
plt.hist <- function (x,main="") {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mar=c(4,4,3,1))
  mu <- mean(x)
  sd <- sdev(x)
  h <- hist (x,xlim=c(-1,8),probability=TRUE,
             breaks=seq(-1,8,by=0.5),
             col="lightblue",
             main=main)
  ym <- max(h$density)
  lines(c(mu, mu), c(0, 1),col="red",lwd=2)
  lines(c(mu-2*sd, mu-2*sd), c(0, ym/2),col="blue",lwd=1)
  lines(c(mu+2*sd, mu+2*sd), c(0, ym/2),col="blue",lwd=1)
  text (-1,0.90*ym,paste("mu = ",format(mu,digits=2)),pos=4) 
  text (-1,0.85*ym,paste("sd = ",format(sd,digits=2)),pos=4) 
}


###################################################
### code chunk number 22: v-simbp.Rnw:415-416 (eval = FALSE)
###################################################
## plt.hist (g, main="Genotypic value, base population") 


###################################################
### code chunk number 23: v-simbp.Rnw:431-434 (eval = FALSE)
###################################################
## phenotype.population("C0", hsq=0.8)
## p <- get.population.pvalue("C0")$pvalue
## plt.hist (p,main="Phenotypic values, base population, hsq = 0.8") 


###################################################
### code chunk number 24: v-simbp.Rnw:445-448 (eval = FALSE)
###################################################
## population.sort("C0", decreasing=TRUE, selection.criterion="P")  
## population.divide("C0sel", "C0", 30)     
## list.populations()


###################################################
### code chunk number 25: v-simbp.Rnw:458-460 (eval = FALSE)
###################################################
## p.C0sel <- get.population.pvalue("C0sel")$pvalue
## plt.hist (p.C0sel, main="Phenotypic values, selected fraction")


###################################################
### code chunk number 26: v-simbp.Rnw:471-473 (eval = FALSE)
###################################################
## g.C0sel <- get.population.gvalue("C0sel")$gvalue   # Gen. values
## plt.hist (g.C0sel, main="Genotypic values, selected fraction")


###################################################
### code chunk number 27: v-simbp.Rnw:484-498
###################################################
sim01 <- function (h.sq, nsel, nrep) {
  R <- rep(0,nrep)
  for (i in 1:nrep) {
    st.set.simpop ( pop.name="C0", data.set="C0" ) 
    genotype.population("C0")
    evaluate.population("C0")
    phenotype.population("C0",h.sq)
    population.sort("C0", decreasing=TRUE, selection.criterion="P") 
    population.divide("C0sel", "C0", nsel)
    g.C0sel <- get.population.gvalue("C0sel")$gvalue
    R[i] <- mean(g.C0sel) - mean(g)
  }
  return( mean(R) )
}


###################################################
### code chunk number 28: v-simbp.Rnw:501-503 (eval = FALSE)
###################################################
## set.info.level(-1) # Only warnings and errors are printed, no messages
## sim01 ( h.sq=0.8, nsel=30, nrep=50 )


###################################################
### code chunk number 29: v-simbp.Rnw:513-527 (eval = FALSE)
###################################################
## nsel <- c(100, 50, 30, 10 )
## h.sq <- c(0.6, 0.7, 0.8, 0.9, 0.99)
## nrep <- 100
## 
## SG <- matrix(nrow=length(nsel),ncol=length(h.sq))
## rownames(SG) <- nsel
## colnames(SG) <- h.sq
## 
## for (i in 1:length(nsel))
##   for (j in 1:length(h.sq)) {
##     cat (sprintf("nsel = %i, h.sq = %4.2f     \r",nsel[i],h.sq[j]))
##     SG[i,j] <- sim01 (nsel=nsel[i], h.sq=h.sq[j],  nrep=nrep )
## }
## SG


###################################################
### code chunk number 30: v-simbp.Rnw:562-571 (eval = FALSE)
###################################################
## nsel <- 30
## h.sq <- 0.8
## 
## st.set.simpop ( pop.name="C0", data.set="C0" ) 
## genotype.population("C0")
## evaluate.population("C0")
## phenotype.population("C0",hsq=0.8)
## population.sort("C0", decreasing=TRUE, selection.criterion="P") 
## population.divide("C0sel", "C0", nsel)


###################################################
### code chunk number 31: v-simbp.Rnw:577-581 (eval = FALSE)
###################################################
## population.copy("t1","C0sel");
## population.copy("t2","C0")
## population.concat("t1","t2")     # Complete population, sorted
## st.get.simpop("t1","t1")         # Make a marker data set out of it    


###################################################
### code chunk number 32: v-simbp.Rnw:589-592 (eval = FALSE)
###################################################
## dist.mat <- st.genetic.distances ( measure="mrd",
##                                    format="m" ,  
##                                    data.set="t1")


###################################################
### code chunk number 33: v-simbp.Rnw:595-596 (eval = FALSE)
###################################################
## PC <- cmdscale(dist.mat,4)      # Principlal coordinate analysis


###################################################
### code chunk number 34: v-simbp.Rnw:599-611 (eval = FALSE)
###################################################
## color  <- c( rep ("red",30),       # selected
##              rep ("darkgrey",201)) # not selected
## symbol <- c( rep (19,30),          # selected
##              rep (22,201))         # not selected
## 
## oldpar <- par(no.readonly = TRUE)
## par(mfrow=c(2,2),mar=c(4,4,2,2))
## plot ( PC[,1], PC[,2], pch=symbol, col=color )
## plot ( PC[,1], PC[,3], pch=symbol, col=color )
## plot ( PC[,2], PC[,3], pch=symbol, col=color )
## plot ( PC[,1], PC[,4], pch=symbol, col=color )
## par(oldpar)


###################################################
### code chunk number 35: v-simbp.Rnw:629-631 (eval = FALSE)
###################################################
## set.info.level(0)
## st.copy.marker.data ( "t1h", source.data.set="default" )


###################################################
### code chunk number 36: v-simbp.Rnw:638-639 (eval = FALSE)
###################################################
## st.get.simpop ( pop.name="C0sel", data.set="t3h" )


###################################################
### code chunk number 37: v-simbp.Rnw:659-662 (eval = FALSE)
###################################################
## st.def.hblocks ( hap = 3,            # Combine Markers in a 3 cM
##                  hap.unit=2,         # window to a haplotype block 
##                  data.set="t1h" )


###################################################
### code chunk number 38: v-simbp.Rnw:675-676 (eval = FALSE)
###################################################
## st.recode.hil  (data.set="t1h")


###################################################
### code chunk number 39: v-simbp.Rnw:681-682 (eval = FALSE)
###################################################
## st.plot.gene.diversity (data.set="t1h") # Base population


###################################################
### code chunk number 40: v-simbp.Rnw:686-689 (eval = FALSE)
###################################################
## st.def.hblocks ( hap = 3, hap.unit=2, data.set="t3h" )
## st.recode.hil ( data.set="t3h" )
## st.plot.gene.diversity (data.set="t3h") # Selected fraction 


###################################################
### code chunk number 41: v-simbp.Rnw:696-697 (eval = FALSE)
###################################################
## st.plot.ggt ( data.set = "t3h" )


###################################################
### code chunk number 42: v-simbp.Rnw:721-730
###################################################
st.set.info.level(-2)  # Errors only
st.copy.marker.data ( "C0", source.data.set="default" )
gs.esteff.rr ( method="BLUP", hsq=0.85, maxiter=0, data.set="C0") 
st.set.simpop ( pop.name="C0", data.set="C0" ) 
genotype.population("C0")
evaluate.population("C0")
phenotype.population("C0",1)
population.sort("C0", decreasing=TRUE,selection.criterion = "P") 
population.divide("C0sel", "C0", 30)


###################################################
### code chunk number 43: v-simbp.Rnw:735-741
###################################################
population.copy("tmp","C0sel")
for (ii in 1:30) {
  nme <- sprintf("l%02i",ii)
  population.divide (nme,"tmp",1)
}
list.populations()[c(1:5,25:33),]


###################################################
### code chunk number 44: v-simbp.Rnw:750-757
###################################################
for (ii in seq(1,29,2)) {         
  nme1 <- sprintf("l%02i",ii)     
  nme2 <- sprintf("l%02i",ii+1)
  cross ("F1",nme1,nme2,1)    
  dh ("DH","F1",20)               
  population.append("C1dh","DH")
}


###################################################
### code chunk number 45: v-simbp.Rnw:760-764
###################################################
genotype.population("C1dh")
evaluate.population("C1dh")
g.C1dh <- get.population.gvalue("C1dh")$gvalue
plt.hist(g.C1dh,main="SSD-lines")


###################################################
### code chunk number 46: v-simbp.Rnw:772-779
###################################################
for (ii in seq(1,29,2)) {          
  nme1 <- sprintf("l%02i",ii)      
  nme2 <- sprintf("l%02i",ii+1)                                                    
  cross ("F1",nme1,nme2,1)                                                         
  ssd.mating ("SSD","F1",20)       
  population.append("C1ssd","SSD")
}


###################################################
### code chunk number 47: v-simbp.Rnw:782-786
###################################################
genotype.population("C1ssd")
evaluate.population("C1ssd")
g.C1ssd <- get.population.gvalue("C1ssd")$gvalue
plt.hist(g.C1ssd)


###################################################
### code chunk number 48: v-simbp.Rnw:802-806
###################################################
st.set.info.level(-2)  # Errors only
st.copy.marker.data ( "C0", source.data.set="default" )
gs.esteff.rr ( method="BLUP", hsq=0.85, maxiter=0, data.set="C0") 
st.set.simpop ( pop.name="C0", data.set="C0" ) 


###################################################
### code chunk number 49: v-simbp.Rnw:809-817
###################################################
hsq             <- 0.85
n.selected      <- 30     #  n.selected/2 crosses
lines.per.cross <- 20     #  15 crosses with 20 lines = 300 lines 

copy.population("base","C0")
genotype.population("base")
evaluate.population("base")
mean(get.population.gvalue("base")$gvalue)


###################################################
### code chunk number 50: v-simbp.Rnw:820-824
###################################################
population.sort("base", decreasing=TRUE)          
population.divide("selected", "base", n.selected)
remove.population("base")
mean(get.population.gvalue("selected")$gvalue)


###################################################
### code chunk number 51: v-simbp.Rnw:829-833
###################################################
for (ii in 1:n.selected) {             
  nme <- sprintf("l%02i",ii)            
  population.divide (nme,"selected",1)  
}


###################################################
### code chunk number 52: v-simbp.Rnw:837-844
###################################################
for (ii in seq(1,n.selected-1,2)) {          
  nme1 <- sprintf("l%02i",ii)                
  nme2 <- sprintf("l%02i",ii+1)              
  cross ("Fx",nme1,nme2,1)                   
  dh ("DH","Fx",lines.per.cross)             
  population.append("base","DH")     
}


###################################################
### code chunk number 53: v-simbp.Rnw:847-850
###################################################
genotype.population("base")              
evaluate.population("base")       
mean(get.population.gvalue("base")$gvalue)  


###################################################
### code chunk number 54: v-simbp.Rnw:857-865 (eval = FALSE)
###################################################
## st.set.simpop ( pop.name="C0", data.set="C0" ) 
## 
## n.selected      <- 30
## lines.per.cross <- 20
## cycles          <- 10
## replications    <-  5
## 
## perf <- matrix(0,nrow=cycles,ncol=replications)


###################################################
### code chunk number 55: v-simbp.Rnw:868-900 (eval = FALSE)
###################################################
## for (r in 1:replications)
##   {
##     copy.population("base","C0")
##     genotype.population("base")
##     evaluate.population("base")
##     phenotype.population("base",hsq)
##     #
##     for (c in 1:cycles)  
##       {
##         cat (sprintf("rep = %i, cycle = %i     \r",r,c))
##         population.sort("base", decreasing=TRUE, selection.criterion="P") 
##         population.divide("selected", "base", n.selected)
##         remove.population("base")
##         #
##         for (ii in 1:n.selected) {
##           nme <- sprintf("l%02i",ii)
##           population.divide (nme,"selected",1)
##         }
##         #
##         for (ii in seq(1,n.selected-1,2)){
##           nme1 <- sprintf("l%02i",ii)
##           nme2 <- sprintf("l%02i",ii+1)
##           cross ("Fx",nme1,nme2,1)    
##           dh ("DH","Fx",lines.per.cross)  
##           population.append("base","DH")
##         }
##         genotype.population("base")
##         evaluate.population("base")
##         perf[c,r] <- mean(get.population.gvalue("base")$gvalue)
##         phenotype.population("base",hsq)
##       }
##   }


###################################################
### code chunk number 56: v-simbp.Rnw:903-904 (eval = FALSE)
###################################################
## perf


###################################################
### code chunk number 57: v-simbp.Rnw:921-922 (eval = FALSE)
###################################################
## apply(perf,1,mean)


###################################################
### code chunk number 58: v-simbp.Rnw:937-997 (eval = FALSE)
###################################################
## lt.res <- function( n.selected      ,
##                     lines.per.cross ,
##                     cycles          ,
##                     replications    ,
##                     crossing.scheme  )
## {
##   perf <- matrix(0,nrow=cycles,ncol=replications)
##   for (r in 1:replications)
##   {
##     copy.population("base","C0")
##     genotype.population("base")
##     evaluate.population("base")
##     phenotype.population("base",hsq)
##     #
##     for (c in 1:cycles)  
##       {
##         cat (sprintf("\r rep = %i, cycle = %i\r",r,c))
##         population.sort("base", decreasing=TRUE, selection.criterion="P") 
##         population.divide("selected", "base", n.selected)
##         remove.population("base")
##         #
##         for (ii in 1:n.selected) {
##           nme <- sprintf("l%02i",ii)
##           population.divide (nme,"selected",1)
##         }
##         #
##         if ("simple-dh"==crossing.scheme)
##           {
##             for (ii in seq(1,n.selected-1,2)){
##               nme1 <- sprintf("l%02i",ii)
##               nme2 <- sprintf("l%02i",ii+1)
##               cross ("Fx",nme1,nme2,1)    
##               dh ("DH","Fx",lines.per.cross)  
##               population.append("base","DH")
##             }
##           }
##          else if ("simple-ssd"==crossing.scheme)
##            {
##              for (ii in seq(1,n.selected-1,2)){
##                nme1 <- sprintf("l%02i",ii)
##                nme2 <- sprintf("l%02i",ii+1)
##                cross ("Fx",nme1,nme2,1)    
##                ssd.mating ("SSD","Fx",lines.per.cross)
##                for (i in 1:1){
##                  copy.population("Parents","SSD")
##                  ssd.mating ("SSD","Parents",1)
##                }
##                population.append("base","SSD")
##              }
##            }
##         else {print("Error");return (NULL);}
##         genotype.population("base")
##         evaluate.population("base")
##         perf[c,r] <- mean(get.population.gvalue("base")$gvalue)
##         phenotype.population("base",hsq)
##       }
##   }
##   P <- apply(perf,1,mean)
##   return (P)
## }


###################################################
### code chunk number 59: v-simbp.Rnw:1000-1027 (eval = FALSE)
###################################################
## P1 <-
## lt.res ( n.selected      =  10,
##          lines.per.cross =  50,
##          cycles          =  20,
##          replications    = 100,
##          crossing.scheme = "simple-dh")
## 
## P2 <-
## lt.res ( n.selected      =  10,
##          lines.per.cross =  50,
##          cycles          =  20,
##          replications    = 100,
##          crossing.scheme = "simple-ssd")
## 
## P3<-
## lt.res ( n.selected      =  30,
##          lines.per.cross =  10,
##          cycles          =  20,
##          replications    = 100,
##          crossing.scheme = "simple-dh")
## 
## P4 <-
## lt.res ( n.selected      =  30,
##          lines.per.cross =  10,
##          cycles          =  20,
##          replications    = 100,
##          crossing.scheme = "simple-ssd")


###################################################
### code chunk number 60: v-simbp.Rnw:1030-1045 (eval = FALSE)
###################################################
## plot ( 1:length(P1), P1, type="n", ylim=c(1,10),
##       xlab="Cylcle of selection", ylab="yield",
##       main="Long term selection response in S2 vs DH schemes")
## 
## points ( 1:length(P1), P1, type="l",lty=1,lwd=2,col="blue" )
## points ( 1:length(P2), P2, type="l",lty=2,lwd=2,col="blue" )
## points ( 1:length(P3), P3, type="l",lty=1,lwd=2,col="red" )
## points ( 1:length(P4), P4, type="l",lty=2,lwd=2,col="red" )
## 
## legend(10,3,legend=c("DH, high selection intensity",
##                      "S2, high selection intensity",
##                      "DH, low selection intensity",
##                      "S2, low selection intensity"),
##        col=c("blue","blue","red","red"),
##        lty=c(1,2,1,2))


###################################################
### code chunk number 61: v-simbp.Rnw:1253-1253
###################################################



###################################################
### code chunk number 62: v-simbp.Rnw:1256-1256 (eval = FALSE)
###################################################
## 


