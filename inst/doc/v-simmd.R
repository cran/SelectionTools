## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(comment = '' )
library(SelectionTools)

## ----readin-------------------------------------------------------------------
map <- data.frame(chrom = 1,
                  pos = 0.0,
                  name = "form")
define.genome(map)

## ----popdef1------------------------------------------------------------------

init.population("P1", homozygote(1))
init.population("P2", homozygote(2))


## ----cross--------------------------------------------------------------------
cross("F1","P1", "P2", 563)

## ----eval1--------------------------------------------------------------------
evaluate.genotype2("P1","form")

## ----eval2--------------------------------------------------------------------
evaluate.genotype2("P2","form")

## ----eval3--------------------------------------------------------------------
evaluate.genotype("F1","form")

## ----echo=FALSE---------------------------------------------------------------
sample.population ("SampleFromF1","F1",252)
cross("F2","SampleFromF1","SampleFromF1",7323,self=1)


## -----------------------------------------------------------------------------
evaluate.genotype("F2","form")

## -----------------------------------------------------------------------------
reset.all()

map <- data.frame(chrom = 1:2, 
                  pos = c(0.0, 0.0), 
                  name = c("form","color")) 
define.genome(map)

init.population("P1",homozygote(1))
init.population("P2",homozygote(2))
cross("F1","P1","P2",10)
cross("F2","F1","F1",1000)

evaluate.genotype2("F2", map$name)

## -----------------------------------------------------------------------------
evaluate.genotype2("F2", map$name, mode=0)

## ----results='hide', message=FALSE, warning=FALSE-----------------------------
reset.all()

map <- data.frame(chrom = 1:2, 
                  pos = c(0.0, 0.0), 
                  name = c("form","color"))
define.genome(map)

init.population("P1",homozygote(1))
init.population("P2",homozygote(2))
cross ("F1","P1","P2",563)
cross ("BC1","F1","P2",1000)
evaluate.genotype2("BC1", map$name)

## ----echo=FALSE---------------------------------------------------------------
evaluate.genotype2("BC1", map$name)

## ----results='hide', message=FALSE, warning=FALSE-----------------------------
reset.all()

map <- data.frame(chrom = c(1,1), 
                  pos = c(0.0, 0.0), 
                  name = c("form","color"))
define.genome(map)

init.population("P1",homozygote(1))
init.population("P2",homozygote(2))
cross ("F1","P1","P2",563)
cross ("BC1","F1","P2",1000)
evaluate.genotype2("BC1", map$name)

## ----echo=FALSE---------------------------------------------------------------
evaluate.genotype2("BC1", map$name)

## -----------------------------------------------------------------------------
map <- data.frame(1,0.01,"form")
define.genome(map)
init.population("P1",homozygote(1))
init.population("P2",homozygote(2))
cross("F1","P1","P2",1)
cross("F2","F1","F1",1000)
evaluate.genotype2 ("F2","form")
evaluate.allele.freq2("F2","form")

## -----------------------------------------------------------------------------

cross("SYN1","F2","F2",1000)
cross("SYN2","SYN1","SYN1",1000)
cross("SYN3","SYN2","SYN2",1000)
evaluate.genotype2("SYN3","form")
evaluate.allele.freq2("SYN3","form")

## -----------------------------------------------------------------------------
remove.all.populations()
init.population("P11",homozygote(1,500))
init.population("P22",homozygote(2,500))
append.population("ADM","P11")
append.population("ADM","P22")
evaluate.genotype2 ("ADM","form")
evaluate.allele.freq2("ADM","form")

## -----------------------------------------------------------------------------
cross("SYN1","ADM","ADM",1000)
evaluate.genotype("SYN1","form")
evaluate.allele.freq("SYN1","form")

## -----------------------------------------------------------------------------
remove.all.populations()
init.population("P11",homozygote(1,100))
init.population("P22",homozygote(2,900))
append.population("ADM","P11")
append.population("ADM","P22")
evaluate.genotype2 ("ADM","form")
evaluate.allele.freq2("ADM","form")

## -----------------------------------------------------------------------------
cross("SYN1","ADM","ADM",1000)
evaluate.genotype("SYN1","form")
evaluate.allele.freq("SYN1","form")

## -----------------------------------------------------------------------------
cross("SYN1","ADM","ADM",1000)
cross("SYN2","SYN1","SYN1",1000)
cross("SYN3","SYN2","SYN2",1000)
evaluate.genotype("SYN3","form")
evaluate.allele.freq("SYN3","form")

## -----------------------------------------------------------------------------
reset.all()
map <- data.frame(1,0.0,"form")
define.genome(map)

init.population("P11",homozygote(1,900))
init.population("P22",homozygote(2,100))
append.population("ADM","P11")
append.population("ADM","P22")
cross("SYN1","ADM","ADM",1000)
evaluate.genotype2("SYN1","form")
evaluate.allele.freq2("SYN1","form")

## -----------------------------------------------------------------------------
sel.crit.form <- data.frame(
  locus    =  c("form", "form") ,
  allele1  = c(1,1),
  allele2  = c(1,2)     
)

sel.crit.form
select.genotypes("allrounds","SYN1",
                 sel.crit.form)

## -----------------------------------------------------------------------------
evaluate.genotype2("allrounds","form")

## -----------------------------------------------------------------------------
cross ("SYNA","allrounds","allrounds",1000)
evaluate.genotype2 ("SYNA","form")
evaluate.allele.freq2("SYNA","form")

## ----message=FALSE, warning=FALSE---------------------------------------------
reset.all()
map <- data.frame(1,0.0,"form")
define.genome(map)

for (i in 1:100){
  init.population("P11",homozygote(1,900))
  init.population("P22",homozygote(2,100))
  append.population("ADM","P11")
  append.population("ADM","P22")
  cross("SYN1","ADM","ADM",1000)
  select.genotypes("allrounds","SYN1",
                 sel.crit.form)  
  cross ("SYNA","allrounds","allrounds",1000)
  append.population("STORE","SYNA")
}


## -----------------------------------------------------------------------------
evaluate.genotype("STORE","form")
evaluate.allele.freq("STORE","form")

## -----------------------------------------------------------------------------
reset.all()
map <- data.frame(1,0.0,"SZA","trait1")
define.genome(map)
init.population("P11",homozygote(1,90))
init.population("P22",homozygote(2,10))
append.population("ADM","P11")
append.population("ADM","P22")
cross("SYN1","ADM","ADM",1000)
evaluate.genotype("SYN1","SZA")
evaluate.allele.freq("SYN1","SZA")

## -----------------------------------------------------------------------------
select.genotypes("R2","SYN1",data.frame("SZA",1,1))
evaluate.genotype("R2","SZA")
select.genotypes("R0","SYN1",data.frame("SZA",2,2))
evaluate.genotype("R0","SZA")
rename.population ("SYN1","R1")
evaluate.genotype("R0","SZA")

## -----------------------------------------------------------------------------
w2 <- 0.5 # Homozygous normal, allele 1
w1 <- 1.0 # Heterozygous
w0 <- 0.2 # Homozygous sickle cell anemia, allele 2

## ----eval=FALSE---------------------------------------------------------------
# ( nR2 <- n.inds$count[n.inds$PopName=="R2"])
# ( nR1 <- n.inds$count[n.inds$PopName=="R1"])
# ( nR0 <- n.inds$count[n.inds$PopName=="R0"])
# ( sR2 <- rbinom(n=1,size=nR2,prob=w2) )
# ( sR1 <- rbinom(n=1,size=nR1,prob=w1) )
# ( sR0 <- rbinom(n=1,size=nR0,prob=w0) )
# sample.population("S2","R2",sR2)
# sample.population("S1","R1",sR1)
# sample.population("S0","R0",sR0)
# remove.population("Sx")
# append.population("Sx","S2")
# append.population("Sx","S1")
# append.population("Sx","S0")
# evaluate.genotype2 ("Sx","SZA")
# evaluate.allele.freq2("Sx","SZA")

## ----eval=FALSE---------------------------------------------------------------
# cross ("SYNA","Sx","Sx",1000)
# evaluate.genotype2 ("SYNA","SZA")
# evaluate.allele.freq2("SYNA","SZA")

## -----------------------------------------------------------------------------
reset.all()
sel.crit.wrinkeled <- data.frame(
  locus    = "form" ,
  allele1  = 2       ,
  allele2  = 2     
)

sel.crit.yellow <- data.frame(
  locus    = c ("color" , "color"),
  allele1  = c (1       ,  1     ),
  allele2  = c( 1       ,  2     )
)

sel.crit.pink <- data.frame(
  locus    = "flower" ,
  allele1  = 1       ,
  allele2  = 2     
)
criteria <- rbind (sel.crit.yellow, sel.crit.wrinkeled, sel.crit.pink)

## -----------------------------------------------------------------------------
map <- data.frame(chrom = c(1,2,3),
                  pos = c(0.01,0.05,0.1),
                  name = c("form","color","flower"))
define.genome(map)

## -----------------------------------------------------------------------------
init.population("P1",homozygote(1))
init.population("P2",homozygote(2))
cross("F1","P1","P2",1)
cross("F2","F1","F1",1000)
cross("F3","F2","F2",1000)

## -----------------------------------------------------------------------------
select.genotypes("Sel","F3", criteria)

## -----------------------------------------------------------------------------
evaluate.genotype2("F3", map$name)
evaluate.genotype2("Sel", map$name)
criteria

