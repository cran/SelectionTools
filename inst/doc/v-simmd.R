### R code from vignette source 'v-simmd.Rnw'

###################################################
### code chunk number 1: setup
###################################################
library(SelectionTools)
options(warn = -1, digits = 5)
set.info.level(-2)
st.set.info.level(-2)


###################################################
### code chunk number 2: readin
###################################################
map <- data.frame(chrom = 1,
                  pos = 0.0,
                  name = "form")
define.genome(map)


###################################################
### code chunk number 3: readin
###################################################
set.genome.par ( no.chrom = 1,       # Only one chromosome
                 no.hom= 2,          # Diploid species
                 chrom.len = 0.01 )  # Chromosome length 1cM

define.map("form trait1 1/0")        # Chromosome1 / Position 0


###################################################
### code chunk number 4: popdef1
###################################################

init.population("P1", homozygote(1))
init.population("P2", homozygote(2))



###################################################
### code chunk number 5: cross
###################################################
cross("F1","P1", "P2", 563)


###################################################
### code chunk number 6: eval1
###################################################
evaluate.genotype("P1","form")


###################################################
### code chunk number 7: eval2
###################################################
evaluate.genotype("P2","form")


###################################################
### code chunk number 8: eval3
###################################################
evaluate.genotype("F1","form")


###################################################
### code chunk number 9: v-simmd.Rnw:277-280
###################################################
sample.population ("SampleFromF1","F1",252)
cross("F2","SampleFromF1","SampleFromF1",7323,self=1)



###################################################
### code chunk number 10: v-simmd.Rnw:283-284
###################################################
evaluate.genotype("F2","form")


###################################################
### code chunk number 11: v-simmd.Rnw:297-310
###################################################
reset.all()

map <- data.frame(chrom = 1:2,              # Two loci placed
                  pos = c(0.0, 0.0),        # at the beginning of
                  name = c("form","color")) # chromosome 1 and 2
define.genome(map)

init.population("P1",homozygote(1))
init.population("P2",homozygote(2))
cross("F1","P1","P2",10)
cross("F2","F1","F1",1000)

evaluate.genotype("F2", map$name)           # for all loci


###################################################
### code chunk number 12: v-simmd.Rnw:313-314
###################################################
evaluate.genotype("F2", map$name, mode=0)


###################################################
### code chunk number 13: v-simmd.Rnw:324-336
###################################################
reset.all()

map <- data.frame(chrom = 1:2, 
                  pos = c(0.0, 0.0), 
                  name = c("form","color"))
define.genome(map)

init.population("P1",homozygote(1))
init.population("P2",homozygote(2))
cross ("F1","P1","P2",563)             
cross ("BC1","F1","P2",1000)            # Backcrossing to parent 2
evaluate.genotype("BC1", map$name)


###################################################
### code chunk number 14: v-simmd.Rnw:342-343
###################################################
evaluate.genotype("BC1", map$name)


###################################################
### code chunk number 15: v-simmd.Rnw:357-369
###################################################
reset.all()

map <- data.frame(chrom = c(1,1), 
                  pos = c(0.0, 0.1), 
                  name = c("form","color"))
define.genome(map)

init.population("P1",homozygote(1))
init.population("P2",homozygote(2))
cross ("F1","P1","P2",563)
cross ("BC1","F1","P2",1000)
evaluate.genotype("BC1", map$name)


###################################################
### code chunk number 16: v-simmd.Rnw:372-373
###################################################
evaluate.genotype("BC1", map$name)


###################################################
### code chunk number 17: v-simmd.Rnw:386-394
###################################################
map <- data.frame(1,0.01,"form")
define.genome(map)
init.population("P1",homozygote(1))
init.population("P2",homozygote(2))
cross("F1","P1","P2",1)
cross("F2","F1","F1",1000)
evaluate.genotype ("F2","form")
evaluate.allele.freq("F2","form")


###################################################
### code chunk number 18: v-simmd.Rnw:400-406
###################################################

cross("SYN1","F2","F2",1000)
cross("SYN2","SYN1","SYN1",1000)
cross("SYN3","SYN2","SYN2",1000)
evaluate.genotype("SYN3","form")
evaluate.allele.freq("SYN3","form")


###################################################
### code chunk number 19: v-simmd.Rnw:416-423
###################################################
remove.all.populations()
init.population("P11",homozygote(1,500))
init.population("P22",homozygote(2,500))
append.population("ADM","P11")
append.population("ADM","P22")
evaluate.genotype ("ADM","form")
evaluate.allele.freq2("ADM","form")


###################################################
### code chunk number 20: v-simmd.Rnw:426-429
###################################################
cross("SYN1","ADM","ADM",1000)
evaluate.genotype("SYN1","form")
evaluate.allele.freq("SYN1","form")


###################################################
### code chunk number 21: v-simmd.Rnw:437-444
###################################################
remove.all.populations()
init.population("P11",homozygote(1,100))
init.population("P22",homozygote(2,900))
append.population("ADM","P11")
append.population("ADM","P22")
evaluate.genotype ("ADM","form")
evaluate.allele.freq2("ADM","form")


###################################################
### code chunk number 22: v-simmd.Rnw:447-450
###################################################
cross("SYN1","ADM","ADM",1000)
evaluate.genotype("SYN1","form")
evaluate.allele.freq("SYN1","form")


###################################################
### code chunk number 23: v-simmd.Rnw:453-458
###################################################
cross("SYN1","ADM","ADM",1000)
cross("SYN2","SYN1","SYN1",1000)
cross("SYN3","SYN2","SYN2",1000)
evaluate.genotype("SYN3","form")
evaluate.allele.freq("SYN3","form")


###################################################
### code chunk number 24: v-simmd.Rnw:470-481
###################################################
reset.all()
map <- data.frame(1,0.0,"form")
define.genome(map)

init.population("P11",homozygote(1,900))
init.population("P22",homozygote(2,100))
append.population("ADM","P11")
append.population("ADM","P22")
cross("SYN1","ADM","ADM",1000)
evaluate.genotype("SYN1","form")
evaluate.allele.freq2("SYN1","form")


###################################################
### code chunk number 25: v-simmd.Rnw:497-506
###################################################
sel.crit.form <- data.frame(
  locus    =  c("form", "form") ,
  allele1  = c(1,1),
  allele2  = c(1,2)     
)

sel.crit.form
select.genotypes("allrounds","SYN1",
                 sel.crit.form)


###################################################
### code chunk number 26: v-simmd.Rnw:511-512
###################################################
evaluate.genotype("allrounds","form")


###################################################
### code chunk number 27: v-simmd.Rnw:518-521
###################################################
cross ("SYNA","allrounds","allrounds",1000)
evaluate.genotype ("SYNA","form")
evaluate.allele.freq("SYNA","form")


###################################################
### code chunk number 28: v-simmd.Rnw:528-544
###################################################
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



###################################################
### code chunk number 29: v-simmd.Rnw:549-551
###################################################
evaluate.genotype("STORE","form")
evaluate.allele.freq("STORE","form")


###################################################
### code chunk number 30: v-simmd.Rnw:562-572
###################################################
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


###################################################
### code chunk number 31: v-simmd.Rnw:579-585
###################################################
select.genotypes("R2","SYN1",data.frame("SZA",1,1))
evaluate.genotype("R2","SZA")
select.genotypes("R0","SYN1",data.frame("SZA",2,2))
evaluate.genotype("R0","SZA")
rename.population ("SYN1","R1")
evaluate.genotype("R0","SZA")


###################################################
### code chunk number 32: v-simmd.Rnw:592-595
###################################################
w2 <- 0.5 # Homozygous normal, allele 1
w1 <- 1.0 # Heterozygous
w0 <- 0.2 # Homozygous sickle cell anemia, allele 2


###################################################
### code chunk number 33: v-simmd.Rnw:601-618
###################################################
n.inds <- population.list()

( nR2 <- n.inds$count[n.inds$PopName=="R2"])
( nR1 <- n.inds$count[n.inds$PopName=="R1"])
( nR0 <- n.inds$count[n.inds$PopName=="R0"])
( sR2 <- rbinom(n=1,size=nR2,prob=w2) )
( sR1 <- rbinom(n=1,size=nR1,prob=w1) )
( sR0 <- rbinom(n=1,size=nR0,prob=w0) )
sample.population("S2","R2",sR2)
sample.population("S1","R1",sR1)
sample.population("S0","R0",sR0)
remove.population("Sx")
append.population("Sx","S2")
append.population("Sx","S1")
append.population("Sx","S0")
evaluate.genotype ("Sx","SZA")
evaluate.allele.freq("Sx","SZA")


###################################################
### code chunk number 34: v-simmd.Rnw:623-626
###################################################
cross ("SYNA","Sx","Sx",1000)
evaluate.genotype ("SYNA","SZA")
evaluate.allele.freq2("SYNA","SZA")


###################################################
### code chunk number 35: v-simmd.Rnw:637-643
###################################################
reset.all()

map <- data.frame(chrom = c(1,2,3),
                  pos = c(0.01,0.05,0.1),
                  name = c("form","col","flow"))
define.genome(map)


###################################################
### code chunk number 36: v-simmd.Rnw:648-666
###################################################
sel.crit.wrinkeled <- data.frame(
  locus    = "form" ,
  allele1  = 2       ,
  allele2  = 2     
)

sel.crit.green <- data.frame(
  locus    = c ("col" , "col"),
  allele1  = c (1       ,  1     ),
  allele2  = c( 1       ,  2     )
)

sel.crit.pink <- data.frame(
  locus    = "flow" ,
  allele1  = 1       ,
  allele2  = 2     
)
criteria <- rbind (sel.crit.green, sel.crit.wrinkeled, sel.crit.pink)


###################################################
### code chunk number 37: v-simmd.Rnw:671-676
###################################################
init.population("P1",homozygote(1))
init.population("P2",homozygote(2))
cross("F1","P1","P2",1)
cross("F2","F1","F1",1000)
cross("F3","F2","F2",1000)


###################################################
### code chunk number 38: v-simmd.Rnw:681-682
###################################################
select.genotypes("Sel","F3", criteria)


###################################################
### code chunk number 39: v-simmd.Rnw:685-687
###################################################
evaluate.genotype("F3", map$name)
evaluate.genotype("Sel", map$name)


###################################################
### code chunk number 40: v-simmd.Rnw:693-694
###################################################
evaluate.allele.freq2("Sel", c("form", "col", "flow"))


