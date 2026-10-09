### R code from vignette source 'v-masbc.Rnw'

###################################################
### code chunk number 1: low_baseline_01 (eval = FALSE)
###################################################
## library(SelectionTools)
## options(width=72, digits=5)
## st.set.num.threads(1)
## 
## st.output.dir <- file.path(tempdir(), "SelectionTools-masbc")
## st.data.dir <- st.output.dir
## dir.create(st.output.dir, recursive=TRUE, showWarnings=FALSE)
## 
## NREP <- 10000
## 
## rng.init(seed=1001)
## reset.all()


###################################################
### code chunk number 2: low_baseline_02 (eval = FALSE)
###################################################
## data("lib00.map1")
## data("lib00.map1a")
## write.table(lib00.map1, file=file.path(st.output.dir,"mab-1.map"),
##             quote=FALSE, row.names=FALSE, col.names=FALSE)
## write.table(lib00.map1a, file=file.path(st.output.dir,"mab-1a.map"),
##             quote=FALSE, row.names=FALSE, col.names=FALSE)
## 
## set.genome.par(
##     no.chrom=10,
##     no.hom=2,
##     chrom.len=c(2.68,1.81,2.08,1.40,1.53,
##                 1.56,0.93,1.44,1.55,1.14)
## )
## load.linkage.map(file.path(st.output.dir,"mab-1.map"))


###################################################
### code chunk number 3: low_baseline_03 (eval = FALSE)
###################################################
## init.population("P1",homozygote(1))
## init.population("P2",homozygote(8))
## cross("F1","P1","P2",NREP)
## 
## gen <- c("P1","F1","BC1","BC2","BC3","BC4","BC5","BC6")
## for (i in 1:7)
##     cross(gen[i+1],gen[i],"P2",NREP)
## 
## rpg <- genome.contribution("BC1 BC2 BC3 BC4 BC5 BC6",8)
## 
## define.effects("markers80","0, mc 8 uniform 0.00625")
## estimated80 <- get.score("BC2","markers80")
## true80 <- evaluate.genome("BC2",8)$ratio


###################################################
### code chunk number 4: low_baseline_04 (eval = FALSE)
###################################################
## load.linkage.map(file.path(st.output.dir,"mab-1a.map"))
## define.effects("markers40","0, mc 8 uniform 0.0125")
## estimated40 <- get.score("BC2","markers40")
## true40 <- evaluate.genome("BC2",8)$ratio
## 
## out <- list(
##     recurrent.parent.genome=round(rpg,4),
##     marker.estimate.correlation=c(
##         markers80=cor(estimated80,true80),
##         markers40=cor(estimated40,true40)
##     )
## )
## print(out)


###################################################
### code chunk number 5: low_foreground_01 (eval = FALSE)
###################################################
## library(SelectionTools)
## options(width=72, digits=5)
## st.set.num.threads(1)
## 
## st.output.dir <- file.path(tempdir(), "SelectionTools-masbc")
## st.data.dir <- st.output.dir
## dir.create(st.output.dir, recursive=TRUE, showWarnings=FALSE)
## 
## NREP <- 10000
## 
## rng.init(seed=1002)
## reset.all()


###################################################
### code chunk number 6: low_foreground_02 (eval = FALSE)
###################################################
## data("lib00.map1")
## write.table(lib00.map1, file=file.path(st.output.dir,"mab-1.map"),
##             quote=FALSE, row.names=FALSE, col.names=FALSE)
## 
## set.genome.par(
##     no.chrom=10,
##     no.hom=2,
##     chrom.len=c(2.68,1.81,2.08,1.40,1.53,
##                 1.56,0.93,1.44,1.55,1.14)
## )
## load.linkage.map(file.path(st.output.dir,"mab-1.map"))
## define.effects("target","0, target 1 uniform 1")


###################################################
### code chunk number 7: low_foreground_03 (eval = FALSE)
###################################################
## init.population("P1",homozygote(1))
## init.population("P2",homozygote(8))
## cross("F1","P1","P2",1)
## 
## pop <- c("F1","BC1","BC2","BC3","BC4","BC5","BC6")
## se1 <- paste(pop,"s1",sep="")
## sel <- paste(pop,"se",sep="")
## store <- paste(pop,"st",sep="")
## sel[1] <- "F1"


###################################################
### code chunk number 8: low_foreground_04 (eval = FALSE)
###################################################
## for (j in 1:NREP) {
##     for (i in 1:6) {
##         cross(pop[i+1],sel[i],"P2",50)
##         select.all.best(se1[i+1],pop[i+1],"target")
##         sample.population(sel[i+1],se1[i+1],1)
##         append.population(store[i+1],sel[i+1])
##     }
## }


###################################################
### code chunk number 9: low_foreground_05 (eval = FALSE)
###################################################
## out <- round(
##     genome.contribution(
##         "BC1st BC2st BC3st BC4st BC5st BC6st",8
##     ),
##     4
## )
## print(out)


###################################################
### code chunk number 10: low_two_stage_01 (eval = FALSE)
###################################################
## library(SelectionTools)
## options(width=72, digits=5)
## st.set.num.threads(1)
## 
## st.output.dir <- file.path(tempdir(), "SelectionTools-masbc")
## st.data.dir <- st.output.dir
## dir.create(st.output.dir, recursive=TRUE, showWarnings=FALSE)
## 
## NREP <- 10000
## 
## rng.init(seed=1003)
## reset.all()


###################################################
### code chunk number 11: low_two_stage_02 (eval = FALSE)
###################################################
## data("lib00.map1")
## write.table(lib00.map1, file=file.path(st.output.dir,"mab-1.map"),
##             quote=FALSE, row.names=FALSE, col.names=FALSE)
## 
## set.genome.par(
##     no.chrom=10,
##     no.hom=2,
##     chrom.len=c(2.68,1.81,2.08,1.40,1.53,
##                 1.56,0.93,1.44,1.55,1.14)
## )
## load.linkage.map(file.path(st.output.dir,"mab-1.map"))
## define.effects("target","0, target 1 uniform 1")
## define.effects("markers","0, mc 8 uniform 0.00625")


###################################################
### code chunk number 12: low_two_stage_03 (eval = FALSE)
###################################################
## init.population("P1",homozygote(1))
## init.population("P2",homozygote(8))
## cross("F1","P1","P2",1)
## 
## pop <- c("F1","BC1","BC2","BC3","BC4","BC5","BC6",
##          "BC7","BC8","BC9","BC10")
## se1 <- paste(pop,"s1",sep="")
## se2 <- paste(pop,"s2",sep="")
## sel <- paste(pop,"se",sep="")
## store <- paste(pop,"st",sep="")
## sel[1] <- "F1"
## 
## mdp <- matrix(NA_real_,nrow=10,ncol=NREP)


###################################################
### code chunk number 13: low_two_stage_04 (eval = FALSE)
###################################################
## for (j in 1:NREP) {
##     reset.mdp()
##     for (i in 1:10) {
##         cross(pop[i+1],sel[i],"P2",100)
##         select.all.best(se1[i+1],pop[i+1],"target")
##         evaluate.mdp(se1[i+1],sel[i],"P2","markers")
##         mdp[i,j] <- get.mdp()
##         select.all.best(se2[i+1],se1[i+1],"markers")
##         sample.population(sel[i+1],se2[i+1],1)
##         append.population(store[i+1],sel[i+1])
##     }
## }


###################################################
### code chunk number 14: low_two_stage_05 (eval = FALSE)
###################################################
## out <- list(
##     recurrent.parent.genome=round(
##         genome.contribution(
##             paste0("BC",1:10,"st",collapse=" "),8
##         ),
##         4
##     ),
##     mean.marker.data.points=data.frame(
##         generation=paste0("BC",1:10),
##         MDP=round(rowMeans(mdp),1)
##     )
## )
## print(out)


###################################################
### code chunk number 15: low_three_stage_01 (eval = FALSE)
###################################################
## library(SelectionTools)
## options(width=72, digits=5)
## st.set.num.threads(1)
## 
## st.output.dir <- file.path(tempdir(), "SelectionTools-masbc")
## st.data.dir <- st.output.dir
## dir.create(st.output.dir, recursive=TRUE, showWarnings=FALSE)
## 
## NREP <- 10000
## 
## rng.init(seed=1004)
## reset.all()


###################################################
### code chunk number 16: low_three_stage_02 (eval = FALSE)
###################################################
## data("lib00.map2")
## write.table(lib00.map2, file=file.path(st.output.dir,"mab-2.map"),
##             quote=FALSE, row.names=FALSE, col.names=FALSE)
## 
## set.genome.par(
##     no.chrom=10,
##     no.hom=2,
##     chrom.len=c(2.68,1.81,2.08,1.40,1.53,
##                 1.56,0.93,1.44,1.55,1.14)
## )
## load.linkage.map(file.path(st.output.dir,"mab-2.map"))
## define.effects("target","0, target 1 uniform 1")
## define.effects("flanking","0, flanking 8 uniform 1")
## define.effects("markers","0, mc 8 uniform 1")


###################################################
### code chunk number 17: low_three_stage_03 (eval = FALSE)
###################################################
## init.population("P1",homozygote(1))
## init.population("P2",homozygote(8))
## cross("F1","P1","P2",1)
## 
## pop <- c("F1","BC1","BC2","BC3","BC4","BC5","BC6")
## s1 <- paste(pop,"s1",sep="")
## s2 <- paste(pop,"s2",sep="")
## sel <- paste(pop,"se",sep="")
## store <- paste(pop,"st",sep="")
## sel[1] <- "F1"
## 
## mdp <- matrix(NA_real_,nrow=6,ncol=NREP)


###################################################
### code chunk number 18: low_three_stage_04 (eval = FALSE)
###################################################
## for (j in 1:NREP) {
##     reset.mdp()
##     for (i in 1:6) {
##         cross(pop[i+1],sel[i],"P2",100)
##         select.all.best(s1[i+1],pop[i+1],"target")
##         evaluate.mdp(s1[i+1],sel[i],"P2","flanking")
##         select.all.best(s2[i+1],s1[i+1],"flanking")
##         evaluate.mdp(s2[i+1],sel[i],"P2","markers")
##         select.n.best(sel[i+1],s2[i+1],"markers",1)
##         append.population(store[i+1],sel[i+1])
##         mdp[i,j] <- get.mdp()
##     }
## }


###################################################
### code chunk number 19: low_three_stage_05 (eval = FALSE)
###################################################
## out <- list(
##     recurrent.parent.genome=round(
##         genome.contribution(
##             "BC1st BC2st BC3st BC4st BC5st BC6st",8
##         ),
##         4
##     ),
##     mean.marker.data.points=data.frame(
##         generation=paste0("BC",1:6),
##         MDP=round(rowMeans(mdp),1)
##     )
## )
## print(out)


###################################################
### code chunk number 20: input_overview_01 (eval = FALSE)
###################################################
## # Input files have been prepared by the shared setup.


###################################################
### code chunk number 21: input_overview_02 (eval = FALSE)
###################################################
## map <- read.table(
##     file.path(st.input.dir,"ex-mab1.map"),
##     stringsAsFactors=FALSE
## )
## names(map) <- c("chromosome","position.cM","locus")


###################################################
### code chunk number 22: input_overview_03 (eval = FALSE)
###################################################
## parent <- read.table(
##     file.path(st.input.dir,"ex-mabp1_m.mda"),
##     stringsAsFactors=FALSE,
##     skip=1
## )
## parent.name <- scan(
##     file.path(st.input.dir,"ex-mabp1_m.mda"),
##     what="character",
##     n=1,
##     quiet=TRUE
## )


###################################################
### code chunk number 23: input_overview_04 (eval = FALSE)
###################################################
## out <- list(
##     map=head(map,6),
##     parent.name=parent.name,
##     parent.markers=head(parent,6)
## )
## print(out)


###################################################
### code chunk number 24: high_core_01 (eval = FALSE)
###################################################
## rng.init(seed=2000)


###################################################
### code chunk number 25: high_core_02 (eval = FALSE)
###################################################
## mab.simulate(
##     simulation.name="core",
##     simulation.run="x",
##     repetitions=NREP,
##     linkage.map="ex-mab1.map",
##     target.loci="target",
##     flanking.loci=c("umc90","umc27"),
##     gen.type=c("f1","bc","bc","bc"),
##     population.size=c(1,100,100,100),
##     sel.strategy=c("n","tfb","tfb","tfb"),
##     result.file="core.RData"
## )


###################################################
### code chunk number 26: high_core_03 (eval = FALSE)
###################################################
## loaded <- mab.load.data(
##     simulation.name="core",
##     simulation.run="x",
##     result.file="core.RData"
## )


###################################################
### code chunk number 27: high_core_04 (eval = FALSE)
###################################################
## tab <- mab.tabulate(
##     simulation.names.runs=c("core","x"),
##     ext.summary=TRUE,
##     TargetAlleles=TRUE,
##     RPG.gw.mean=TRUE,
##     RPG.gw.Q10=TRUE,
##     DonorSegments.gw=TRUE,
##     LinkageDrag=TRUE,
##     MarkerAnalyses=TRUE,
##     result.files="core.RData"
## )


###################################################
### code chunk number 28: high_core_05 (eval = FALSE)
###################################################
## out <- list(
##     result.components=names(loaded),
##     selected.summary=tab
## )
## print(out)


###################################################
### code chunk number 29: ex01_01 (eval = FALSE)
###################################################
## rng.init(seed=2101)


###################################################
### code chunk number 30: ex01_02 (eval = FALSE)
###################################################
## mab.simulate(
##     simulation.name="fbm99-01",
##     simulation.run="x",
##     repetitions=NREP,
##     linkage.map="ex-mab1.map",
##     target.loci="target",
##     flanking.loci=c("umc90","umc27"),
##     gen.type=c("f1","bc","bc","bc","bc","bc","bc","s","s"),
##     population.size=c(1,100,100,100,100,100,100,100,100),
##     sel.strategy=c("n","t","t","t","t","t","t","t","t"),
##     result.file="fbm99-01.RData"
## )


###################################################
### code chunk number 31: ex01_03 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=c("fbm99-01","x"),
##     TargetAlleles=TRUE,
##     RPG.gw.mean=TRUE,
##     LinkageDrag=TRUE,
##     result.files="fbm99-01.RData"
## )
## print(out)


###################################################
### code chunk number 32: ex02_01 (eval = FALSE)
###################################################
## rng.init(seed=2102)


###################################################
### code chunk number 33: ex02_02 (eval = FALSE)
###################################################
## mab.simulate(
##     simulation.name="fbm99-01",
##     simulation.run="x",
##     repetitions=NREP,
##     linkage.map="ex-mab1.map",
##     target.loci="target",
##     flanking.loci=c("umc90","umc27"),
##     gen.type=c("f1","bc","bc","bc","bc","bc","bc","s","s"),
##     population.size=c(1,100,100,100,100,100,100,100,100),
##     sel.strategy=c("n","t","t","t","t","t","t","t","t"),
##     result.file="fbm99-01.RData"
## )


###################################################
### code chunk number 34: ex02_03 (eval = FALSE)
###################################################
## mab.simulate(
##     simulation.name="fbm99-02",
##     simulation.run="x",
##     repetitions=NREP,
##     linkage.map="ex-mab1.map",
##     target.loci="target",
##     flanking.loci=c("umc90","umc27"),
##     gen.type=c("f1","bc","bc","bc","bc","bc"),
##     population.size=c(1,100,100,100,100,100),
##     sel.strategy=c("n","tb","tb","tb","tb","tb"),
##     result.file="fbm99-02.RData"
## )


###################################################
### code chunk number 35: ex02_04 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=c(
##         "fbm99-01","x",
##         "fbm99-02","x"
##     ),
##     TargetAlleles=TRUE,
##     RPG.gw.mean=TRUE,
##     RPG.gw.Q10=TRUE,
##     LinkageDrag=TRUE,
##     result.files=c("fbm99-01.RData","fbm99-02.RData")
## )
## print(out)


###################################################
### code chunk number 36: ex03_01 (eval = FALSE)
###################################################
## rng.init(seed=2103)


###################################################
### code chunk number 37: ex03_02 (eval = FALSE)
###################################################
## mab.simulate(
##     simulation.name="fbm99-03",
##     simulation.run="x",
##     repetitions=NREP,
##     linkage.map="ex-mab1.map",
##     target.loci="target",
##     flanking.loci=c("umc90","umc27"),
##     gen.type=c("f1","bc","bc","bc"),
##     population.size=c(1,33,100,166),
##     sel.strategy=c("n","tfb","tfb","tfb"),
##     result.file="fbm99-03.RData"
## )


###################################################
### code chunk number 38: ex03_03 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=c("fbm99-03","x"),
##     RPG.gw.Q10=TRUE,
##     LinkageDrag=TRUE,
##     MarkerAnalyses=TRUE,
##     result.files="fbm99-03.RData"
## )
## print(out)


###################################################
### code chunk number 39: ex04_01 (eval = FALSE)
###################################################
## rng.init(seed=2104)


###################################################
### code chunk number 40: ex04_02 (eval = FALSE)
###################################################
## mab.compare(
##     a.simulation.name="s01-060",
##     a.simulation.run="x",
##     a.repetitions=NREP,
##     a.linkage.map="ex-mab1.map",
##     a.target.loci="target",
##     a.flanking.loci=c("umc90","umc27"),
##     a.gen.type=c("f1","bc","bc","bc"),
##     a.population.size=c(1,60,60,60),
##     a.sel.strategy=c("n","tb","tb","tb"),
##     b.simulation.name="s01-080",
##     b.population.size=c(1,80,80,80),
##     c.simulation.name="s01-100",
##     c.population.size=c(1,100,100,100),
##     d.simulation.name="s01-125",
##     d.population.size=c(1,125,125,125),
##     e.simulation.name="s01-150",
##     e.population.size=c(1,150,150,150),
##     f.simulation.name="s01-200",
##     f.population.size=c(1,200,200,200),
##     result.files=c(
##         "s01-060.RData",
##         "s01-080.RData",
##         "s01-100.RData",
##         "s01-125.RData",
##         "s01-150.RData",
##         "s01-200.RData"
##     )
## )


###################################################
### code chunk number 41: ex04_03 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=c(
##         "s01-060","x",
##         "s01-080","x",
##         "s01-100","x",
##         "s01-125","x",
##         "s01-150","x",
##         "s01-200","x"
##     ),
##     RPG.gw.Q10=TRUE,
##     LinkageDrag=TRUE,
##     MarkerAnalyses=TRUE,
##     result.files=c(
##         "s01-060.RData",
##         "s01-080.RData",
##         "s01-100.RData",
##         "s01-125.RData",
##         "s01-150.RData",
##         "s01-200.RData"
##     )
## )
## print(out)


###################################################
### code chunk number 42: ex05_01 (eval = FALSE)
###################################################
## rng.init(seed=2105)


###################################################
### code chunk number 43: ex05_02 (eval = FALSE)
###################################################
## mab.simulate(
##     simulation.name="fm01u-01",
##     simulation.run="x",
##     repetitions=NREP,
##     linkage.map="ex-mab2.map",
##     target.loci=c("target1","target2"),
##     flanking.loci=c("umc90","umc27","umc111","bnl1407"),
##     gen.type=c("f1","bc","bc","bc"),
##     population.size=c(1,60,120,180),
##     sel.strategy=c("n","tfb","tfb","tfb"),
##     result.file="fm01u-01.RData"
## )


###################################################
### code chunk number 44: ex05_03 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=c("fm01u-01","x"),
##     TargetAlleles=TRUE,
##     RPG.gw.mean=TRUE,
##     LinkageDrag=TRUE,
##     MarkerAnalyses=TRUE,
##     result.files="fm01u-01.RData"
## )
## print(out)


###################################################
### code chunk number 45: ex06_01 (eval = FALSE)
###################################################
## rng.init(seed=2106)


###################################################
### code chunk number 46: ex06_02 (eval = FALSE)
###################################################
## mab.compare(
##     a.simulation.name="ld-t",
##     a.simulation.run="x",
##     a.repetitions=NREP,
##     a.linkage.map="ex-mab1.map",
##     a.target.loci="target",
##     a.flanking.loci=c("umc90","umc27"),
##     a.gen.type=c("f1","bc","bc","bc"),
##     a.population.size=c(1,100,100,100),
##     a.sel.strategy=c("n","t","t","t"),
##     a.no.selected=c(1,1,1,1),
##     b.simulation.name="ld-tb",
##     b.sel.strategy=c("n","tb","tb","tb"),
##     c.simulation.name="ld-tfb",
##     c.sel.strategy=c("n","tfb","tfb","tfb"),
##     result.files=c("ld-t.RData","ld-tb.RData","ld-tfb.RData")
## )


###################################################
### code chunk number 47: ex06_03 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=c(
##         "ld-t","x",
##         "ld-tb","x",
##         "ld-tfb","x"
##     ),
##     RPG.gw.Q10=TRUE,
##     LinkageDrag=TRUE,
##     MarkerAnalyses=TRUE,
##     result.files=c("ld-t.RData","ld-tb.RData","ld-tfb.RData")
## )
## print(out)


###################################################
### code chunk number 48: ex07_01 (eval = FALSE)
###################################################
## rng.init(seed=2107)


###################################################
### code chunk number 49: ex07_02 (eval = FALSE)
###################################################
## mab.compare(
##     a.simulation.name="fm013g-01",
##     a.simulation.run="x",
##     a.repetitions=NREP,
##     a.linkage.map="ex-mab3.map",
##     a.target.loci=c("target1","target2","target3"),
##     a.flanking.loci=c(
##         "umc90","umc27",
##         "umc111","bnl1407",
##         "umc37","bnl1518"
##     ),
##     a.gen.type=c("f1","bc","bc","bc"),
##     a.population.size=c(1,100,100,100),
##     a.sel.strategy=c("n","t","t","t"),
##     b.simulation.name="fm013g-02",
##     b.sel.strategy=c("n","tb","tb","tb"),
##     c.simulation.name="fm013g-03",
##     c.sel.strategy=c("n","tfb","tfb","tfb"),
##     d.simulation.name="fm013g-04",
##     d.population.size=c(1,200,100,100),
##     d.sel.strategy=c("n","tfb","tb","tb"),
##     result.files=c(
##         "fm013g-01.RData",
##         "fm013g-02.RData",
##         "fm013g-03.RData",
##         "fm013g-04.RData"
##     )
## )


###################################################
### code chunk number 50: ex07_03 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=c(
##         "fm013g-01","x",
##         "fm013g-02","x",
##         "fm013g-03","x",
##         "fm013g-04","x"
##     ),
##     RPG.gw.Q10=TRUE,
##     TargetAlleles=TRUE,
##     LinkageDrag=TRUE,
##     MarkerAnalyses=TRUE,
##     result.files=c(
##         "fm013g-01.RData",
##         "fm013g-02.RData",
##         "fm013g-03.RData",
##         "fm013g-04.RData"
##     )
## )
## print(out)


###################################################
### code chunk number 51: ex08_01 (eval = FALSE)
###################################################
## rng.init(seed=2108)


###################################################
### code chunk number 52: ex08_02 (eval = FALSE)
###################################################
## mab.compare(
##     a.simulation.name="marker-01",
##     a.simulation.run="x",
##     a.repetitions=NREP,
##     a.donor.parent="ex-mabp1_m.mda",
##     a.recurrent.parent="ex-mabp2_m.mda",
##     a.linkage.map="ex-mab4.map",
##     a.target.loci="s7e2190xxx",
##     a.flanking.loci=c("s7e5163xxx","s7e4815xxx"),
##     a.gen.type=c("f1","bc","bc"),
##     a.population.size=c(1,25,25),
##     a.sel.strategy=c("n","t","t"),
##     a.missing.allele="9",
##     b.simulation.name="marker-02",
##     b.population.size=c(1,50,50),
##     b.sel.strategy=c("n","tb","tb"),
##     result.files=c("marker-01.RData","marker-02.RData")
## )


###################################################
### code chunk number 53: ex08_03 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=c(
##         "marker-01","x",
##         "marker-02","x"
##     ),
##     RPG.gw.Q10=TRUE,
##     DonorSegments.gw=TRUE,
##     LinkageDrag=TRUE,
##     result.files=c("marker-01.RData","marker-02.RData")
## )
## print(out)


###################################################
### code chunk number 54: ex09_01 (eval = FALSE)
###################################################
## rng.init(seed=2109)


###################################################
### code chunk number 55: ex09_02 (eval = FALSE)
###################################################
## mab.compare(
##     a.simulation.name="bc1i1tfb",
##     a.simulation.run="x",
##     a.repetitions=NREP,
##     a.donor.parent="ex-mabi1.mda",
##     a.recurrent.parent="ex-mabp2.mda",
##     a.linkage.map="ex-mab1.map",
##     a.target.loci="target",
##     a.flanking.loci=c("umc90","umc27"),
##     a.gen.type=c("bc","bc"),
##     a.population.size=c(100,100),
##     a.sel.strategy=c("tfb","tfb"),
##     b.simulation.name="bc1i2tfb",
##     b.donor.parent="ex-mabi2.mda",
##     c.simulation.name="bc1i3tfb",
##     c.donor.parent="ex-mabi3.mda",
##     d.simulation.name="bc1i4tfb",
##     d.donor.parent="ex-mabi4.mda",
##     e.simulation.name="bc1i5tfb",
##     e.donor.parent="ex-mabi5.mda",
##     f.simulation.name="bc1i6tfb",
##     f.donor.parent="ex-mabi6.mda",
##     result.files=paste0("bc1i",1:6,"tfb.RData")
## )


###################################################
### code chunk number 56: ex09_03 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=as.vector(rbind(
##         paste0("bc1i",1:6,"tfb"),
##         rep("x",6)
##     )),
##     RPG.gw.Q10=TRUE,
##     DonorSegments.gw=TRUE,
##     LinkageDrag=TRUE,
##     result.files=paste0("bc1i",1:6,"tfb.RData")
## )
## print(out)


###################################################
### code chunk number 57: ex09a_01 (eval = FALSE)
###################################################
## rng.init(seed=2110)


###################################################
### code chunk number 58: ex09a_02 (eval = FALSE)
###################################################
## mab.compare(
##     a.simulation.name="bc1i1",
##     a.simulation.run="x",
##     a.donor.parent="ex-mabi1.mda",
##     a.recurrent.parent="ex-mabp2.mda",
##     a.linkage.map="ex-mab1.map",
##     a.target.loci="target",
##     a.flanking.loci=c("umc90","umc27"),
##     b.simulation.name="bc1i2",
##     b.donor.parent="ex-mabi2.mda",
##     c.simulation.name="bc1i3",
##     c.donor.parent="ex-mabi3.mda",
##     d.simulation.name="bc1i4",
##     d.donor.parent="ex-mabi4.mda",
##     e.simulation.name="bc1i5",
##     e.donor.parent="ex-mabi5.mda",
##     f.simulation.name="bc1i6",
##     f.donor.parent="ex-mabi6.mda",
##     result.files=paste0("bc1i",1:6,".RData")
## )


###################################################
### code chunk number 59: ex09a_03 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=as.vector(rbind(
##         paste0("bc1i",1:6),
##         rep("x",6)
##     )),
##     RPG.gw.Q10=TRUE,
##     DonorSegments.gw=TRUE,
##     LinkageDrag=TRUE,
##     result.files=paste0("bc1i",1:6,".RData")
## )
## print(out)


###################################################
### code chunk number 60: ex10_01 (eval = FALSE)
###################################################
## rng.init(seed=2111)


###################################################
### code chunk number 61: ex10_02 (eval = FALSE)
###################################################
## mab.compare(
##     a.simulation.name="bc1i1tfs",
##     a.simulation.run="x",
##     a.repetitions=NREP,
##     a.donor.parent="ex-mabi1.mda",
##     a.recurrent.parent="ex-mabp2.mda",
##     a.linkage.map="ex-mab1.map",
##     a.target.loci="target",
##     a.flanking.loci=c("umc90","umc27"),
##     a.gen.type=c("bc","bc"),
##     a.population.size=c(100,100),
##     a.sel.strategy=c("tfs","tfs"),
##     a.no.selected=c(1,1),
##     b.simulation.name="bc1i2tfs",
##     b.donor.parent="ex-mabi2.mda",
##     c.simulation.name="bc1i3tfs",
##     c.donor.parent="ex-mabi3.mda",
##     d.simulation.name="bc1i4tfs",
##     d.donor.parent="ex-mabi4.mda",
##     e.simulation.name="bc1i5tfs",
##     e.donor.parent="ex-mabi5.mda",
##     f.simulation.name="bc1i6tfs",
##     f.donor.parent="ex-mabi6.mda",
##     result.files=paste0("bc1i",1:6,"tfs.RData")
## )


###################################################
### code chunk number 62: ex10_03 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=as.vector(rbind(
##         paste0("bc1i",1:6,"tfs"),
##         rep("x",6)
##     )),
##     RPG.gw.Q10=TRUE,
##     DonorSegments.gw=TRUE,
##     LinkageDrag=TRUE,
##     result.files=paste0("bc1i",1:6,"tfs.RData")
## )
## print(out)


###################################################
### code chunk number 63: ex11_01 (eval = FALSE)
###################################################
## rng.init(seed=2112)


###################################################
### code chunk number 64: ex11_02 (eval = FALSE)
###################################################
## mab.simulate(
##     simulation.name="regions-1",
##     simulation.run="x",
##     repetitions=NREP,
##     linkage.map="ex-mab1.map",
##     target.loci="target",
##     flanking.loci=c("umc90","umc27"),
##     gen.type=c("f1","bc","bc","bc"),
##     population.size=c(1,500,500,500),
##     sel.strategy=c("n","t","t","t"),
##     no.selected=c(1,1,1,1),
##     reg.chr=c(1,5,5),
##     reg.begin=c(0,0,0),
##     reg.end=c(0,0,60),
##     result.file="regions-1.RData"
## )


###################################################
### code chunk number 65: ex11_03 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=c("regions-1","x"),
##     ext.summary=TRUE,
##     RPG.gw.mean=TRUE,
##     RPG.reg.mean=TRUE,
##     DonorSegments.gw=TRUE,
##     DonorSegments.reg=TRUE,
##     result.files="regions-1.RData"
## )
## print(out)


###################################################
### code chunk number 66: ex12_01 (eval = FALSE)
###################################################
## recipient <- c(
##     "bnl845","umc53","umc6","umc34","umc131",
##     "umc98a","umc4","umc10","bnl844b","umc36"
## )


###################################################
### code chunk number 67: ex12_02 (eval = FALSE)
###################################################
## run.recipient <- function(name,strategy,seed) {
##     rng.init(seed=seed)
##     mab.simulate(
##         simulation.name=name,
##         simulation.run="x",
##         repetitions=NREP,
##         linkage.map="ex-mab1.map",
##         target.loci="target",
##         flanking.loci=c("umc90","umc27"),
##         recipient.loci=recipient,
##         gen.type=c("f1","bc","bc","bc"),
##         population.size=c(1,1000,1000,1000),
##         sel.strategy=c("n",strategy,strategy,strategy),
##         no.selected=c(1,1,1,1),
##         reg.chr=c(1,2,5,5),
##         reg.begin=c(0,0,0,0),
##         reg.end=c(0,0,0,60),
##         result.file=paste0(name,".RData")
##     )
## }


###################################################
### code chunk number 68: ex12_03 (eval = FALSE)
###################################################
## run.recipient("ss-t",   "t",   21131)
## run.recipient("ss-tf",  "tf",  21132)
## run.recipient("ss-tr",  "tr",  21133)
## run.recipient("ss-tfb", "tfb", 21134)
## run.recipient("ss-trb", "trb", 21135)
## run.recipient("ss-trfb","trfb",21136)
## run.recipient("ss-tfrb","tfrb",21137)


###################################################
### code chunk number 69: ex12_04 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=c(
##         "ss-t","x",
##         "ss-tf","x",
##         "ss-tr","x",
##         "ss-tfb","x",
##         "ss-trb","x",
##         "ss-trfb","x",
##         "ss-tfrb","x"
##     ),
##     ext.summary=TRUE,
##     RPG.gw.mean=TRUE,
##     RPG.reg.mean=TRUE,
##     LinkageDrag=TRUE,
##     MarkerAnalyses=TRUE,
##     result.files=c(
##         "ss-t.RData","ss-tf.RData","ss-tr.RData",
##         "ss-tfb.RData","ss-trb.RData",
##         "ss-trfb.RData","ss-tfrb.RData"
##     )
## )
## print(out)


###################################################
### code chunk number 70: ex13_01 (eval = FALSE)
###################################################
## rng.init(seed=2114)


###################################################
### code chunk number 71: ex13_02 (eval = FALSE)
###################################################
## mab.compare(
##     a.simulation.name="success-1",
##     a.simulation.run="x",
##     a.repetitions=NREP,
##     a.linkage.map="ex-mab1.map",
##     a.target.loci="target",
##     a.flanking.loci=c("umc90","umc27"),
##     a.gen.type=c("f1","bc","bc","bc","bc","bc"),
##     a.population.size=c(1,100,100,100,100,100),
##     a.sel.strategy=c("n","tb","tb","tb","tb","tb"),
##     a.no.selected=c(1,10,10,10,10,10),
##     a.success.factor=1,
##     b.simulation.name="success-2",
##     b.no.selected=c(1,10,10,10,10,10),
##     b.success.factor=0.5,
##     c.simulation.name="success-3",
##     c.no.selected=c(1,5,5,5,5,5),
##     c.success.factor=1,
##     result.files=c(
##         "success-1.RData",
##         "success-2.RData",
##         "success-3.RData"
##     )
## )


###################################################
### code chunk number 72: ex13_03 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=c(
##         "success-1","x",
##         "success-2","x",
##         "success-3","x"
##     ),
##     ext.summary=TRUE,
##     RPG.gw.mean=TRUE,
##     RPG.gw.Q10=TRUE,
##     result.files=c(
##         "success-1.RData",
##         "success-2.RData",
##         "success-3.RData"
##     )
## )
## print(out)


###################################################
### code chunk number 73: ex14_01 (eval = FALSE)
###################################################
## rng.init(seed=2115)


###################################################
### code chunk number 74: ex14_02 (eval = FALSE)
###################################################
## mab.simulate(
##     simulation.name="rec-01",
##     simulation.run="x",
##     donor.parent="ex-mabp1_m.mda",
##     recurrent.parent="ex-mabp2_m.mda",
##     linkage.map="ex-mab4.map",
##     result.file="rec-01.RData",
##     recode.infiles=TRUE,
##     recode.file="recoded-parents.mda"
## )


###################################################
### code chunk number 75: ex14_03 (eval = FALSE)
###################################################
## out <- list(
##     result=names(mab.load.data(
##         simulation.name="rec-01",
##         simulation.run="x",
##         result.file="rec-01.RData"
##     )),
##     recoded.file=head(
##         readLines(file.path(st.output.dir,"recoded-parents.mda")),
##         12
##     )
## )
## print(out)


###################################################
### code chunk number 76: ex15_01 (eval = FALSE)
###################################################
## rng.init(seed=2116)


###################################################
### code chunk number 77: ex15_02 (eval = FALSE)
###################################################
## mab.compare(
##     a.simulation.name="sel-50-0a-01",
##     a.simulation.run="x",
##     a.repetitions=NREP,
##     a.linkage.map="ex-mab1.map",
##     a.target.loci="target",
##     a.flanking.loci=c("umc90","umc27"),
##     a.gen.type=c("f1","bc","bc","bc"),
##     a.population.size=c(1,50,50,50),
##     a.sel.strategy=c("n","tfb","tfb","tfb"),
##     b.simulation.name="sel-50-05-01",
##     b.no.preselected=c(1,5,5,5),
##     result.files=c("sel-50-0a-01.RData","sel-50-05-01.RData")
## )


###################################################
### code chunk number 78: ex15_03 (eval = FALSE)
###################################################
## out <- mab.tabulate(
##     simulation.names.runs=c(
##         "sel-50-0a-01","x",
##         "sel-50-05-01","x"
##     ),
##     ext.summary=TRUE,
##     RPG.gw.Q10=TRUE,
##     LinkageDrag=TRUE,
##     MarkerAnalyses=TRUE,
##     result.files=c("sel-50-0a-01.RData","sel-50-05-01.RData")
## )
## print(out)


