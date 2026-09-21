### R code from vignette source 'v-hapld.Rnw'

###################################################
### code chunk number 1: v-hapld.Rnw:191-196 (eval = FALSE)
###################################################
## library("SelectionTools")
## 
## vignette.data <- new.env(parent=emptyenv())
## data("v-tropmaize-vcf", package="SelectionTools", envir=vignette.data)
## st.load.vcf.data(vignette.data$v.tropmaize.vcf, data.set="base")


###################################################
### code chunk number 2: v-hapld.Rnw:204-208 (eval = FALSE)
###################################################
## st.copy.marker.data("clustering", "base")
## st.restrict.marker.data(NoAll.MAX=2, MaMis.MAX=0.1, data.set="clustering")
## st.restrict.marker.data(InMis.MAX=0.1, data.set="clustering")
## st.restrict.marker.data(ExHet.MIN=0.1, data.set="clustering")


###################################################
### code chunk number 3: v-hapld.Rnw:211-214 (eval = FALSE)
###################################################
## dist.mat <- st.genetic.distances( measure="rd",
##                                   format="m",
##                                   data.set="clustering")


###################################################
### code chunk number 4: v-hapld.Rnw:219-224 (eval = FALSE)
###################################################
## pco <- cmdscale( dist.mat,
##                  k  = 4,
## 		 eig= TRUE)
## PC <- pco$points
## var.exp <- round(pco$eig / sum(pco$eig) * 100, 2)


###################################################
### code chunk number 5: v-hapld.Rnw:228-232 (eval = FALSE)
###################################################
## set.seed(1)
## clus <- kmeans(as.matrix(dist.mat), centers=5, nstart=25)
## clusters <- clus$cluster
## cl.c <- as.character(clusters)


###################################################
### code chunk number 6: v-hapld.Rnw:235-243 (eval = FALSE)
###################################################
## cl <- function(i)
##   paste0("Coordinate ", i, " (", sprintf("%.2f", var.exp[i]), "%)")
## 
## par(mfrow=c(2,2), mar=c(4,4,2,2))
## plot(PC[,1], PC[,2], pch=cl.c, col=clusters, xlab=cl(1), ylab=cl(2))
## plot(PC[,1], PC[,3], pch=cl.c, col=clusters, xlab=cl(1), ylab=cl(3))
## plot(PC[,2], PC[,3], pch=cl.c, col=clusters, xlab=cl(2), ylab=cl(3))
## plot(PC[,1], PC[,4], pch=cl.c, col=clusters, xlab=cl(1), ylab=cl(4))


###################################################
### code chunk number 7: v-hapld.Rnw:249-251 (eval = FALSE)
###################################################
## target.cluster <- unname(clusters["142"])
## target.ind <- names(clusters)[clusters == target.cluster]


###################################################
### code chunk number 8: v-hapld.Rnw:255-260 (eval = FALSE)
###################################################
## st.copy.marker.data("q3", "base")
## st.restrict.marker.data(ind.list=target.ind, data.set="q3")
## st.restrict.marker.data(MaMis.MAX=0.1, data.set="q3")
## st.restrict.marker.data(InMis.MAX=0.1, data.set="q3")
## st.restrict.marker.data(ExHet.MIN=0.1, data.set="q3")


###################################################
### code chunk number 9: v-hapld.Rnw:272-275 (eval = FALSE)
###################################################
## ld <- st.calc.ld( ld.measure="r2",
##                   data.set="q3")
## ld[16042:16052,]


###################################################
### code chunk number 10: v-hapld.Rnw:300-308 (eval = FALSE)
###################################################
## chromosome <- 8
## l <- st.LDplot.ld(ld, chromosome)
## m <- st.LDplot.map(chromosome, data.set="q3")
## st.LDheatmap(l, m,
##              distances="genetical",
##              LDmeasure="r2",
##              title=paste("r2, Chrom.", chromosome),
##              col=heat.colors(20))


###################################################
### code chunk number 11: v-hapld.Rnw:316-324 (eval = FALSE)
###################################################
## ld <- st.calc.ld(ld.measure="Dp", data.set="q3")
## l <- st.LDplot.ld(ld, chrom=chromosome)
## m <- st.LDplot.map(chrom=chromosome, data.set="q3")
## st.LDheatmap(l, m,
##              distances="genetical",
##              LDmeasure="Dprime",
##              title=paste("Dprime, Chrom.", chromosome),
##              col=heat.colors(20))


###################################################
### code chunk number 12: v-hapld.Rnw:332-344 (eval = FALSE)
###################################################
## LDplot <- function(ld.list, chrom, data.set, ld.measure="r2",
##                    title="", no.colors=20)
## {
##   l <- st.LDplot.ld(ld.list, chrom)
##   m <- st.LDplot.map(chrom, data.set)
##   if ("Dp" == ld.measure) LDmeasure <- "Dprime" else LDmeasure <- "r2"
##   st.LDheatmap(l, m,
##                distances="genetical",
##                LDmeasure=LDmeasure,
##                title=paste("Chrom.", chrom, title),
##                col=heat.colors(no.colors))
## }


###################################################
### code chunk number 13: v-hapld.Rnw:372-383 (eval = FALSE)
###################################################
## st.copy.marker.data("t3", "q3")
## ld <- st.calc.ld(ld.measure="r2", data.set="t3")
## LDplot(ld, chrom=8, data.set="t3", no.colors=5, title="before building blocks")
## 
## st.def.hblocks(ld.threshold=0.8,
##                ld.criterion="flanking",
##                data.set="t3")
## st.recode.hil(data.set="t3")
## 
## ld <- st.calc.ld(ld.measure="r2", data.set="t3")
## LDplot(ld, chrom=8, data.set="t3", no.colors=5, title="after building blocks, r2=0.8")


###################################################
### code chunk number 14: v-hapld.Rnw:396-405 (eval = FALSE)
###################################################
## st.copy.marker.data("t3r06", "q3")
## st.calc.ld(ld.measure="r2", data.set="t3r06")
## st.def.hblocks(ld.threshold=0.6,
##                ld.criterion="flanking",
##                data.set="t3r06")
## st.recode.hil(data.set="t3r06")
## ld06 <- st.calc.ld(ld.measure="r2", data.set="t3r06")
## LDplot(ld06, chrom=8, data.set="t3r06",
##        no.colors=5, title="after building blocks, r2=0.6")


###################################################
### code chunk number 15: v-hapld.Rnw:408-417 (eval = FALSE)
###################################################
## st.copy.marker.data("t3r04", "q3")
## st.calc.ld(ld.measure="r2", data.set="t3r04")
## st.def.hblocks(ld.threshold=0.4,
##                ld.criterion="flanking",
##                data.set="t3r04")
## st.recode.hil(data.set="t3r04")
## ld04 <- st.calc.ld(ld.measure="r2", data.set="t3r04")
## LDplot(ld04, chrom=8, data.set="t3r04",
##        no.colors=5, title="after building blocks, r2=0.4")


###################################################
### code chunk number 16: v-hapld.Rnw:429-439 (eval = FALSE)
###################################################
## st.copy.marker.data("t3tol", "q3")
## st.calc.ld(ld.measure="r2", data.set="t3tol")
## st.def.hblocks(ld.threshold=0.8,
##                tolerance=1,
##                ld.criterion="flanking",
##                data.set="t3tol")
## st.recode.hil(data.set="t3tol")
## ldtol <- st.calc.ld(ld.measure="r2", data.set="t3tol")
## LDplot(ldtol, chrom=8, data.set="t3tol",
##        no.colors=5, title="after building blocks, r2=0.8, tolerance=1")


###################################################
### code chunk number 17: v-hapld.Rnw:448-457 (eval = FALSE)
###################################################
## st.copy.marker.data("t3avg", "q3")
## st.calc.ld(ld.measure="r2", data.set="t3avg")
## st.def.hblocks(ld.threshold=0.6,
##                ld.criterion="average",
##                data.set="t3avg")
## st.recode.hil(data.set="t3avg")
## ldavg <- st.calc.ld(ld.measure="r2", data.set="t3avg")
## LDplot(ldavg, chrom=8, data.set="t3avg",
##        no.colors=5, title="after building blocks, r2=0.6 average LD")


###################################################
### code chunk number 18: v-hapld.Rnw:469-472 (eval = FALSE)
###################################################
## options(width=120)
## st.copy.marker.data("t3detail", "q3")
## ld <- st.calc.ld(ld.measure="r2", data.set="t3detail")


###################################################
### code chunk number 19: v-hapld.Rnw:488-489 (eval = FALSE)
###################################################
## x1 <- st.marker.data.statistics(data.set="t3detail")


###################################################
### code chunk number 20: v-hapld.Rnw:510-513 (eval = FALSE)
###################################################
## h <- st.def.hblocks(ld.threshold=0.8,
##                     ld.criterion="flanking",
##                     data.set="t3detail")


###################################################
### code chunk number 21: v-hapld.Rnw:532-533 (eval = FALSE)
###################################################
## r <- st.recode.hil(data.set="t3detail")


###################################################
### code chunk number 22: v-hapld.Rnw:551-552 (eval = FALSE)
###################################################
## x2 <- st.marker.data.statistics(data.set="t3detail")


###################################################
### code chunk number 23: v-hapld.Rnw:575-600 (eval = FALSE)
###################################################
## # Use the globally filtered clustering data so both populations have the
## # same marker set and marker order. The comparison population is the
## # largest of the four clusters other than the target cluster.
## other.sizes <- table(clusters[clusters != target.cluster])
## other.cluster <- as.integer(names(other.sizes)[which.max(other.sizes)])
## other.ind <- names(clusters)[clusters == other.cluster]
## 
## st.copy.marker.data("qtarget", "clustering")
## st.copy.marker.data("qother", "clustering")
## st.restrict.marker.data(ind.list=target.ind, data.set="qtarget")
## st.restrict.marker.data(ind.list=other.ind, data.set="qother")
## 
## st.copy.marker.data("ttarget", "qtarget")
## st.copy.marker.data("tother", "qother")
## st.calc.ld(ld.measure="r2", data.set="ttarget")
## h <- st.def.hblocks(ld.threshold=0.8,
##                     tolerance=3,
##                     ld.criterion="flanking",
##                     data.set="ttarget")
## 
## h2 <- st.set.hblocks(haplotype.list=h,
##                      hap.symbol="c",
##                      data.set="tother")
## 
## st.recode.hil(data.set="tother")


###################################################
### code chunk number 24: v-hapld.Rnw:693-701 (eval = FALSE)
###################################################
## st.copy.marker.data("q3phase", "q3")
## st.restrict.marker.data(NoAll.MAX=2, data.set="q3phase")
## st.restrict.marker.data(ExHet.MIN=0.001, data.set="q3phase")
## 
## ld.phase.r2 <- st.calc.ld(ld.measure="r2-estimate-phases",
##                           data.set="q3phase")
## ld.phase.Dp <- st.calc.ld(ld.measure="Dp-estimate-phases",
##                           data.set="q3phase")


