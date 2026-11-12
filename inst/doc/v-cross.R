### R code from vignette source 'v-cross.Rnw'

###################################################
### code chunk number 1: v-cross.Rnw:205-213 (eval = FALSE)
###################################################
## library("SelectionTools")
## st.set.info.level(0)
## 
## data("v-tropmaize-vcf")
## data("v-tropmaize-phe")
## 
## st.load.vcf.data(v.tropmaize.vcf)
## st.load.performance.data(v.tropmaize.phe)


###################################################
### code chunk number 2: v-cross.Rnw:216-220 (eval = FALSE)
###################################################
## st.restrict.marker.data( NoAll.MAX = 2    )
## st.restrict.marker.data( MaMis.MAX = 0.05 )
## st.restrict.marker.data( ExHet.MIN = 0.1  )
## st.restrict.marker.data( InMis.MAX = 0.1  )


###################################################
### code chunk number 3: v-cross.Rnw:223-224 (eval = FALSE)
###################################################
## st.copy.marker.data("preprocessed","default")


###################################################
### code chunk number 4: v-cross.Rnw:232-236 (eval = FALSE)
###################################################
## st.copy.marker.data( "hblocks01", "preprocessed" )
## h <- st.def.hblocks( hap      = 5,     
##                      hap.unit = 1,  
##                      data.set = "hblocks01" )


###################################################
### code chunk number 5: v-cross.Rnw:240-241 (eval = FALSE)
###################################################
## st.recode.hil ( data.set="hblocks01" )      


###################################################
### code chunk number 6: v-cross.Rnw:244-247 (eval = FALSE)
###################################################
## st.plot.ggt( f.ind    = 1,
##              l.ind    = 20,
##              data.set = "hblocks01")    


###################################################
### code chunk number 7: v-cross.Rnw:258-260 (eval = FALSE)
###################################################
## h.first <- h[1:10, ]
## h.first


###################################################
### code chunk number 8: v-cross.Rnw:289-297 (eval = FALSE)
###################################################
## st.copy.marker.data( "hblocks02", "preprocessed" )
## 
## h2 <- st.def.hblocks( hap      = 5,     
##                       hap.unit = 2,  
##                       data.set = "hblocks02" )
## 
## h2.first <- h2[1:10, ]
## h2.first


###################################################
### code chunk number 9: v-cross.Rnw:335-341 (eval = FALSE)
###################################################
## st.plot.ggt(f.ind                  = 1,
##             l.ind                  = 30,
##             sort.performance       = TRUE,
##             show.performance       = TRUE,
##             performance.decreasing = TRUE,
##             data.set               = "hblocks01")


###################################################
### code chunk number 10: v-cross.Rnw:354-355 (eval = FALSE)
###################################################
## st.copy.marker.data( "hblocks03", "hblocks01" )


###################################################
### code chunk number 11: v-cross.Rnw:359-361 (eval = FALSE)
###################################################
## st.recode.ref( reference = "142",
##                data.set  = "hblocks03")


###################################################
### code chunk number 12: v-cross.Rnw:365-371 (eval = FALSE)
###################################################
## st.plot.ggt( f.ind                  = 1,
##              l.ind                  = 30,
##              sort.performance       = TRUE,
##              show.performance       = TRUE,
##              performance.decreasing = TRUE,
##              data.set               = "hblocks03")


###################################################
### code chunk number 13: v-cross.Rnw:381-382 (eval = FALSE)
###################################################
## st.copy.marker.data( "gd01", "preprocessed" )


###################################################
### code chunk number 14: v-cross.Rnw:386-387 (eval = FALSE)
###################################################
## best <- st.select.phen( n=50, data.set="gd01" )


###################################################
### code chunk number 15: v-cross.Rnw:390-393 (eval = FALSE)
###################################################
## 
## st.restrict.marker.data ( ind.list=best$i,
##                           data.set="gd01" )


###################################################
### code chunk number 16: v-cross.Rnw:402-408 (eval = FALSE)
###################################################
## dist.mat <- st.genetic.distances( measure  = "rd",
##                                   format   = "m",
##                                   data.set = "gd01" )
## 
## dm <- as.matrix(dist.mat)
## heatmap( dm, scale="none", col=heat.colors(12) )


###################################################
### code chunk number 17: v-cross.Rnw:419-423 (eval = FALSE)
###################################################
## gs.cross.eval.gd ( dist="mrd", data.set="gd01" )
## crosses <- gs.cross.info ( sortby="gd", bestn=50, data.set="gd01" )
## distance.crosses <- crosses[1:10, ]
## distance.crosses


###################################################
### code chunk number 18: v-cross.Rnw:441-443 (eval = FALSE)
###################################################
## st.def.hblocks( hap=5, data.set="gd01" )
## st.recode.hil ( data.set="gd01" )


###################################################
### code chunk number 19: v-cross.Rnw:446-461 (eval = FALSE)
###################################################
## st.set.info.level(-2)
## for (i in seq_len(min(10,nrow(crosses)))) {
##     a <- crosses$P1Name[i]
##     b <- crosses$P2Name[i]
##     f.nme <- sprintf("ex-crs-07-%03i.pdf",100+i)
##     pdf(f.nme,width=8,height=0.8,pointsize=10)
##     st.plot.ggt(i.list=c(a,b),
##                 show.performance=TRUE,
##                 cex.chrom=0.8,
##                 p.t=FALSE,
##                 p.s=FALSE,
##                 data.set="gd01")
##     dev.off()
## }
## st.set.info.level(0)


###################################################
### code chunk number 20: v-cross.Rnw:487-488 (eval = FALSE)
###################################################
## st.copy.marker.data( "effects01", "preprocessed" )


###################################################
### code chunk number 21: v-cross.Rnw:494-495 (eval = FALSE)
###################################################
## gs.esteff.rr( method="BLUP", data.set="effects01" )


###################################################
### code chunk number 22: v-cross.Rnw:506-509 (eval = FALSE)
###################################################
## gs.cross.eval.mu( data.set="effects01" )
## gs.cross.eval.mi( data.set="effects01" )
## gs.cross.eval.ma( data.set="effects01" )


###################################################
### code chunk number 23: v-cross.Rnw:517-518 (eval = FALSE)
###################################################
## gs.cross.eval.va( pop.type="DH", data.set="effects01" )


###################################################
### code chunk number 24: v-cross.Rnw:528-529 (eval = FALSE)
###################################################
## gs.cross.eval.es( alpha=0.25, data.set="effects01" )


###################################################
### code chunk number 25: v-cross.Rnw:539-540 (eval = FALSE)
###################################################
## gs.cross.eval.es( N=1, G=4, data.set="effects01" )


###################################################
### code chunk number 26: v-cross.Rnw:547-548 (eval = FALSE)
###################################################
## gs.cross.eval.gd(dist="rd",data.set="effects01")


###################################################
### code chunk number 27: v-cross.Rnw:559-568 (eval = FALSE)
###################################################
## gs.cross.eval.gd(dist="rd",data.set="effects01")
## gs.cross.eval.mu(data.set="effects01")
## gs.cross.eval.ma(data.set="effects01")
## gs.cross.eval.mi(data.set="effects01")
## gs.cross.eval.va(pop.type="DH",data.set="effects01")
## gs.cross.eval.es(N=1,G=4,data.set="effects01")
## crosses <- gs.cross.info(sortby="es",bestn=20,data.set="effects01")
## formatted.crosses <- format(crosses,digits=2,nsmall=2)
## formatted.crosses


###################################################
### code chunk number 28: v-cross.Rnw:619-625 (eval = FALSE)
###################################################
## values <- gs.estimate.gv(estimation.set="effects01",
##                          prediction.set="effects01")
## merit <- setNames(values$yhat, values$i)
## A <- st.get.A("effects01")
## merit <- merit[rownames(A)]
## K <- A/2


###################################################
### code chunk number 29: v-cross.Rnw:636-644 (eval = FALSE)
###################################################
## contributions <- st.ocs(
##     merit = merit,                             # Parental genetic values
##     K = K,                                     # Kinship matrix
##     theta = 0.02,                              # Group relatedness limit
##     lower = 0,                                 # Minimum parental contribution
##     upper = 0.1                                # Maximum parental contribution
## )
## contributions


###################################################
### code chunk number 30: v-cross.Rnw:732-738 (eval = FALSE)
###################################################
## marker.effects <- gs.return.effects(data.set="effects01")
## st.copy.marker.data("cross.parents", "effects01")
## st.restrict.marker.data(ind.list=names(merit)[1:30],
##                         data.set="cross.parents")
## st.restrict.marker.data(MaMis.MAX=0, data.set="cross.parents")
## cross.map <- st.get.map(data.set="cross.parents")


###################################################
### code chunk number 31: v-cross.Rnw:745-750 (eval = FALSE)
###################################################
## cross.effects <- marker.effects[
##   marker.effects$marker %in% c("mu", cross.map$Name), ]
## gs.set.effects(cross.effects, data.set="cross.parents")
## cross.values <- gs.estimate.gv(estimation.set="cross.parents",
##                                prediction.set="cross.parents")


###################################################
### code chunk number 32: v-cross.Rnw:757-761 (eval = FALSE)
###################################################
## cross.merit <- setNames(cross.values$yhat, cross.values$i)
## cross.A <- st.get.A("cross.parents")
## cross.merit <- cross.merit[rownames(cross.A)]
## cross.K <- cross.A/2


###################################################
### code chunk number 33: v-cross.Rnw:771-780 (eval = FALSE)
###################################################
## cross.plan <- st.cross.plan(
##     merit = cross.merit,                       # Parental genetic values
##     K = cross.K,                               # Kinship matrix
##     n.crosses = 15,                            # Number of distinct parent pairs
##     n.progeny = 100,                           # Progeny per family
##     theta = 0.02,                              # Group relatedness limit
##     upper = 0.1                                # Maximum parental contribution
## )
## cross.plan


###################################################
### code chunk number 34: v-cross.Rnw:854-868 (eval = FALSE)
###################################################
## fixed.plan <- st.cross.plan(
##     merit = cross.merit,                       # Parental genetic values
##     K = cross.K,                               # Kinship matrix
##     n.crosses = 15,                            # Number of distinct parent pairs
##     n.progeny = 100,                           # Progeny per family
##     theta = 0.02,                              # Group relatedness limit
##     upper = 0.1,                               # Maximum parental contribution
##     contributions = cross.plan$contributions,  # Fixed parental gene shares
##     mode = "fixed",                            # Pairing at fixed contributions
##     control = list(                            # Search settings
##         start = cross.plan                     # Starting crossing plan
##     )
## )
## fixed.plan


###################################################
### code chunk number 35: v-cross.Rnw:910-929 (eval = FALSE)
###################################################
## st.set.simpop(pop.name="Founders", data.set="cross.parents")
## rng.init(seed=1001)
## founder.names <- population.pedigree("Founders")$name
## planned.crosses <- cross.plan$crosses
## 
## for (i in seq_len(nrow(planned.crosses))) {
##     population.copy("P1", "Founders",
##                     start=match(planned.crosses$P1Name[i], founder.names),
##                     n=1)
##     population.copy("P2", "Founders",
##                     start=match(planned.crosses$P2Name[i], founder.names),
##                     n=1)
##     cross("F1", "P1", "P2", NoPg=1, crossing.scheme=3)
##     dh(sprintf("DH%03d", i), "F1",
##        NoPg=planned.crosses$n.progeny[i])
## }
## population.remove(c("P1", "P2", "F1"))
## families <- population.list()
## families


###################################################
### code chunk number 36: v-cross.Rnw:960-965 (eval = FALSE)
###################################################
## for (i in seq_len(nrow(planned.crosses))) {
##     population.concat("DH", sprintf("DH%03d", i))
## }
## merged.populations <- population.list()
## merged.populations


###################################################
### code chunk number 37: v-cross.Rnw:977-981 (eval = FALSE)
###################################################
## genotype.population("Founders DH")
## evaluate.population("Founders DH", effName="effect")
## genetic.summary <- summarize.gvalue("Founders DH", effectfile="effect")
## genetic.summary


###################################################
### code chunk number 38: v-cross.Rnw:994-997 (eval = FALSE)
###################################################
## selection.gain <- genetic.summary["DH", "Mean"] -
##                   genetic.summary["Founders", "Mean"]
## selection.gain


###################################################
### code chunk number 39: v-cross.Rnw:1007-1016 (eval = FALSE)
###################################################
## prediction.comparison <- data.frame(
##     Predicted = c(cross.plan$gain,
##                   cross.plan$gain - genetic.summary["Founders", "Mean"]),
##     Simulated = c(genetic.summary["DH", "Mean"], selection.gain),
##     row.names = c("Population mean", "Selection gain")
## )
## prediction.comparison$Difference <- prediction.comparison$Simulated -
##                                     prediction.comparison$Predicted
## prediction.comparison


###################################################
### code chunk number 40: v-cross.Rnw:1046-1050 (eval = FALSE)
###################################################
## gs.cross.eval.mu(data.set="cross.parents")
## gs.cross.eval.va(pop.type="DH", data.set="cross.parents")
## gs.cross.eval.es(N=10, G=100, data.set="cross.parents")
## cross.scores <- gs.cross.info(data.set="cross.parents")


###################################################
### code chunk number 41: v-cross.Rnw:1059-1073 (eval = FALSE)
###################################################
## mean.plan <- st.cross.plan(
##     merit = cross.merit,                       # Parental genetic values
##     K = cross.K,                               # Kinship matrix
##     n.crosses = 15,                            # Number of distinct parent pairs
##     n.progeny = 100,                           # Progeny per family
##     theta = 0.02,                              # Group relatedness limit
##     upper = 0.1,                               # Maximum parental contribution
##     scores = cross.scores,                     # Precomputed cross scores
##     criterion = "mean",                        # Objective: maximum cross mean
##     control = list(                            # Search settings
##         start = cross.plan                     # Starting crossing plan
##     )
## )
## mean.plan


###################################################
### code chunk number 42: v-cross.Rnw:1115-1130 (eval = FALSE)
###################################################
## usefulness.plan <- st.cross.plan(
##     merit = cross.merit,                       # Parental genetic values
##     K = cross.K,                               # Kinship matrix
##     n.crosses = 15,                            # Number of distinct parent pairs
##     n.progeny = 100,                           # Progeny per family
##     theta = 0.02,                              # Group relatedness limit
##     upper = 0.1,                               # Maximum parental contribution
##     scores = cross.scores,                     # Precomputed cross scores
##     criterion = "usefulness",                  # Objective: maximum usefulness
##     accuracy = 1,                              # Within-family ranking accuracy
##     control = list(                            # Search settings
##         start = cross.plan                     # Starting crossing plan
##     )
## )
## usefulness.plan


###################################################
### code chunk number 43: v-cross.Rnw:1192-1196 (eval = FALSE)
###################################################
## oma.parents <- data.frame(id=rownames(cross.K), min=0, max=0.1)
## oma.mean.matings <- data.frame(
##   parent1=cross.scores$P1Name, parent2=cross.scores$P2Name,
##   merit=cross.scores$mu, min=0, max=0.1)


###################################################
### code chunk number 44: v-cross.Rnw:1211-1219 (eval = FALSE)
###################################################
## oma.mean <- st.oma(
##     parents = oma.parents,                     # Parent contribution bounds
##     matings = oma.mean.matings,                # Allowed crosses and scores
##     K = cross.K,                               # Kinship matrix
##     dF = c(0, 0.01),                           # F1 change range; F2 change cap
##     base = "RM"                                # Random-mating reference
## )
## oma.mean


###################################################
### code chunk number 45: v-cross.Rnw:1306-1308 (eval = FALSE)
###################################################
## oma.usefulness.matings <- oma.mean.matings
## oma.usefulness.matings$merit <- cross.scores$es


###################################################
### code chunk number 46: v-cross.Rnw:1315-1323 (eval = FALSE)
###################################################
## oma.usefulness <- st.oma(
##     parents = oma.parents,                     # Parent contribution bounds
##     matings = oma.usefulness.matings,          # Allowed crosses and scores
##     K = cross.K,                               # Kinship matrix
##     dF = c(0, 0.01),                           # F1 change range; F2 change cap
##     base = "RM"                                # Random-mating reference
## )
## oma.usefulness


###################################################
### code chunk number 47: v-cross.Rnw:1379-1386 (eval = FALSE)
###################################################
## parent.slots <- c(usefulness.plan$crosses$P1Name,
##                   usefulness.plan$crosses$P2Name)
## diversity.loss <- propsd(
##     G = cross.K,                               # Complete reference covariance
##     parents = parent.slots                     # Parental occurrences in families
## )
## diversity.loss


###################################################
### code chunk number 48: v-cross.Rnw:1411-1432 (eval = FALSE)
###################################################
## diversity.plan <- st.cross.plan(
##     merit = cross.merit,                       # Parental genetic values
##     K = cross.K,                               # Kinship matrix
##     n.crosses = 15,                            # Number of distinct parent pairs
##     n.progeny = 100,                           # Progeny per family
##     theta = 0.02,                              # Group relatedness limit
##     upper = 0.1,                               # Maximum parental contribution
##     scores = cross.scores,                     # Precomputed cross scores
##     criterion = "usefulness",                  # Objective: maximum usefulness
##     accuracy = 1,                              # Within-family ranking accuracy
##     diversity = list(                          # Additional diversity constraint
##         G = cross.K,                           # Complete reference covariance
##         max = 0.025,                           # Strict upper limit on PropSD
##         sample = "family.slots"                # Two parental occurrences per family
##     ),
##     control = list(                            # Search settings
##         starts = 4,                            # Number of search starts
##         maxit = 5000                           # Proposals per search start
##     )
## )
## diversity.plan


###################################################
### code chunk number 49: v-cross.Rnw:1484-1494 (eval = FALSE)
###################################################
## progeny.allocation <- st.cross.allocate(
##     plan = diversity.plan,                     # Fixed parent pairs and constraints
##     scores = cross.scores,                     # Cross means and segregation variances
##     total.progeny = 1500,                      # Total produced progeny budget
##     n.selected = 10,                           # Retained descendants per family
##     min.progeny = 10,                          # Minimum produced count per family
##     max.progeny = 1500,                        # Maximum produced count per family
##     accuracy = 1                               # Within-family ranking accuracy
## )
## progeny.allocation


