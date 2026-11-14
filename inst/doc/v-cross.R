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
### code chunk number 19: v-cross.Rnw:446-460 (eval = FALSE)
###################################################
## st.set.info.level(-2)
## for (i in seq_len(min(10,nrow(crosses)))) {
##     a <- crosses$P1Name[i]
##     b <- crosses$P2Name[i]
##     f.nme <- sprintf("ex-crs-07-%03i.pdf",100+i)
##     pdf(f.nme,width=8,height=0.8,pointsize=10)
##     st.plot.ggt(i.list=c(a,b),
##                 show.performance=TRUE,
##                 cex.chrom=1,
##                 p.t=FALSE,
##                 p.s=FALSE,
##                 data.set="gd01")
##     dev.off()
## }


###################################################
### code chunk number 20: v-cross.Rnw:486-487 (eval = FALSE)
###################################################
## st.copy.marker.data( "effects01", "preprocessed" )


###################################################
### code chunk number 21: v-cross.Rnw:493-494 (eval = FALSE)
###################################################
## gs.esteff.rr( method="BLUP", data.set="effects01" )


###################################################
### code chunk number 22: v-cross.Rnw:505-508 (eval = FALSE)
###################################################
## gs.cross.eval.mu( data.set="effects01" )
## gs.cross.eval.mi( data.set="effects01" )
## gs.cross.eval.ma( data.set="effects01" )


###################################################
### code chunk number 23: v-cross.Rnw:516-517 (eval = FALSE)
###################################################
## gs.cross.eval.va( pop.type="DH", data.set="effects01" )


###################################################
### code chunk number 24: v-cross.Rnw:527-528 (eval = FALSE)
###################################################
## gs.cross.eval.es( alpha=0.25, data.set="effects01" )


###################################################
### code chunk number 25: v-cross.Rnw:538-539 (eval = FALSE)
###################################################
## gs.cross.eval.es( N=1, G=4, data.set="effects01" )


###################################################
### code chunk number 26: v-cross.Rnw:546-547 (eval = FALSE)
###################################################
## gs.cross.eval.gd(dist="rd",data.set="effects01")


###################################################
### code chunk number 27: v-cross.Rnw:558-567 (eval = FALSE)
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
### code chunk number 28: v-cross.Rnw:618-624 (eval = FALSE)
###################################################
## values <- gs.estimate.gv(estimation.set="effects01",
##                          prediction.set="effects01")
## merit <- setNames(values$yhat, values$i)
## A <- st.get.A("effects01")
## merit <- merit[rownames(A)]
## K <- A/2


###################################################
### code chunk number 29: v-cross.Rnw:635-643 (eval = FALSE)
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
### code chunk number 30: v-cross.Rnw:731-737 (eval = FALSE)
###################################################
## marker.effects <- gs.return.effects(data.set="effects01")
## st.copy.marker.data("cross.parents", "effects01")
## st.restrict.marker.data(ind.list=names(merit)[1:30],
##                         data.set="cross.parents")
## st.restrict.marker.data(MaMis.MAX=0, data.set="cross.parents")
## cross.map <- st.get.map(data.set="cross.parents")


###################################################
### code chunk number 31: v-cross.Rnw:744-749 (eval = FALSE)
###################################################
## cross.effects <- marker.effects[
##   marker.effects$marker %in% c("mu", cross.map$Name), ]
## gs.set.effects(cross.effects, data.set="cross.parents")
## cross.values <- gs.estimate.gv(estimation.set="cross.parents",
##                                prediction.set="cross.parents")


###################################################
### code chunk number 32: v-cross.Rnw:756-760 (eval = FALSE)
###################################################
## cross.merit <- setNames(cross.values$yhat, cross.values$i)
## cross.A <- st.get.A("cross.parents")
## cross.merit <- cross.merit[rownames(cross.A)]
## cross.K <- cross.A/2


###################################################
### code chunk number 33: v-cross.Rnw:770-779 (eval = FALSE)
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
### code chunk number 34: v-cross.Rnw:853-867 (eval = FALSE)
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
### code chunk number 35: v-cross.Rnw:909-928 (eval = FALSE)
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
### code chunk number 36: v-cross.Rnw:959-964 (eval = FALSE)
###################################################
## for (i in seq_len(nrow(planned.crosses))) {
##     population.concat("DH", sprintf("DH%03d", i))
## }
## merged.populations <- population.list()
## merged.populations


###################################################
### code chunk number 37: v-cross.Rnw:976-980 (eval = FALSE)
###################################################
## genotype.population("Founders DH")
## evaluate.population("Founders DH", effName="effect")
## genetic.summary <- summarize.gvalue("Founders DH", effectfile="effect")
## genetic.summary


###################################################
### code chunk number 38: v-cross.Rnw:993-996 (eval = FALSE)
###################################################
## selection.gain <- genetic.summary["DH", "Mean"] -
##                   genetic.summary["Founders", "Mean"]
## selection.gain


###################################################
### code chunk number 39: v-cross.Rnw:1006-1015 (eval = FALSE)
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
### code chunk number 40: v-cross.Rnw:1045-1049 (eval = FALSE)
###################################################
## gs.cross.eval.mu(data.set="cross.parents")
## gs.cross.eval.va(pop.type="DH", data.set="cross.parents")
## gs.cross.eval.es(N=10, G=100, data.set="cross.parents")
## cross.scores <- gs.cross.info(data.set="cross.parents")


###################################################
### code chunk number 41: v-cross.Rnw:1058-1072 (eval = FALSE)
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
### code chunk number 42: v-cross.Rnw:1114-1129 (eval = FALSE)
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
### code chunk number 43: v-cross.Rnw:1191-1195 (eval = FALSE)
###################################################
## oma.parents <- data.frame(id=rownames(cross.K), min=0, max=0.1)
## oma.mean.matings <- data.frame(
##   parent1=cross.scores$P1Name, parent2=cross.scores$P2Name,
##   merit=cross.scores$mu, min=0, max=0.1)


###################################################
### code chunk number 44: v-cross.Rnw:1210-1218 (eval = FALSE)
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
### code chunk number 45: v-cross.Rnw:1305-1307 (eval = FALSE)
###################################################
## oma.usefulness.matings <- oma.mean.matings
## oma.usefulness.matings$merit <- cross.scores$es


###################################################
### code chunk number 46: v-cross.Rnw:1314-1322 (eval = FALSE)
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
### code chunk number 47: v-cross.Rnw:1378-1385 (eval = FALSE)
###################################################
## parent.slots <- c(usefulness.plan$crosses$P1Name,
##                   usefulness.plan$crosses$P2Name)
## diversity.loss <- propsd(
##     G = cross.K,                               # Complete reference covariance
##     parents = parent.slots                     # Parental occurrences in families
## )
## diversity.loss


###################################################
### code chunk number 48: v-cross.Rnw:1410-1431 (eval = FALSE)
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
### code chunk number 49: v-cross.Rnw:1485-1495 (eval = FALSE)
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


###################################################
### code chunk number 50: v-cross.Rnw:1548-1562 (eval = FALSE)
###################################################
## progeny.order.allocation <- st.cross.allocate(
##     plan = diversity.plan,                     # Fixed parent pairs and constraints
##     scores = cross.scores,                     # Cross means and segregation variances
##     total.progeny = 1500,                      # Total produced progeny budget
##     n.selected = 10,                           # Retained descendants per family
##     min.progeny = 10,                          # Minimum produced count per family
##     max.progeny = 1500,                        # Maximum produced count per family
##     accuracy = 1,                              # Within-family ranking accuracy
##     selection.method = "order.statistics",     # Finite-normal selection expectation
##     control = list(                            # Allocation settings
##         algorithm = "auto"                     # Automatic discrete allocator
##     )
## )
## progeny.order.allocation


###################################################
### code chunk number 51: v-cross.Rnw:1609-1614 (eval = FALSE)
###################################################
## cross.population <- st.cross.evaluate(
##     plan = progeny.order.allocation,           # Fixed families and produced counts
##     K = cross.K                                # Original named kinship reference
## )
## cross.population


###################################################
### code chunk number 52: v-cross.Rnw:1680-1681 (eval = FALSE)
###################################################
## allocation.validation


###################################################
### code chunk number 53: v-cross.Rnw:1718-1747 (eval = FALSE)
###################################################
## practical.ids <- names(cross.merit)[1:6]
## practical.K <- cross.K[practical.ids, practical.ids]
## st.copy.marker.data(
##     "practical.parents",                       # Six-parent copy for this policy
##     "cross.parents"                            # Existing real marker panel
## )
## st.restrict.marker.data(
##     ind.list = practical.ids,                  # Six parents from the real panel
##     data.set = "practical.parents"             # Six-parent marker data set
## )
## gs.set.effects(
##     cross.effects,                             # Fitted effects for the same markers
##     data.set = "practical.parents"             # Six-parent marker data set
## )
## gs.cross.eval.mu(
##     data.set = "practical.parents"             # Recomputed cross means
## )
## gs.cross.eval.va(
##     pop.type = "DH",                           # Segregation among derived DH lines
##     data.set = "practical.parents"             # Six-parent marker data set
## )
## gs.cross.eval.es(
##     N = 2,                                     # Retained individuals per family
##     G = 20,                                    # Produced individuals per family
##     data.set = "practical.parents"             # Six-parent marker data set
## )
## practical.scores <- gs.cross.info(
##     data.set = "practical.parents"             # Scores for exactly the six parents
## )


###################################################
### code chunk number 54: v-cross.Rnw:1757-1777 (eval = FALSE)
###################################################
## practical.parents <- data.frame(
##     parent = practical.ids,                    # Actual parent identifiers
##     pool = rep(c("seed", "pollen"), each=3),   # Operational parent groups
##     seed = rep(c(TRUE, FALSE), each=3),        # First group eligible as seed parent
##     pollen = rep(c(FALSE, TRUE), each=3)       # Second group eligible as pollen parent
## )
## practical.limits <- data.frame(
##     parent = NA_character_,                    # Group-level parent counts
##     pool = c(NA, "seed", "pollen"),            # All parents or the named group
##     role = c("all", "seed", "pollen"),         # Count parents in the stated role
##     basis = "family.slots",                    # Distinct-family participation
##     measure = "parents",                       # Number of distinct used parents
##     lower = c(5, 2, 3),                        # Exact total, seed and pollen counts
##     upper = c(5, 2, 3)                         # Equal lower and upper bounds
## )
## practical.equal.use <- data.frame(
##     pool = "seed",                             # Seed-parent group
##     role = "seed",                             # Seed use only
##     basis = "family.slots"                     # Equal integer family use
## )


###################################################
### code chunk number 55: v-cross.Rnw:1785-1803 (eval = FALSE)
###################################################
## practical.plan <- st.cross.plan(
##     merit = cross.merit[practical.ids],        # Merits of the six eligible parents
##     K = practical.K,                           # Original-scale kinship submatrix
##     n.crosses = 4,                             # Four distinct biological families
##     n.progeny = 20,                            # Equal produced count at this stage
##     theta = max(diag(practical.K)),            # Nonrestrictive relatedness ceiling
##     scores = practical.scores,                 # Scores for the six-parent panel
##     criterion = "usefulness",                  # Expected selected-family performance
##     constraints = list(                        # Intersecting practical restrictions
##         parents = practical.parents,           # Pool labels and role eligibility
##         limits = practical.limits,             # Exact used-parent counts
##         equal.use = practical.equal.use        # Equal use among selected seed parents
##     ),
##     control = list(                            # Search settings
##         method = "exact"                       # Bounded exhaustive cross search
##     )
## )
## practical.plan


###################################################
### code chunk number 56: v-cross.Rnw:1844-1867 (eval = FALSE)
###################################################
## practical.capacity <- data.frame(
##     parent = practical.ids[4:6],               # Individual pollen parents
##     pool = NA_character_,                      # Individual selectors need no pool filter
##     role = "pollen",                           # Pollen contribution only
##     basis = "produced",                        # Produced descendant slots
##     measure = "uses",                          # Integer incidences before division
##     lower = 0,                                 # No minimum produced pollen use
##     upper = 30                                 # At most thirty slots per pollen parent
## )
## practical.allocation <- st.cross.allocate(
##     plan = practical.plan,                     # Fixed families and declared directions
##     scores = practical.scores,                 # Cross means and segregation variances
##     total.progeny = 80,                        # Total produced descendant budget
##     n.selected = 2,                            # Retained descendants per family
##     min.progeny = 2,                           # Minimum produced count per family
##     max.progeny = 40,                          # Maximum produced count per family
##     accuracy = 1,                              # Within-family ranking accuracy
##     selection.method = "order.statistics",     # Finite-normal selection expectation
##     constraints = list(                        # Additional produced-capacity restriction
##         limits = practical.capacity            # Individual pollen-parent capacities
##     )
## )
## practical.allocation


###################################################
### code chunk number 57: v-cross.Rnw:1908-1912 (eval = FALSE)
###################################################
## practical.population <- st.cross.evaluate(
##     plan = practical.allocation                # Independent inherited-constraint check
## )
## practical.population


###################################################
### code chunk number 58: v-cross.Rnw:1979-2005 (eval = FALSE)
###################################################
## joint.equal.time <- system.time({
##     joint.equal <- st.cross.plan(
##         merit = cross.merit[practical.ids],    # Merits of the same six founders
##         K = practical.K,                       # Original-scale named kinship reference
##         n.crosses = 3,                         # Three distinct biological families
##         theta = max(diag(practical.K)),        # Nonrestrictive family-slot ceiling
##         scores = practical.scores,             # Means and variances for all fifteen pairs
##         criterion = "usefulness",              # Expected selected G1 mean
##         accuracy = 1,                          # Within-family ranking accuracy
##         progeny = list(                        # Explicit common production contract
##             total.progeny = 90,                # Total produced DH descendants
##             n.selected = 4,                    # Retained descendants per active family
##             min.progeny = 30,                  # Equal families in this baseline
##             max.progeny = 30,                  # Exactly thirty produced per family
##             selection.method =                 # Finite-normal selection expectation
##                 "order.statistics",
##             control = list(                    # Shared allocation settings
##                 algorithm = "auto"             # Qualified automatic integer allocator
##             )
##         ),
##         control = list(                        # Family-set search settings
##             method = "exact"                   # Bounded exhaustive family-set search
##         )
##     )
## })
## joint.equal


###################################################
### code chunk number 59: v-cross.Rnw:2042-2058 (eval = FALSE)
###################################################
## joint.fixed.time <- system.time({
##     joint.fixed <- st.cross.allocate(
##         plan = joint.equal,                    # The three baseline families remain fixed
##         scores = practical.scores,             # Identical means and segregation variances
##         total.progeny = 90,                    # Identical total production budget
##         n.selected = 4,                        # Identical retention per family
##         min.progeny = 4,                       # Every retained line has a produced line
##         max.progeny = 50,                      # At most fifty produced per family
##         accuracy = 1,                          # Identical ranking accuracy
##         selection.method = "order.statistics", # Identical finite-normal scoring
##         control = list(                        # Fixed-family allocation settings
##             algorithm = "auto"                 # Qualified automatic integer allocator
##         )
##     )
## })
## joint.fixed


###################################################
### code chunk number 60: v-cross.Rnw:2088-2115 (eval = FALSE)
###################################################
## joint.plan.time <- system.time({
##     joint.plan <- st.cross.plan(
##         merit = cross.merit[practical.ids],    # The same six eligible founders
##         K = practical.K,                       # The same kinship reference and scale
##         n.crosses = 3,                         # Identical number of active families
##         theta = max(diag(practical.K)),        # Identical family-slot relatedness ceiling
##         scores = practical.scores,             # The same fifteen candidate families
##         criterion = "usefulness",              # Expected selected G1 mean
##         accuracy = 1,                          # Identical ranking accuracy
##         progeny = list(                        # Joint pair and production decisions
##             total.progeny = 90,                # Identical total production budget
##             n.selected = 4,                    # Identical retention per active family
##             min.progeny = 4,                   # Lower produced count per family
##             max.progeny = 50,                  # Upper produced count per family
##             selection.method =                 # Identical finite-normal scoring
##                 "order.statistics",
##             control = list(                    # Shared allocation settings
##                 algorithm = "auto"             # Qualified automatic integer allocator
##             )
##         ),
##         control = list(                        # Family-set search settings
##             method = "exact",                  # Bounded exhaustive family-set search
##             start = joint.fixed                # Checked feasible allocation as incumbent
##         )
##     )
## })
## joint.plan


###################################################
### code chunk number 61: v-cross.Rnw:2153-2175 (eval = FALSE)
###################################################
## joint.equal.eval <- st.cross.evaluate(
##     plan = joint.equal                         # Recomputed equal-family population
## )
## joint.fixed.eval <- st.cross.evaluate(
##     plan = joint.fixed                         # Recomputed fixed-pair population
## )
## joint.plan.eval <- st.cross.evaluate(
##     plan = joint.plan                          # Recomputed joint population
## )
## joint.metrics <- c("retained.baseline", "response", "retained.mean")
## joint.comparison <- rbind(
##     Equal = unlist(joint.equal.eval$population[joint.metrics]),
##     Fixed = unlist(joint.fixed.eval$population[joint.metrics]),
##     Joint = unlist(joint.plan.eval$population[joint.metrics])
## )
## colnames(joint.comparison) <- c("Baseline", "Response", "Selected")
## joint.comparison <- cbind(
##     joint.comparison,
##     Seconds = c(joint.equal.time["elapsed"], joint.fixed.time["elapsed"],
##                 joint.plan.time["elapsed"])
## )
## joint.comparison


###################################################
### code chunk number 62: v-cross.Rnw:2230-2264 (eval = FALSE)
###################################################
## retention.time <- system.time({
##     retention.plan <- st.cross.plan(
##         merit = cross.merit[practical.ids],    # The same six eligible maize founders
##         K = practical.K,                       # Original-scale named kinship reference
##         n.crosses = 3,                         # Three distinct produced families
##         theta = max(diag(practical.K)),        # The same family-slot relatedness ceiling
##         upper = 0.5,                           # The same family-slot parental bound
##         scores = practical.scores,             # The same fifteen cross means and variances
##         criterion = "usefulness",              # Retained-individual G1 mean
##         accuracy = 1,                          # Identical within-family ranking accuracy
##         progeny = list(                        # Fixed production and retention totals
##             total.progeny = 90,                # Ninety produced DH descendants
##             total.selected = 12,               # Twelve retained DH descendants
##             min.progeny = 4,                   # Lower produced count for each active pair
##             max.progeny = 50,                  # Upper produced count for each active pair
##             min.selected = 1,                  # At least one retained from each family
##             max.selected = 12,                 # Retained upper bound before r <= n
##             selection.method =                 # Identical finite-normal scoring
##                 "order.statistics",
##             control = list(                    # Cumulative allocation settings
##                 algorithm = "auto",            # Supported bounded integer allocation
##                 max.transitions = 1e7          # Shared work ceiling for the whole call
##             )
##         ),
##         control = list(                        # Reproducible bounded family search
##             method = "heuristic",              # Best-found search over family choices
##             start = joint.plan,                # Actual produced and retained counts
##             seed = 1,                          # Local search seed
##             starts = 1,                        # One start from the checked incumbent
##             maxit = 2000                       # Ceiling on outer search iterations
##         )
##     )
## })
## retention.plan


###################################################
### code chunk number 63: v-cross.Rnw:2313-2317 (eval = FALSE)
###################################################
## retention.population <- st.cross.evaluate(
##     plan = retention.plan                      # Independent actual-count population check
## )
## retention.population


###################################################
### code chunk number 64: v-cross.Rnw:2353-2371 (eval = FALSE)
###################################################
## retention.comparison <- rbind(
##     Fixed = c(
##         unlist(joint.plan.eval$population[joint.metrics]),
##         joint.plan.eval$relatedness$family.slots,
##         joint.plan.eval$relatedness$retained,
##         unname(joint.plan.time["elapsed"])
##     ),
##     Variable = c(
##         unlist(retention.population$population[joint.metrics]),
##         retention.population$relatedness$family.slots,
##         retention.population$relatedness$retained,
##         unname(retention.time["elapsed"])
##     )
## )
## colnames(retention.comparison) <- c(
##     "Baseline", "Response", "Selected", "Family.q", "Retained.q", "Seconds"
## )
## retention.comparison


###################################################
### code chunk number 65: v-cross.Rnw:2432-2484 (eval = FALSE)
###################################################
## genetic.time <- system.time({
##     genetic.plan <- st.cross.plan(
##         merit = cross.merit[practical.ids],    # The same six candidate founders
##         K = practical.K,                       # The original named kinship reference
##         n.crosses = 3,                         # Three produced biological families
##         theta = max(diag(practical.K)),        # The original family-slot ceiling
##         upper = 0.5,                           # The original family-slot parent bound
##         scores = practical.scores,             # The same precise cross means and variances
##         criterion = "usefulness",              # Expected selected G1 mean
##         accuracy = 1,                          # Ranking by the fitted genetic value
##         progeny = list(                        # The same production and retention budgets
##             total.progeny = 90,                # Produced DH descendants
##             total.selected = 12,               # Retained DH descendants
##             min.progeny = 4,                   # Minimum produced from each active pair
##             max.progeny = 50,                  # Maximum produced from each active pair
##             min.selected = 1,                  # At least one retained from each family
##             max.selected = 12,                 # Retention also remains below production
##             selection.method =                 # Finite-normal selection prediction
##                 "order.statistics",
##             control = list(                    # Cumulative count-search settings
##                 algorithm = "auto",            # Supported bounded integer search
##                 max.transitions = 1e7          # Shared count-work ceiling
##             )
##         ),
##         genetic.model = list(                  # Optional additive genotype model
##             data.set = "cross.parents",        # Thirty-row fixed genomic reference
##             population = "DH",                 # Immediate F1-derived DH descendants
##             t = 0,                             # No additional intermating generation
##             map.function = "Haldane",          # cM map with independent crossovers
##             index.model =                      # Declared independent-noise index law
##                 "trait.independent.noise",
##             marker.weights = NULL,             # Uniform weights on all model markers
##             selected.founder.theta = NULL,     # Contributions reported without a new cap
##             scenarios = list(                  # Final genotype assessment only
##                 replicates = 250,              # Independent retained-population replicates
##                 seed = 12052026                # Reproducible private scenario bank
##             ),
##             control = list(                    # Whole-call genetic resource settings
##                 max.cells = 1e7,               # Genetic workspace and snapshot budget
##                 max.work = 1e8                 # Cumulative genetic-work ceiling
##             )
##         ),
##         control = list(                        # The same reproducible family search
##             method = "heuristic",              # Bounded best-found family search
##             start = joint.plan,                # The same actual starting counts
##             seed = 1,                          # The same local optimizer seed
##             starts = 1,                        # One checked starting plan
##             maxit = 2000                       # The same outer iteration ceiling
##         )
##     )
## })
## genetic.plan


###################################################
### code chunk number 66: v-cross.Rnw:2527-2531 (eval = FALSE)
###################################################
## genetic.population <- st.cross.evaluate(
##     plan = genetic.plan                        # Frozen core inherited; no scenario request
## )
## genetic.population


###################################################
### code chunk number 67: v-cross.Rnw:2561-2565 (eval = FALSE)
###################################################
## genetic.actual <- genetic.plan$scenarios$population[
##     , c("metric", "estimate", "SE")
## ]
## genetic.actual


###################################################
### code chunk number 68: v-cross.Rnw:2598-2616 (eval = FALSE)
###################################################
## origin.prediction <- genetic.plan$selected.contributions$families[
##     , c("P1Name", "P2Name", "P1")
## ]
## origin.scenario <- genetic.plan$scenarios$families[
##     genetic.plan$scenarios$families$metric == "P1.origin",
##     c("P1Name", "P2Name", "estimate", "SE")
## ]
## origin.comparison <- merge(
##     origin.prediction,
##     origin.scenario,
##     by = c("P1Name", "P2Name")
## )
## names(origin.comparison) <- c(
##     "Parent1", "Parent2", "Predicted", "Simulated", "SE"
## )
## origin.comparison$Difference <- origin.comparison$Simulated -
##                                 origin.comparison$Predicted
## origin.comparison


###################################################
### code chunk number 69: v-cross.Rnw:2642-2705 (eval = FALSE)
###################################################
## genetic.fixed <- st.cross.allocate(
##     plan = genetic.plan,                       # Inherited frozen model for these pairs
##     scores = practical.scores,                 # The same precise family moments
##     total.progeny = 90,                        # The same total production budget
##     n.selected = 4,                            # A fixed-retention feasible comparator
##     min.progeny = 4,                           # The same produced lower bound
##     max.progeny = 50,                          # The same produced upper bound
##     accuracy = 1,                              # The same ranking accuracy
##     selection.method = "order.statistics"      # The same finite-normal prediction
## )
## genetic.cap <- 0.33
## genetic.feasible.start <- genetic.fixed$crosses
## 
## genetic.cap.time <- system.time({
##     genetic.capped.plan <- st.cross.plan(
##         merit = cross.merit[practical.ids],    # The same six candidate founders
##         K = practical.K,                       # The original named kinship reference
##         n.crosses = 3,                         # Three produced biological families
##         theta = max(diag(practical.K)),        # The original family-slot ceiling
##         upper = 0.5,                           # The original family-slot parent bound
##         scores = practical.scores,             # The same precise cross means and variances
##         criterion = "usefulness",              # Expected selected G1 mean
##         accuracy = 1,                          # Ranking by the fitted genetic value
##         progeny = list(                        # The same production and retention budgets
##             total.progeny = 90,                # Produced DH descendants
##             total.selected = 12,               # Retained DH descendants
##             min.progeny = 4,                   # Minimum produced from each active pair
##             max.progeny = 50,                  # Maximum produced from each active pair
##             min.selected = 1,                  # At least one retained from each family
##             max.selected = 12,                 # Retention also remains below production
##             selection.method =                 # Finite-normal selection prediction
##                 "order.statistics",
##             control = list(                    # Cumulative count-search settings
##                 algorithm = "auto",            # Supported bounded integer search
##                 max.transitions = 1e7          # Shared count-work ceiling
##             )
##         ),
##         genetic.model = list(                  # Optional additive genotype model
##             data.set = "cross.parents",        # Thirty-row fixed genomic reference
##             population = "DH",                 # Immediate F1-derived DH descendants
##             t = 0,                             # No additional intermating generation
##             map.function = "Haldane",          # cM map with independent crossovers
##             index.model =                      # Declared independent-noise index law
##                 "trait.independent.noise",
##             marker.weights = NULL,             # Uniform weights on all model markers
##             selected.founder.theta =           # Verified real-data proxy ceiling
##                 genetic.cap,
##             scenarios = NULL,                  # No genotype sampling during this call
##             control = list(                    # Whole-call genetic resource settings
##                 max.cells = 1e7,               # Genetic workspace and snapshot budget
##                 max.work = 1e8                 # Cumulative genetic-work ceiling
##             )
##         ),
##         control = list(                        # The same reproducible family search
##             method = "heuristic",              # Bounded best-found family search
##             start = genetic.feasible.start,    # Independently checked feasible cap incumbent
##             seed = 1,                          # The same local optimizer seed
##             starts = 1,                        # One checked starting plan
##             maxit = 2000                       # The same outer iteration ceiling
##         )
##     )
## })
## genetic.capped.plan


###################################################
### code chunk number 70: v-cross.Rnw:2776-2855 (eval = FALSE)
###################################################
## future.equal <- st.cross.allocate(
##     plan = joint.equal,                        # Actual equal-production baseline pairs
##     scores = practical.scores,                 # The same precise cross means and variances
##     total.progeny = 90,                        # Identical total produced descendants
##     n.selected = 4,                            # Explicit retained count in each family
##     min.progeny = 30,                          # Equal production in all three families
##     max.progeny = 30,                          # Exactly thirty candidates per family
##     accuracy = 1,                              # The same current ranking accuracy
##     selection.method = "order.statistics"      # The same analytic G1 selection method
## )
## future.progeny <- list(
##     total.progeny = 90,                        # Current produced DH descendants
##     total.selected = 12,                       # Current retained DH descendants
##     min.progeny = 4,                           # Minimum produced per active family
##     max.progeny = 50,                          # Maximum produced per active family
##     min.selected = 1,                          # At least one current retained line
##     max.selected = 12,                         # Retention also stays below production
##     selection.method = "order.statistics",     # Analytic G1 selection expectation
##     control = list(                            # Current integer-count settings
##         algorithm = "auto",                    # Supported bounded count search
##         max.transitions = 1e7                  # Cumulative count-work ceiling
##     )
## )
## future.policy <- list(
##     policy = "random.disjoint",                # Every retained G1 line is mated once
##     n.progeny = 20,                            # Produced DHs from each new G2 F1
##     n.selected = 2,                            # Retained DHs from each G2 family
##     accuracy = 1,                              # G2 ranking by fitted genetic value
##     training = list(                           # Fixed bank used by future search
##         replicates = 32,                       # Complete replicates per candidate
##         seed = 0                               # Reproducible private training seed
##     ),
##     evaluation = list(                         # Separate held-out assessment bank
##         replicates = 256,                      # Requested whole common replicates
##         seed = 1                               # Role-separated evaluation seed
##     ),
##     comparators = list(                        # Actual equal-resource designs
##         Equal = future.equal,                  # Evaluated actual n=30/r=4 reference
##         G1 = genetic.plan                      # Existing G1-focused actual plan
##     )
## )
## future.model <- list(
##     data.set = "cross.parents",                # Original thirty-row genomic reference
##     population = "DH",                         # Immediate F1-derived DH candidates
##     t = 0,                                     # No extra intermating generation
##     map.function = "Haldane",                  # cM map and independent crossovers
##     index.model =                              # Independent-noise ranking law
##         "trait.independent.noise",
##     marker.weights = NULL,                     # Uniform weights on all 243 markers
##     selected.founder.theta = NULL,             # No additional G1 founder-proxy cap
##     scenarios = NULL,                          # No separate single-cycle scenario bank
##     future = future.policy,                    # Actual two-cycle request declared above
##     control = list(                            # Whole-call genetic resource settings
##         max.cells = 1e7,                       # Shared genetic workspace ceiling
##         max.work = 2e8                         # Search and final assessment together
##     )
## )
## future.time <- system.time({
##     future.plan <- st.cross.plan(
##         merit = cross.merit[practical.ids],    # The same six candidate founders
##         K = practical.K,                       # The original named kinship reference
##         n.crosses = 3,                         # Current produced biological families
##         theta = max(diag(practical.K)),        # The original family-slot ceiling
##         upper = 0.5,                           # The original parental family-slot bound
##         scores = practical.scores,             # Precise scores for all fifteen pairs
##         criterion = "future",                  # Actual selected G2 mean on training bank
##         accuracy = 1,                          # Current within-family ranking accuracy
##         progeny = future.progeny,              # Current production and retention budgets
##         genetic.model = future.model,          # Explicit core and next-cycle request
##         control = list(                        # Bounded joint family and count search
##             method = "heuristic",              # Best-found finite-bank search
##             start = genetic.plan,              # Actual checked G1 pairs and n/r counts
##             seed = 1,                          # Reproducible search seed
##             starts = 1,                        # Requested starting schedule
##             maxit = 2000                       # Current outer iteration ceiling
##         )
##     )
## })
## future.plan


###################################################
### code chunk number 71: v-cross.Rnw:2913-2919 (eval = FALSE)
###################################################
## future.comparison <- future.plan$future$population[
##     future.plan$future$population$design %in% c("Equal", "G1", "selected") &
##         future.plan$future$population$metric == "selected.mean",
##     c("design", "generation", "estimate", "SE")
## ]
## future.comparison


###################################################
### code chunk number 72: v-cross.Rnw:2943-2948 (eval = FALSE)
###################################################
## future.differences <- future.plan$future$contrasts[
##     future.plan$future$contrasts$comparator %in% c("Equal", "G1"),
##     c("comparator", "metric", "estimate", "SE")
## ]
## future.differences


###################################################
### code chunk number 73: v-cross.Rnw:2986-2991 (eval = FALSE)
###################################################
## future.pair.moments <- future.plan$future$pair.moments[
##     future.plan$future$pair.moments$design %in% c("Equal", "G1", "selected"),
##     c("design", "metric", "estimate", "SE")
## ]
## future.pair.moments


###################################################
### code chunk number 74: v-cross.Rnw:3016-3023 (eval = FALSE)
###################################################
## future.genomics <- future.plan$future$population[
##     future.plan$future$population$design %in% c("Equal", "G1", "selected") &
##         future.plan$future$population$generation == "G2" &
##         future.plan$future$population$metric %in% c("genomic.q", "gene.diversity"),
##     c("design", "metric", "estimate", "SE")
## ]
## future.genomics


###################################################
### code chunk number 75: v-cross.Rnw:3054-3081 (eval = FALSE)
###################################################
## future.theta.limits <- c(0.25, 0.30, 0.35)
## future.frontier.model <- future.model
## future.frontier.model$future$comparators <- NULL
## future.frontier.time <- system.time({
##     future.frontier <- st.cross.frontier(
##         merit = cross.merit[practical.ids],    # The same six candidate founders
##         K = practical.K,                       # Original-scale fixed kinship reference
##         n.crosses = 3,                         # Three current produced families
##         theta = max(diag(practical.K)),        # Replaced by each declared theta limit
##         upper = 0.5,                           # The same parental family-slot bound
##         scores = practical.scores,             # The same fifteen actual candidate pairs
##         criterion = "future",                  # Actual selected G2 training objective
##         accuracy = 1,                          # The same current ranking accuracy
##         progeny = future.progeny,              # The same current count domain and totals
##         genetic.model = future.frontier.model, # Same core and banks; no extra comparators
##         limits = future.theta.limits,          # Ordered G1 family-slot theta ceilings
##         metric = "theta",                      # Axis concerns nominal G1 family slots
##         warmstart = FALSE,                     # Each point keeps its standalone start
##         control = list(                        # The same declared search schedule
##             method = "heuristic",              # Bounded finite-bank search at each cap
##             seed = 1,                          # The same recorded search seed
##             starts = 1,                        # The same requested starting schedule
##             maxit = 2000                       # Per-point outer iteration ceiling
##         )
##     )
## })
## future.frontier


###################################################
### code chunk number 76: v-cross.Rnw:4433-4493 (eval = FALSE)
###################################################
## input <- readRDS("crossing-example-inputs.rds")
## cross.example.prepare <- function(input) {
##     X <- input$markers
##     stopifnot(is.matrix(X), is.numeric(X), nrow(X) >= 10,
##               ncol(X) >= 100)
##     X <- X[seq_len(10), seq_len(100), drop = FALSE]
##     stopifnot(!anyNA(X), all(X %in% c(0, 2)),
##               !is.null(rownames(X)), !is.null(colnames(X)),
##               !anyDuplicated(rownames(X)), !anyDuplicated(colnames(X)))
##     loci <- colnames(X)
##     effects <- function(x) {
##         stopifnot(is.numeric(x), !is.null(names(x)),
##                   !anyDuplicated(names(x)), all(loci %in% names(x)))
##         out <- unname(x[loci])
##         stopifnot(all(is.finite(out)))
##         out
##     }
##     a.yield <- effects(input$additive$YLD)
##     stopifnot(is.numeric(input$weights),
##               !anyDuplicated(names(input$weights)),
##               all(c("YLD", "MAT") %in% names(input$weights)))
##     weights <- input$weights[c("YLD", "MAT")]
##     stopifnot(all(is.finite(weights)))
##     a.index <- weights[1] * a.yield +
##                weights[2] * effects(input$additive$MAT)
##     d.index <- weights[1] * effects(input$dominance$YLD) +
##                weights[2] * effects(input$dominance$MAT)
##     centered <- sweep(X, 2, colMeans(X), "-")
##     p <- colMeans(X) / 2
##     scale <- 2 * sum(p * (1 - p))
##     stopifnot(is.finite(scale), scale > 0)
##     G <- tcrossprod(centered) / scale
##     dimnames(G) <- list(rownames(X), rownames(X))
##     pairs <- t(combn(rownames(X), 2))
##     h <- 2^(-3)
##     intensity <- dnorm(qnorm(1 - 5 / 10)) / (5 / 10)
##     score <- function(a, d = rep(0, ncol(X))) {
##         mu <- va <- numeric(nrow(pairs))
##         for (j in seq_len(nrow(pairs))) {
##             x1 <- X[pairs[j, 1], ]
##             x2 <- X[pairs[j, 2], ]
##             segregating <- x1 != x2
##             mu[j] <- sum((x1 + x2) * a / 2) +
##                      h * sum(d[segregating])
##             va[j] <- sum(((1 - h) * a^2 +
##                           h * (1 - h) * d^2)[segregating])
##         }
##         data.frame(P1Name = pairs[, 1], P2Name = pairs[, 2],
##                    mu = mu, va = va,
##                    es = mu + sqrt(va) * intensity)
##     }
##     list(markers = X, G = G, K = G / 2,
##          effects = list(yield = a.yield, index = a.index,
##                         dominance = a.index),
##          scores = list(yield = score(a.yield), index = score(a.index),
##                        dominance = score(a.index, d.index)),
##          accuracy = 0.5, loss.limit = 0.05,
##          n.crosses = 5L, n.progeny = 10L, n.selected = 5L)
## }
## example <- cross.example.prepare(input)


###################################################
### code chunk number 77: v-cross.Rnw:4515-4532 (eval = FALSE)
###################################################
## cross.example.plan <- function(example, scenario) {
##     merit <- setNames(drop(example$markers %*% example$effects[[scenario]]),
##                       rownames(example$markers))
##     SelectionTools::st.cross.plan(
##         merit = merit, K = example$K,
##         n.crosses = example$n.crosses, n.progeny = example$n.progeny,
##         theta = max(diag(example$K)), upper = 0.5,
##         scores = example$scores[[scenario]],
##         criterion = "usefulness", accuracy = example$accuracy,
##         diversity = list(G = example$G, max = example$loss.limit,
##                          sample = "family.slots"),
##         control = list(method = "exact", max.nodes = 3000000,
##                        starts = 4, maxit = 2500, seed = 1L)
##     )
## }
## yield.plan <- cross.example.plan(example, "yield")
## yield.plan


###################################################
### code chunk number 78: v-cross.Rnw:4557-4559 (eval = FALSE)
###################################################
## index.plan <- cross.example.plan(example, "index")
## index.plan


###################################################
### code chunk number 79: v-cross.Rnw:4584-4586 (eval = FALSE)
###################################################
## dominance.plan <- cross.example.plan(example, "dominance")
## dominance.plan


