### R code from vignette source 'v-phase.Rnw'

###################################################
### code chunk number 1: startup (eval = FALSE)
###################################################
## library(SelectionTools)
## st.set.num.threads(4)


###################################################
### code chunk number 2: f2-load (eval = FALSE)
###################################################
## data("v-tropmaize-vcf")
## st.load.vcf.data(v.tropmaize.vcf,data.set="tropical")
## 
## st.copy.marker.data( "P1","tropical" )
## st.restrict.marker.data( ind.list="1", data.set="P1")
## st.restrict.marker.data( NoAll.MAX=1, MaMis.MAX=0, data.set="P1")
## 
## st.copy.marker.data("P2","tropical")
## st.restrict.marker.data(ind.list="2",data.set="P2")
## st.restrict.marker.data( NoAll.MAX=1, MaMis.MAX=0, data.set="P2")


###################################################
### code chunk number 3: f2-parents (eval = FALSE)
###################################################
## retained.markers <- intersect( st.get.map(data.set="P1")$Name,
##                                st.get.map(data.set="P2")$Name)
## 
## st.copy.marker.data("parents","tropical")
## st.restrict.marker.data( ind.list = c("1","2"),
##                          mar.list = retained.markers,
##                          data.set = "parents")
## st.recode.ref.2( r1=1,
##                  r2=2,
## 		 data.set="parents")


###################################################
### code chunk number 4: f2-sim (eval = FALSE)
###################################################
## st.set.simpop( pop.name="Parents", data.set="parents" )
## population.divide( "P1", "Parents", 1 )
## population.rename( "Parents", "P2" )
## cross( "F1", "P1", "P2", NoPg=1 )
## st.get.simpop( "F1", data.set="f1.reference")
## cross( "F2", "F1", "F1", NoPg=200 )
## st.get.simpop( "F2", data.set="f2.true" )


###################################################
### code chunk number 5: f2-unphase (eval = FALSE)
###################################################
## st.copy.marker.data( "f2", "f2.true")
## set.seed(314159)
## st.destroy.phase( data.set="f2")
## st.copy.marker.data( "f2.unphased", "f2" )
## st.copy.marker.data( "f2.hetero",   "f2" )


###################################################
### code chunk number 6: f2-phase (eval = FALSE)
###################################################
## st.phase( reference.data.set="parents",
##           data.set="f2")


###################################################
### code chunk number 7: f2-switch-error (eval = FALSE)
###################################################
## st.switch.error( true.data.set="f2.true",
##                  data.set="f2")


###################################################
### code chunk number 8: plot-data (eval = FALSE)
###################################################
## st.copy.marker.data( "informative", "parents")
## st.restrict.marker.data( ExHet.MIN=0.5,
##                          data.set="informative")
## informative.markers <- st.get.map(data.set="informative")$Name
## 
## st.copy.marker.data( "f2.plot.before", "f2.unphased")
## st.restrict.marker.data( mar.list=informative.markers,
##                          data.set="f2.plot.before")
## 
## st.copy.marker.data( "f2.plot.after", "f2" )
## st.restrict.marker.data( mar.list=informative.markers,
##                          data.set="f2.plot.after")
## 
## st.copy.marker.data("f2.plot.true","f2.true")
## st.restrict.marker.data( mar.list=informative.markers,
##                          data.set="f2.plot.true")


###################################################
### code chunk number 9: plot-three (eval = FALSE)
###################################################
## st.plot.ggt(data.set="f2.plot.before",f.ind=1,l.ind=20,
##             color=c("yellow","blue"),z.min=1,z.max=2)
## st.plot.ggt(data.set="f2.plot.after",f.ind=1,l.ind=20,
##             color=c("yellow","blue"),z.min=1,z.max=2)
## st.plot.ggt(data.set="f2.plot.true",f.ind=1,l.ind=20,
##             color=c("yellow","blue"),z.min=1,z.max=2)


###################################################
### code chunk number 10: hetero-phase (eval = FALSE)
###################################################
## st.phase( reference.data.set="f1.reference",
##           population.type="F2",
##           data.set="f2.hetero")


###################################################
### code chunk number 11: hetero-switch-error (eval = FALSE)
###################################################
## st.switch.error( true.data.set="f2.true",
##                  data.set="f2.hetero")


###################################################
### code chunk number 12: cohort-sim (eval = FALSE)
###################################################
## st.set.simpop( pop.name="Base", data.set="tropical" )
## m.s <- st.marker.data.statistics( data.set="tropical",
## 				  mar     = FALSE,
##                                   gen     = FALSE)$individual.list
## n.base <- nrow(m.s)
## cross( "popSYN1", "Base", "Base", NoPg=n.base )
## st.get.simpop( "popSYN1", data.set="popSYN1.true" )


###################################################
### code chunk number 13: cohort-unphase (eval = FALSE)
###################################################
## st.copy.marker.data( "popSYN1", "popSYN1.true")
## set.seed(271828)
## st.destroy.phase( data.set="popSYN1")


###################################################
### code chunk number 14: cohort-phase (eval = FALSE)
###################################################
## st.phase( impute   = FALSE,
##           data.set = "popSYN1")


###################################################
### code chunk number 15: cohort-switch-error (eval = FALSE)
###################################################
## st.switch.error( true.data.set = "popSYN1.true",
##                  data.set      = "popSYN1")


