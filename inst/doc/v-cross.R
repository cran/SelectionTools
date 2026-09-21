### R code from vignette source 'v-cross.Rnw'

###################################################
### code chunk number 1: v-cross.Rnw:196-203
###################################################
library("SelectionTools"); st.set.info.level(-2) # Only error messages

data("v-tropmaize-vcf")
data("v-tropmaize-phe")

st.load.vcf.data(v.tropmaize.vcf)
st.load.performance.data(v.tropmaize.phe)


###################################################
### code chunk number 2: v-cross.Rnw:206-210
###################################################
st.restrict.marker.data( NoAll.MAX = 2    )
st.restrict.marker.data( MaMis.MAX = 0.05 )
st.restrict.marker.data( ExHet.MIN = 0.1  )
st.restrict.marker.data( InMis.MAX = 0.1  )


###################################################
### code chunk number 3: v-cross.Rnw:213-214
###################################################
st.load.performance.data(v.tropmaize.phe)


###################################################
### code chunk number 4: v-cross.Rnw:217-218
###################################################
st.copy.marker.data("preprocessed","default")


###################################################
### code chunk number 5: v-cross.Rnw:226-230
###################################################
st.copy.marker.data( "hblocks01", "preprocessed" )
h <- st.def.hblocks( hap      = 5,     
                     hap.unit = 1,  
                     data.set = "hblocks01" )


###################################################
### code chunk number 6: v-cross.Rnw:234-235
###################################################
st.recode.hil ( data.set="hblocks01" )      


###################################################
### code chunk number 7: v-cross.Rnw:238-241 (eval = FALSE)
###################################################
## st.plot.ggt( f.ind    = 1,
##              l.ind    = 20,
##              data.set = "hblocks01")    


###################################################
### code chunk number 8: v-cross.Rnw:252-253
###################################################
h[1:10,]     


###################################################
### code chunk number 9: v-cross.Rnw:258-265 (eval = FALSE)
###################################################
## st.copy.marker.data( "hblocks02", "preprocessed" )
## 
## h2 <- st.def.hblocks( hap      = 5,     
##                       hap.unit = 2,  
##                       data.set = "hblocks02" )
## 
## h2[1:10,]   


###################################################
### code chunk number 10: v-cross.Rnw:303-309 (eval = FALSE)
###################################################
## st.plot.ggt(f.ind                  = 1,
##             l.ind                  = 30,
##             sort.performance       = TRUE,
##             show.performance       = TRUE,
##             performance.decreasing = TRUE,
##             data.set               = "hblocks01")


###################################################
### code chunk number 11: v-cross.Rnw:322-323
###################################################
st.copy.marker.data( "hblocks03", "hblocks01" )


###################################################
### code chunk number 12: v-cross.Rnw:327-329
###################################################
st.recode.ref( reference = "142",
               data.set  = "hblocks03")


###################################################
### code chunk number 13: v-cross.Rnw:333-339 (eval = FALSE)
###################################################
## st.plot.ggt( f.ind                  = 1,
##              l.ind                  = 30,
##              sort.performance       = TRUE,
##              show.performance       = TRUE,
##              performance.decreasing = TRUE,
##              data.set               = "hblocks03")


###################################################
### code chunk number 14: v-cross.Rnw:349-350 (eval = FALSE)
###################################################
## st.copy.marker.data( "gd01", "preprocessed" )


###################################################
### code chunk number 15: v-cross.Rnw:354-355 (eval = FALSE)
###################################################
## best <- st.select.phen( n=50, data.set="gd01" )


###################################################
### code chunk number 16: v-cross.Rnw:358-361 (eval = FALSE)
###################################################
## 
## st.restrict.marker.data ( ind.list=best$i,
##                           data.set="gd01" )


###################################################
### code chunk number 17: v-cross.Rnw:370-376 (eval = FALSE)
###################################################
## dist.mat <- st.genetic.distances( measure  = "rd",
##                                   format   = "m",
##                                   data.set = "gd01" )
## 
## dm <- as.matrix(dist.mat)
## heatmap( dm, scale="none", col=rev(heat.colors(12) ) )


###################################################
### code chunk number 18: v-cross.Rnw:387-390 (eval = FALSE)
###################################################
## gs.cross.eval.gd ( dist="mrd", data.set="gd01" )
## crosses <- gs.cross.info ( sortby="gd", bestn=50, data.set="gd01" )
## crosses[1:10,]


###################################################
### code chunk number 19: v-cross.Rnw:409-411 (eval = FALSE)
###################################################
## st.def.hblocks( hap=5, data.set="gd01" )
## st.recode.hil ( data.set="gd01" )


###################################################
### code chunk number 20: v-cross.Rnw:414-426 (eval = FALSE)
###################################################
## st.set.info.level(-2)
## for (i in seq_len(min(10,nrow(crosses)))) {
##     a <- crosses$P1Name[i]
##     b <- crosses$P2Name[i]
##     f.nme <- sprintf("ex-crs-07-%03i.pdf",100+i)
##     pdf(f.nme,width=8,height=0.8,pointsize=10)
##     st.plot.ggt(i.list=c(a,b),
##                 show.performance=TRUE,
##                 data.set="gd01")
##     dev.off()
## }
## st.set.info.level(0)


###################################################
### code chunk number 21: v-cross.Rnw:452-453 (eval = FALSE)
###################################################
## st.copy.marker.data( "effects01", "preprocessed" )


###################################################
### code chunk number 22: v-cross.Rnw:459-460 (eval = FALSE)
###################################################
## gs.esteff.rr( method="BLUP", data.set="effects01" )


###################################################
### code chunk number 23: v-cross.Rnw:471-474 (eval = FALSE)
###################################################
## gs.cross.eval.mu( data.set="effects01" )
## gs.cross.eval.mi( data.set="effects01" )
## gs.cross.eval.ma( data.set="effects01" )


###################################################
### code chunk number 24: v-cross.Rnw:482-483 (eval = FALSE)
###################################################
## gs.cross.eval.va( pop.type="DH", data.set="effects01" )


###################################################
### code chunk number 25: v-cross.Rnw:492-493 (eval = FALSE)
###################################################
## gs.cross.eval.es( alpha=0.25, data.set="effects01" )


###################################################
### code chunk number 26: v-cross.Rnw:498-499 (eval = FALSE)
###################################################
## gs.cross.eval.es( N=1, G=4, data.set="effects01" )


###################################################
### code chunk number 27: v-cross.Rnw:506-507 (eval = FALSE)
###################################################
## gs.cross.eval.gd(dist="rd",data.set="effects01")


###################################################
### code chunk number 28: v-cross.Rnw:514-522 (eval = FALSE)
###################################################
## gs.cross.eval.gd(dist="rd",data.set="effects01")
## gs.cross.eval.mu(data.set="effects01")
## gs.cross.eval.ma(data.set="effects01")
## gs.cross.eval.mi(data.set="effects01")
## gs.cross.eval.va(pop.type="DH",data.set="effects01")
## gs.cross.eval.es(N=1,G=4,data.set="effects01")
## crosses <- gs.cross.info(sortby="es",bestn=20,data.set="effects01")
## format(crosses,digits=2,nsmall=2)


