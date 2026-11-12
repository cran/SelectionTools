### R code from vignette source 'v-enrich.Rnw'

###################################################
### code chunk number 1: setup-display (eval = FALSE)
###################################################
## library(SelectionTools)
## options(
##     width = 68,
##     digits = 5,
##     bitmapType = "cairo"
## )
## st.set.num.threads(1)
## st.set.info.level(-1)


###################################################
### code chunk number 2: setup-figures (eval = FALSE)
###################################################
## output.dir <- getOption("SelectionTools.enrichment.dir", getwd())
## figure.dir <- file.path(output.dir, "figures")
## dir.create(
##     figure.dir,
##     recursive = TRUE,
##     showWarnings = FALSE
## )


###################################################
### code chunk number 3: load (eval = FALSE)
###################################################
## data("v-tropmaize-vcf")
## st.load.vcf.data(
##     v.tropmaize.vcf,
##     data.set = "base"
## )
## st.copy.marker.data("clustering", "base")


###################################################
### code chunk number 4: filter-markers (eval = FALSE)
###################################################
## st.restrict.marker.data(
##     NoAll.MAX = 2,  # At most two alleles per marker
##     MaMis.MAX = 0.1, # Missing-call limit per marker
##     data.set = "clustering"
## )


###################################################
### code chunk number 5: filter-lines (eval = FALSE)
###################################################
## st.restrict.marker.data(
##     InMis.MAX = 0.1, # Missing-call limit per line
##     data.set = "clustering"
## )
## st.restrict.marker.data(
##     ExHet.MIN = 0.1, # Retain polymorphic, informative markers
##     data.set = "clustering"
## )


###################################################
### code chunk number 6: filter-summary (eval = FALSE)
###################################################
## filtered <- st.marker.data.statistics(
##     data.set = "clustering",
##     gen = FALSE # Omit the large genotype table
## )
## st.data.set.list()


###################################################
### code chunk number 7: cluster-distances (eval = FALSE)
###################################################
## dist.mat <- st.genetic.distances(
##     measure = "rd", # Rogers distance between lines
##     format = "m",   # R 'dist' object for cmdscale
##     data.set = "clustering"
## )


###################################################
### code chunk number 8: cluster-groups (eval = FALSE)
###################################################
## set.seed(1) # Reproduce the random starts and cluster labels
## clus <- kmeans(
##     as.matrix(dist.mat), # Rows are distance profiles, not SNPs
##     centers = 5,
##     nstart = 25          # Keep the best of 25 random starts
## )
## clusters <- clus$cluster


###################################################
### code chunk number 9: cluster-coordinates (eval = FALSE)
###################################################
## pco <- cmdscale( # Principal coordinates from Rogers distances
##     dist.mat,
##     k = 4,       # Retain four axes; plot the first two below
##     eig = TRUE   # Keep eigenvalues to assess the embedding
## )
## PC <- pco$points
## print(table(clusters))


###################################################
### code chunk number 10: pcoa-initial-device (eval = FALSE)
###################################################
## cluster.colours <- c("black", "#DF536B", "#61D04F",
##                      "#2297E6", "#28E2E5")
## png(
##     file.path(figure.dir, "enrichment-clusters-initial.png"),
##     width = 1000,
##     height = 750,
##     pointsize = 18,
##     res = 120
## )
## par(mar = c(4.5, 4.5, 2, 1))


###################################################
### code chunk number 11: pcoa-initial-plot (eval = FALSE)
###################################################
## plot(
##     PC[, 1],
##     PC[, 2],
##     col = cluster.colours[clusters], # Stable cluster colours
##     pch = 16,
##     xlab = "Principal coordinate 1",
##     ylab = "Principal coordinate 2"
## )


###################################################
### code chunk number 12: pcoa-initial-legend (eval = FALSE)
###################################################
## legend(
##     "topright",
##     pch = 16,
##     bty = "n",
##     col = cluster.colours,
##     legend = c("1: excluded", "2: excluded", "3: main",
##                 "4: pool A", "5: pool B")
## )
## invisible(dev.off())


###################################################
### code chunk number 13: pools (eval = FALSE)
###################################################
## main.ind <- names(clusters)[clusters == 3]
## a.pool <- names(clusters)[clusters == 4]
## b.pool <- names(clusters)[clusters == 5]
## pools <- list(A = sort(a.pool), B = sort(b.pool))
## joint.ind <- c(main.ind, a.pool, b.pool)


###################################################
### code chunk number 14: pool-sizes (eval = FALSE)
###################################################
## populations <- data.frame(
##     Population = c("Main", "A", "B"),
##     Cluster = c(3, 4, 5),
##     Lines = c(length(main.ind), length(a.pool), length(b.pool))
## )
## print(populations, row.names = FALSE)


###################################################
### code chunk number 15: copy-base (eval = FALSE)
###################################################
## st.copy.marker.data("base.phased", "clustering")
## st.restrict.marker.data(ind.list=main.ind, data.set="base.phased")


###################################################
### code chunk number 16: copy-a (eval = FALSE)
###################################################
## st.copy.marker.data("A.phased", "clustering")
## st.restrict.marker.data(ind.list=a.pool, data.set="A.phased")


###################################################
### code chunk number 17: copy-b (eval = FALSE)
###################################################
## st.copy.marker.data("B.phased", "clustering")
## st.restrict.marker.data(ind.list=b.pool, data.set="B.phased")


###################################################
### code chunk number 18: phase-base (eval = FALSE)
###################################################
## st.phase(
##     impute = TRUE, # Retain inferred calls at missing genotypes
##     phasing.options = list(seed = 20261003), # Reproducible calls
##     data.set = "base.phased"
## )


###################################################
### code chunk number 19: phase-a (eval = FALSE)
###################################################
## st.phase(
##     impute = TRUE, # Retain inferred calls at missing genotypes
##     phasing.options = list(seed = 20261004),
##     data.set = "A.phased"
## )


###################################################
### code chunk number 20: phase-b (eval = FALSE)
###################################################
## st.phase(
##     impute = TRUE, # Retain inferred calls at missing genotypes
##     phasing.options = list(seed = 20261005),
##     data.set = "B.phased"
## )


###################################################
### code chunk number 21: options (eval = FALSE)
###################################################
## haplotype.options <- list(
##     hap = 5,
##     hap.unit = 1 # Block size is measured in markers
## )


###################################################
### code chunk number 22: select-a (eval = FALSE)
###################################################
## # Shared donor variants count once in the combined gain.
## forward.A <- st.select.enrichment(
##     base.pool          = "base.phased",
##     enrichment.pool    = "A.phased",
##     n                  = 5,
##     haplotype.options  = haplotype.options,
##     exclude.compatible = TRUE,
##     method             = "forward"
## )
## print(forward.A)


###################################################
### code chunk number 23: restrict-a (eval = FALSE)
###################################################
## st.copy.marker.data("selected.A", "A.phased")
## st.restrict.marker.data(ind.list=forward.A, data.set="selected.A")


###################################################
### code chunk number 24: calculate-a (eval = FALSE)
###################################################
## gain.A <- st.calculate.enrichment(
##     base.pool = "base.phased",
##     enrichment.pool = "selected.A", # Only the five chosen donors
##     haplotype.options = haplotype.options,
##     exclude.compatible = TRUE # Exclude uncertain additions
## )
## print(gain.A$summary, row.names = FALSE)


###################################################
### code chunk number 25: select-b (eval = FALSE)
###################################################
## # Keep the same base and blocks for a comparable pool-B gain.
## forward.B <- st.select.enrichment(
##     base.pool          = "base.phased",
##     enrichment.pool    = "B.phased",
##     n                  = 5,
##     haplotype.options  = haplotype.options,
##     exclude.compatible = TRUE,
##     method             = "forward"
## )
## print(forward.B)


###################################################
### code chunk number 26: restrict-b (eval = FALSE)
###################################################
## st.copy.marker.data("selected.B", "B.phased")
## st.restrict.marker.data(ind.list=forward.B, data.set="selected.B")


###################################################
### code chunk number 27: calculate-b (eval = FALSE)
###################################################
## gain.B <- st.calculate.enrichment(
##     base.pool = "base.phased",
##     enrichment.pool = "selected.B", # Only the five chosen donors
##     haplotype.options = haplotype.options,
##     exclude.compatible = TRUE # Exclude uncertain additions
## )
## print(gain.B$summary, row.names = FALSE)


###################################################
### code chunk number 28: details (eval = FALSE)
###################################################
## print(head(gain.B$by.block), row.names = FALSE)
## print(head(gain.B$alleles), row.names = FALSE) # Added sequences


###################################################
### code chunk number 29: fixed-blocks (eval = FALSE)
###################################################
## blocks <- st.def.hblocks(
##     hap = 5,
##     hap.unit = 1, # Block size is measured in markers
##     data.set = "base.phased"
## )


###################################################
### code chunk number 30: fixed-block-display (eval = FALSE)
###################################################
## print(head(blocks, 5), row.names = FALSE)


###################################################
### code chunk number 31: fixed-select-a (eval = FALSE)
###################################################
## forward.fixed.A <- st.select.enrichment(
##     base.pool          = "base.phased",
##     enrichment.pool    = "A.phased",
##     n                  = 5,
##     blocks             = blocks,
##     exclude.compatible = TRUE,
##     method             = "forward"
## )
## print(forward.fixed.A)


###################################################
### code chunk number 32: plot-results (eval = FALSE)
###################################################
## a.ind <- forward.A
## b.ind <- forward.B
## result <- rbind(gain.A$summary, gain.B$summary)
## result$Pool <- c("A", "B")


###################################################
### code chunk number 33: pcoa-device (eval = FALSE)
###################################################
## png(
##     file.path(figure.dir, "enrichment-clusters.png"),
##     width = 1000,
##     height = 750,
##     pointsize = 18,
##     res = 120
## )
## par(mar = c(4.5, 4.5, 2, 1))
## plot(
##     PC[, 1],
##     PC[, 2],
##     col = cluster.colours[clusters], # Stable cluster colours
##     pch = 16,
##     xlab = "Principal coordinate 1",
##     ylab = "Principal coordinate 2"
## )


###################################################
### code chunk number 34: pcoa-labels (eval = FALSE)
###################################################
## points(
##     PC[c(a.ind, b.ind), 1],
##     PC[c(a.ind, b.ind), 2],
##     pch = 1,
##     cex = 1.6,
##     lwd = 1.5
## )
## legend(
##     "topright",
##     pch = 16,
##     bty = "n",
##     col = cluster.colours,
##     legend = c("1: excluded", "2: excluded", "3: main",
##                 "4: pool A", "5: pool B")
## )
## invisible(dev.off())


###################################################
### code chunk number 35: gain-plot (eval = FALSE)
###################################################
## png(
##     file.path(figure.dir, "enrichment-gain.png"),
##     width = 900,
##     height = 650,
##     pointsize = 18,
##     res = 120
## )
## barplot(
##     rbind(result$No.match, result$Hidden),
##     names.arg = paste("Main +", result$Pool),
##     ylim = c(0, 1.8 * max(result$Additional)),
##     col = c("cornflowerblue", "coral"),
##     ylab = "Additional block-specific alleles",
##     legend.text = c("No compatible missing call", "Potentially hidden"),
##     args.legend = list(x = "topright", bty = "n")
## )
## invisible(dev.off())


###################################################
### code chunk number 36: plot-pools (eval = FALSE)
###################################################
## pool.names <- c("base.phased", "A.phased", "B.phased")
## phased <- lapply(pool.names, function(x) {
##     st.marker.data.statistics(data.set = x)
## })


###################################################
### code chunk number 37: plot-input (eval = FALSE)
###################################################
## snp <- do.call(cbind, lapply(phased, function(x) {
##     as.matrix(x$genotypes[, -1, drop = FALSE])
## }))
## rownames(snp) <- phased[[1]]$genotypes[[1]]
## missing.codes <- c("-1", "--", "-", ".", "NA", "255")
## block.markers <- strsplit(gain.A$blocks$Markers, ";", fixed = TRUE)


###################################################
### code chunk number 38: split-homologues (eval = FALSE)
###################################################
## homologue <- function(x, which) {
##     parts <- strsplit(as.vector(x), "/", fixed = TRUE)
##     values <- vapply(parts, function(z) {
##         if (length(z) == 2L) z[which] else NA_character_
##     }, character(1))
##     # Never treat missing calls as an observed haplotype allele.
##     values[values %in% missing.codes] <- NA_character_
##     matrix(values, nrow(x), dimnames = dimnames(x))
## }


###################################################
### code chunk number 39: make-homologues (eval = FALSE)
###################################################
## h1 <- homologue(snp, 1)
## h2 <- homologue(snp, 2)
## stopifnot(identical(dimnames(h1), dimnames(h2)))


###################################################
### code chunk number 40: plot-sequences (eval = FALSE)
###################################################
## block.sequence <- function(markers) {
##     cbind(h1[markers, , drop = FALSE], h2[markers, , drop = FALSE])
## }
## sequences <- lapply(block.markers, block.sequence)
## complete.calls <- function(x) {
##     apply(x, 2, function(z) {
##         if (anyNA(z)) NA_character_ else paste(z, collapse = ";")
##     })
## }


###################################################
### code chunk number 41: plot-codes (eval = FALSE)
###################################################
## plot.block <- function(x) {
##     z <- complete.calls(x)
##     # Rank alleles across all three pools for shared figure colours.
##     counts <- sort(table(z), decreasing = TRUE)
##     ranks <- match(z, names(counts))
##     ranks[is.na(ranks)] <- -1L
##     n <- ncol(x) / 2L
##     paste(ranks[seq_len(n)], ranks[n + seq_len(n)], sep = "/")
## }
## plot.genotypes <- t(vapply(sequences, plot.block,
##                           character(ncol(snp))))


###################################################
### code chunk number 42: plot-map (eval = FALSE)
###################################################
## colnames(plot.genotypes) <- colnames(snp)
## rownames(plot.genotypes) <- gain.A$blocks$Name
## plot.map <- gain.A$blocks[, c("Chrom", "Pos", "Name")]
## max.rank <- max(vapply(sequences, function(x) {
##     length(unique(na.omit(complete.calls(x))))
## }, integer(1)))
## base.colours <- eval(formals(st.plot.ggt)$color)
## ggt.colours <- c(base.colours,
##     hcl.colors(max(1, max.rank - length(base.colours)), "Dark 3"))
## ggt.colours <- ggt.colours[seq_len(max.rank)]


###################################################
### code chunk number 43: plot-choice (eval = FALSE)
###################################################
## # Read base IDs from the actual enrichment dataset, not a new pool.
## base.ids <- sort(colnames(phased[[1]]$genotypes)[-1])
## selected <- list(A = list(Lines = forward.A), B = list(Lines = forward.B))


###################################################
### code chunk number 44: ggt-function (eval = FALSE)
###################################################
## draw.ggt <- function(group) {
##     ids <- c(selected[[group]]$Lines, base.ids)
##     labels <- paste(c(rep(group, 5), rep("M", length(base.ids))), ids)
##     st.plot.ggt.src(
##         map = plot.map,
##         pop = as.data.frame(plot.genotypes[, ids]),
##         f.lr = 0.7,
##         f.tb = 2,
##         d.h = 0.001,
##         d.v = 0.0003,
##         p.t = FALSE,
##         p.s = FALSE,
##         d.t = 50,
##         d.map = 100,
##         c.nme = as.character(sort(unique(plot.map$Chrom))),
##         i.nme = labels,
##         color = ggt.colours,
##         z.min = 1,
##         z.max = max.rank,
##         cex.chrom = 0.6,
##         cex.ind = 0.53
##     )
## }


###################################################
### code chunk number 45: ggt-a (eval = FALSE)
###################################################
## png(
##     file.path(figure.dir, "enrichment-ggt-A.png"),
##     width = 1000,
##     height = 1414,
##     pointsize = 20,
##     res = 120,
##     antialias = "none"
## )
## draw.ggt("A")
## invisible(dev.off())


###################################################
### code chunk number 46: ggt-b (eval = FALSE)
###################################################
## png(
##     file.path(figure.dir, "enrichment-ggt-B.png"),
##     width = 1000,
##     height = 1414,
##     pointsize = 20,
##     res = 120,
##     antialias = "none"
## )
## draw.ggt("B")
## invisible(dev.off())


###################################################
### code chunk number 47: core-select (eval = FALSE)
###################################################
## core <- st.detect.core(
##     pool               = "base.phased",
##     n.max              = length(base.ids),
##     blocks             = blocks,
##     exclude.compatible = TRUE,
##     method             = "forward"
## )


###################################################
### code chunk number 48: core-summary (eval = FALSE)
###################################################
## core.curve <- data.frame(
##     Lines = lengths(core$genotypes), # Each entry is the whole core
##     Haplotypes = core$segments      # Cumulative, not marginal gains
## )
## print(rbind(head(core.curve, 3), tail(core.curve, 1)),
##       row.names = FALSE)


