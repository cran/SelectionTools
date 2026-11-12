### R code from vignette source 'v-lmm.Rnw'

###################################################
### code chunk number 1: I01-fit (eval = FALSE)
###################################################
## library(SelectionTools)
## data("v-lmm-rapeseed", package = "SelectionTools")
## 
## m.rapeseed <- st.lmm(
##     fixed = hgt ~ geno * env,
##     random = ~env:rep,
##     data = rapeseed
## )
## m.rapeseed


###################################################
### code chunk number 2: I01-lsm (eval = FALSE)
###################################################
## st.lsm(m.rapeseed, "geno")


###################################################
### code chunk number 3: I01-blup (eval = FALSE)
###################################################
## data("v-lmm-tropmaize", package = "SelectionTools")
## data("v-tropmaize-vcf", package = "SelectionTools")
## st.load.vcf.data(v.tropmaize.vcf, data.set = "maize.markers")
## A <- st.get.A("maize.markers")
## maize.model <- st.lmm(
##     fixed  = y ~ 1,
##     random = ~ gen + block,
##     ucov   = list(gen = A),
##     data   = tropmaize
## )


###################################################
### code chunk number 4: I01-predictions (eval = FALSE)
###################################################
## maize.predictions <- st.predict(maize.model, "gen")
## maize.predictions


###################################################
### code chunk number 5: I01-chicken (eval = FALSE)
###################################################
## data("v-lmm-chickens", package = "SelectionTools")
## chickens$weight.s <- sqrt(chickens$weight)
## chicken.model <- st.lmm(
##     fixed  = weight.s ~ Time:Diet,
##     random = ~ (0 + Time) | Chick,
##     cov    = ~ at(Time):units,
##     data   = chickens
## )


###################################################
### code chunk number 6: I01-slopes (eval = FALSE)
###################################################
## st.lsm(
##     chicken.model,
##     "Diet",
##     trend = "Time"
## )


###################################################
### code chunk number 7: I01-minimum (eval = FALSE)
###################################################
## st.lsm(
##     chicken.model,
##     "Diet",
##     trend = "Time",
##     method = "trt.vs.ctrl",
##     ref = "d1",
##     rhs = .05,
##     side = ">",
##     adjust = "single-step"
## )


###################################################
### code chunk number 8: I01-means (eval = FALSE)
###################################################
## st.lsm(
##     chicken.model,
##     "Diet",
##     at = list(Time = 21)
## )


###################################################
### code chunk number 9: I02-selection (eval = FALSE)
###################################################
## data("v-lmm-mtcars", package = "SelectionTools")
## model.selection <- st.lmm(
##     mpg ~ disp + hp + wt + qsec + am,
##     data = mtcars,
##     method = "moments",
##     selection = list(
##         search = "exhaustive",
##         nvmax = 5,
##         criterion = "BICm"
##     )
## )
## summary(model.selection)


###################################################
### code chunk number 10: PRINT-options (eval = FALSE)
###################################################
## st.options()
## st.options(list(
##     print.digits = 3,
##     print.char.perline = 60,
##     print.top = 5,
##     print.bottom = 2
## ))
## coef(m.rapeseed)


###################################################
### code chunk number 11: OPT-optimizers (eval = FALSE)
###################################################
## data("v-lmm-rapeseed", package = "SelectionTools")
## model.auto <- st.lmm(
##     hgt ~ geno * env,
##     random = ~env:rep,
##     data = rapeseed,
##     control = list(optimizer = "auto", backend = "dense")
## )


###################################################
### code chunk number 12: CM01-1 (eval = FALSE)
###################################################
## data("v-lmm-climate", package = "SelectionTools")
## m.climate <- st.lmm(
##     yield ~ climate * variety,
##     data = climate,
##     method = "moments"
## )
## m.climate


###################################################
### code chunk number 13: CM01-2 (eval = FALSE)
###################################################
## st.lsm(m.climate, ~climate)


###################################################
### code chunk number 14: CM01-2-b (eval = FALSE)
###################################################
## st.lsm(
##     m.climate,
##     ~climate,
##     method = "pairwise",
##     adjust = "tukey"
## )


###################################################
### code chunk number 15: CM01-3 (eval = FALSE)
###################################################
## st.lsm(m.climate, ~climate * variety)


###################################################
### code chunk number 16: CM02-1 (eval = FALSE)
###################################################
## data("v-lmm-soils", package = "SelectionTools")
## m.soils <- st.lmm(
##     fixed = yield ~ region,
##     random = ~region:village,
##     data = soils,
##     method = "moments"
## )
## m.soils


###################################################
### code chunk number 17: CM02-details (eval = FALSE)
###################################################
## st.options(list(
##     print.ems = TRUE,
##     print.test.construction = TRUE
## ))
## m.soils


###################################################
### code chunk number 18: CM02-region-means (eval = FALSE)
###################################################
## st.lsm(m.soils, ~region)


###################################################
### code chunk number 19: CM02-region-effects (eval = FALSE)
###################################################
## st.lsm(
##     m.soils,
##     ~region,
##     method = "eff"
## )


###################################################
### code chunk number 20: CM02-villages (eval = FALSE)
###################################################
## m.villages <- st.lmm(
##     fixed = yield ~ village,
##     data = soils,
##     method = "moments"
## )
## m.villages


###################################################
### code chunk number 21: CM02-villages-cld (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(m.villages, ~village),
##     sort = TRUE,
##     adjust = "tukey"
## )


###################################################
### code chunk number 22: CM15-1 (eval = FALSE)
###################################################
## data("v-lmm-rapeseed", package = "SelectionTools")
## st.options(list(
##     print.ems = TRUE,
##     print.test.construction = TRUE
## ))


###################################################
### code chunk number 23: CM15-2 (eval = FALSE)
###################################################
## m.rape <- st.lmm(
##     fixed = hgt ~ geno * env,
##     random = ~rep:env,
##     data = rapeseed,
##     method = "moments"
## )
## m.rape


###################################################
### code chunk number 24: CM15-8 (eval = FALSE)
###################################################
## st.lsm(m.rape, ~geno)


###################################################
### code chunk number 25: CM15-combinations (eval = FALSE)
###################################################
## st.lsm(m.rape, ~geno * env)


###################################################
### code chunk number 26: CM15-e2-cld (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(
##         m.rape,
##         ~geno | env,
##         at = list(env = "e2")
##     ),
##     alpha = 0.10,
##     adjust = "tukey"
## )


###################################################
### code chunk number 27: CM16-1 (eval = FALSE)
###################################################
## data("v-lmm-barley-met", package = "SelectionTools")
## st.options(list(
##     print.ems = TRUE,
##     print.test.construction = TRUE
## ))
## m.m23.moments <- st.lmm(
##     fixed = TKW ~ 1,
##     random = ~gen + loc + year + gen:loc + gen:year + loc:year +
##         gen:loc:year + rep:loc:year,
##     data = barley.met,
##     method = "moments"
## )
## m.m23.moments


###################################################
### code chunk number 28: M01-1 (eval = FALSE)
###################################################
## data("v-lmm-soils", package = "SelectionTools")
## st.lmm(
##     fixed = yield ~ region,
##     random = ~region:village,
##     data = soils
## )


###################################################
### code chunk number 29: M02-1 (eval = FALSE)
###################################################
## data("v-lmm-rapeseed", package = "SelectionTools")
## m <- st.lmm(
##     fixed = hgt ~ geno * env,
##     random = ~env:rep,
##     data = rapeseed
## )
## m


###################################################
### code chunk number 30: M03_NA_NA-2 (eval = FALSE)
###################################################
## lsmeans.G <- st.lsm(m, "geno")
## lsmeans.G


###################################################
### code chunk number 31: M03_NA_NA-3 (eval = FALSE)
###################################################
## data.frame(lsmeans.G)$SE[1]


###################################################
### code chunk number 32: M03_NA_NA-4 (eval = FALSE)
###################################################
## data.frame(lsmeans.G)$DF[1]


###################################################
### code chunk number 33: M03_NA_NA-5 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     "env",
##     by = "geno"
## )


###################################################
### code chunk number 34: M23-1 (eval = FALSE)
###################################################
## data("v-lmm-barley-met", package = "SelectionTools")
## m.m23.reml <- st.lmm(
##     fixed = TKW ~ 1,
##     random = ~gen + loc + year + gen:loc + gen:year + loc:year +
##         gen:loc:year + rep:loc:year,
##     data = barley.met
## )
## m.m23.reml


###################################################
### code chunk number 35: M30-1 (eval = FALSE)
###################################################
## data("v-lmm-material-groups", package = "SelectionTools")
## st.lmm(
##     fixed = yield ~ 1,
##     random = ~at(grp):gen + env,
##     data = material.groups
## )


###################################################
### code chunk number 36: M37_NA_NA-data (eval = FALSE)
###################################################
## data("v-lmm-barley-met", package = "SelectionTools")
## xtabs(~loc + year, data = barley.met, subset = !is.na(TKW))


###################################################
### code chunk number 37: M37_NA_NA-1 (eval = FALSE)
###################################################
## m <- st.lmm(
##     fixed = TKW ~ gen + loc + year + loc:year + year:loc:rep,
##     random = ~gen:loc + gen:year + gen:loc:year,
##     data = barley.met
## )
## m


###################################################
### code chunk number 38: M37_NA_NA-details (eval = FALSE)
###################################################
## st.options(list(print.test.construction = TRUE))
## print(m)


###################################################
### code chunk number 39: M37_NA_NA-2 (eval = FALSE)
###################################################
## st.lsm(m, "gen")


###################################################
### code chunk number 40: M39_NA_NA-data (eval = FALSE)
###################################################
## data("v-lmm-met-unbalanced", package = "SelectionTools")
## xtabs(~loc + year, data = met.unbalanced, subset = !is.na(yld))


###################################################
### code chunk number 41: M39_NA_NA-1 (eval = FALSE)
###################################################
## m <- st.lmm(
##     fixed = yld ~ gen + env + rep:env,
##     random = ~gen:env,
##     data = met.unbalanced
## )
## m


###################################################
### code chunk number 42: M39_NA_NA-2 (eval = FALSE)
###################################################
## st.lsm(m, "gen")


###################################################
### code chunk number 43: M41_NA_NA-data (eval = FALSE)
###################################################
## data("v-lmm-met-missing", package = "SelectionTools")
## xtabs(~loc + year, data = met.missing, subset = !is.na(YLD))


###################################################
### code chunk number 44: M41_NA_NA-1 (eval = FALSE)
###################################################
## m <- st.lmm(
##     fixed = YLD ~ gen + loc + year + loc:year + year:loc:rep,
##     random = ~gen:loc + gen:year + gen:loc:year,
##     data = met.missing
## )
## m


###################################################
### code chunk number 45: M41_NA_NA-2 (eval = FALSE)
###################################################
## st.lsm(m, "gen")


###################################################
### code chunk number 46: M42_NA_NA-data (eval = FALSE)
###################################################
## data("v-lmm-met-missing", package = "SelectionTools")
## xtabs(~loc + year, data = met.missing, subset = !is.na(YLD))


###################################################
### code chunk number 47: M42_NA_NA-1 (eval = FALSE)
###################################################
## met.missing$REP <- met.missing$rep:met.missing$year:met.missing$loc
## m <- st.lmm(
##     fixed = YLD ~ gen + REP,
##     random = ~gen:loc + gen:year + gen:loc:year,
##     data = met.missing
## )
## m


###################################################
### code chunk number 48: M42_NA_NA-2 (eval = FALSE)
###################################################
## st.lsm(m, "gen")


###################################################
### code chunk number 49: MCH-inference-fit (eval = FALSE)
###################################################
## data("v-lmm-machines", package = "SelectionTools")
## m.machines <- st.lmm(
##     fixed = score ~ machine,
##     random = ~operator + machine:operator,
##     data = machines
## )


###################################################
### code chunk number 50: MCH-inference-broad (eval = FALSE)
###################################################
## st.lsm(m.machines, ~machine)


###################################################
### code chunk number 51: MCH-inference-intermediate-operator (eval = FALSE)
###################################################
## st.lsm(
##     m.machines,
##     ~machine,
##     inference.space = ~operator
## )


###################################################
### code chunk number 52: MCH-inference-intermediate-interaction (eval = FALSE)
###################################################
## st.lsm(
##     m.machines,
##     ~machine,
##     inference.space = ~machine:operator
## )


###################################################
### code chunk number 53: MCH-inference-narrow (eval = FALSE)
###################################################
## st.lsm(
##     m.machines,
##     ~machine,
##     inference.space = "narrow"
## )


###################################################
### code chunk number 54: CM04-1 (eval = FALSE)
###################################################
## data("v-lmm-fertilizer", package = "SelectionTools")
## m.fertilizer <- st.lmm(
##     yield ~ fert,
##     data = fertilizer,
##     method = "moments"
## )
## m.fertilizer


###################################################
### code chunk number 55: CM04-2 (eval = FALSE)
###################################################
## st.lsm(m.fertilizer, ~fert)


###################################################
### code chunk number 56: CM04-2-b (eval = FALSE)
###################################################
## st.lsm(
##     m.fertilizer,
##     ~fert,
##     method = "pairwise",
##     adjust = "none"
## )


###################################################
### code chunk number 57: CM04-3 (eval = FALSE)
###################################################
## st.lsm(
##     m.fertilizer,
##     ~fert,
##     method = "pairwise",
##     adjust = "tukey"
## )


###################################################
### code chunk number 58: CM05-1 (eval = FALSE)
###################################################
## data("v-lmm-cotton", package = "SelectionTools")
## m.cotton <- st.lmm(
##     Strength ~ Blk + Fert,
##     data = cotton,
##     method = "moments"
## )
## m.cotton


###################################################
### code chunk number 59: CM05-2 (eval = FALSE)
###################################################
## st.lsm(m.cotton, ~Fert)


###################################################
### code chunk number 60: CM05-2-b (eval = FALSE)
###################################################
## st.lsm(
##     m.cotton,
##     ~Fert,
##     method = "pairwise",
##     adjust = "tukey"
## )


###################################################
### code chunk number 61: CM06-1 (eval = FALSE)
###################################################
## data("v-lmm-sugar", package = "SelectionTools")
## m.sugar <- st.lmm(
##     Yield ~ Row + Col + Variety,
##     data = sugar,
##     method = "moments"
## )
## m.sugar


###################################################
### code chunk number 62: CM06-2 (eval = FALSE)
###################################################
## st.lsm(m.sugar, ~Variety)


###################################################
### code chunk number 63: CM06-2-b (eval = FALSE)
###################################################
## st.lsm(
##     m.sugar,
##     ~Variety,
##     method = "pairwise",
##     adjust = "none"
## )


###################################################
### code chunk number 64: CM07-1 (eval = FALSE)
###################################################
## data("v-lmm-peas", package = "SelectionTools")
## m.peas <- st.lmm(
##     Yield ~ Block + Spacing + Variety + Spacing:Variety,
##     data = peas,
##     method = "moments"
## )
## m.peas


###################################################
### code chunk number 65: CM07-2 (eval = FALSE)
###################################################
## st.lsm(m.peas, ~Spacing | Variety)


###################################################
### code chunk number 66: CM07-2-b (eval = FALSE)
###################################################
## st.lsm(
##     m.peas,
##     ~Spacing | Variety,
##     method = "pairwise",
##     adjust = "none"
## )


###################################################
### code chunk number 67: CM08-1 (eval = FALSE)
###################################################
## data("v-lmm-alfalfa", package = "SelectionTools")
## m.alfalfa <- st.lmm(
##     fixed = Yield ~ Variety * Date,
##     random = ~Block + Block:Variety,
##     data = alfalfa,
##     method = "moments"
## )
## m.alfalfa


###################################################
### code chunk number 68: CM08-2 (eval = FALSE)
###################################################
## st.lsm(
##     m.alfalfa,
##     ~Variety,
##     method = "pairwise",
##     adjust = "none"
## )


###################################################
### code chunk number 69: CM08-2-b (eval = FALSE)
###################################################
## st.lsm(
##     m.alfalfa,
##     ~Date,
##     method = "pairwise",
##     adjust = "none"
## )


###################################################
### code chunk number 70: CM08-2-c (eval = FALSE)
###################################################
## st.lsm(
##     m.alfalfa,
##     ~Date | Variety,
##     method = "pairwise",
##     adjust = "none"
## )


###################################################
### code chunk number 71: CM08-3 (eval = FALSE)
###################################################
## st.lsm(m.alfalfa, ~Variety | Date)


###################################################
### code chunk number 72: CM08-4 (eval = FALSE)
###################################################
## c.variety.date <- st.lsm(
##     m.alfalfa,
##     ~Variety | Date,
##     method = "pairwise",
##     adjust = "none"
## )
## c.variety.date


###################################################
### code chunk number 73: CM08-5 (eval = FALSE)
###################################################
## st.lsm.test(c.variety.date)


###################################################
### code chunk number 74: CM08-6 (eval = FALSE)
###################################################
## # Mean order: Cossack/A, Ladack/A, Cossack/B, Ladack/B.
## C.alf <- rbind(
##     within.plot = c(1, 0, -1, 0),
##     across.plots = c(1, 0, 0, -1)
## )
## # Within Cossack: A - B. Across varieties: Cossack/A - Ladack/B.
## h.alf.cells <- st.lsm(
##     m.alfalfa,
##     ~Variety * Date,
##     at = list(
##         Variety = c("Cossack", "Ladack"),
##         Date = c("A", "B")
##     ),
##     method = C.alf,
##     adjust = "none"
## )
## h.alf.cells


###################################################
### code chunk number 75: CM08-7 (eval = FALSE)
###################################################
## # Joint test of the two contrasts, allowing their different error strata.
## st.lsm(
##     h.alf.cells,
##     joint = TRUE
## )


###################################################
### code chunk number 76: L01-1 (eval = FALSE)
###################################################
## data("v-lmm-fertilizer", package = "SelectionTools")
## m <- st.lmm(
##     yield ~ fert,
##     data = fertilizer
## )
## m


###################################################
### code chunk number 77: L01-2 (eval = FALSE)
###################################################
## st.lsm(m, ~fert)


###################################################
### code chunk number 78: L01-3 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~fert,
##     method = "pairwise",
##     adjust = "none"
## )


###################################################
### code chunk number 79: L01-4 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~fert,
##     method = "pairwise",
##     adjust = "tukey"
## )


###################################################
### code chunk number 80: L01-5 (eval = FALSE)
###################################################
## st.lsm.cld(st.lsm(m, ~fert))


###################################################
### code chunk number 81: L01-6 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~fert,
##     method = "trt.vs.ctrl",
##     ref = 3,
##     side = "<",
##     adjust = "mvt"
## )


###################################################
### code chunk number 82: L01-7 (eval = FALSE)
###################################################
## C1 <- list(
##     CCr1_CCr2 = c(1, -1, 0, 0, 0),
##     MinA_MinB = c(0, 0, 0, 1, -1),
##     Org_Min = c(1/3, 1/3, 1/3, -1/2, -1/2),
##     CCr_Man = c(1/2, 1/2, -1, 0, 0)
## )
## st.lsm(
##     m,
##     ~fert,
##     method = C1
## )


###################################################
### code chunk number 83: L01-8 (eval = FALSE)
###################################################
## C2 <- list(
##     Man_CCr1 = c(-1, 0, 1, 0, 0),
##     Man_CCr2 = c(0, -1, 1, 0, 0),
##     Man_MinA = c(0, 0, 1, -1, 0),
##     Man_MinB = c(0, 0, 1, 0, -1)
## )
## st.lsm(
##     m,
##     ~fert,
##     method = C2,
##     adjust = "hochberg"
## )


###################################################
### code chunk number 84: L01-9 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~fert,
##     method = C2,
##     side = ">",
##     adjust = "hochberg"
## )


###################################################
### code chunk number 85: L01-10 (eval = FALSE)
###################################################
## C3 <- list(
##     MinA_CCr1 = c(-1, 0, 0, 1, 0),
##     MinA_CCr2 = c(0, -1, 0, 1, 0),
##     MinA_Man = c(0, 0, -1, 1, 0),
##     MinB_CCr1 = c(-1, 0, 0, 0, 1),
##     MinB_CCr2 = c(0, -1, 0, 0, 1),
##     MinB_Man = c(0, 0, -1, 0, 1)
## )
## st.lsm(
##     m,
##     ~fert,
##     method = C3,
##     adjust = "fdr"
## )


###################################################
### code chunk number 86: L01-11 (eval = FALSE)
###################################################
## # Numeric contrast rows act on the five displayed fertilizer means.
## C.fert <- do.call(rbind, C1)
## h.fert <- st.lsm(
##     m,
##     ~fert,
##     method = C.fert,
##     adjust = "single-step"
## )
## h.fert


###################################################
### code chunk number 87: L01-12 (eval = FALSE)
###################################################
## # The four independent contrasts jointly test equality of all five means.
## st.lsm(
##     h.fert,
##     joint = TRUE
## )


###################################################
### code chunk number 88: L02-1 (eval = FALSE)
###################################################
## data("v-lmm-cotton", package = "SelectionTools")
## m <- st.lmm(
##     Strength ~ Blk + Fert,
##     data = cotton
## )
## m


###################################################
### code chunk number 89: L02-2 (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(m, ~Fert),
##     adjust = "none"
## )


###################################################
### code chunk number 90: L02-3 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~Fert,
##     method = "pairwise",
##     adjust = "none"
## )


###################################################
### code chunk number 91: L02-4 (eval = FALSE)
###################################################
## st.lsm.cld(st.lsm(m, ~Fert))


###################################################
### code chunk number 92: L02-5 (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(m, ~Fert),
##     alpha = 0.1
## )


###################################################
### code chunk number 93: L02-6 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~Fert,
##     method = "pairwise",
##     adjust = "tukey"
## )


###################################################
### code chunk number 94: L03-1 (eval = FALSE)
###################################################
## data("v-lmm-sugar", package = "SelectionTools")
## m <- st.lmm(
##     Yield ~ Row + Col + Variety,
##     data = sugar
## )
## m


###################################################
### code chunk number 95: L03-2 (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(m, ~Variety),
##     adjust = "none"
## )


###################################################
### code chunk number 96: L03-3 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~Variety,
##     method = "pairwise",
##     adjust = "none"
## )


###################################################
### code chunk number 97: L04-1 (eval = FALSE)
###################################################
## data("v-lmm-peas", package = "SelectionTools")
## m <- st.lmm(
##     Yield ~ Block + Spacing * Variety,
##     data = peas
## )
## m


###################################################
### code chunk number 98: L04-2 (eval = FALSE)
###################################################
## st.lsm(m, ~Spacing:Variety)


###################################################
### code chunk number 99: L04-3 (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(m, ~Spacing | Variety),
##     adjust = "none",
##     sort = FALSE
## )


###################################################
### code chunk number 100: L04-4 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~Spacing | Variety,
##     method = "pairwise",
##     adjust = "none"
## )


###################################################
### code chunk number 101: L04-5 (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(m, ~Variety | Spacing),
##     adjust = "none",
##     sort = FALSE
## )


###################################################
### code chunk number 102: L04-6 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~Variety | Spacing,
##     method = "pairwise",
##     adjust = "none"
## )


###################################################
### code chunk number 103: M04_NA_NA-1 (eval = FALSE)
###################################################
## data("v-lmm-soybean-blocks", package = "SelectionTools")
## soybean.blocks$Block.in.Rep <- soybean.blocks$Block:soybean.blocks$Rep
## m <- st.lmm(
##     fixed = Yield ~ Variety,
##     random = ~Rep + Block.in.Rep,
##     data = soybean.blocks
## )
## m


###################################################
### code chunk number 104: M04_NA_NA-2 (eval = FALSE)
###################################################
## st.lsm(m, ~Variety)


###################################################
### code chunk number 105: M06_NA_NA-1 (eval = FALSE)
###################################################
## data("v-lmm-oats", package = "SelectionTools")
## m <- st.lmm(
##     fixed = Yield ~ Variety + Nitro + Variety:Nitro,
##     random = ~Block + Block:Variety,
##     data = oats
## )
## m


###################################################
### code chunk number 106: M06_NA_NA-2 (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(m, ~Nitro),
##     adjust = "hochberg"
## )


###################################################
### code chunk number 107: M06_NA_NA-3 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~Nitro,
##     method = "pairwise",
##     adjust = "hochberg"
## )


###################################################
### code chunk number 108: M06_NA_NA-4 (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(m, ~Variety),
##     adjust = "hochberg"
## )


###################################################
### code chunk number 109: M06_NA_NA-5 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~Variety,
##     method = "pairwise",
##     adjust = "hochberg"
## )


###################################################
### code chunk number 110: M06_NA_NA-6 (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(m, ~Variety:Nitro),
##     adjust = "hochberg",
##     alpha = 0.1,
##     sort = FALSE
## )


###################################################
### code chunk number 111: M06_NA_NA-7 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~Variety:Nitro,
##     method = "pairwise",
##     adjust = "hochberg"
## )


###################################################
### code chunk number 112: M07_NA_NA-1 (eval = FALSE)
###################################################
## data("v-lmm-barley-lattice", package = "SelectionTools")
## m.barley.lattice <- st.lmm(
##     fixed = yield ~ Variety,
##     random = ~Rep + RowBlk,
##     data = barley.lattice
## )
## m.barley.lattice


###################################################
### code chunk number 113: M08_NA_NA-1 (eval = FALSE)
###################################################
## data("v-lmm-barley-lattice", package = "SelectionTools")
## m <- st.lmm(
##     fixed = yield ~ Variety,
##     random = ~Rep + RowBlk + ColBlk,
##     data = barley.lattice
## )
## m


###################################################
### code chunk number 114: M08_NA_NA-2 (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(m, ~Variety),
##     adjust = "hochberg"
## )


###################################################
### code chunk number 115: M08_NA_NA-3 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~Variety,
##     method = "pairwise",
##     adjust = "hochberg"
## )


###################################################
### code chunk number 116: M08_NA_NA-4 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~Variety,
##     method = "trt.vs.ctrl",
##     ref = 4,
##     side = ">",
##     adjust = "none"
## )


###################################################
### code chunk number 117: LMM-barley-ar1xar1 (eval = FALSE)
###################################################
## data("v-lmm-barley-lattice", package = "SelectionTools")
## m.barley.ar1xar1 <- st.lmm(
##     fixed = yield ~ Variety,
##     cov = ~ar1(Row):ar1(Column),
##     data = barley.lattice
## )
## m.barley.ar1xar1


###################################################
### code chunk number 118: M09_NA_NA-1 (eval = FALSE)
###################################################
## data("v-lmm-rapeseed-blocks", package = "SelectionTools")
## m <- st.lmm(
##     fixed = hgt ~ geno + env + geno:env,
##     random = ~rep:env + blk:rep:env,
##     data = rapeseed.blocks
## )
## m


###################################################
### code chunk number 119: M09_NA_NA-2 (eval = FALSE)
###################################################
## st.lsm(m, ~geno)


###################################################
### code chunk number 120: M09_NA_NA-3 (eval = FALSE)
###################################################
## st.lsm(m, ~geno | env)


###################################################
### code chunk number 121: M09_NA_NA-4 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~geno:env,
##     at = list(env = "s")
## )


###################################################
### code chunk number 122: M09_NA_NA-5 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~geno:env,
##     at = list(env = "s"),
##     method = "pairwise",
##     adjust = "none",
##     level = 0.9
## )


###################################################
### code chunk number 123: M09_NA_NA-6 (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(
##         m,
##         ~geno:env,
##         at = list(env = "s")
##     ),
##     adjust = "none",
##     alpha = 0.1
## )


###################################################
### code chunk number 124: M09_NA_NA-7 (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(
##         m,
##         ~geno:env,
##         at = list(env = "s")
##     ),
##     adjust = "FDR",
##     alpha = 0.1
## )


###################################################
### code chunk number 125: L05-TKW-1 (eval = FALSE)
###################################################
## data("v-lmm-barley-met", package = "SelectionTools")
## m <- st.lmm(
##     TKW ~ gen * loc * year + year:loc:rep,
##     data = barley.met
## )
## m


###################################################
### code chunk number 126: L05-TKW-2 (eval = FALSE)
###################################################
## st.lsm(m, "gen")


###################################################
### code chunk number 127: L05-TKW-3 (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(m, "gen"),
##     reversed = TRUE,
##     alpha = 0.1,
##     adjust = "none"
## )


###################################################
### code chunk number 128: M21_y2016_HOH_TKW-1 (eval = FALSE)
###################################################
## data("v-lmm-barley-met", package = "SelectionTools")
## barley.met$env <- barley.met$year:barley.met$loc
## barley.met <- subset(barley.met, env == "y2016:HOH")
## st.lmm(
##     fixed = TKW ~ 1,
##     random = ~gen + rep,
##     data = barley.met
## )


###################################################
### code chunk number 129: M15-1 (eval = FALSE)
###################################################
## data("v-lmm-maize-counties", package = "SelectionTools")
## m <- st.lmm(
##     fixed = yield ~ gen + county + gen:county,
##     random = ~county:rep + county:rep:block,
##     cov = ~at(county):units,
##     data = maize.counties
## )
## m


###################################################
### code chunk number 130: M14-2 (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(m, ~gen),
##     adjust = "none",
##     alpha = 0.1
## )


###################################################
### code chunk number 131: M16-2 (eval = FALSE)
###################################################
## st.lsm(m, ~gen)


###################################################
### code chunk number 132: M16-3 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~gen,
##     method = "pairwise",
##     adjust = "none",
##     level = 0.9
## )


###################################################
### code chunk number 133: M16-4 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~gen,
##     at = list(county = "C5")
## )


###################################################
### code chunk number 134: M16-5 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     ~gen,
##     at = list(county = "C5"),
##     method = "pairwise",
##     adjust = "none",
##     level = 0.9
## )


###################################################
### code chunk number 135: M16-6 (eval = FALSE)
###################################################
## st.lsm.cld(
##     st.lsm(
##         m,
##         ~gen,
##         at = list(county = "C5")
##     ),
##     adjust = "none",
##     alpha = 0.1
## )


###################################################
### code chunk number 136: M26-1 (eval = FALSE)
###################################################
## data("v-lmm-soybean", package = "SelectionTools")
## st.lmm(
##     fixed = yield ~ 1,
##     random = ~gen + env + gen:env,
##     cov = ~at(env):units,
##     data = soybean
## )


###################################################
### code chunk number 137: M48-1 (eval = FALSE)
###################################################
## m <- st.lmm(
##     fixed = HLW ~ gen,
##     random = ~env + gen:env + env:rep,
##     cov = ~at(env):units,
##     data = barley.seven.traits
## )
## m


###################################################
### code chunk number 138: M48-2 (eval = FALSE)
###################################################
## st.lsm(m, "gen")


###################################################
### code chunk number 139: HER-barley-data (eval = FALSE)
###################################################
## data("v-lmm-barley-met", package = "SelectionTools")
## xtabs(~loc + year, data = barley.met)


###################################################
### code chunk number 140: HER-barley-variance (eval = FALSE)
###################################################
## m.barley.reml <- st.lmm(
##     fixed = TKW ~ 1,
##     random = ~gen + loc + year + gen:loc + gen:year + loc:year +
##         gen:loc:year + rep:loc:year,
##     data = barley.met
## )
## v <- st.lmm.components(m.barley.reml, "variance")[, "TKW"]
## v.g <- as.numeric(v[["gen"]])
## v.g


###################################################
### code chunk number 141: HER-barley-components (eval = FALSE)
###################################################
## L <- nlevels(barley.met$loc)
## Y <- nlevels(barley.met$year)
## R <- nlevels(barley.met$rep)
## v.masking <- v[["gen:loc"]] / L +
##     v[["gen:year"]] / Y +
##     v[["gen:loc:year"]] / (L * Y) +
##     v[["units"]] / (L * Y * R)
## v.masking


###################################################
### code chunk number 142: HER-barley-heritability (eval = FALSE)
###################################################
## heritability <- v.g / (v.g + v.masking)
## heritability


###################################################
### code chunk number 143: HER-barley-means-fit (eval = FALSE)
###################################################
## m.barley.means <- st.lmm(
##     fixed = TKW ~ gen + loc + year + loc:year + year:loc:rep,
##     random = ~gen:loc + gen:year + gen:loc:year,
##     data = barley.met
## )


###################################################
### code chunk number 144: HER-barley-means (eval = FALSE)
###################################################
## barley.comparisons <- as.data.frame(st.lsm(
##     m.barley.means,
##     ~gen,
##     method = "pairwise",
##     adjust = "none"
## ))
## v.masking <- mean(barley.comparisons$SE^2) / 2
## v.masking


###################################################
### code chunk number 145: HER-barley-means-heritability (eval = FALSE)
###################################################
## heritability <- v.g / (v.g + v.masking)
## heritability


###################################################
### code chunk number 146: HER-barley-blup-errors (eval = FALSE)
###################################################
## PEV <- vcov(st.predict(m.barley.reml, "0 + gen"))
## PEV.differences <- outer(diag(PEV), diag(PEV), "+") - 2 * PEV
## v.diff <- mean(PEV.differences[upper.tri(PEV.differences)])
## v.diff


###################################################
### code chunk number 147: HER-barley-blup-heritability (eval = FALSE)
###################################################
## heritability <- 1 - v.diff / (2 * v.g)
## heritability


###################################################
### code chunk number 148: HER-missing-data (eval = FALSE)
###################################################
## data("v-lmm-met-missing", package = "SelectionTools")
## xtabs(~loc + year, data = met.missing, subset = !is.na(YLD))


###################################################
### code chunk number 149: HER-missing-variance (eval = FALSE)
###################################################
## met.vc <- st.lmm(
##     fixed = YLD ~ 1,
##     random = ~gen + loc + year + gen:loc + gen:year + loc:year +
##         gen:loc:year + year:loc:rep,
##     data = met.missing
## )
## v <- st.lmm.components(met.vc, "variance")
## v.g <- as.numeric(v["gen", "YLD"])
## v.g


###################################################
### code chunk number 150: HER-missing-means-fit (eval = FALSE)
###################################################
## met.missing$REP <- with(met.missing, rep:year:loc)
## met.means <- st.lmm(
##     fixed = YLD ~ gen + REP,
##     random = ~gen:loc + gen:year + gen:loc:year,
##     data = met.missing
## )


###################################################
### code chunk number 151: HER-missing-means (eval = FALSE)
###################################################
## met.comparisons <- as.data.frame(st.lsm(
##     met.means,
##     ~gen,
##     method = "pairwise",
##     adjust = "none"
## ))
## v.masking <- mean(met.comparisons$SE^2) / 2
## v.masking


###################################################
### code chunk number 152: HER-missing-heritability (eval = FALSE)
###################################################
## heritability <- v.g / (v.g + v.masking)
## heritability


###################################################
### code chunk number 153: HER-missing-blup-errors (eval = FALSE)
###################################################
## PEV <- vcov(st.predict(met.vc, "0 + gen"))
## PEV.differences <- outer(diag(PEV), diag(PEV), "+") - 2 * PEV
## v.diff <- mean(PEV.differences[upper.tri(PEV.differences)])
## v.diff


###################################################
### code chunk number 154: HER-missing-blup-heritability (eval = FALSE)
###################################################
## heritability <- 1 - v.diff / (2 * v.g)
## heritability


###################################################
### code chunk number 155: CM09-1 (eval = FALSE)
###################################################
## data("v-lmm-iris-two", package = "SelectionTools")
## m.hotelling <- st.lmm(
##     cbind(SepalLength, SepalWidth) ~ Species,
##     data = iris.two,
##     method = "moments"
## )
## anova(m.hotelling, effects = "manova", test = "Hotelling-Lawley")


###################################################
### code chunk number 156: CM10-1 (eval = FALSE)
###################################################
## data("v-lmm-iris", package = "SelectionTools")
## m.iris <- st.lmm(
##     cbind(SepalLength, SepalWidth, PetalLength, PetalWidth) ~ Species,
##     data = iris.manova,
##     method = "moments"
## )
## anova(m.iris, effects = "manova")


###################################################
### code chunk number 157: CM10-2 (eval = FALSE)
###################################################
## st.lmm.components(m.iris, "sscp.H", term = "Species")


###################################################
### code chunk number 158: CM10-2-b (eval = FALSE)
###################################################
## st.lmm.components(m.iris, "sscp.E")


###################################################
### code chunk number 159: M60-4 (eval = FALSE)
###################################################
## data("v-lmm-maize-met", package = "SelectionTools")
## m.maize.three <- st.lmm(
##     fixed = cbind(yld, dmc, pht) ~ 1,
##     random = ~gen + env + gen:env + rep:env,
##     data = maize.met,
##     method = "moments",
##     random.test = "F",
##     control = list(design = "general")
## )
## m.maize.three


###################################################
### code chunk number 160: M60-5 (eval = FALSE)
###################################################
## st.lmm.components(m.maize.three, "covariance", block = "gen")


###################################################
### code chunk number 161: M60-6 (eval = FALSE)
###################################################
## st.lmm.components(m.maize.three, "correlation", block = "gen")


###################################################
### code chunk number 162: maize-moments-means (eval = FALSE)
###################################################
## data("v-lmm-maize-met", package = "SelectionTools")
## m.maize.moments.means <- st.lmm(
##     fixed = cbind(yld, dmc, pht) ~ gen,
##     random = ~env + gen:env + rep:env,
##     data = maize.met,
##     method = "moments",
##     control = list(design = "general")
## )
## maize.moments.lsmeans <- as.data.frame(st.lsm(m.maize.moments.means, ~gen))
## maize.moments.means <- xtabs(lsmeans ~ gen + trait,
##     data = maize.moments.lsmeans)[, c("yld", "dmc", "pht")]
## maize.moments.means


###################################################
### code chunk number 163: maize-moments-phenotypic-cov (eval = FALSE)
###################################################
## round(cov(maize.moments.means), 4)


###################################################
### code chunk number 164: maize-moments-phenotypic-cor (eval = FALSE)
###################################################
## round(cor(maize.moments.means), 4)


###################################################
### code chunk number 165: maize-moments-means-two (eval = FALSE)
###################################################
## st.lsm(
##     m.maize.moments.means,
##     ~gen,
##     at = list(gen = c("G01", "G02"))
## )


###################################################
### code chunk number 166: maize-moments-pair (eval = FALSE)
###################################################
## st.lsm(
##     m.maize.moments.means,
##     ~gen,
##     at = list(gen = c("G01", "G02")),
##     method = "pairwise",
##     adjust = "none"
## )


###################################################
### code chunk number 167: maize-joint-means (eval = FALSE)
###################################################
## data("v-lmm-maize-met", package = "SelectionTools")
## m.maize.means <- st.lmm(
##     fixed = cbind(yld, dmc, pht) ~ gen,
##     random = ~env + gen:env + rep:env,
##     data = maize.met
## )
## maize.lsmeans <- as.data.frame(st.lsm(m.maize.means, ~gen))
## maize.means <- xtabs(lsmeans ~ gen + trait,
##     data = maize.lsmeans)[, c("yld", "dmc", "pht")]
## maize.means


###################################################
### code chunk number 168: maize-phenotypic-cov (eval = FALSE)
###################################################
## round(cov(maize.means), 4)


###################################################
### code chunk number 169: maize-phenotypic-cor (eval = FALSE)
###################################################
## round(cor(maize.means), 4)


###################################################
### code chunk number 170: maize-genetic-reml (eval = FALSE)
###################################################
## data("v-lmm-maize-met", package = "SelectionTools")
## m.maize.genetic <- st.lmm(
##     fixed = cbind(yld, dmc, pht) ~ 1,
##     random = ~gen + env + gen:env + rep:env,
##     data = maize.met
## )
## st.lmm.components(m.maize.genetic, "covariance", block = "gen")


###################################################
### code chunk number 171: maize-genetic-cor (eval = FALSE)
###################################################
## st.lmm.components(m.maize.genetic, "correlation", block = "gen")


###################################################
### code chunk number 172: M22-1 (eval = FALSE)
###################################################
## data("v-lmm-barley-met", package = "SelectionTools")
## barley.met$env <- factor(paste(barley.met$year, barley.met$loc))
## m.barley <- st.lmm(
##     fixed = cbind(GYLD, TKW, PLH) ~ 1,
##     random = ~gen + loc + year + gen:loc + gen:year + loc:year +
##         gen:loc:year + year:loc:rep,
##     cov = ~at(env):units,
##     data = barley.met,
##     control = list(maxit = 1000)
## )
## m.barley


###################################################
### code chunk number 173: M22-variance (eval = FALSE)
###################################################
## st.lmm.components(m.barley, "variance", traits = c("GYLD", "TKW", "PLH"))


###################################################
### code chunk number 174: M22-gen-cov (eval = FALSE)
###################################################
## st.lmm.components(m.barley, "covariance", block = "gen")


###################################################
### code chunk number 175: M47-plots (eval = FALSE)
###################################################
## addmargins(xtabs(~ env + rep, data = barley.seven.traits))


###################################################
### code chunk number 176: M47-genotypes (eval = FALSE)
###################################################
## addmargins(xtabs(~ gen + env, data = barley.seven.traits))


###################################################
### code chunk number 177: M47-1 (eval = FALSE)
###################################################
## m.met.seven <- st.lmm(
##     fixed = cbind(YLD, HLW, TKW, PLH, DH, DM, SC) ~ env + rep:env,
##     random = ~gen + gen:env,
##     cov = ~at(env):units,
##     data = barley.seven.traits,
##     random.test = "none",
##     control = list(backend = "block")
## )
## st.lmm.components(m.met.seven, "covariance", block = "gen")


###################################################
### code chunk number 178: M47-cor (eval = FALSE)
###################################################
## st.lmm.components(m.met.seven, "correlation", block = "gen")


###################################################
### code chunk number 179: M47-lsmeans (eval = FALSE)
###################################################
## m.met.seven.means <- st.lmm(
##     fixed = cbind(YLD, HLW, TKW, PLH, DH, DM, SC) ~ gen + env + rep:env,
##     random = ~gen:env,
##     cov = ~at(env):units,
##     data = barley.seven.traits,
##     random.test = "none",
##     control = list(backend = "block")
## )
## barley.lsmeans <- as.data.frame(st.lsm(m.met.seven.means, ~gen))
## round(xtabs(lsmeans ~ gen + trait, data = barley.lsmeans)[
##     , c("YLD", "HLW", "TKW", "PLH", "DH", "DM", "SC")
## ], 2)


###################################################
### code chunk number 180: CM11-1 (eval = FALSE)
###################################################
## data("v-lmm-fertilizer-regression", package = "SelectionTools")
## m.regression <- st.lmm(
##     ertrag ~ dueng,
##     data = fertilizer.regression,
##     method = "moments"
## )
## m.regression


###################################################
### code chunk number 181: CM14-1 (eval = FALSE)
###################################################
## data("v-lmm-fertilizer-replicated", package = "SelectionTools")
## m.replicated <- st.lmm(
##     ertrag ~ dueng + factor(dueng),
##     data = fertilizer.replicated,
##     method = "moments"
## )
## m.replicated


###################################################
### code chunk number 182: CM12-1 (eval = FALSE)
###################################################
## data("v-lmm-rats", package = "SelectionTools")
## m.rats <- st.lmm(
##     y ~ x.1 + x.2,
##     data = rats,
##     method = "moments"
## )
## m.rats


###################################################
### code chunk number 183: CM12-2 (eval = FALSE)
###################################################
## m.rats.reverse <- st.lmm(
##     y ~ x.2 + x.1,
##     data = rats,
##     method = "moments"
## )
## m.rats.reverse


###################################################
### code chunk number 184: CM12-4 (eval = FALSE)
###################################################
## L.rats <- matrix(
##         c(0, 1, 0), nrow = 1,
##         dimnames = list("initial weight", names(coef(m.rats)))
##     )
## st.lsm(
##     m.rats,
##     L = L.rats,
##     joint = TRUE
## )


###################################################
### code chunk number 185: CM12-4-b (eval = FALSE)
###################################################
## m.rats.reduced <- st.lmm(
##     y ~ x.2,
##     data = rats,
##     method = "moments"
## )
## anova(m.rats.reduced, m.rats)


###################################################
### code chunk number 186: CM17-search (eval = FALSE)
###################################################
## data("v-lmm-mtcars", package = "SelectionTools")
## m.mtcars.subsets <- st.lmm(
##     mpg ~ .,
##     data = mtcars,
##     method = "moments",
##     selection = list(
##         search = "exhaustive",
##         nvmax = 7,
##         nbest = 1,
##         criterion = "adjr2"
##     )
## )


###################################################
### code chunk number 187: CM17-summary (eval = FALSE)
###################################################
## summary(m.mtcars.subsets)


###################################################
### code chunk number 188: CM17-plot (eval = FALSE)
###################################################
## blues <- c(
##     "#bef7ff", "#b4ecfd", "#aae1fb", "#a0d6f9", "#96cbf7",
##     "#8cc1f5", "#82b6f3", "#78abf1", "#6ea0ef", "#6495ed"
## )
## blues <- blues[length(blues):1]
## par(mar = c(6, 5.5, 3, 1), mgp = c(4, 1, 0))
## plot(m.mtcars.subsets, scale = "adjr2", col = blues,
##      main = "Fuel economy: best subsets")


###################################################
### code chunk number 189: CM17-refit (eval = FALSE)
###################################################
## m.mtcars <- st.lmm(
##     mpg ~ disp + hp + wt + qsec + am,
##     data = mtcars,
##     method = "moments"
## )


###################################################
### code chunk number 190: CM17-fit (eval = FALSE)
###################################################
## m.mtcars


###################################################
### code chunk number 191: CM13-1 (eval = FALSE)
###################################################
## data("v-lmm-tomato", package = "SelectionTools")
## m.tomato <- st.lmm(
##     y ~ 0 + x.v,
##     data = tomato,
##     method = "moments"
## )
## m.tomato


###################################################
### code chunk number 192: CM13-2 (eval = FALSE)
###################################################
## L.tomato <- matrix(
##         c(1, -1, 0), nrow = 1,
##         dimnames = list("group 1 - group 2", names(coef(m.tomato)))
##     )
## st.lsm(
##     m.tomato,
##     L = L.tomato,
##     joint = TRUE
## )


###################################################
### code chunk number 193: CM13-2-b (eval = FALSE)
###################################################
## st.lsm(m.tomato, ~x.v)


###################################################
### code chunk number 194: CM13-3 (eval = FALSE)
###################################################
## # Illustrative null: each difference from group 1 equals five response units.
## C.tomato <- rbind(
##     group2_group1 = c(-1, 1, 0),
##     group3_group1 = c(-1, 0, 1)
## )
## st.lsm(
##     m.tomato,
##     ~x.v,
##     method = C.tomato,
##     rhs = c(5, 5),
##     adjust = "none"
## )


###################################################
### code chunk number 195: CM13-4 (eval = FALSE)
###################################################
## # L columns refer to fitted coefficients; method columns refer to means.
## # This no-intercept cell-means model makes the two matrices identical.
## L.tomato <- C.tomato
## colnames(L.tomato) <- names(coef(m.tomato))
## h.tomato.L <- st.lsm(
##     m.tomato,
##     L = L.tomato,
##     rhs = c(5, 5),
##     adjust = "none"
## )
## h.tomato.L


###################################################
### code chunk number 196: CM13-5 (eval = FALSE)
###################################################
## # Jointly test both five-unit null differences.
## st.lsm(
##     h.tomato.L,
##     joint = TRUE
## )


###################################################
### code chunk number 197: M52-1 (eval = FALSE)
###################################################
## data("v-lmm-chickens", package = "SelectionTools")
## chickens <- subset(chickens, Time == 21)
## m <- st.lmm(
##     fixed = weight ~ Diet,
##     data = chickens
## )
## m


###################################################
### code chunk number 198: M52-2 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     "Diet",
##     method = "pairwise",
##     side = "two.sided",
##     adjust = "tukey"
## )


###################################################
### code chunk number 199: M53-1 (eval = FALSE)
###################################################
## data("v-lmm-chickens", package = "SelectionTools")
## chickens$weight.s <- sqrt(chickens$weight)
## st.lmm(
##     fixed = weight.s ~ Time:Diet,
##     data = chickens
## )


###################################################
### code chunk number 200: M54-1 (eval = FALSE)
###################################################
## data("v-lmm-chickens", package = "SelectionTools")
## chickens$weight.s <- sqrt(chickens$weight)
## st.lmm(
##     fixed = weight.s ~ Diet:Time,
##     random = ~Chick,
##     data = chickens
## )


###################################################
### code chunk number 201: M55-1 (eval = FALSE)
###################################################
## data("v-lmm-chickens", package = "SelectionTools")
## chickens$weight.s <- sqrt(chickens$weight)
## st.lmm(
##     fixed = weight.s ~ 0 + Diet + Diet:Time,
##     random = ~Chick,
##     data = chickens
## )


###################################################
### code chunk number 202: M56-1 (eval = FALSE)
###################################################
## data("v-lmm-chickens", package = "SelectionTools")
## chickens$weight.s <- sqrt(chickens$weight)
## st.lmm(
##     fixed = weight.s ~ Diet:Time,
##     random = ~Time | Chick,
##     data = chickens
## )


###################################################
### code chunk number 203: M57-1 (eval = FALSE)
###################################################
## data("v-lmm-chickens", package = "SelectionTools")
## chickens$weight.s <- sqrt(chickens$weight)
## st.lmm(
##     fixed = weight.s ~ 0 + Diet + Diet:Time,
##     random = ~Time | Chick,
##     data = chickens
## )


###################################################
### code chunk number 204: M58-1 (eval = FALSE)
###################################################
## data("v-lmm-chickens", package = "SelectionTools")
## chickens$weight.s <- sqrt(chickens$weight)
## st.lmm(
##     fixed = weight.s ~ Diet:Time,
##     random = ~Time | Chick,
##     cov = ~at(Time):units,
##     data = chickens
## )


###################################################
### code chunk number 205: M59-1 (eval = FALSE)
###################################################
## data("v-lmm-chickens", package = "SelectionTools")
## chickens$weight.s <- sqrt(chickens$weight)
## m <- st.lmm(
##     fixed = weight.s ~ Diet:Time,
##     random = ~0 + Time | Chick,
##     cov = ~at(Time):units,
##     data = chickens
## )
## m


###################################################
### code chunk number 206: M59-2 (eval = FALSE)
###################################################
## st.lsm(
##     m,
##     "Diet",
##     trend = "Time",
##     method = "trt.vs.ctrl",
##     ref = "d1",
##     side = ">",
##     adjust = "fdr"
## )


###################################################
### code chunk number 207: M59-3 (eval = FALSE)
###################################################
## # Illustrative minimum improvement: 0.05 sqrt(weight) units per day.
## # The object retains the null, direction, adjustment and confidence level.
## chicken.minimum <- st.lsm(
##     m,
##     "Diet",
##     trend = "Time",
##     method = "trt.vs.ctrl",
##     ref = "d1",
##     rhs = .05,
##     side = ">",
##     adjust = "single-step",
##     level = .90
## )
## chicken.minimum


###################################################
### code chunk number 208: M59-4 (eval = FALSE)
###################################################
## # Null values retained for the three slope comparisons.
## coef(
##     chicken.minimum,
##     rhs = TRUE
## )


###################################################
### code chunk number 209: E01-data (eval = FALSE)
###################################################
## data("v-lmm-halfdiallel", package = "SelectionTools")
## head(halfdiallel, 5)


###################################################
### code chunk number 210: E01-fit (eval = FALSE)
###################################################
## m.diallel <- st.lmm(
##     fixed = sugar ~ 1,
##     random = ~overlay(female, male) + cross,
##     data = halfdiallel
## )
## m.diallel


###################################################
### code chunk number 211: E01-parent (eval = FALSE)
###################################################
## diallel.parents <- st.predict(
##     m.diallel,
##     "0 + overlay(female, male)"
## )
## diallel.parents


###################################################
### code chunk number 212: E01-sca (eval = FALSE)
###################################################
## st.predict(m.diallel, "0 + cross")


###################################################
### code chunk number 213: E01-cross (eval = FALSE)
###################################################
## st.predict(m.diallel, "male + female + cross")


###################################################
### code chunk number 214: tropmaize-field (eval = FALSE)
###################################################
## data("v-lmm-tropmaize", package = "SelectionTools")
## maize.field <- st.lmm(
##     fixed = y ~ gen,
##     random = ~block,
##     data = tropmaize
## )
## maize.field


###################################################
### code chunk number 215: tropmaize-blue (eval = FALSE)
###################################################
## st.lsm(maize.field, "gen")


###################################################
### code chunk number 216: tropmaize-blupfit (eval = FALSE)
###################################################
## maize.random <- st.lmm(
##     fixed = y ~ 1,
##     random = ~gen + block,
##     data = tropmaize
## )
## maize.random


###################################################
### code chunk number 217: tropmaize-blup (eval = FALSE)
###################################################
## maize.blup <- st.predict(maize.random, "gen")
## maize.blup


###################################################
### code chunk number 218: tropmaize-matrix (eval = FALSE)
###################################################
## data("v-tropmaize-vcf", package = "SelectionTools")
## st.load.vcf.data(v.tropmaize.vcf, data.set = "maize.lmm")


###################################################
### code chunk number 219: tropmaize-matrix-b (eval = FALSE)
###################################################
## maize.Z <- gs.get.Z("maize.lmm", names = TRUE)
## maize.A <- gs.make.A(maize.Z)


###################################################
### code chunk number 220: tropmaize-genomic (eval = FALSE)
###################################################
## maize.genomic <- st.lmm(
##     fixed = y ~ 1,
##     random = ~gen + block,
##     ucov = list(gen = maize.A),
##     data = tropmaize
## )
## maize.genomic


###################################################
### code chunk number 221: tropmaize-gblup (eval = FALSE)
###################################################
## maize.gblup <- st.predict(maize.genomic, "gen")
## maize.gblup


###################################################
### code chunk number 222: tropmaize-series-setup (eval = FALSE)
###################################################
## library(SelectionTools)
## st.set.num.threads(1)
## data("v-lmm-tropmaize-lattice", package = "SelectionTools")
## data("v-tropmaize-vcf", package = "SelectionTools")
## st.load.vcf.data(v.tropmaize.vcf, data.set = "maize.series")
## envs <- unique(tropmaize.lattice$env)
## G <- st.get.A("maize.series")


###################################################
### code chunk number 223: tropmaize-series-env1 (eval = FALSE)
###################################################
## fit.env1 <- st.lmm(y ~ gen, random = ~ rep + rep:block,
##     data = droplevels(subset(tropmaize.lattice, env == envs[1])))


###################################################
### code chunk number 224: tropmaize-series-env2 (eval = FALSE)
###################################################
## fit.env2 <- st.lmm(y ~ gen, random = ~ rep + rep:block,
##     data = droplevels(subset(tropmaize.lattice, env == envs[2])))


###################################################
### code chunk number 225: tropmaize-series-env3 (eval = FALSE)
###################################################
## fit.env3 <- st.lmm(y ~ gen, random = ~ rep + rep:block,
##     data = droplevels(subset(tropmaize.lattice, env == envs[3])))


###################################################
### code chunk number 226: tropmaize-series-fit (eval = FALSE)
###################################################
## fit.series.G <- st.lmm(y ~ env,
##     random = ~ gen + gen:env,
##     ucov = list(gen = G),
##     series = list("env", fit.env1, fit.env2, fit.env3))
## fit.series.G


###################################################
### code chunk number 227: barley-series-setup (eval = FALSE)
###################################################
## library(SelectionTools)
## data("v-lmm-barley-met", package = "SelectionTools")
## barley.met$env <- interaction(barley.met$loc, barley.met$year,
##     drop = TRUE)


###################################################
### code chunk number 228: barley-series-first (eval = FALSE)
###################################################
## barley.trials <- split(barley.met, barley.met$env)
## barley.first <- lapply(barley.trials, function(trial) {
##     st.lmm(TKW ~ gen, random = ~ rep,
##         data = droplevels(trial))
## })


###################################################
### code chunk number 229: barley-series-fit (eval = FALSE)
###################################################
## barley.series <- st.lmm(TKW ~ env,
##     random = ~ gen + gen:loc + gen:year + gen:env,
##     series = c(list("env"), barley.first))
## barley.series


###################################################
### code chunk number 230: barley-series-predict (eval = FALSE)
###################################################
## barley.prediction <- st.predict(barley.series, "gen")
## barley.prediction


###################################################
### code chunk number 231: E03-load (eval = FALSE)
###################################################
## data("v-rapeseed", package = "SelectionTools")
## st.set.num.threads(1)
## st.load.vcf.data(v.rapeseed.vcf,
##                  data.set = "rapeseed")


###################################################
### code chunk number 232: E03-phenotypes (eval = FALSE)
###################################################
## head(v.rapeseed.phe)


###################################################
### code chunk number 233: E03-k-p3d (eval = FALSE)
###################################################
## rapeseed.k <- st.gwas(
##     fixed = flt ~ 1,
##     random = ~gen,
##     data = v.rapeseed.phe,
##     marker.data.set = "rapeseed",
##     K = "A",
##     min.MAF = 0.05,
##     P3D = TRUE)
## print(rapeseed.k, top = 3, bottom = 0)


###################################################
### code chunk number 234: E03-figure-k-manhattan (eval = FALSE)
###################################################
## plot(rapeseed.k, which = "manhattan",
##      main = "Flowering time: K model", threshold = FALSE,
##      cex = 0.42)


###################################################
### code chunk number 235: E03-k-reml (eval = FALSE)
###################################################
## rapeseed.k.reml <- st.gwas(
##     fixed = flt ~ 1,
##     random = ~gen,
##     data = v.rapeseed.phe,
##     marker.data.set = "rapeseed",
##     K = "A",
##     P3D = FALSE)
## print(rapeseed.k.reml, top = 3, bottom = 0)


###################################################
### code chunk number 236: E03-pk-three (eval = FALSE)
###################################################
## rapeseed.pk3 <- st.gwas(
##     fixed = flt ~ 1,
##     random = ~gen,
##     data = v.rapeseed.phe,
##     marker.data.set = "rapeseed",
##     K = "A",
##     n.PC = 3,
##     P3D = TRUE)


###################################################
### code chunk number 237: E03-model-summary (eval = FALSE)
###################################################
## compare.gwas.fits(
##     rapeseed.k,
##     rapeseed.pk3,
##     adjust = "FDR",
##     alpha = 0.05
## )


###################################################
### code chunk number 238: E03-figure-pk-qq (eval = FALSE)
###################################################
## par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.2, 0.7),
##     mgp = c(2.7, 0.8, 0), cex.axis = 0.85)
## plot(rapeseed.k, which = "qq", main = "K", level = 0.95,
##      col = "cornflowerblue", cex = 0.38)
## plot(rapeseed.pk3, which = "qq", main = "PK: 3 PCs", level = 0.95,
##      col = "cornflowerblue", cex = 0.38)


###################################################
### code chunk number 239: E03-hits (eval = FALSE)
###################################################
## rapeseed.summary <- summarize.gwas.fit(
##     rapeseed.pk3,
##     adjust = "FDR",
##     alpha = 0.05
## )
## rapeseed.summary


###################################################
### code chunk number 240: E03-simulated-threshold (eval = FALSE)
###################################################
## rapeseed.pk3.sim <- st.gwas(
##     fixed = flt ~ 1,
##     random = ~gen,
##     data = v.rapeseed.phe,
##     marker.data.set = "rapeseed",
##     K = "A",
##     n.PC = 3,
##     P3D = TRUE,
##     control = list(batch.size = 1024L),
##     sim.threshold = list(method = "parametricFixedVC",
##                          nsim = 999L, alpha = 0.05, seed = 1L))
## rapeseed.pk3.sim


###################################################
### code chunk number 241: E03-simulation-hits (eval = FALSE)
###################################################
## summarize.gwas.simulation(rapeseed.pk3.sim)


###################################################
### code chunk number 242: E03-figure-pk-manhattan (eval = FALSE)
###################################################
## plot(rapeseed.pk3.sim, which = "manhattan", adjust = "simulation",
##      main = "Flowering time: K + 3 PCs",
##      cex = 0.42)


###################################################
### code chunk number 243: E03-region (eval = FALSE)
###################################################
## rapeseed.lead <- rapeseed.summary$markers[1, ]


###################################################
### code chunk number 244: E03-figure-region (eval = FALSE)
###################################################
## plot(rapeseed.pk3, which = "region",
##      chromosome = rapeseed.lead$chromosome,
##      xlim = rapeseed.lead$position + c(-2e6, 2e6),
##      position.scale = 1e6, position.unit = "Mb",
##      shade = rapeseed.lead$position + c(-250000, 250000),
##      highlight = rapeseed.lead$marker, labels = TRUE,
##      threshold = 1e-5, main = "Flowering time: K + 3 PCs",
##      col = "cornflowerblue", cex = 0.65)


###################################################
### code chunk number 245: tropmaize-prediction-load (eval = FALSE)
###################################################
## library(SelectionTools)
## data("v-lmm-tropmaize", package = "SelectionTools")
## data("v-tropmaize-vcf", package = "SelectionTools")
## st.load.vcf.data(v.tropmaize.vcf, data.set = "maize.lmm")


###################################################
### code chunk number 246: tropmaize-prediction-matrix (eval = FALSE)
###################################################
## maize.A <- st.get.A("maize.lmm")


###################################################
### code chunk number 247: tropmaize-mask (eval = FALSE)
###################################################
## maize.genotypes <- levels(tropmaize$gen)
## maize.masked <- tropmaize
## maize.masked$y[maize.masked$gen %in% maize.genotypes[1:10]] <- NA


###################################################
### code chunk number 248: tropmaize-two-step-field (eval = FALSE)
###################################################
## maize.train.field <- st.lmm(
##     fixed = y ~ gen,
##     random = ~block,
##     data = maize.masked
## )


###################################################
### code chunk number 249: tropmaize-two-step-means (eval = FALSE)
###################################################
## maize.means <- as.data.frame(st.lsm(maize.train.field, "gen"))
## maize.means <- maize.means[, c("gen", "lsmeans")]


###################################################
### code chunk number 250: tropmaize-two-step-fit (eval = FALSE)
###################################################
## maize.two.step <- st.lmm(
##     fixed = lsmeans ~ 1,
##     random = ~gen,
##     ucov = list(gen = maize.A),
##     data = maize.means
## )
## maize.two.step


###################################################
### code chunk number 251: tropmaize-two-step-predict (eval = FALSE)
###################################################
## maize.two.pred <- st.predict(maize.two.step, "gen")
## maize.two.pred


###################################################
### code chunk number 252: tropmaize-split-load (eval = FALSE)
###################################################
## library(SelectionTools)
## data("v-tropmaize-vcf", package = "SelectionTools")
## st.load.vcf.data(v.tropmaize.vcf, data.set = "maize.split")


###################################################
### code chunk number 253: tropmaize-split-statistics (eval = FALSE)
###################################################
## maize.marker.statistics <- st.marker.data.statistics(
##     data.set = "maize.split",
##     mar = FALSE,
##     gen = FALSE
## )


###################################################
### code chunk number 254: tropmaize-split-ids (eval = FALSE)
###################################################
## maize.individuals <- maize.marker.statistics$individual.list$Name
## maize.prediction.ids <- maize.individuals[1:10]
## maize.training.ids <- maize.individuals[-(1:10)]


###################################################
### code chunk number 255: tropmaize-split-training-markers (eval = FALSE)
###################################################
## st.copy.marker.data("maize.training", "maize.split")
## st.restrict.marker.data(
##     ind.list = maize.training.ids,
##     data.set = "maize.training"
## )


###################################################
### code chunk number 256: tropmaize-split-prediction-markers (eval = FALSE)
###################################################
## st.copy.marker.data("maize.prediction", "maize.split")
## st.restrict.marker.data(
##     ind.list = maize.prediction.ids,
##     data.set = "maize.prediction"
## )


###################################################
### code chunk number 257: tropmaize-split-performance (eval = FALSE)
###################################################
## data("v-tropmaize-phe", package = "SelectionTools")
## st.load.performance.data(
##     v.tropmaize.phe,
##     data.set = "maize.training"
## )


###################################################
### code chunk number 258: tropmaize-split-fit (eval = FALSE)
###################################################
## maize.training.model <- st.lmm(
##     fixed = y ~ 1,
##     random = ~gen,
##     ucov = list(gen = "A"),
##     data = "maize.training"
## )
## maize.training.model


###################################################
### code chunk number 259: tropmaize-split-predict (eval = FALSE)
###################################################
## maize.predictions <- st.predict(
##     maize.training.model,
##     newdata = "maize.prediction"
## )
## maize.predictions


###################################################
### code chunk number 260: tropmaize-lattice-mask (eval = FALSE)
###################################################
## data("v-lmm-tropmaize-lattice", package = "SelectionTools")
## data("v-tropmaize-vcf", package = "SelectionTools")
## st.load.vcf.data(v.tropmaize.vcf, data.set = "maize.lattice")


###################################################
### code chunk number 261: tropmaize-lattice-matrix (eval = FALSE)
###################################################
## maize.lattice.A <- st.get.A("maize.lattice")
## maize.lattice.missing.gen <- levels(tropmaize.lattice$gen)[seq_len(10)]


###################################################
### code chunk number 262: tropmaize-lattice-mask-b (eval = FALSE)
###################################################
## maize.lattice.masked <- tropmaize.lattice
## maize.lattice.masked$y[
##     maize.lattice.masked$gen %in% maize.lattice.missing.gen
## ] <- NA_real_


###################################################
### code chunk number 263: tropmaize-lattice-one-step-fit (eval = FALSE)
###################################################
## maize.lattice.one.step <- st.lmm(
##     fixed = y ~ env,
##     random = ~gen + env:rep + env:rep:block + gen:env,
##     ucov = list(gen = maize.lattice.A),
##     data = maize.lattice.masked
## )
## maize.lattice.one.step


###################################################
### code chunk number 264: tropmaize-lattice-one-step-predict (eval = FALSE)
###################################################
## maize.lattice.one.pred <- st.predict(maize.lattice.one.step, "gen")
## maize.lattice.one.pred


###################################################
### code chunk number 265: H01-load (eval = FALSE)
###################################################
## data("v-flintmaize-vcf", package = "SelectionTools")
## data("v-dentmaize-vcf", package = "SelectionTools")
## data("v-hybmaize-phe", package = "SelectionTools")
## st.load.vcf.data(v.flintmaize.vcf, data.set = "flintmaize")


###################################################
### code chunk number 266: H01-load-b (eval = FALSE)
###################################################
## st.load.vcf.data(v.dentmaize.vcf, data.set = "dentmaize")


###################################################
### code chunk number 267: H01-A (eval = FALSE)
###################################################
## st.restrict.marker.data(ExHet.MIN = 1e-8, NoAll.MAX = 2,
##                         data.set = "flintmaize")


###################################################
### code chunk number 268: H01-A-b (eval = FALSE)
###################################################
## st.restrict.marker.data(ExHet.MIN = 1e-8, NoAll.MAX = 2,
##                         data.set = "dentmaize")


###################################################
### code chunk number 269: H01-A-c (eval = FALSE)
###################################################
## Aflint <- st.get.A("flintmaize")
## Adent <- st.get.A("dentmaize")


###################################################
### code chunk number 270: H01-SCA (eval = FALSE)
###################################################
## hybmaize <- v.hybmaize.phe
## hybmaize$flint <- factor(hybmaize$flint, levels = rownames(Aflint))
## hybmaize$dent <- factor(hybmaize$dent, levels = rownames(Adent))
## flint <- as.character(hybmaize$flint)
## dent <- as.character(hybmaize$dent)


###################################################
### code chunk number 271: H01-SCA-matrix (eval = FALSE)
###################################################
## Asca <- Aflint[flint, flint] * Adent[dent, dent]
## hybrid.names <- paste(flint, dent, sep = ":")
## dimnames(Asca) <- list(hybrid.names, hybrid.names)
## hybrid.covariance.dimensions <- dim(Asca)
## hybrid.covariance.dimensions


###################################################
### code chunk number 272: H01-mask (eval = FALSE)
###################################################
## hybmaize[seq_len(10), c("GY", "GM")] <- NA_real_
## hybrids.without.responses <- hybmaize[seq_len(10), c("hybrid", "GY")]
## hybrids.without.responses


###################################################
### code chunk number 273: H01-fit (eval = FALSE)
###################################################
## hybrid.fit <- st.lmm(
##     fixed = GY ~ 1,
##     random = ~flint + dent + flint:dent,
##     ucov = list(
##         flint = Aflint,
##         dent = Adent,
##         "flint:dent" = Asca
##     ),
##     data = hybmaize
## )
## hybrid.fit


###################################################
### code chunk number 274: H01-predict (eval = FALSE)
###################################################
## hybrid.predicted <- predict(
##     hybrid.fit,
##     newdata = hybmaize,
##     type = "conditional"
## )$fit
## hybrid.predicted


###################################################
### code chunk number 275: H01-predict-b (eval = FALSE)
###################################################
## hybrid.missing.phenotypes <- sum(is.na(hybmaize$GY))
## hybrid.missing.phenotypes


###################################################
### code chunk number 276: H01-missing-predictions (eval = FALSE)
###################################################
## hybrid.missing.predictions <- sum(is.na(hybrid.predicted))
## hybrid.missing.predictions


###################################################
### code chunk number 277: tropmaize-dataset-load (eval = FALSE)
###################################################
## data("v-tropmaize-vcf", package = "SelectionTools")
## data("v-tropmaize-phe", package = "SelectionTools")
## st.load.vcf.data(v.tropmaize.vcf, data.set = "maize")


###################################################
### code chunk number 278: tropmaize-dataset-load-b (eval = FALSE)
###################################################
## v.tropmaize.phe$yy[
##     v.tropmaize.phe$ii %in% v.tropmaize.vcf$individual[seq_len(10)]
## ] <- NA_real_
## st.load.performance.data(v.tropmaize.phe, data.set = "maize")


###################################################
### code chunk number 279: tropmaize-dataset-mask (eval = FALSE)
###################################################
## maize.y <- gs.get.y("maize", names = TRUE)
## maize.missing.responses <- head(maize.y, 10)
## maize.missing.responses


###################################################
### code chunk number 280: tropmaize-dataset-fit (eval = FALSE)
###################################################
## maize.two.step <- st.lmm(
##     fixed = y ~ 1,
##     random = ~gen,
##     ucov = list(gen = "A"),
##     data = "maize"
## )
## maize.two.step


###################################################
### code chunk number 281: tropmaize-dataset-predict (eval = FALSE)
###################################################
## maize.two.pred <- st.predict(maize.two.step, "gen")
## maize.two.pred


###################################################
### code chunk number 282: tropmaize-dataset-complete (eval = FALSE)
###################################################
## maize.marker.genotypes <- length(maize.y)
## maize.marker.genotypes


###################################################
### code chunk number 283: tropmaize-dataset-missing-phenotypes (eval = FALSE)
###################################################
## maize.missing.phenotypes <- sum(is.na(maize.y))
## maize.missing.phenotypes


###################################################
### code chunk number 284: tropmaize-dataset-missing-predictions (eval = FALSE)
###################################################
## maize.missing.predictions <- sum(is.na(coef(maize.two.pred)))
## maize.missing.predictions


###################################################
### code chunk number 285: H02-load (eval = FALSE)
###################################################
## data("v-flintmaize-vcf", package = "SelectionTools")
## data("v-dentmaize-vcf", package = "SelectionTools")
## data("v-hybmaize-phe", package = "SelectionTools")
## st.load.vcf.data(v.flintmaize.vcf, data.set = "flint")


###################################################
### code chunk number 286: H02-load-b (eval = FALSE)
###################################################
## st.load.vcf.data(v.dentmaize.vcf, data.set = "dent")


###################################################
### code chunk number 287: H02-load-c (eval = FALSE)
###################################################
## st.restrict.marker.data(ExHet.MIN = 1e-8, NoAll.MAX = 2,
##                         data.set = "flint")


###################################################
### code chunk number 288: H02-load-d (eval = FALSE)
###################################################
## st.restrict.marker.data(ExHet.MIN = 1e-8, NoAll.MAX = 2,
##                         data.set = "dent")


###################################################
### code chunk number 289: H02-load-e (eval = FALSE)
###################################################
## hybrid.performance <- v.hybmaize.phe
## hybrid.performance[seq_len(10), c("GY", "GM")] <- NA_real_


###################################################
### code chunk number 290: H02-store (eval = FALSE)
###################################################
## st.hybrid.data.set(
##     markerdata1 = "flint",
##     markerdata2 = "dent",
##     hybrid.performance = hybrid.performance,
##     hybrid.data.set = "maizehybrids"
## )


###################################################
### code chunk number 291: H02-missing-rows (eval = FALSE)
###################################################
## hybrid.missing.rows <- hybrid.performance[seq_len(10), c("hybrid", "GY")]
## hybrid.missing.rows


###################################################
### code chunk number 292: H02-fit (eval = FALSE)
###################################################
## hybrid.gca <- st.lmm(
##     fixed = GY ~ 1,
##     random = ~flint + dent,
##     ucov = list("flint" = "A", "dent" = "A"),
##     data = "maizehybrids",
##     prediction = "Full"
## )
## hybrid.gca


###################################################
### code chunk number 293: H02-predict (eval = FALSE)
###################################################
## gca.predictions <- predict(
##     hybrid.gca,
##     newdata = hybrid.performance,
##     type = "conditional"
## )$fit
## gca.predictions


###################################################
### code chunk number 294: H02-predict-complete (eval = FALSE)
###################################################
## gca.crosses <- st.lmm.components(hybrid.gca, "predictions")
## gca.crosses.preview <- head(gca.crosses[, c("flint", "dent", "prediction")])
## gca.crosses.preview


###################################################
### code chunk number 295: H02-cross-count (eval = FALSE)
###################################################
## gca.cross.count <- nrow(gca.crosses)
## gca.cross.count


###################################################
### code chunk number 296: H02-missing-phenotypes (eval = FALSE)
###################################################
## gca.missing.phenotypes <- sum(is.na(hybrid.performance$GY))
## gca.missing.phenotypes


###################################################
### code chunk number 297: H02-missing-predictions (eval = FALSE)
###################################################
## gca.missing.predictions <- sum(is.na(gca.predictions))
## gca.missing.predictions


###################################################
### code chunk number 298: H03-fit (eval = FALSE)
###################################################
## hybrid.gca.sca <- st.lmm(
##     fixed = GY ~ 1,
##     random = ~flint + dent + flint:dent,
##     ucov = list(
##         "flint" = "A",
##         "dent" = "A",
##         "flint:dent" = "kro(a,a)"
##     ),
##     data = "maizehybrids",
##     prediction = "NAonly"
## )
## hybrid.gca.sca


###################################################
### code chunk number 299: H03-predict (eval = FALSE)
###################################################
## sca.predictions <- st.lmm.components(hybrid.gca.sca, "predictions")
## sca.predictions


###################################################
### code chunk number 300: H03-predict-b (eval = FALSE)
###################################################
## sca.missing.phenotypes <- sum(is.na(hybrid.performance$GY))
## sca.missing.phenotypes


###################################################
### code chunk number 301: H03-missing-predictions (eval = FALSE)
###################################################
## sca.missing.predictions <- sum(is.na(sca.predictions$prediction))
## sca.missing.predictions


###################################################
### code chunk number 302: H03-list (eval = FALSE)
###################################################
## hybrid.status <- st.hybrid.data.set.list()
## hybrid.status.summary <- hybrid.status[
##     , c("hybrid.data.set", "markerdata1", "markerdata2",
##         "no.records", "no.hybrids", "SCA.matrix")
## ]
## hybrid.status.summary


###################################################
### code chunk number 303: H04-load (eval = FALSE)
###################################################
## data("v-flintmaize2-vcf")
## data("v-dentmaize2-vcf")
## data("v-hybmaize2-phe")


###################################################
### code chunk number 304: H04-load-remove (eval = FALSE)
###################################################
## st.hybrid.data.set.remove("maizehybrids")
## st.hybrid.data.set.remove("maizeTRUE")
## st.hybrid.data.set.remove("maizeTS")


###################################################
### code chunk number 305: H04-load-markers (eval = FALSE)
###################################################
## st.load.vcf.data(v.flintmaize2.vcf, data.set = "flint")
## st.load.vcf.data(v.dentmaize2.vcf, data.set = "dent")


###################################################
### code chunk number 306: H04-hybrid-data (eval = FALSE)
###################################################
## st.hybrid.data.set("flint", "dent", v.hybmaize2.phe, "maizeTRUE")


###################################################
### code chunk number 307: H04-fit-real (eval = FALSE)
###################################################
## maize.founder.model <- st.lmm(
##     GY ~ 1,
##     random = ~flint + dent + flint:dent,
##     ucov = list(flint = "A", dent = "A", "flint:dent" = "kro(a,a)"),
##     data = "maizeTRUE"
## )
## maize.founder.model


###################################################
### code chunk number 308: H04-simulate (eval = FALSE)
###################################################
## n.markers <- 9000L
## n.dh      <- 100L
## n.hybrids <- 200L
## h2        <- 0.8


###################################################
### code chunk number 309: H04-simulate-initialize (eval = FALSE)
###################################################
## st.set.num.threads(1)
## rng.init(seed = 20261024L)
## reset.all()


###################################################
### code chunk number 310: H04-simulate-map (eval = FALSE)
###################################################
## st.set.sim.gp("flint")
## st.set.sim.mp("flint")
## st.set.sim.pp("simFlint", "flint")
## st.set.sim.pp("simDent", "dent")
## extend.marker.map(n.markers)


###################################################
### code chunk number 311: H04-simulate-dh (eval = FALSE)
###################################################
## cross("simFlintF1", "simFlint", "simFlint", n.dh, mode = 0)
## cross("simDentF1", "simDent", "simDent", n.dh, mode = 0)
## dh("simflintDH", "simFlintF1", n.dh)
## dh("simdentDH", "simDentF1", n.dh)


###################################################
### code chunk number 312: H04-simulate-export (eval = FALSE)
###################################################
## st.get.simpop("simflintDH", data.set = "TSflint")
## st.get.simpop("simdentDH", data.set = "TSdent")


###################################################
### code chunk number 313: H04-pedigree (eval = FALSE)
###################################################
## flint.pedigree <- head(
##     population.pedigree("TSflint", source = "marker"),
##     3
## )
## flint.pedigree


###################################################
### code chunk number 314: H04-crosses (eval = FALSE)
###################################################
## hybrid.design <- population.cross(
##     "simHybrids", "simflintDH", "simdentDH", n.hybrids
## )
## names(hybrid.design) <- c("hybrid", "TSflint", "TSdent")
## st.hybrid.data.set("TSflint", "TSdent", hybrid.design, "maizeTS")


###################################################
### code chunk number 315: H04-crosses-phenotypes (eval = FALSE)
###################################################
## hybrid.predictions <- as.data.frame(
##     st.predict("maizeTRUE", newdata = "maizeTS")
## )
## st.generate.phenotypes("maizeTS", hybrid.predictions, h2 = h2, trait = "GY")


###################################################
### code chunk number 316: H04-fit-sim (eval = FALSE)
###################################################
## maize.simulated.model <- st.lmm(
##     GY ~ 1,
##     random = ~TSflint + TSdent + TSflint:TSdent,
##     ucov = list(
##         TSflint = "A", TSdent = "A",
##         "TSflint:TSdent" = "kro(a,a)"
##     ),
##     data = "maizeTS"
## )
## maize.simulated.model


###################################################
### code chunk number 317: H04-predict (eval = FALSE)
###################################################
## genotypic.values <- coef(st.predict("maizeTS"))
## first.genotypic.values <- head(genotypic.values)
## first.genotypic.values


