\name{v-lmm-data}
\alias{v-lmm-data}
\alias{v-lmm-peas}
\alias{peas}
\alias{v-lmm-sugar}
\alias{sugar}
\alias{v-lmm-cotton}
\alias{cotton}
\alias{v-lmm-fertilizer}
\alias{fertilizer}
\alias{v-lmm-soils}
\alias{v-lmm-rapeseed}
\alias{v-lmm-soybean-blocks}
\alias{v-lmm-alfalfa}
\alias{v-lmm-oats}
\alias{v-lmm-barley-met}
\alias{barley.met}
\alias{v-lmm-barley-lattice}
\alias{v-lmm-rapeseed-blocks}
\alias{v-lmm-maize-counties}
\alias{v-lmm-soybean}
\alias{v-lmm-soybean-met}
\alias{v-lmm-material-groups}
\alias{v-lmm-maize-met}
\alias{v-lmm-met-unbalanced}
\alias{v-lmm-met-missing}
\alias{v-lmm-met-seven-traits}
\alias{v-lmm-barley-seven-traits}
\alias{v-lmm-chickens}
\alias{soils}
\alias{rapeseed}
\alias{soybean.blocks}
\alias{alfalfa}
\alias{oats}
\alias{barley.lattice}
\alias{rapeseed.blocks}
\alias{maize.counties}
\alias{soybean}
\alias{soybean.met}
\alias{material.groups}
\alias{maize.met}
\alias{met.unbalanced}
\alias{met.missing}
\alias{met.seven.traits}
\alias{barley.seven.traits}
\alias{chickens}
\docType{data}
\title{Datasets for Linear Mixed Models}
\description{Agricultural and repeated-measures datasets for linear mixed
models. Dataset entries describe the loaded objects, variables, storage types,
factor levels and missing values.
The additional datasets for moment estimators are documented separately in
\code{\link{v-lmm-climate}}, \code{\link{v-lmm-corn}},
\code{\link{v-lmm-iris-two}}, \code{\link{v-lmm-iris}},
\code{\link{v-lmm-fertilizer-regression}},
\code{\link{v-lmm-fertilizer-replicated}}, \code{\link{v-lmm-rats}}
and \code{\link{v-lmm-tomato}}.
The machine and operator experiment is documented in \code{\link{v-lmm-machines}}.
The additional regression dataset \code{v-lmm-mtcars} is documented in its
own help topic and loads the numeric data frame \code{mtcars}.}
\usage{data("v-lmm-oats", package = "SelectionTools")}
\format{Each file contains one \code{data.frame}. Row counts, storage types
and factor levels describe the stored data. Column names are case-sensitive.
\describe{
\item{\code{v-lmm-soils}}{Yield observations from fields within villages and regions.
File \file{v-lmm-soils.rda} loads \code{soils},
with 24 rows and 4 columns.
The stored data contain no missing values.
\describe{
\item{\code{region}}{Region identifier. Factor with 3 levels: \code{I}, \code{II}, \code{III}.}
\item{\code{village}}{Village identifier. Factor with 6 levels: \code{A} to \code{F}. Villages \code{A} and \code{B} belong to region \code{I}, \code{C} and \code{D} to region \code{II}, and \code{E} and \code{F} to region \code{III}.}
\item{\code{field}}{Field identifier. Factor with 4 levels: \code{f1}, \code{f2}, \code{f3}, \code{f4}.}
\item{\code{yield}}{Yield response. Integer.}
}}
\item{\code{v-lmm-rapeseed}}{Rapeseed height observations across environments.
File \file{v-lmm-rapeseed.rda} loads \code{rapeseed},
with 120 rows and 4 columns.
The stored data contain no missing values.
\describe{
\item{\code{geno}}{Genotype identifier. Factor with 20 levels.}
\item{\code{env}}{Environment identifier. Factor with 2 levels: \code{e1}, \code{e2}.}
\item{\code{rep}}{Replicate identifier. Factor with 3 levels: \code{r1}, \code{r2}, \code{r3}.}
\item{\code{hgt}}{Height response. Integer.}
}}
\item{\code{v-lmm-soybean-blocks}}{Soybean yield observations in an incomplete-block layout.
File \file{v-lmm-soybean-blocks.rda} loads \code{soybean.blocks},
with 50 rows and 5 columns.
The stored data contain no missing values.
\describe{
\item{\code{Rep}}{Replicate identifier. Factor with 2 levels: \code{r1}, \code{r2}.}
\item{\code{Block}}{Block identifier. Factor with 10 levels.}
\item{\code{Plot}}{Plot identifier. Factor with 5 levels: \code{p1}, \code{p2}, \code{p3}, \code{p4}, \code{p5}.}
\item{\code{Variety}}{Variety identifier. Factor with 25 levels.}
\item{\code{Yield}}{Yield response. Integer.}
}}
\item{\code{v-lmm-alfalfa}}{Balanced alfalfa split-plot experiment: varieties are whole-plot treatments within blocks and dates are subplot treatments.
File \file{v-lmm-alfalfa.rda} loads \code{alfalfa},
with 72 rows and 4 columns.
The stored data contain no missing values.
\describe{
\item{\code{Variety}}{Variety identifier. Factor with 3 levels: \code{Cossack}, \code{Ladack}, \code{Ranger}.}
\item{\code{Date}}{Date treatment. Factor with 4 levels: \code{A}, \code{B}, \code{C}, \code{D}.}
\item{\code{Block}}{Block identifier. Factor with 6 levels: \code{b1}, \code{b2}, \code{b3}, \code{b4}, \code{b5}, \code{b6}.}
\item{\code{Yield}}{Yield response. Numeric (double).}
}}
\item{\code{v-lmm-oats}}{Balanced oats split-plot experiment: varieties are whole-plot treatments within blocks and nitrogen levels are subplot treatments.
File \file{v-lmm-oats.rda} loads \code{oats},
with 72 rows and 4 columns.
The stored data contain no missing values.
\describe{
\item{\code{Block}}{Block identifier. Factor with 6 levels: \code{I}, \code{II}, \code{III}, \code{IV}, \code{V}, \code{VI}.}
\item{\code{Variety}}{Variety identifier. Factor with 3 levels: \code{Golden Rain}, \code{Marvellous}, \code{Victory}.}
\item{\code{Nitro}}{Nitrogen treatment. Factor with 4 levels: \code{n0.0}, \code{n0.2}, \code{n0.4}, \code{n0.6}.}
\item{\code{Yield}}{Yield response. Integer.}
}}
\item{\code{v-lmm-barley-lattice}}{Barley yield observations with row and column blocks.
File \file{v-lmm-barley-lattice.rda} loads \code{barley.lattice},
with 150 rows and 7 columns. These observations are from the 1976
balanced lattice square trial at Slate Hall Farm, UK, described by
Gilmour, Thompson and Cullis (1995). The six design and treatment columns
are stored as factors; the zero-padded field coordinates require no
releveling before spatial analysis.
The stored data contain no missing values.
\describe{
\item{\code{Rep}}{Replicate identifier. Factor with 6 levels: \code{re1}, \code{re2}, \code{re3}, \code{re4}, \code{re5}, \code{re6}.}
\item{\code{RowBlk}}{Row-block identifier. Factor with 30 levels.}
\item{\code{ColBlk}}{Column-block identifier. Factor with 30 levels.}
\item{\code{Row}}{Physical field-row identifier. Factor with levels
\code{row01}, \ldots, \code{row10}, in physical and lexical order.}
\item{\code{Column}}{Physical field-column identifier. Factor with levels
\code{col01}, \ldots, \code{col15}, in physical and lexical order.}
\item{\code{Variety}}{Variety identifier. Factor with 25 levels.}
\item{\code{yield}}{Yield response. Integer.}
}}
\item{\code{v-lmm-rapeseed-blocks}}{Rapeseed height observations with alternative block layouts.
File \file{v-lmm-rapeseed-blocks.rda} loads \code{rapeseed.blocks},
with 120 rows and 7 columns.
The stored data contain no missing values.
\describe{
\item{\code{rep}}{Replicate identifier. Factor with 3 levels: \code{r1}, \code{r2}, \code{r3}.}
\item{\code{blk}}{Block identifier. Factor with 10 levels.}
\item{\code{geno}}{Genotype identifier. Factor with 20 levels.}
\item{\code{RowBlk}}{Row-block identifier. Factor with 6 levels: \code{rb01}, \code{rb02}, \code{rb03}, \code{rb04}, \code{rb05}, \code{rb06}.}
\item{\code{ColBlk}}{Column-block identifier. Factor with 20 levels.}
\item{\code{env}}{Environment identifier. Factor with 2 levels: \code{k}, \code{s}.}
\item{\code{hgt}}{Height response. Integer.}
}}
\item{\code{v-lmm-maize-counties}}{Maize yield observations from trials in several counties.
File \file{v-lmm-maize-counties.rda} loads \code{maize.counties},
with 1152 rows and 7 columns.
The stored data contain no missing values.
\describe{
\item{\code{county}}{County identifier. Factor with 6 levels: \code{C1}, \code{C2}, \code{C3}, \code{C4}, \code{C5}, \code{C6}.}
\item{\code{row}}{Row identifier. Factor with 18 levels.}
\item{\code{col}}{Column identifier. Factor with 11 levels.}
\item{\code{rep}}{Replicate identifier. Factor with 3 levels: \code{R1}, \code{R2}, \code{R3}.}
\item{\code{block}}{Block identifier. Factor with 8 levels: \code{B1}, \code{B2}, \code{B3}, \code{B4}, \code{B5}, \code{B6}, \code{B7}, \code{B8}.}
\item{\code{yield}}{Yield response. Numeric (double).}
\item{\code{gen}}{Genotype identifier. Factor with 64 levels.}
}}
\item{\code{v-lmm-soybean}}{Soybean yield observations across locations and years.
File \file{v-lmm-soybean.rda} loads \code{soybean},
with 306 rows and 6 columns.
The stored data contain no missing values.
\describe{
\item{\code{yield}}{Yield response. Numeric (double).}
\item{\code{rep}}{Replicate identifier. Factor with 4 levels: \code{R1}, \code{R2}, \code{R3}, \code{R4}.}
\item{\code{gen}}{Genotype identifier. Factor with 7 levels: \code{Chip}, \code{Cors}, \code{Evan}, \code{Hodg}, \code{S200}, \code{Well}, \code{Wilk}.}
\item{\code{env}}{Environment identifier. Factor with 12 levels.}
\item{\code{loc}}{Location identifier. Factor with 4 levels: \code{A}, \code{G}, \code{L}, \code{V}.}
\item{\code{year}}{Year. Factor with 3 levels: \code{y80}, \code{y81}, \code{y82}.}
}}
\item{\code{v-lmm-soybean-met}}{An extended soybean multi-environment yield dataset.
File \file{v-lmm-soybean-met.rda} loads \code{soybean.met},
with 1454 rows and 6 columns.
The stored data contain no missing values.
\describe{
\item{\code{yield}}{Yield response. Numeric (double).}
\item{\code{rep}}{Replicate identifier. Factor with 4 levels: \code{R1}, \code{R2}, \code{R3}, \code{R4}.}
\item{\code{gen}}{Genotype identifier. Factor with 7 levels: \code{Chip}, \code{Cors}, \code{Evan}, \code{Hodg}, \code{S200}, \code{Well}, \code{Wilk}.}
\item{\code{env}}{Environment identifier. Factor with 55 levels.}
\item{\code{loc}}{Location identifier. Factor with 10 levels.}
\item{\code{year}}{Year. Factor with 12 levels.}
}}
\item{\code{v-lmm-material-groups}}{Multi-environment observations with local and reference material groups.
File \file{v-lmm-material-groups.rda} loads \code{material.groups},
with 464 rows and 11 columns.
The stored data contain no missing values.
\describe{
\item{\code{env}}{Environment identifier. Factor with 8 levels: \code{B70}, \code{B71}, \code{L70}, \code{L71}, \code{N70}, \code{N71}, \code{R70}, \code{R71}.}
\item{\code{loc}}{Location identifier. Factor with 4 levels: \code{Brookstead}, \code{Lawes}, \code{Nambour}, \code{RedlandBay}.}
\item{\code{year}}{Year. Factor with 2 levels: \code{1970}, \code{1971}.}
\item{\code{gen}}{Genotype identifier. Factor with 58 levels.}
\item{\code{yield}}{Yield response. Numeric (double).}
\item{\code{height}}{Height response. Numeric (double).}
\item{\code{lodging}}{Lodging response. Numeric (double).}
\item{\code{size}}{Recorded size trait. Numeric (double).}
\item{\code{protein}}{Protein response. Numeric (double).}
\item{\code{oil}}{Oil response. Numeric (double).}
\item{\code{grp}}{Material group (local or reference). Factor with 2 levels: \code{local}, \code{ref}.}
}}
\item{\code{v-lmm-maize-met}}{Maize plot records for yield, dry-matter content and plant height across
environments. The analyses estimate adjusted genotype means and genetic
covariances among traits.
File \file{v-lmm-maize-met.rda} loads \code{maize.met},
with 1200 rows and 6 columns.
The stored data contain no missing values.
\describe{
\item{\code{env}}{Environment identifier. Factor with 24 levels.}
\item{\code{gen}}{Genotype identifier. Factor with 25 levels.}
\item{\code{rep}}{Replicate identifier. Factor with 2 levels: \code{R1}, \code{R2}.}
\item{\code{pht}}{Plant-height response. Integer.}
\item{\code{dmc}}{Dry-matter-content response. Numeric (double).}
\item{\code{yld}}{Yield response. Numeric (double).}
}}
\item{\code{v-lmm-barley-met}}{Anonymized barley multi-environment observations.
File \file{v-lmm-barley-met.rda} loads \code{barley.met},
with 900 rows and 9 columns: 15 lines in 5 locations, 3 years and 4 replicates.
The stored data contain no missing values. The record for \code{HFH},
\code{y2018}, \code{lne04}, \code{r3} contains trait means from the other
three replicates in that cell to complete the balanced layout.
Location and line identifiers are coded.
\describe{
\item{\code{loc}}{Location identifier. Factor with levels \code{HOH}, \code{EWE}, \code{EIN}, \code{OLI}, \code{HFH}, in that order.}
\item{\code{year}}{Year. Factor with levels \code{y2016}, \code{y2017}, \code{y2018}.}
\item{\code{gen}}{Line identifier. Factor with levels \code{lne01} to \code{lne15}.}
\item{\code{rep}}{Replicate identifier. Factor with levels \code{r1} to \code{r4}.}
\item{\code{DHE}}{Recorded DHE trait. Numeric.}
\item{\code{DMA}}{Recorded DMA trait. Numeric.}
\item{\code{PLH}}{Plant-height response. Numeric.}
\item{\code{GYLD}}{Grain-yield response. Numeric.}
\item{\code{TKW}}{Thousand-kernel-weight response. Numeric.}
}}
\item{\code{v-lmm-met-unbalanced}}{Yield observations from unbalanced multi-environment trials.
File \file{v-lmm-met-unbalanced.rda} loads \code{met.unbalanced},
with 450 rows and 6 columns.
The stored data contain no missing values.
\describe{
\item{\code{gen}}{Genotype identifier. Factor with 15 levels and zero-padded
identifiers such as \code{g01}.}
\item{\code{rep}}{Replicate identifier. Factor with 3 levels: \code{r1}, \code{r2}, \code{r3}.}
\item{\code{env}}{Environment identifier. Factor with 10 levels.}
\item{\code{yld}}{Yield response. Numeric (double).}
\item{\code{loc}}{Location identifier. Factor with 6 levels: \code{ASE}, \code{BOV}, \code{MOS}, \code{NIE}, \code{RHH}, \code{ROS}.}
\item{\code{year}}{Year. Factor with 2 levels: \code{2015}, \code{2016}.}
}}
\item{\code{v-lmm-met-missing}}{Multi-trait observations with missing yield values.
File \file{v-lmm-met-missing.rda} loads \code{met.missing},
with 720 rows and 8 columns.
\describe{
\item{\code{loc}}{Location identifier. Factor with 5 levels: \code{loc1}, \code{loc2}, \code{loc3}, \code{loc4}, \code{loc5}.}
\item{\code{rep}}{Replicate identifier. Factor with 4 levels: \code{r1}, \code{r2}, \code{r3}, \code{r4}.}
\item{\code{year}}{Year. Factor with 3 levels: \code{y2017}, \code{y2018}, \code{y2019}.}
\item{\code{gen}}{Genotype identifier. Factor with 20 levels.}
\item{\code{DHE}}{Recorded DHE trait. Integer.}
\item{\code{DMA}}{Recorded DMA trait. Integer.}
\item{\code{PLH}}{Plant-height response. Numeric (double).}
\item{\code{YLD}}{Yield response. Numeric (double). Contains 2 missing values (\code{NA}).}
}}
\item{\code{v-lmm-met-seven-traits}}{Seven recorded traits across locations and years.
File \file{v-lmm-met-seven-traits.rda} loads \code{met.seven.traits},
with 300 rows and 11 columns.
\describe{
\item{\code{year}}{Year. Factor with 2 levels: \code{y2019}, \code{y2020}.}
\item{\code{loc}}{Location identifier. Factor with 3 levels: \code{loc1}, \code{loc2}, \code{loc3}.}
\item{\code{rep}}{Replicate identifier. Factor with 3 levels: \code{r1}, \code{r2}, \code{r3}.}
\item{\code{gen}}{Genotype identifier. Factor with 28 levels.}
\item{\code{YLD}}{Yield response. Numeric (double).}
\item{\code{HLW}}{Hectolitre-weight response. Numeric (double). Contains 9 missing values (\code{NA}).}
\item{\code{TKW}}{Thousand-kernel-weight response. Numeric (double).}
\item{\code{PLH}}{Plant-height response. Numeric (double).}
\item{\code{DH}}{Recorded DH trait. Integer.}
\item{\code{DM}}{Recorded DM trait. Integer.}
\item{\code{SC}}{Recorded SC score. Integer.}
}}
\item{\code{v-lmm-barley-seven-traits}}{Barley data with all seven responses
observed. File \file{v-lmm-barley-seven-traits.rda} loads
\code{barley.seven.traits}, with 291 rows and 14 columns. It contains the
records from \code{met.seven.traits} with all seven responses observed. The
nine records with missing \code{HLW} are excluded. Three factors describe
the groups used in the fit. \code{env} identifies the five year-location
combinations, \code{gen_env} identifies genotype within environment, and
\code{env_rep} identifies replicate within environment. The original
measurement units and row identifiers are retained.}
\item{\code{v-lmm-chickens}}{Repeated chick body-weight measurements under different diets.
File \file{v-lmm-chickens.rda} loads \code{chickens},
with 578 rows and 4 columns.
The stored data contain no missing values.
\describe{
\item{\code{weight}}{Body-weight response. Integer.}
\item{\code{Time}}{Measurement time. Integer.}
\item{\code{Chick}}{Chick identifier. Factor with 50 levels.}
\item{\code{Diet}}{Diet treatment. Factor with 4 levels: \code{d1}, \code{d2}, \code{d3}, \code{d4}.}
}}
}}
\details{The file name and loaded object name are listed for every dataset.
All listed datasets are available through \code{data()}; its \code{envir}
argument specifies the destination environment. Object names follow the
names listed above, with dots separating words in multiword names.

All categorical grouping variables, including years, are stored as factors. Responses and
continuous covariates (including chick measurement time) remain numeric.
Observations and missing values are retained except for constructions
explicitly described above. The barley multi-environment identifiers are
coded, and its balancing record is documented in its dataset entry.
Measurement units are retained; no unit conversion or standardization is
applied when storing these datasets.
Trait abbreviations without a recorded definition are used as variable names.

The stored missing values occur in \code{YLD} of \code{met.missing}
(two values) and \code{HLW} of \code{met.seven.traits} (nine values).
Missing experimental combinations in an unbalanced design are distinct from
explicit \code{NA} values. Repeated group labels can identify different
groups across environments or higher-level factors. The model formula must
specify the required interactions and nesting.}
\section{Constructed tropical maize trial}{
The additional file \file{v-lmm-tropmaize.rda} loads \code{tropmaize}, a
792-row data frame with columns \code{block}, \code{plot}, \code{gen}
and \code{y}. It contains synthetic observations from three complete
blocks in one environment, with all 264 genotypes in every block.
Field-only genotype BLUEs reproduce the values in \code{v-tropmaize-phe}.
The construction and full format are described in \code{\link{v-lmm-tropmaize}}.
The data contain no missing values.
}
\section{Additional experimental datasets}{
\describe{
\item{\code{v-lmm-fertilizer}}{Object \code{fertilizer}: 25 observations of \code{yield}, fertilizer factor \code{fert} and replication factor \code{rep}.}
\item{\code{v-lmm-cotton}}{Object \code{cotton}: 15 observations of \code{Strength}, block factor \code{Blk} and fertilizer factor \code{Fert}.}
\item{\code{v-lmm-sugar}}{Object \code{sugar}: 36 observations of \code{Yield} and factors \code{Row}, \code{Col} and \code{Variety}.}
\item{\code{v-lmm-peas}}{Object \code{peas}: 36 observations of \code{Yield} and factors \code{Variety}, \code{Spacing} and \code{Block}.}
}}
\value{Each file supplies the named data frame described above.
\code{data} loads the data frame into the
environment specified by \code{envir}.}
\references{Gilmour AR, Thompson R, Cullis BR (1995). Average information
REML: An efficient algorithm for variance parameter estimation in linear
mixed models. Biometrics 51:1440--1450. \doi{10.2307/2533274}.}
\seealso{\code{\link{st.lmm}}, \code{\link{st.lsm}}}
\keyword{datasets}
