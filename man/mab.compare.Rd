\name{mab.compare}
\alias{mab.compare}

\title{Runs and Compares a Series of Marker-Assisted Backcrossing Simulations}

\description{
Runs up to nine marker-assisted backcrossing simulation scenarios by calling
\code{mab.simulate()} for each active scenario. Scenario A provides the base
settings; scenarios B through I can override individual settings.
}

\usage{
mab.compare(
    a.simulation.name,
    a.simulation.run,
    a.repetitions = 1000,
    a.recurrent.parent = "H",
    a.donor.parent = "H",
    a.linkage.map,
    a.target.loci = NULL,
    a.flanking.loci = NULL,
    a.recipient.loci = NULL,
    a.gen.type = NULL,
    a.population.size = NULL,
    a.sel.strategy = NULL,
    a.no.selected = NULL,
    a.no.preselected = NULL,
    a.reg.chr = NULL,
    a.reg.begin = NULL,
    a.reg.end = NULL,
    a.missing.allele = "9",
    a.success.factor = 1,
    b.simulation.name = NULL,
    c.simulation.name = NULL,
    d.simulation.name = NULL,
    b.simulation.run = NULL,
    c.simulation.run = NULL,
    d.simulation.run = NULL,
    b.repetitions = NULL,
    c.repetitions = NULL,
    d.repetitions = NULL,
    b.recurrent.parent = NULL,
    c.recurrent.parent = NULL,
    d.recurrent.parent = NULL,
    b.donor.parent = NULL,
    c.donor.parent = NULL,
    d.donor.parent = NULL,
    b.linkage.map = NULL,
    c.linkage.map = NULL,
    d.linkage.map = NULL,
    b.target.loci = NULL,
    c.target.loci = NULL,
    d.target.loci = NULL,
    b.flanking.loci = NULL,
    c.flanking.loci = NULL,
    d.flanking.loci = NULL,
    b.recipient.loci = NULL,
    c.recipient.loci = NULL,
    d.recipient.loci = NULL,
    b.gen.type = NULL,
    c.gen.type = NULL,
    d.gen.type = NULL,
    b.population.size = NULL,
    c.population.size = NULL,
    d.population.size = NULL,
    b.sel.strategy = NULL,
    c.sel.strategy = NULL,
    d.sel.strategy = NULL,
    b.no.selected = NULL,
    c.no.selected = NULL,
    d.no.selected = NULL,
    b.no.preselected = NULL,
    c.no.preselected = NULL,
    d.no.preselected = NULL,
    b.reg.chr = NULL,
    c.reg.chr = NULL,
    d.reg.chr = NULL,
    b.reg.begin = NULL,
    c.reg.begin = NULL,
    d.reg.begin = NULL,
    b.reg.end = NULL,
    c.reg.end = NULL,
    d.reg.end = NULL,
    b.missing.allele = NULL,
    c.missing.allele = NULL,
    d.missing.allele = NULL,
    b.success.factor = NULL,
    c.success.factor = NULL,
    d.success.factor = NULL,
    e.simulation.name = NULL,
    f.simulation.name = NULL,
    g.simulation.name = NULL,
    e.simulation.run = NULL,
    f.simulation.run = NULL,
    g.simulation.run = NULL,
    e.repetitions = NULL,
    f.repetitions = NULL,
    g.repetitions = NULL,
    e.recurrent.parent = NULL,
    f.recurrent.parent = NULL,
    g.recurrent.parent = NULL,
    e.donor.parent = NULL,
    f.donor.parent = NULL,
    g.donor.parent = NULL,
    e.linkage.map = NULL,
    f.linkage.map = NULL,
    g.linkage.map = NULL,
    e.target.loci = NULL,
    f.target.loci = NULL,
    g.target.loci = NULL,
    e.flanking.loci = NULL,
    f.flanking.loci = NULL,
    g.flanking.loci = NULL,
    e.recipient.loci = NULL,
    f.recipient.loci = NULL,
    g.recipient.loci = NULL,
    e.gen.type = NULL,
    f.gen.type = NULL,
    g.gen.type = NULL,
    e.population.size = NULL,
    f.population.size = NULL,
    g.population.size = NULL,
    e.sel.strategy = NULL,
    f.sel.strategy = NULL,
    g.sel.strategy = NULL,
    e.no.selected = NULL,
    f.no.selected = NULL,
    g.no.selected = NULL,
    e.no.preselected = NULL,
    f.no.preselected = NULL,
    g.no.preselected = NULL,
    e.reg.chr = NULL,
    f.reg.chr = NULL,
    g.reg.chr = NULL,
    e.reg.begin = NULL,
    f.reg.begin = NULL,
    g.reg.begin = NULL,
    e.reg.end = NULL,
    f.reg.end = NULL,
    g.reg.end = NULL,
    e.missing.allele = NULL,
    f.missing.allele = NULL,
    g.missing.allele = NULL,
    e.success.factor = NULL,
    f.success.factor = NULL,
    g.success.factor = NULL,
    h.simulation.name = NULL,
    i.simulation.name = NULL,
    h.simulation.run = NULL,
    i.simulation.run = NULL,
    h.repetitions = NULL,
    i.repetitions = NULL,
    h.recurrent.parent = NULL,
    i.recurrent.parent = NULL,
    h.donor.parent = NULL,
    i.donor.parent = NULL,
    h.linkage.map = NULL,
    i.linkage.map = NULL,
    h.target.loci = NULL,
    i.target.loci = NULL,
    h.flanking.loci = NULL,
    i.flanking.loci = NULL,
    h.recipient.loci = NULL,
    i.recipient.loci = NULL,
    h.gen.type = NULL,
    i.gen.type = NULL,
    h.population.size = NULL,
    i.population.size = NULL,
    h.sel.strategy = NULL,
    i.sel.strategy = NULL,
    h.no.selected = NULL,
    i.no.selected = NULL,
    h.no.preselected = NULL,
    i.no.preselected = NULL,
    h.reg.chr = NULL,
    i.reg.chr = NULL,
    h.reg.begin = NULL,
    i.reg.begin = NULL,
    h.reg.end = NULL,
    i.reg.end = NULL,
    h.missing.allele = NULL,
    i.missing.allele = NULL,
    h.success.factor = NULL,
    i.success.factor = NULL,
    result.files
)
}

\arguments{
  \item{a.simulation.name, b.simulation.name, c.simulation.name, d.simulation.name, e.simulation.name, f.simulation.name, g.simulation.name, h.simulation.name, i.simulation.name}{Simulation name for each scenario. Scenario A is required. For scenarios B through I, a non-\code{NULL} simulation name activates that scenario; a \code{NULL} value means that scenario is not run.}
  \item{a.simulation.run, b.simulation.run, c.simulation.run, d.simulation.run, e.simulation.run, f.simulation.run, g.simulation.run, h.simulation.run, i.simulation.run}{Simulation-run identifier passed to \code{mab.simulate()}. For active scenarios B through I, \code{NULL} inherits \code{a.simulation.run}.}
  \item{a.repetitions, b.repetitions, c.repetitions, d.repetitions, e.repetitions, f.repetitions, g.repetitions, h.repetitions, i.repetitions}{Number of simulation repetitions. Scenario A defaults to 1000. For active scenarios B through I, \code{NULL} inherits \code{a.repetitions}.}
  \item{a.recurrent.parent, b.recurrent.parent, c.recurrent.parent, d.recurrent.parent, e.recurrent.parent, f.recurrent.parent, g.recurrent.parent, h.recurrent.parent, i.recurrent.parent}{Description of the recurrent parent passed to \code{mab.simulate()}; \code{"H"} denotes the homozygous-parent mode used by that function. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{a.donor.parent, b.donor.parent, c.donor.parent, d.donor.parent, e.donor.parent, f.donor.parent, g.donor.parent, h.donor.parent, i.donor.parent}{Description of the donor parent passed to \code{mab.simulate()}; it follows the same representation as \code{recurrent.parent}. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{a.linkage.map, b.linkage.map, c.linkage.map, d.linkage.map, e.linkage.map, f.linkage.map, g.linkage.map, h.linkage.map, i.linkage.map}{Name of the linkage-map input file passed to \code{mab.simulate()}. Scenario A has no default. For active scenarios B through I, \code{NULL} inherits \code{a.linkage.map}.}
  \item{a.target.loci, b.target.loci, c.target.loci, d.target.loci, e.target.loci, f.target.loci, g.target.loci, h.target.loci, i.target.loci}{Character vector of target-locus names passed to \code{mab.simulate()}, or \code{NULL}. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{a.flanking.loci, b.flanking.loci, c.flanking.loci, d.flanking.loci, e.flanking.loci, f.flanking.loci, g.flanking.loci, h.flanking.loci, i.flanking.loci}{Character vector of flanking-locus names used for recurrent-parent preselection, or \code{NULL}. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{a.recipient.loci, b.recipient.loci, c.recipient.loci, d.recipient.loci, e.recipient.loci, f.recipient.loci, g.recipient.loci, h.recipient.loci, i.recipient.loci}{Character vector of loci used for preselection of recurrent-parent alleles, or \code{NULL}. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{a.gen.type, b.gen.type, c.gen.type, d.gen.type, e.gen.type, f.gen.type, g.gen.type, h.gen.type, i.gen.type}{Character vector describing the generation types passed to \code{mab.simulate()} (for example \code{"f1"}, \code{"bc"}, \code{"s"}, and \code{"dh"}), or \code{NULL}. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{a.population.size, b.population.size, c.population.size, d.population.size, e.population.size, f.population.size, g.population.size, h.population.size, i.population.size}{Numeric vector of population sizes for the generations of the backcrossing program, or \code{NULL}. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{a.sel.strategy, b.sel.strategy, c.sel.strategy, d.sel.strategy, e.sel.strategy, f.sel.strategy, g.sel.strategy, h.sel.strategy, i.sel.strategy}{Character vector of selection-strategy codes passed to \code{mab.simulate()}, or \code{NULL}. Further details appear in \code{mab.simulate()} for the implemented strategies. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{a.no.selected, b.no.selected, c.no.selected, d.no.selected, e.no.selected, f.no.selected, g.no.selected, h.no.selected, i.no.selected}{Numeric vector giving the number of selected individuals by generation, or \code{NULL}. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{a.no.preselected, b.no.preselected, c.no.preselected, d.no.preselected, e.no.preselected, f.no.preselected, g.no.preselected, h.no.preselected, i.no.preselected}{Numeric vector giving the minimum number of individuals retained by the flanking-marker preselection step, or \code{NULL}. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{a.reg.chr, b.reg.chr, c.reg.chr, d.reg.chr, e.reg.chr, f.reg.chr, g.reg.chr, h.reg.chr, i.reg.chr}{Chromosome numbers for regions evaluated separately for recurrent-parent genome content, or \code{NULL}. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{a.reg.begin, b.reg.begin, c.reg.begin, d.reg.begin, e.reg.begin, f.reg.begin, g.reg.begin, h.reg.begin, i.reg.begin}{Beginning positions of the regions specified by \code{reg.chr}, in cM from the telomere, or \code{NULL}. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{a.reg.end, b.reg.end, c.reg.end, d.reg.end, e.reg.end, f.reg.end, g.reg.end, h.reg.end, i.reg.end}{End positions of the regions specified by \code{reg.chr}, in cM from the telomere; a value of zero is interpreted by \code{mab.simulate()} as the chromosome end. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{a.missing.allele, b.missing.allele, c.missing.allele, d.missing.allele, e.missing.allele, f.missing.allele, g.missing.allele, h.missing.allele, i.missing.allele}{Character string used by \code{mab.simulate()} to identify missing marker data. Scenario A defaults to \code{"9"}. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{a.success.factor, b.success.factor, c.success.factor, d.success.factor, e.success.factor, f.success.factor, g.success.factor, h.success.factor, i.success.factor}{Success factor passed to \code{mab.simulate()}. Scenario A defaults to 1. When it is below 1, \code{mab.simulate()} draws the number of successful selected individuals from a binomial distribution and retains at least one individual. For active scenarios B through I, \code{NULL} inherits the scenario-A setting.}
  \item{result.files}{Character vector containing one explicit result-file name for each active scenario, in scenario order A through I. The names must be non-empty and unique. Relative names are resolved against \code{st.data.dir}; absolute paths are used directly.}
}

\details{
The scenario prefixes \code{a.} through \code{i.} identify up to nine
simulation definitions. Scenario A is always defined by
\code{a.simulation.name}. A scenario B through I is run only when its
\code{simulation.name} argument is non-\code{NULL}.

For every setting other than \code{simulation.name}, an active scenario B
through I inherits the corresponding scenario-A value when its own argument is
\code{NULL}; a non-\code{NULL} scenario-specific value overrides the
scenario-A value. The resulting settings are passed directly to
\code{mab.simulate()}. Further details appear in \code{mab.simulate()} for the detailed meaning of
the generation and selection settings and for the structure of the simulation
results. The file name for every active scenario is taken from
\code{result.files}; no result-file name is generated automatically.
}

\value{
\code{NULL}, returned invisibly as the value of the final \code{for} loop.
The function is called for its side effect of running the active simulation
scenarios; each call to \code{mab.simulate()} saves its result in the
corresponding explicitly supplied file from \code{result.files}.
}

\seealso{
\code{\link{mab.simulate}}, \code{\link{mab.load.data}},
\code{\link{mab.tabulate}}, \code{\link{mab.Input.files}},
\code{\link{mab.Examples}}
}
