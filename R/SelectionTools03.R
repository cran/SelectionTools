###########################################################################
# SelectionTools03.R
# Chromosome-wide statistical phasing
###########################################################################

"st.phase" <- function(data.set="default", ref=NULL, phasing.options=list(),
                         impute=FALSE, commit=TRUE)
{
    invisible(.Call("st_phase", data.set, ref, phasing.options, impute, commit,
                    PACKAGE="SelectionTools"))
}


"st.phase.v1" <- function(
    reference.data.set = NULL,
    reference.individuals = NULL,
    initial.pair.weights = NULL,
    population.type = NULL,
    expected.breaks.per.morgan = NULL,
    max.states = NULL,
    error.rate = 0.001,
    break.rate1.per.morgan = NULL,
    break.rate2.per.morgan = NULL,
    confidence.threshold = 0.50,
    impute = TRUE,
    imputation.threshold = 0.90,
    individuals = NULL,
    chromosomes = NULL,
    burnin.max.iterations = 4,
    burnin.min.iterations = 2,
    burnin.convergence.fraction = 0.001,
    progressive.max.iterations = 2,
    progressive.min.iterations = 1,
    progressive.convergence.fraction = 0.001,
    burnin.states = 12,
    local.states = 8,
    stage1.minor.allele.count = 4,
    provisional.threshold = 0.80,
    commit = TRUE,
    data.set = "default"
)
{

    if (!is.null(reference.data.set)) {
        cohort.arguments.supplied <-
            !missing(burnin.max.iterations) ||
            !missing(burnin.min.iterations) ||
            !missing(burnin.convergence.fraction) ||
            !missing(progressive.max.iterations) ||
            !missing(progressive.min.iterations) ||
            !missing(progressive.convergence.fraction) ||
            !missing(burnin.states) ||
            !missing(local.states) ||
            !missing(stage1.minor.allele.count) ||
            !missing(provisional.threshold)

        if (cohort.arguments.supplied)
            stop("Reference-data-set mode does not accept cohort-specific arguments")

        burnin.max.iterations <- NULL
        burnin.min.iterations <- NULL
        burnin.convergence.fraction <- NULL
        progressive.max.iterations <- NULL
        progressive.min.iterations <- NULL
        progressive.convergence.fraction <- NULL
        burnin.states <- NULL
        local.states <- NULL
        stage1.minor.allele.count <- NULL
        provisional.threshold <- NULL
    }

    .Call(
        "st_phase_03_step1", data.set, reference.data.set,
        reference.individuals, initial.pair.weights, population.type,
        expected.breaks.per.morgan, max.states, individuals, chromosomes,
        error.rate, break.rate1.per.morgan, break.rate2.per.morgan,
        confidence.threshold, impute, imputation.threshold,
        burnin.max.iterations, burnin.min.iterations,
        burnin.convergence.fraction, progressive.max.iterations,
        progressive.min.iterations, progressive.convergence.fraction,
        burnin.states, local.states, stage1.minor.allele.count,
        provisional.threshold, commit
    )
}


"st.destroy.phase" <- function(data.set="default")
{
    invisible(.Call("st_destroy_phase", data.set))
}


"st.switch.error" <- function(true.data.set, data.set="default")
{
    result <- .Call("st_switch_error", true.data.set, data.set)
    result <- data.frame(result=unname(unlist(result,use.names=FALSE)),
                         row.names=names(result))
    class(result) <- c("st.switch.error","data.frame")
    result
}


"print.st.switch.error" <- function(x, complete=FALSE, digits=3, ...)
{
    digits <- as.integer(digits)
    if (length(digits) != 1L || is.na(digits) || digits < 1L || digits > 22L)
        stop("digits must be one integer between 1 and 22")
    if (length(complete) != 1L || is.na(complete))
        stop("complete must be TRUE or FALSE")
    out <- if (isTRUE(complete)) x else
           x[c("switch.error","coverage","genotype.mismatches"),,drop=FALSE]
    class(out) <- "data.frame"
    print(out,digits=digits,...)
    invisible(x)
}


"show.phase" <- function(label,result,score)
{
    if (is.data.frame(score))
        score <- setNames(score[[1L]],row.names(score))

    suffix <- if (isTRUE(result$breaks.estimated)) " (estimated)" else ""
    cat(label,"\n",sep="")
    cat("  population.type: ",result$population.type,"\n",sep="")
    cat(sprintf("  breaks/Morgan: %.3f%s\n",
                result$expected.breaks.per.morgan,suffix))
    cat(sprintf("  heterozygotes: %d; written: %d; imputed: %d\n",
                result$n.heterozygous,result$n.written,result$n.imputed))
    cat(sprintf("  switch error: %.3f%% (%d/%d links); coverage: %.2f%%\n",
                score[["switch.error.percent"]],score[["errors"]],
                score[["links"]],100*score[["coverage"]]))
    invisible(result)
}
