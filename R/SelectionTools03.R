




"st.phase" <- function(data.set="default", ref=NULL, phasing.options=list(),
                         impute=FALSE, commit=TRUE)
{
    diagnostics <- .Call("st_phase", data.set, ref, phasing.options, impute,
                         commit, PACKAGE="SelectionTools")
    if (isTRUE(commit)) invisible(NULL) else diagnostics
}


"st.destroy.phase" <- function(data.set="default")
{
    invisible(.Call("st_destroy_phase", data.set))
}


"st.switch.error" <- function(true.data.set, data.set="default")
{
    .Call("st_switch_error", true.data.set, data.set)

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

    cat(label,"\n",sep="")
    if (!is.null(result))
        print(result[c("actual.seed", "windows", "maximum.worker.team",
                       "filled.genotypes", "committed")])
    print(score)
    invisible(result)
}
