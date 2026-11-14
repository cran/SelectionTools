

st.options <- function(options = NULL) {
    result <- .Call("st_lmm_frontend", list(action = "st.options", options = options),
        PACKAGE = "SelectionTools")
    if (is.null(options)) result else invisible(result)
}

.st_mxd_mv <- function(y, X, random = list(), residual = NULL, method = "REML", control = list()) .Call("st_mxd_mv_interface",
    list(y = y, X = X, random = random, residual = residual, method = method, control = control), PACKAGE = "SelectionTools")

.st_mxd_mv_block <- function(group = NULL, basis = NULL, A = NULL, covariance = "diag", start = NULL,
    kernel = NULL, hom = NULL) .Call("st_mxd_mv_block", list(group = group, basis = basis, A = A, covariance = covariance,
    start = start, kernel = kernel, hom = hom), PACKAGE = "SelectionTools")

.st_mxd_mv_linear <- function(object, L, tol = 1e-07) .Call("st_mxd_mv_linear", list(object = object,
    L = L, tol = tol), PACKAGE = "SelectionTools")

st.lmm <- function(fixed, random = NULL, cov = ~units, data, method = "REML",
    ddf = "Satterthwaite", type = 3L, random.test = "auto", control = list(),
    ucov = NULL, prediction = "NAonly", selection = NULL, weights = NULL,
    contrasts = NULL, subset = NULL, na.action = NULL, offset = NULL,
    rcov = NULL, series = NULL, constraints = NULL) .Call("st_lmm",
        list(call = match.call(), fixed = fixed, random = random, cov = cov,
            data = if (missing(data)) NULL else data, method = method, ddf = ddf, type = type,
            random.test = random.test, control = control,
            ddf.supplied = !missing(ddf),
            random.test.supplied = !missing(random.test),
            ucov = substitute(ucov), prediction = prediction,
            selection = selection, weights = substitute(weights),
            contrasts = contrasts, weights.supplied = !missing(weights),
            subset = substitute(subset), na.action = na.action,
            offset = substitute(offset), rcov = rcov, series = series,
            constraints = constraints,
            data.supplied = !missing(data), cov.supplied = !missing(cov),
            envir = parent.frame()),
        PACKAGE = "SelectionTools")

.st_lmm_method <- function(object, action, ...) .Call("st_lmm_method", list(object = object, action = action,
    ...), PACKAGE = "SelectionTools")

.st_lmm_result <- function(x) .st_lmm_method(x, "result", recursive = TRUE)

.st_lmm_plain <- function(x, recursive = FALSE) .st_lmm_method(x, "plain.result", recursive = recursive)

print.st.lmm.result <- function(x, digits = NULL, top = NULL, bottom = NULL, char.per.line = NULL, ...) invisible(.st_lmm_method(x,
    "print.result", digits = digits, top = top, bottom = bottom, char.per.line = char.per.line, ...))

Ops.st.lmm.result <- function(e1, e2) .Call("st_lmm_frontend", list(action = "result.Ops", generic = .Generic,
    values = if (missing(e2)) list(e1) else list(e1, e2)), PACKAGE = "SelectionTools")

Math.st.lmm.result <- function(x, ...) .Call("st_lmm_frontend", list(action = "result.Math", generic = .Generic,
    values = list(x, ...)), PACKAGE = "SelectionTools")

Summary.st.lmm.result <- function(..., na.rm = FALSE) .Call("st_lmm_frontend", list(action = "result.Summary",
    generic = .Generic, values = list(...), na.rm = na.rm), PACKAGE = "SelectionTools")

as.data.frame.st.lmm.result <- function(x, row.names = NULL, optional = FALSE, ...) .Call("st_lmm_frontend",
    list(action = "result.as.data.frame", object = x, row.names = row.names, optional = optional, dots = list(...),
        expression = substitute(x)), PACKAGE = "SelectionTools")

as.matrix.st.lmm.result <- function(x, ...) .Call("st_lmm_frontend", list(action = "result.as.matrix",
    object = x, dots = list(...)), PACKAGE = "SelectionTools")

as.vector.st.lmm.result <- function(x, mode = "any") .Call("st_lmm_frontend", list(action = "result.as.vector",
    object = x, mode = mode, dots = list()), PACKAGE = "SelectionTools")

as.list.st.lmm.result <- function(x, ...) .Call("st_lmm_frontend", list(action = "result.as.list", object = x,
    dots = list(...)), PACKAGE = "SelectionTools")

`[.st.lmm.result` <- function(x, ..., drop = FALSE) .Call("st_lmm_frontend", list(action = "result.subset",
    object = x, indices = substitute(list(...)), envir = parent.frame(), drop = drop), PACKAGE = "SelectionTools")

`[[.st.lmm.result` <- function(x, ..., exact = TRUE) .Call("st_lmm_frontend", list(action = "result.extract",
    object = x, indices = substitute(list(...)), envir = parent.frame(), exact = exact), PACKAGE = "SelectionTools")

t.st.lmm.result <- function(x) .Call("st_lmm_frontend", list(action = "result.transpose", object = x),
    PACKAGE = "SelectionTools")

print.st.lmm <- function(x, digits = NULL, top = NULL, bottom = NULL, char.per.line = NULL, ...) invisible(.st_lmm_method(x,
    "print", digits = digits, top = top, bottom = bottom, char.per.line = char.per.line, ...))

summary.st.lmm <- function(object, selection = TRUE, ...) .st_lmm_method(object, "summary", selection = selection,
    ...)

print.st.lmm.selection.summary <- function(x, digits = NULL, top = NULL, bottom = NULL, char.per.line = NULL,
    ...) invisible(.st_lmm_method(x, "print.selection", digits = digits, top = top, bottom = bottom,
    char.per.line = char.per.line, ...))

plot.st.lmm <- function(x, scale = NULL, labels = NULL, main = NULL, col = NULL, which = NULL, type = "conditional",
    id.n = 3L, ...) invisible(.Call("st_lmm_frontend", list(action = "plot.dispatch", object = x, scale = scale,
    labels = labels, main = main, col = col, which = which, type = type, id.n = id.n, dots = list(...)),
    PACKAGE = "SelectionTools"))

plot.st.lmm.selection.summary <- function(x, scale = NULL, labels = NULL, main = NULL, col = NULL, ...) invisible(.st_lmm_method(x,
    "plot.selection", scale = scale, labels = labels, main = main, col = col, ...))

anova.st.lmm <- function(object, effects = "fixed", ..., refit = TRUE, test = "Chisq", model.names = NULL) .st_lmm_method(object,
    "anova", effects = effects, extra = list(...), refit = refit, test = test, model.names = model.names,
    model.expressions = substitute(list(object, effects, ...)), test.supplied = !missing(test))

drop1.st.lmm <- function(object, scope, test = "none", k = 2, trace = FALSE, ..., refit = TRUE) .st_lmm_method(object,
    "drop1", scope = if (missing(scope)) NULL else scope, test = test, k = k, trace = trace, extra = list(...),
    refit = refit)

coef.st.lmm <- function(object, id = NULL, vcov = FALSE, ...) .st_lmm_method(object, "coef", id = id,
    vcov = vcov, ...)

vcov.st.lmm <- function(object, adjusted = NULL, id = NULL, ...) .st_lmm_method(object, "vcov", adjusted = adjusted,
    id = id, ...)

fitted.st.lmm <- function(object, ..., type = "response") .st_lmm_method(object, "fitted", type = type)

residuals.st.lmm <- function(object, type = "response", ...) .st_lmm_method(object, "residuals", type = type)

hatvalues.st.lmm <- function(model, type = "conditional", ...) .st_lmm_method(model, "hatvalues", type = type)

rstandard.st.lmm <- function(model, type = "conditional", ...) .st_lmm_method(model, "rstandard", type = type)

st.lmm.diagnostics <- function(object, type = "conditional") .st_lmm_method(object, "diagnostics", type = type)

influence.st.lmm <- function(model, group = NULL, units = NULL, control = list(), keep.fits = FALSE, ...) .Call("st_lmm_frontend", list(action = "influence", object = model,
    group = group, units = units, control = control, keep.fits = keep.fits,
    dots = list(...), call = match.call()), PACKAGE = "SelectionTools")

rstudent.st.lmm <- function(model, influence = NULL, control = list(), ...) .st_lmm_method(model, "rstudent",
    influence = influence, control = control)

.st_lmm_influence_inventory <- function(object, group = NULL) .st_lmm_method(object, "influence.inventory",
    group = group)

.st_lmm_delete_unit <- function(object, inventory, id, control = list(), keep.fit = FALSE) .st_lmm_method(object,
    "delete.unit", inventory = inventory, id = id, control = control, keep.fit = keep.fit)

.st_lmm_influence_collect <- function(object, inventory, results, control = list()) .st_lmm_method(object,
    "influence.collect", inventory = inventory, results = results, control = control)

.st_lmm_diagnostic_capture <- function(object, action, request) .Call("st_lmm_frontend", list(action = "capture.diagnostic",
    object = object, method = action, request = request), PACKAGE = "SelectionTools")

summary.st.lmm.influence <- function(object, ...) .st_lmm_result(object$summary)

print.st.lmm.influence <- function(x, ...) invisible(.Call("st_lmm_frontend", list(action = "print.influence",
    object = x, dots = list(...)), PACKAGE = "SelectionTools"))

plot.st.lmm.influence <- function(x, which = c("cook", "dfbetas"), parm = NULL, labels = NULL, id.n = 3L,
    main = NULL, col = NULL, ...) invisible(.st_lmm_method(x, "plot.influence", which = which, parm = parm,
    labels = labels, id.n = id.n, main = main, col = col, ...))

cooks.distance.st.lmm.influence <- function(model, ...) .st_lmm_result(model$Cook)

dfbeta.st.lmm.influence <- function(model, ...) .st_lmm_result(model$fixed.change)

dfbetas.st.lmm.influence <- function(model, ...) .st_lmm_result(model$dfbetas)

cooks.distance.st.lmm <- function(model, influence = NULL, ...) .Call("st_lmm_frontend", list(action = "measure",
    object = model, influence = influence, measure = "Cook", dots = list(...)), PACKAGE = "SelectionTools")

dfbetas.st.lmm <- function(model, influence = NULL, ...) .Call("st_lmm_frontend", list(action = "measure",
    object = model, influence = influence, measure = "dfbetas", dots = list(...)), PACKAGE = "SelectionTools")

dfbeta.st.lmm <- function(model, influence = NULL, ...) .Call("st_lmm_frontend", list(action = "measure",
    object = model, influence = influence, measure = "fixed.change", dots = list(...)), PACKAGE = "SelectionTools")

logLik.st.lmm <- function(object, ...) .st_lmm_method(object, "logLik")

nobs.st.lmm <- function(object, ...) .st_lmm_method(object, "nobs")

df.residual.st.lmm <- function(object, ...) .st_lmm_method(object, "df.residual")

formula.st.lmm <- function(x, id = NULL, ...) .st_lmm_method(x, "formula", id = id, ...)

model.matrix.st.lmm <- function(object, ...) .st_lmm_method(object, "model.matrix")

predict.st.lmm <- function(object, newdata = NULL, type = "marginal", se.fit = FALSE, allow.new.levels = FALSE,
    new.levels = NULL,
    id = NULL, offset = NULL, re.form = NULL, newparams = NULL, ...) .st_lmm_method(object, "predict",
    newdata = newdata, type = type, se.fit = se.fit, allow.new.levels = allow.new.levels,
    new.levels = new.levels, id = id, offset = offset,
    offset.supplied = !missing(offset), re.form = re.form, re.form.supplied = !missing(re.form), newparams = newparams,
    ...)

confint.st.lmm <- function(object, parm, level = 0.95, method = "Wald", ...) .Call("st_lmm_frontend",
    list(action = "confint", object = object, parm = if (missing(parm)) NULL else parm, level = level,
        method = method, dots = list(...)), PACKAGE = "SelectionTools")

AIC.st.lmm <- function(object, ..., k = 2) .Call("st_lmm_frontend", list(action = "AIC", object = object,
    dots = list(...), k = k, expressions = substitute(list(object, ...))), PACKAGE = "SelectionTools")

BIC.st.lmm <- function(object, ...) .Call("st_lmm_frontend", list(action = "BIC", object = object, dots = list(...),
    expressions = substitute(list(object, ...))), PACKAGE = "SelectionTools")

st.lsm <- function(object, specs, by = NULL, at = list(), weights = "equal", ddf = NULL, level = 0.95,
    adjust = "none", type = "link", transform = NULL, trend = NULL, delta.var = NULL, cov.reduce = mean,
    cov.keep = "2", rg.limit = 10000L, nesting = NULL, bias.adjust = FALSE, sigma = NULL, regrid = NULL,
    max.degree = 1L, L = NULL, rhs = 0, joint = FALSE, side = "two.sided", test = "F", DF = NULL, vcov. = NULL,
    mvt.control = NULL, calpha = NULL, method = NULL, ref = NULL, wts = NULL, scores = NULL, reverse = FALSE,
    infer = NULL, inference.space = "broad") .Call("st_lsm",
    list(object = object, specs = if (missing(specs)) NULL else specs,
    specs.supplied = !missing(specs), by = by, at = at, weights = weights, ddf = ddf, level = level,
    adjust = adjust, type = type, transform = transform, trend = trend, delta.var = delta.var, cov.reduce = cov.reduce,
    cov.keep = cov.keep, rg.limit = rg.limit, nesting = nesting, bias.adjust = bias.adjust, sigma = sigma,
    regrid = regrid, max.degree = max.degree, L = L, rhs = rhs, joint = joint, rhs.supplied = !missing(rhs),
    side = side, test = test, DF = DF, vcov. = vcov., mvt.control = mvt.control, calpha = calpha, method = method,
    ref = ref, wts = wts, scores = scores, reverse = reverse, infer = infer,
    inference.space = inference.space, inference.space.supplied = !missing(inference.space),
    contrast.supplied = list(level = !missing(level),
        adjust = !missing(adjust), type = !missing(type), bias.adjust = !missing(bias.adjust), max.degree = !missing(max.degree))),
    PACKAGE = "SelectionTools")

.st_lsm_method <- function(object, action, ...) .Call("st_lsm_method", list(object = object, action = action,
    ...), PACKAGE = "SelectionTools")

print.st.lsm <- function(x, digits = NULL, top = NULL, bottom = NULL, char.per.line = NULL, tests = NULL, ...) invisible(.st_lsm_method(x,
    "print", digits = digits, top = top, bottom = bottom, char.per.line = char.per.line, tests = tests, ...))

print.summary.st.lsm <- function(x, digits = NULL, top = NULL, bottom = NULL, char.per.line = NULL, tests = NULL, ...) invisible(.st_lsm_method(x,
    "print.summary", digits = digits, top = top, bottom = bottom, char.per.line = char.per.line, tests = tests, ...))

summary.st.lsm <- function(object, level = NULL, adjust = NULL, type = NULL, rhs = 0, side = "two.sided",
    mvt.control = NULL, bias.adjust = NULL, sigma = NULL, DF = NULL, calpha = NULL, ...) .st_lsm_method(object,
    "summary", level = level, adjust = adjust, type = type, null = if (missing(rhs)) NULL else rhs, side = if (missing(side)) NULL else side,
    mvt.control = mvt.control, bias.adjust = bias.adjust, sigma = sigma, DF = DF, calpha = calpha)

as.data.frame.st.lsm <- function(x, ...) summary.st.lsm(x, ...)

coef.st.lsm <- function(object, ..., rhs = FALSE) .st_lsm_method(object, "coef", rhs = rhs)

vcov.st.lsm <- function(object, ...) .st_lsm_method(object, "vcov")

confint.st.lsm <- function(object, parm, level = 0.95, ...) summary.st.lsm(object, level = level, ...)

pairs.st.lsm <- function(x, reverse = FALSE, ...) st.lsm(x, method = if (reverse) "revpairwise" else "pairwise",
    ...)

st.lsm.cld <- function(object, alpha = 0.05, adjust = "tukey", sort = TRUE, reversed = FALSE) .st_lsm_method(object,
    "cld", alpha = alpha, adjust = adjust, sort = sort, reversed = reversed)

st.lsm.test <- function(object, joint = TRUE, by = NULL, rhs = 0, ...) .st_lsm_method(object, "test",
    joint = joint, by = by, null = if (missing(rhs)) NULL else rhs, ...)

st.lsm.mvcontrast <- function(object, method = "eff", mult.name = "trait", by = NULL, adjust = "sidak",
    rhs = 0, show.ests = FALSE, ...) .st_lsm_method(object, "mvcontrast", method = method, mult.name = mult.name,
    by = by, adjust = adjust, null = rhs, show.ests = show.ests, ...)

print.st.lsm.mvcontrast <- function(x, digits = NULL, top = NULL, bottom = NULL, char.per.line = NULL, tests = NULL, ...) invisible(.st_lsm_method(x,
    "print.mvcontrast", digits = digits, top = top, bottom = bottom, char.per.line = char.per.line, tests = tests, ...))

st.predict <- function(object, term = NULL, level = 0.95, newdata = NULL) .Call("st_lmm_frontend", list(action = "st.predict",
    object = object, term = term, level = level, newdata = newdata), PACKAGE = "SelectionTools")

.st_predict_method <- function(object, action, ...) .Call("st_predict_method", list(object = object,
    action = action, ...), PACKAGE = "SelectionTools")

print.st.predict <- function(x, digits = NULL, top = NULL, bottom = NULL, char.per.line = NULL, ...) invisible(.st_predict_method(x,
    "print", digits = digits, top = top, bottom = bottom, char.per.line = char.per.line, ...))

summary.st.predict <- function(object, level = NULL, ...) .st_predict_method(object, "summary", level = level)

print.summary.st.predict <- print.summary.st.lsm

as.data.frame.st.predict <- function(x, ...) summary.st.predict(x, ...)

coef.st.predict <- function(object, ...) .st_predict_method(object, "coef")

vcov.st.predict <- function(object, ...) .st_predict_method(object, "vcov")

confint.st.predict <- function(object, parm, level = 0.95, ...) summary.st.predict(object, level = level,
    ...)

st.kernel.structure <- function(K, factors = NULL, index = NULL, dimensions = NULL, tolerance = 1e-10) {
    if (!is.null(index))
        index <- as.matrix(index)
    .Call("st_kernel_structure", list(K = K, factors = factors, index = index, dimensions = dimensions,
        tolerance = tolerance), PACKAGE = "SelectionTools")
}

st.get.A <- function(data.set = "default", coding = "0 1 2") .Call("st_get_A", data.set, coding, PACKAGE = "SelectionTools")

st.ocs <- function(merit, K, theta, lower = 0, upper = 1, control = list()) .Call("st_lmm_frontend",
    list(action = "ocs", merit = merit, K = K, theta = theta, lower = lower, upper = upper,
        control = control), PACKAGE = "SelectionTools")

print.st.ocs <- function(x, digits = NULL, top = NULL, bottom = NULL, char.per.line = NULL, ...) {
    .Call("st_lmm_frontend", list(action = "ocs.print", object = x, digits = digits, top = top,
        bottom = bottom, char.per.line = char.per.line), PACKAGE = "SelectionTools")
    invisible(x)
}

st.oma <- function(parents, matings, K, dF, base = "current", control = list()) .Call("st_lmm_frontend",
    list(action = "oma", parents = parents, matings = matings, K = K, dF = dF,
        base = base, control = control), PACKAGE = "SelectionTools")

print.st.oma <- function(x, digits = NULL, top = NULL, bottom = NULL, char.per.line = NULL, ...) {
    .Call("st_lmm_frontend", list(action = "oma.print", object = x, digits = digits,
        top = top, bottom = bottom, char.per.line = char.per.line), PACKAGE = "SelectionTools")
    invisible(x)
}

propsd <- function(G, parents = NULL, counts = NULL) .Call("st_lmm_frontend",
    list(action = "propsd", G = G, parents = parents, counts = counts), PACKAGE = "SelectionTools")

print.propsd <- function(x, digits = NULL, ...) {
    .Call("st_lmm_frontend", list(action = "propsd.print", object = x, digits = digits),
        PACKAGE = "SelectionTools")
    invisible(x)
}

st.cross.plan <- function(merit, K, n.crosses, theta, n.progeny = 100, lower = 0, upper = 0.5,
    pairs = NULL, contributions = NULL, mode = "joint", control = list(), scores = NULL,
    criterion = "mean", accuracy = 1, diversity = NULL, constraints = NULL,
    progeny = NULL, genetic.model = NULL) .Call("st_lmm_frontend",
        list(action = "cross.plan", merit = merit, K = K, n.crosses = n.crosses, theta = theta,
            n.progeny = n.progeny, lower = lower, upper = upper, pairs = pairs,
            contributions = contributions, mode = mode, control = control, scores = scores,
            criterion = criterion, accuracy = accuracy, diversity = diversity,
            constraints = constraints, progeny = progeny, genetic.model = genetic.model,
            n.progeny.supplied = !missing(n.progeny)), PACKAGE = "SelectionTools")

st.cross.frontier <- function(merit, K, n.crosses, theta, n.progeny = 100, lower = 0, upper = 0.5,
    pairs = NULL, contributions = NULL, mode = "joint", control = list(), scores = NULL,
    criterion = "mean", accuracy = 1, diversity = NULL, constraints = NULL,
    progeny = NULL, genetic.model = NULL, limits, metric = "theta", warmstart = FALSE) .Call("st_lmm_frontend",
        list(action = "cross.frontier", merit = merit, K = K, n.crosses = n.crosses, theta = theta,
            n.progeny = n.progeny, lower = lower, upper = upper, pairs = pairs,
            contributions = contributions, mode = mode, control = control, scores = scores,
            criterion = criterion, accuracy = accuracy, diversity = diversity,
            constraints = constraints, progeny = progeny, genetic.model = genetic.model,
            n.progeny.supplied = !missing(n.progeny), limits = limits, metric = metric,
            warmstart = warmstart), PACKAGE = "SelectionTools")

print.st.cross.frontier <- function(x, digits = NULL, top = NULL, bottom = NULL, char.per.line = NULL, ...) {
    .Call("st_lmm_frontend", list(action = "cross.frontier.print", object = x, digits = digits, top = top,
        bottom = bottom, char.per.line = char.per.line), PACKAGE = "SelectionTools")
    invisible(x)
}

print.st.cross.plan <- function(x, digits = NULL, top = NULL, bottom = NULL, char.per.line = NULL, ...) {
    .Call("st_lmm_frontend", list(action = "cross.plan.print", object = x, digits = digits, top = top,
        bottom = bottom, char.per.line = char.per.line), PACKAGE = "SelectionTools")
    invisible(x)
}

st.hybrid.data.set <- function(markerdata1, markerdata2, hybrid.performance, hybrid.data.set) invisible(.Call("st_hybrid_data_set",
    list(action = "create", markerdata1 = markerdata1, markerdata2 = markerdata2, hybrid.performance = hybrid.performance,
        hybrid.data.set = hybrid.data.set), PACKAGE = "SelectionTools"))

st.hybrid.data.set.list <- function(digits = getOption("digits")) .Call("st_lmm_frontend", list(action = "hybrid.list",
    digits = digits), PACKAGE = "SelectionTools")

print.st.hybrid.data.set.list <- function(x, digits = NULL, ...) invisible(.Call("st_lmm_frontend",
    list(action = "hybrid.print", object = x, digits = digits, dots = list(...)), PACKAGE = "SelectionTools"))

st.hybrid.data.set.remove <- function(hybrid.data.set) invisible(.Call("st_hybrid_data_set", list(action = "remove",
    hybrid.data.set = hybrid.data.set), PACKAGE = "SelectionTools"))

st.hybrid.data.set.remove.all <- function() invisible(.Call("st_hybrid_data_set", list(action = "remove.all"),
    PACKAGE = "SelectionTools"))

st.set.hybrid.performance.data <- function(hybrid.performance, hybrid.data.set) invisible(.Call("st_hybrid_data_set",
    list(action = "update", hybrid.performance = hybrid.performance, hybrid.data.set = hybrid.data.set),
    PACKAGE = "SelectionTools"))

st.get.hybrid.performance.data <- function(hybrid.data.set) .Call("st_hybrid_data_set", list(action = "get",
    hybrid.data.set = hybrid.data.set), PACKAGE = "SelectionTools")

st.repack.data.set <- function(data.set = "default", commit = TRUE) .Call("st_allele_repack", data.set,
    commit, PACKAGE = "SelectionTools")

model.frame.st.lmm <- function(formula, ...) .st_lmm_method(formula, "model.frame")

terms.st.lmm <- function(x, ...) .st_lmm_method(x, "terms")

weights.st.lmm <- function(object, ...) .st_lmm_method(object, "weights")

update.st.lmm <- function(object, formula., ..., evaluate = TRUE) .st_lmm_method(object, "update", formula = if (missing(formula.)) NULL else formula.,
    extra = substitute(list(...)), evaluate = evaluate, envir = parent.frame())

st.lmm.components <- function(object, name = NULL, traits = NULL, block = NULL,
    term = NULL) .st_lmm_method(object, "components", name = name, traits = traits,
    block = block, term = term)

st.lmm.singular <- function(object, tol = 1e-04) .st_lmm_method(object, "singular", tol = tol)

st.lmm.pca <- function(object, tol = 1e-07) .st_lmm_method(object, "pca", tol = tol)

st.lmm.refit <- function(object, response, weights, method = NULL, control = list(), warm = TRUE) .st_lmm_method(object,
    "refit", response = if (missing(response)) NULL else response, response.supplied = !missing(response),
    weights = if (missing(weights)) NULL else weights, weights.supplied = !missing(weights), method = method,
    control = control, warm = warm, call = match.call())

st.lmm.objective <- function(object, parameters = NULL, method = NULL, gradient = TRUE) .st_lmm_method(object,
    "objective", parameters = parameters, method = method, gradient = gradient)

profile.st.lmm <- function(fitted, ...) st.lmm.profile(fitted, ...)

confint.st.lmm.profile <- function(object, parm, level = object$level, ...) .st_lmm_method(object, "profile.confint",
    parm = if (missing(parm)) NULL else parm, level = level)

print.st.lmm.profile <- function(x, ...) invisible(.st_lmm_method(x, "profile.print", ...))

plot.st.lmm.profile <- function(x, which = NULL, ...) invisible(.st_lmm_method(x, "profile.plot", which = which,
    dots = list(...)))

.st_lmm_profile_prepare <- function(object, level = 0.95, control = list()) .st_lmm_method(object, "profile.prepare",
    parm = character(), level = level, control = control, call = match.call())

.st_lmm_profile_task <- function(prepared, parm) .st_lmm_method(prepared$object, "profile.task",
    prepared = prepared, parm = parm, level = prepared$level, control = prepared$control,
    call = match.call())

simulate.st.lmm <- function(object, nsim = 1, seed = NULL, re.form = NA, conditional = NULL,
    newdata = NULL, newparams = NULL,
    allow.new.levels = FALSE, offset = NULL, weights = NULL, ...) .Call("st_lmm_frontend", list(action = "simulate",
    object = object, nsim = nsim, seed = seed, re.form = re.form, conditional = conditional,
    newdata = newdata, newparams = newparams,
    allow.new.levels = allow.new.levels, offset = offset, offset.supplied = !missing(offset), weights = weights,
    weights.supplied = !missing(weights), dots = list(...)), PACKAGE = "SelectionTools")

st.lmm.simulate <- function(fixed, random = NULL, cov = ~units, data, parameters, nsim = 1, seed = NULL,
    weights = NULL, offset = NULL, ucov = NULL, contrasts = NULL, traits = NULL, ...) .Call("st_lmm_frontend",
    list(action = "simulate.specification", fixed = fixed, random = random, cov = cov, data = substitute(data),
        parameters = parameters, nsim = nsim, seed = seed, weights = substitute(weights), offset = substitute(offset),
        ucov = substitute(ucov), contrasts = contrasts, traits = traits, envir = parent.frame(), dots = list(...)),
    PACKAGE = "SelectionTools")

st.lmm.bootstrap <- function(object, statistic = stats::coef, nsim = 999L, seed = NULL, type = c("parametric",
    "residual"), re.form = NA, statistic.args = list(), responses = NULL, control = list(), series = NULL) .Call("st_lmm_frontend", list(action = "bootstrap", object = object, statistic = statistic,
    nsim = nsim, seed = seed, type = type, re.form = re.form, statistic.args = statistic.args, responses = responses,
    control = control, series = series, call = match.call()), PACKAGE = "SelectionTools")

summary.st.lmm.bootstrap <- function(object, ...) .st_lmm_method(object, "bootstrap.summary", parm = NULL)

print.st.lmm.bootstrap <- function(x, ...) invisible(.Call("st_lmm_frontend", list(action = "print.bootstrap",
    object = x, dots = list(...)), PACKAGE = "SelectionTools"))

confint.st.lmm.bootstrap <- function(object, parm, level = 0.95, type = c("percentile", "basic", "normal"),
    allow.failed = FALSE, ...) .Call("st_lmm_frontend", list(action = "bootstrap.intervals", object = object,
    parm = if (missing(parm)) NULL else parm, level = level, type = type, allow.failed = allow.failed,
    dots = list(...)), PACKAGE = "SelectionTools")

st.lmm.allfit <- function(object, optimizers = NULL, starts = c("fitted", "default"), control = list(), keep.fits = TRUE) .Call("st_lmm_frontend", list(action = "allfit",
    object = object, optimizers = optimizers, starts = starts, control = control, keep.fits = keep.fits, call = match.call()), PACKAGE = "SelectionTools")

summary.st.lmm.allfit <- function(object, ...) .st_lmm_method(object, "repeat.allfit.summary")

print.st.lmm.allfit <- function(x, ...) invisible(.Call("st_lmm_frontend", list(action = "print.allfit",
    object = x, dots = list(...)), PACKAGE = "SelectionTools"))

st.lmm.profile <- function(object, parm = NULL, level = 0.95, control = list()) .Call("st_lmm_frontend",
    list(action = "profile", object = object, parm = parm, level = level, control = control, call = match.call()), PACKAGE = "SelectionTools")

st.lmm.by <- function(fixed, group, data, random = NULL, cov = ~units, method = "REML", ddf = "Satterthwaite",
    type = 3L, random.test = "none", control = list(), ucov = NULL, weights = NULL, contrasts = NULL,
    subset = NULL, na.action = NULL, offset = NULL, pool = FALSE) .Call("st_lmm_frontend",
    list(action = "by", request = list(call = match.call(),
        fixed = fixed, group = substitute(group), data = substitute(data), random = random, cov = cov,
        method = method, ddf = ddf, type = type, random.test = random.test, control = control, ucov = substitute(ucov),
        weights = substitute(weights), contrasts = contrasts, subset = substitute(subset), na.action = na.action,
        offset = substitute(offset), pool = pool, envir = parent.frame(), ddf.supplied = !missing(ddf),
        random.test.supplied = !missing(random.test), weights.supplied = !missing(weights))), PACKAGE = "SelectionTools")

.st_lmm_by_prepare <- function(fixed, group, data, random = NULL, cov = ~units, method = "REML", ddf = "Satterthwaite",
    type = 3L, random.test = "none", control = list(), ucov = NULL, weights = NULL, contrasts = NULL,
    subset = NULL, na.action = NULL, offset = NULL, pool = FALSE) .st_lmm_method(NULL, "by.prepare",
    request = list(call = match.call(), fixed = fixed, group = substitute(group), data = substitute(data),
        random = random, cov = cov, method = method, ddf = ddf, type = type, random.test = random.test,
        control = control, ucov = substitute(ucov), weights = substitute(weights), contrasts = contrasts,
        subset = substitute(subset), na.action = na.action, offset = substitute(offset), pool = pool,
        envir = parent.frame(), ddf.supplied = !missing(ddf), random.test.supplied = !missing(random.test),
        weights.supplied = !missing(weights)))

.st_lmm_by_task <- function(prepared, id) .st_lmm_method(prepared, "by.task", id = id)

.st_lmm_by_collect <- function(prepared, results) .st_lmm_method(prepared, "by.collect", results = results)

.st_lmm_by_capture <- function(args) .Call("st_lmm_frontend", list(action = "capture.method", request = args),
    PACKAGE = "SelectionTools")

summary.st.lmm.by <- function(object, pool = NULL, level = 0.95, ...) .st_lmm_method(object, "by.summary",
    pool = pool, level = level)

coef.st.lmm.by <- function(object, ...) .st_lmm_method(object, "by.coef")

confint.st.lmm.by <- function(object, parm, level = 0.95, pool = NULL, ...) .st_lmm_method(object, "by.confint",
    parm = if (missing(parm)) NULL else parm, level = level, pool = pool)

print.st.lmm.by <- function(x, digits = NULL, top = NULL, bottom = NULL, char.per.line = NULL, ...) invisible(.st_lmm_method(x,
    "by.print", digits = digits, top = top, bottom = bottom, char.per.line = char.per.line, ...))

print.summary.st.lmm.by <- print.st.lmm.by

st.marker.results <- function(data.set = "default") .Call("st_marker_results", as.character(data.set),
    PACKAGE = "SelectionTools")

st.generate.phenotypes <- function(data.set, values, h2 = 0.8, trait = "y") invisible(.Call("st_marker_phenotypes",
    list(data.set = as.character(data.set), values = values, h2 = h2, trait = as.character(trait)), PACKAGE = "SelectionTools"))


st.diallel.cross <- function(mother, father) .Call("st_lmm_frontend",
    list(action = "diallel.cross", mother = mother, father = father),
    PACKAGE = "SelectionTools")


st.gwas <- function(fixed, random, data, marker.data.set = "default",
    K = NULL, n.PC = 0L, min.MAF = 0.05, P3D = TRUE, max.missing = 1,
    missing = "omit", map = NULL, rcov = NULL, control = list(),
    sim.threshold = list())
    .Call("st_lmm_frontend", list(action = "gwas", fixed = fixed,
        random = random, data = data, marker.data.set = marker.data.set,
        K = substitute(K), envir = parent.frame(),
        n.PC = n.PC, min.MAF = min.MAF, P3D = P3D,
        max.missing = max.missing, missing = missing, map = map, rcov = rcov,
        control = control, sim.threshold = sim.threshold),
        PACKAGE = "SelectionTools")

as.data.frame.st.gwas <- function(x, ...) .Call("st_lmm_frontend",
    list(action = "gwas.data.frame", object = x), PACKAGE = "SelectionTools")

print.st.gwas <- function(x, digits = NULL, top = NULL, bottom = NULL,
    char.per.line = NULL, ...) invisible(.Call("st_lmm_frontend",
    list(action = "gwas.print", object = x, digits = digits, top = top,
        bottom = bottom, char.per.line = char.per.line, ...),
    PACKAGE = "SelectionTools"))



plot.st.gwas <- function(x, which = "manhattan", trait = NULL, alpha = 0.05,
    adjust = "BH", threshold = NULL, chromosome = NULL, xlim = NULL,
    position.scale = 1e6, position.unit = "Mb", shade = NULL,
    highlight = NULL, labels = FALSE, level = NULL, main = NULL,
    col = c("grey35", "cornflowerblue"), cex = 0.5, pch = 16, ...)
{
    dots <- list(...)
    kinds <- match.arg(which, c("manhattan", "qq", "region"),
                       several.ok = TRUE)
    adjustment <- match.arg(adjust, c("BH", "bonferroni", "none", "simulation"))
    simulation.explicit <- is.null(threshold) && identical(adjustment, "simulation")
    simulation.use <- simulation.explicit ||
        (any(kinds %in% c("manhattan", "region")) && is.null(threshold) &&
         missing(adjust) && is.list(x) && !is.null(x$simulation))
    alpha.explicit <- !missing(alpha)
    legacy <- !simulation.use && all(which %in% c("manhattan", "qq")) &&
        is.null(threshold) && is.null(chromosome) && is.null(xlim) &&
        is.null(shade) && is.null(highlight) && identical(labels, FALSE) &&
        is.null(level) && is.null(main) && is.null(col) &&
        missing(position.scale) && missing(position.unit) &&
        missing(cex) && missing(pch) && !length(dots)
    if (legacy)
        return(invisible(.Call("st_lmm_frontend",
            list(action = "gwas.plot", object = x, which = which,
                 trait = trait, alpha = alpha, adjust = adjustment),
            PACKAGE = "SelectionTools")))
    if (!inherits(x, "st.gwas") || !is.data.frame(x$results))
        stop("x must be an st.gwas object.")
    which <- kinds
    adjust <- adjustment
    if (!is.numeric(alpha) || length(alpha) != 1L ||
        !is.finite(alpha) || alpha <= 0 || alpha >= 1)
        stop("alpha must be a finite number between zero and one.")
    if (!is.null(threshold) && !identical(threshold, FALSE) &&
        (!is.numeric(threshold) || !length(threshold) ||
         any(!is.finite(threshold) | threshold <= 0 | threshold > 1)))
        stop("threshold must be NULL, FALSE, or raw p-values in (0, 1].")
    if (!is.null(level) && (!is.numeric(level) || length(level) != 1L ||
                           !is.finite(level) || level <= 0 || level >= 1))
        stop("level must be NULL or a finite number between zero and one.")
    if (!is.numeric(position.scale) || length(position.scale) != 1L ||
        !is.finite(position.scale) || position.scale <= 0)
        stop("position.scale must be positive and finite.")
    if (!is.character(position.unit) || length(position.unit) != 1L ||
        is.na(position.unit))
        stop("position.unit must be one character string.")
    interval <- function(z, name) {
        if (!is.null(z) && (!is.numeric(z) || length(z) != 2L ||
                           any(!is.finite(z)) || z[1L] > z[2L]))
            stop(name, " must contain two increasing finite positions.")
    }
    interval(xlim, "xlim")
    interval(shade, "shade")
    if (!is.null(main) && (!is.character(main) || length(main) != 1L ||
                          is.na(main)))
        stop("main must be NULL or one character string.")
    if (!is.null(highlight) && (!is.character(highlight) || anyNA(highlight)))
        stop("highlight must be NULL or a character vector of marker names.")
    if (!(is.logical(labels) && length(labels) == 1L && !is.na(labels)) &&
        !(is.character(labels) && !anyNA(labels)))
        stop("labels must be TRUE, FALSE, or a character vector of marker names.")
    if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) ||
                        anyDuplicated(names(dots)) ||
                        any(names(dots) %in% c("x", "y", "type"))))
        stop("Additional graphical arguments must be uniquely named; x, y and type are reserved.")
    table <- x$results
    needed <- c("trait", "marker", "chromosome", "position", "status",
                "score", "p", "p.BH")
    if (!all(needed %in% names(table)))
        stop("The GWAS results table is incomplete.")
    traits <- as.character(x$traits)
    if (!is.null(trait)) {
        if (!is.character(trait) || length(trait) != 1L || is.na(trait) ||
            !trait %in% traits)
            stop("trait must name one trait in x.")
        traits <- trait
    }
    if ("region" %in% which &&
        (is.null(chromosome) || length(chromosome) != 1L || is.na(chromosome)))
        stop("A regional plot requires one chromosome.")
    palette <- if (is.null(col)) c("grey35", "cornflowerblue") else col
    if (!length(palette) || anyNA(palette))
        stop("col must contain at least one nonmissing colour.")
    draw <- function(arguments) {
        arguments[names(dots)] <- dots
        do.call(graphics::plot, arguments)
    }
    simulated.cutoff <- function(trait.name) {
        simulation <- x$simulation
        thresholds <- if (is.list(simulation)) simulation$thresholds else NULL
        valid <- is.data.frame(thresholds) &&
            all(c("trait", "alpha", "score", "status") %in% names(thresholds)) &&
            "p.sim" %in% names(table)
        if (valid) {
            row <- which(!is.na(thresholds$trait) &
                         thresholds$trait == trait.name)
            valid <- length(row) == 1L
            if (valid) {
                cutoff <- thresholds[row, , drop = FALSE]
                valid <- !is.na(cutoff$status) && cutoff$status == "ok" &&
                    is.numeric(cutoff$score) && is.finite(cutoff$score) &&
                    cutoff$score >= 0 && is.numeric(cutoff$alpha) &&
                    is.finite(cutoff$alpha) && cutoff$alpha > 0 && cutoff$alpha < 1
            }
        }
        if (!valid) {
            message <- paste0("A complete simulated threshold is unavailable for trait '",
                              trait.name, "'.")
            if (simulation.explicit) stop(message, call. = FALSE)
            warning(message, " No significance line is drawn.", call. = FALSE)
            return(NULL)
        }
        if (alpha.explicit && alpha != cutoff$alpha)
            stop("alpha differs from the stored simulated level for trait '",
                 trait.name, "'; use that level or refit with the requested sim.threshold$alpha.",
                 call. = FALSE)
        cutoff
    }
    for (trait.name in traits) {
        cutoff <- if (simulation.use) simulated.cutoff(trait.name) else NULL
        keep <- !is.na(table$trait) & table$trait == trait.name &
            !is.na(table$status) & table$status == "ok" &
            is.finite(table$p) & table$p >= 0 &
            table$p <= 1 & !is.na(table$score)
        result <- table[keep, , drop = FALSE]
        if (!nrow(result)) {
            warning("No successful marker tests for trait ", trait.name, ".")
            next
        }
        cut <- threshold
        cut.alpha <- alpha
        cut.label <- NULL
        if (simulation.use) {
            cut <- FALSE
            cut.score <- numeric()
            if (!is.null(cutoff)) {
                cut.score <- cutoff$score
                cut.alpha <- cutoff$alpha
                cut.label <- paste0("Simulation: alpha = ",
                    format(cut.alpha, digits = 3L, trim = TRUE),
                    ", cutoff = ", format(cut.score, digits = 4L, trim = TRUE))
            }
        } else if (is.null(cut)) {
            cut <- switch(adjust,
                none = alpha,
                bonferroni = alpha / nrow(result),
                BH = {
                    significant <- !is.na(result$p.BH) & result$p.BH <= alpha
                    if (any(significant)) max(result$p[significant]) else FALSE
                })
        }
        if (!simulation.use)
            cut.score <- if (identical(cut, FALSE)) numeric() else
                -log10(pmax(cut, .Machine$double.xmin))
        for (kind in which) {
            title <- if (is.null(main)) trait.name else main
            if (kind == "qq") {
                values <- sort(result$p)
                n <- length(values)
                expected <- -log10((seq_len(n) - 0.5) / n)
                observed <- -log10(pmax(values, .Machine$double.xmin))
                lower <- upper <- numeric()
                if (!is.null(level)) {
                    tail <- (1 - level) / 2
                    lower <- -log10(pmax(stats::qbeta(1 - tail, seq_len(n),
                        n + 1 - seq_len(n)), .Machine$double.xmin))
                    upper <- -log10(pmax(stats::qbeta(tail, seq_len(n),
                        n + 1 - seq_len(n)), .Machine$double.xmin))
                }
                draw(list(x = expected, y = observed, type = "n",
                    xlim = range(0, expected), ylim = range(0, observed, upper),
                    xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
                    main = title, bty = "l"))
                if (!is.null(level))
                    graphics::polygon(c(expected, rev(expected)),
                        c(lower, rev(upper)), border = NA, col = "grey90")
                graphics::abline(0, 1, col = "coral", lty = 2)
                graphics::points(expected, observed, pch = pch,
                                 cex = cex, col = palette[1L])
                next
            }
            current <- result
            axis.at <- axis.labels <- NULL
            if (kind == "region") {
                current <- current[!is.na(current$chromosome) &
                    as.character(current$chromosome) == as.character(chromosome) &
                    is.finite(current$position), , drop = FALSE]
                span <- xlim
                if (is.null(span) && nrow(current)) span <- range(current$position)
                if (!is.null(span))
                    current <- current[current$position >= span[1L] &
                                       current$position <= span[2L], , drop = FALSE]
                if (!nrow(current))
                    stop("No successful mapped marker tests in the requested region.")
                horizontal <- current$position / position.scale
                horizontal.limit <- span / position.scale
                colour <- rep(palette[1L], nrow(current))
                if (simulation.use && length(cut.score))
                    colour[!is.na(current$p.sim) & current$p.sim <= cut.alpha] <- "coral"
                else if (!identical(cut, FALSE))
                    colour[current$p <= max(cut)] <- "coral"
                xlab <- paste("Chromosome", chromosome, "position",
                              paste0("(", position.unit, ")"))
            } else {
                mapped <- !anyNA(current$chromosome) &&
                    all(is.finite(current$position))
                horizontal <- seq_len(nrow(current))
                colour <- rep(palette[1L], nrow(current))
                xlab <- "Marker"
                if (mapped) {
                    chromosomes <- unique(as.character(current$chromosome))
                    numeric.chr <- suppressWarnings(as.numeric(chromosomes))
                    if (all(is.finite(numeric.chr)))
                        chromosomes <- chromosomes[order(numeric.chr)]
                    group <- match(as.character(current$chromosome), chromosomes)
                    low <- vapply(seq_along(chromosomes), function(j)
                        min(current$position[group == j]), numeric(1))
                    high <- vapply(seq_along(chromosomes), function(j)
                        max(current$position[group == j]), numeric(1))
                    span <- high - low
                    gap <- if (sum(span) > 0) 0.02 * sum(span) / length(span) else 1
                    offset <- c(0, utils::head(cumsum(span + gap), -1L))
                    horizontal <- offset[group] + current$position - low[group]
                    colour <- rep(palette, length.out = length(chromosomes))[group]
                    axis.at <- offset + span / 2
                    axis.labels <- chromosomes
                    xlab <- "Chromosome"
                }
                horizontal.limit <- range(horizontal)
            }
            score <- current$score
            nonfinite <- !is.finite(score)
            score[nonfinite] <- -log10(pmax(current$p[nonfinite],
                                           .Machine$double.xmin))
            vertical.max <- max(c(1, score, cut.score))
            arguments <- list(x = horizontal, y = score, type = "n",
                xlim = horizontal.limit, ylim = c(0, 1.12 * vertical.max),
                xlab = xlab, ylab = "-log10(p)", main = title, bty = "l")
            if (!is.null(axis.at)) arguments$xaxt <- "n"
            draw(arguments)
            if (kind == "region" && !is.null(shade)) {
                bounds <- graphics::par("usr")
                graphics::rect(shade[1L] / position.scale, bounds[3L],
                               shade[2L] / position.scale, bounds[4L],
                               col = "grey94", border = NA)
            }
            if (length(cut.score))
                graphics::abline(h = cut.score, col = "coral", lty = 2)
            if (!is.null(cut.label)) {
                bounds <- graphics::par("usr")
                graphics::text(bounds[2L] - 0.02 * (bounds[2L] - bounds[1L]),
                    cut.score, labels = cut.label, adj = c(1, -0.35),
                    col = "coral", cex = 0.75)
            }
            if (!is.null(highlight))
                colour[current$marker %in% highlight] <- "coral"
            graphics::points(horizontal, score, pch = pch, cex = cex, col = colour)
            if (!is.null(axis.at) && !identical(dots$axes, FALSE))
                graphics::axis(1, at = axis.at, labels = axis.labels,
                               tick = FALSE, cex.axis = 0.8)
            marked <- if (isTRUE(labels)) highlight else if (is.character(labels))
                labels else character()
            annotate <- current$marker %in% marked
            if (any(annotate))
                graphics::text(horizontal[annotate], score[annotate],
                               labels = current$marker[annotate], pos = 3, cex = 0.85)
        }
    }
    invisible(x)
}



st.plot.genetic.distances <- function(x, labels = FALSE, cluster = TRUE,
    col = grDevices::heat.colors(200), hclust.method = "complete", ...)
{
    distance <- if (inherits(x, "dist")) as.matrix(x) else x
    if (!is.matrix(distance) || !is.numeric(distance) ||
        nrow(distance) != ncol(distance) || nrow(distance) < 2L ||
        any(!is.finite(distance)) || any(distance < 0) ||
        !isTRUE(all.equal(distance, t(distance), check.attributes = FALSE,
                          tolerance = 1e-10)) ||
        any(abs(diag(distance)) > 1e-10))
        stop("x must be a finite symmetric nonnegative distance matrix with a zero diagonal.")
    if (!is.logical(cluster) || length(cluster) != 1L || is.na(cluster))
        stop("cluster must be TRUE or FALSE.")
    if (is.logical(labels) && length(labels) == 1L && !is.na(labels)) {
        label <- if (labels) rownames(distance) else rep(NA_character_, nrow(distance))
        if (is.null(label)) label <- as.character(seq_len(nrow(distance)))
    } else if (is.character(labels) && length(labels) == nrow(distance)) {
        label <- labels
    } else stop("labels must be TRUE, FALSE, or one character label per individual.")
    dots <- list(...)
    reserved <- c("x", "Rowv", "Colv", "symm", "scale", "col", "labRow", "labCol")
    if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) ||
                        anyDuplicated(names(dots)) || any(names(dots) %in% reserved)))
        stop("Additional heatmap arguments must be uniquely named and must not override distance or clustering arguments.")
    rows <- if (cluster)
        stats::as.dendrogram(stats::hclust(stats::as.dist(distance),
                                         method = hclust.method)) else NA
    arguments <- c(list(x = distance, Rowv = rows,
        Colv = if (cluster) "Rowv" else NA, symm = TRUE, scale = "none",
        col = col, labRow = label, labCol = label), dots)
    invisible(do.call(stats::heatmap, arguments))
}

st.plot.pcoa <- function(x, which = c("scores", "variance", "cumulative"),
    axes = c(1L, 2L), groups = NULL, col = NULL, pch = NULL,
    legend = TRUE, cex = 0.8, ...)
{
    which <- match.arg(which)
    if (!is.list(x) || !is.matrix(x$points) || !is.numeric(x$points) ||
        !is.numeric(x$eig) || !length(x$eig) || any(!is.finite(x$eig)))
        stop("x must be a cmdscale result with points and finite eigenvalues (eig = TRUE).")
    positive <- pmax(x$eig, 0)
    if (sum(positive) <= 0)
        stop("The principal-coordinate result has no positive eigenvalues.")
    percent <- 100 * positive / sum(positive)
    dots <- list(...)
    if (length(dots) && (is.null(names(dots)) || any(!nzchar(names(dots))) ||
                        anyDuplicated(names(dots)) ||
                        any(names(dots) %in% c("x", "y", "col", "pch", "cex"))))
        stop("Additional graphical arguments must be uniquely named; x, y, col, pch and cex are reserved.")
    if (!is.logical(legend) || length(legend) != 1L || is.na(legend))
        stop("legend must be TRUE or FALSE.")
    if (which == "scores") {
        if (!is.numeric(axes) || length(axes) != 2L ||
            any(!is.finite(axes) | axes != as.integer(axes)) ||
            any(axes < 1 | axes > ncol(x$points)) || anyDuplicated(axes))
            stop("axes must name two different available principal coordinates.")
        if (any(!is.finite(x$points[, axes, drop = FALSE])))
            stop("The requested principal coordinates contain nonfinite scores.")
        has.groups <- !is.null(groups)
        if (has.groups) {
            if (length(groups) != nrow(x$points) || anyNA(groups))
                stop("groups must contain one nonmissing group per individual.")
            if (!is.null(names(groups)) && !is.null(rownames(x$points))) {
                if (anyDuplicated(names(groups)) ||
                    !all(rownames(x$points) %in% names(groups)))
                    stop("Named groups must cover each principal-coordinate row exactly once.")
                groups <- groups[rownames(x$points)]
            }
            groups <- if (is.factor(groups)) droplevels(groups) else factor(groups)
        } else groups <- factor(rep("All", nrow(x$points)))
        colours <- if (is.null(col)) c("cornflowerblue", "coral", "darkgreen",
                                      "orchid4", "goldenrod3") else col
        symbols <- if (is.null(pch)) c(16, 17, 15, 18, 8) else pch
        if (!length(colours) || !length(symbols))
            stop("col and pch must each contain at least one value.")
        colours <- rep(colours, length.out = nlevels(groups))
        symbols <- rep(symbols, length.out = nlevels(groups))
        label <- function(axis) sprintf("PCoA %d (%.2f%%)", axis, percent[axis])
        arguments <- list(x = x$points[, axes[1L]], y = x$points[, axes[2L]],
            xlab = label(axes[1L]), ylab = label(axes[2L]),
            col = colours[as.integer(groups)], pch = symbols[as.integer(groups)],
            cex = cex, bty = "l")
        arguments[names(dots)] <- dots
        do.call(graphics::plot, arguments)
        if (legend && has.groups)
            graphics::legend("topright", legend = levels(groups), col = colours,
                             pch = symbols, title = "Group", bty = "n", cex = 0.85)
    } else {
        values <- if (which == "variance") percent else cumsum(percent)
        arguments <- list(x = seq_along(values), y = values,
            xlab = if (which == "variance") "Principal coordinate" else "Number of coordinates",
            ylab = if (which == "variance") "Explained variance (%)" else "Cumulative variance (%)",
            col = if (is.null(col)) "cornflowerblue" else col[1L],
            pch = if (is.null(pch)) 16 else pch[1L], cex = cex, bty = "l")
        if (which == "cumulative") arguments$ylim <- c(0, 100)
        arguments[names(dots)] <- dots
        do.call(graphics::plot, arguments)
    }
    invisible(list(points = x$points, variance.percent = percent, groups = groups))
}

recover_data.st.lmm <- function(object, ...) .Call("st_lmm_frontend",
    list(action = "adapter.recover", object = object, dots = list(...)), PACKAGE = "SelectionTools")

emm_basis.st.lmm <- function(object, trms, xlev, grid, ddf = NULL, ...) .Call("st_lmm_frontend",
    list(action = "adapter.basis", object = object, grid = grid, ddf = ddf,
        DF.fun = .st_lmm_emm_DF, dots = list(...)), PACKAGE = "SelectionTools")

.st_lmm_emm_DF <- function(k, dfargs) .Call("st_lmm_frontend",
    list(action = "adapter.DF", object = dfargs$object, ddf = dfargs$ddf, k = k), PACKAGE = "SelectionTools")

modelparm.st.lmm <- function(model, coef. = NULL, vcov. = NULL, df = NULL, ...) .Call("st_lmm_frontend",
    list(action = "adapter.modelparm", object = model, coef. = coef., vcov. = vcov., df = df,
        dots = list(...)), PACKAGE = "SelectionTools")

compare.gwas.fits <- function(..., adjust = "FDR", alpha = 0.05,
    trait = NULL) .Call("st_lmm_frontend",
    list(action = "gwas.summary", fits = list(...),
        expressions = substitute(list(...)), comparison = TRUE,
        adjust = adjust, alpha = alpha, trait = trait),
    PACKAGE = "SelectionTools")

summarize.gwas.fit <- function(fit, adjust = "FDR", alpha = 0.05,
    top = 6L, trait = NULL) .Call("st_lmm_frontend",
    list(action = "gwas.summary", fits = list(fit),
        expressions = substitute(list(fit)), comparison = FALSE,
        adjust = adjust, alpha = alpha, top = top, trait = trait),
    PACKAGE = "SelectionTools")

summarize.gwas.simulation <- function(fit, alpha = 0.05, top = 6L,
    trait = NULL) .Call("st_lmm_frontend",
    list(action = "gwas.summary", fits = list(fit),
        expressions = substitute(list(fit)), comparison = FALSE,
        adjust = "simulation", alpha = alpha, top = top, trait = trait),
    PACKAGE = "SelectionTools")

print.st.gwas.summary <- function(x, digits = NULL, char.per.line = NULL, ...)
    invisible(.Call("st_lmm_frontend",
        list(action = "gwas.summary.print", object = x, digits = digits,
            char.per.line = char.per.line, ...), PACKAGE = "SelectionTools"))
