options(stringsAsFactors = FALSE, warn = 1)
set.seed(123)
suppressPackageStartupMessages(library(data.table))
if (!requireNamespace("brglm2", quietly = TRUE)) stop("Missing brglm2")
setDTthreads(percent = 75)
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 5L) stop("Usage: Rscript 05_canada_external_assessment.R RAW_DIR WORK_DIR LOCK_DIR CLASS_MAP FAERS_REFERENCE")
RAW_DIR <- normalizePath(a[1], winslash = "/", mustWork = TRUE)
WORK_DIR <- normalizePath(a[2], winslash = "/", mustWork = FALSE)
LOCK_DIR <- normalizePath(a[3], winslash = "/", mustWork = TRUE)
CLASS_MAP <- normalizePath(a[4], winslash = "/", mustWork = TRUE)
FAERS_REFERENCE <- normalizePath(a[5], winslash = "/", mustWork = TRUE)
PRE <- file.path(WORK_DIR, "canada_readiness")
readiness <- file.path(WORK_DIR, "canada_readiness_results")
dir.create(PRE, recursive = TRUE, showWarnings = FALSE)
dir.create(readiness, recursive = TRUE, showWarnings = FALSE)
logmsg <- function(...) invisible(NULL)
stop0 <- function(...) stop(paste(..., collapse = " "), call. = FALSE)
clean <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("^\\ufeff", "", x)
    x[x %in% c("", "NULL")] <- NA_character_
    x
}
node_norm <- function(x) tolower(gsub("[[:space:]]+", " ", trimws(as.character(x))))
date_parse <- function(x) {
    z <- clean(x)
    o <- as.Date(rep(NA_character_, length(z)))
    i <- !is.na(z) & grepl("^[0-9]{1,2}-[A-Za-z]{3}-[0-9]{2,4}$", z)
    if (any(i)) {
        d <- as.integer(sub("^([0-9]{1,2})-.*$", "\\1", z[i]))
        m <- match(toupper(sub("^[0-9]{1,2}-([A-Za-z]{3})-.*$", "\\1", z[i])), toupper(month.abb))
        yy <- sub("^.*-([0-9]{2,4})$", "\\1", z[i])
        y <- as.integer(yy)
        two <- nchar(yy) == 2L
        y[two & y <= 26] <- y[two & y <= 26] + 2000L
        y[two & y >= 65 & y <= 99] <- y[two & y >= 65 & y <= 99] + 1900L
        o[i] <- as.Date(sprintf("%04d-%02d-%02d", y, m, d))
    }
    o
}
qperiod <- function(d) paste0(format(d, "%Y"), "Q", ((as.integer(format(d, "%m")) - 
    1L)%/%3L) + 1L)
bgroup <- function(x) {
    q <- unique(as.numeric(quantile(as.numeric(x), probs = c(0, 0.5, 0.8, 0.95, 1), 
        na.rm = TRUE, type = 7)))
    if (length(q) < 3L) 
        return(factor(ifelse(is.na(x), "missing", "one_group")))
    z <- as.character(cut(x, breaks = q, include.lowest = TRUE, ordered_result = TRUE))
    z[is.na(z)] <- "missing"
    factor(z)
}
onefile <- function(n) {
    x <- list.files(RAW_DIR, recursive = TRUE, full.names = TRUE)
    x <- x[tolower(basename(x)) == tolower(n)]
    if (length(x) != 1L) 
        stop0("Expected exactly one", n, "found", length(x))
    x
}
readpos <- function(p, cols, nms) {
    z <- fread(p, sep = "$", quote = "\"", header = FALSE, select = cols, fill = TRUE, 
        showProgress = TRUE, encoding = "unknown")
    setnames(z, nms)
    for (n in nms) set(z, j = n, value = clean(z[[n]]))
    z
}
read_report_drug_logical <- function(p) {
    tmp <- tempfile(pattern = "canada_report_drug_reconstructed_", fileext = ".txt")
    ic <- file(p, open = "r", encoding = "UTF-8")
    oc <- file(tmp, open = "w", encoding = "UTF-8")
    on.exit({
        try(close(ic), silent = TRUE)
        try(close(oc), silent = TRUE)
        unlink(tmp, force = TRUE)
    }, add = TRUE)
    carry <- NULL
    nrec <- 0L
    rx <- "^\\s*\"?[0-9]+\"?\\$\"?[0-9]*\"?\\$\"?[0-9]*\"?\\$"
    emit <- function(x) {
        if (length(x)) {
            x <- gsub("[\\r\\n]+", " ", x, perl = TRUE)
            writeLines(x, oc, useBytes = TRUE)
            nrec <<- nrec + length(x)
        }
    }
    repeat {
        lines <- readLines(ic, n = 100000L, warn = FALSE)
        if (!length(lines)) 
            break
        st <- which(grepl(rx, lines, perl = TRUE))
        if (!length(st)) {
            if (is.null(carry)) 
                stop0("Report_Drug logical parser found non-record prefix")
            carry <- paste(c(carry, lines), collapse = " ")
            next
        }
        if (st[1] > 1L) {
            if (is.null(carry)) 
                stop0("Report_Drug logical parser found unassigned prefix")
            carry <- paste(c(carry, lines[seq_len(st[1] - 1L)]), collapse = " ")
        }
        if (!is.null(carry)) {
            emit(carry)
            carry <- NULL
        }
        if (length(st) > 1L) 
            emit(vapply(seq_len(length(st) - 1L), function(k) paste(lines[st[k]:(st[k + 
                1L] - 1L)], collapse = " "), character(1)))
        carry <- paste(lines[st[length(st)]:length(lines)], collapse = " ")
    }
    if (!is.null(carry)) 
        emit(carry)
    close(ic)
    close(oc)
    if (nrec != 5669687L) 
        stop0("Report_Drug reconstructed logical-record anchor failed: observed", 
            nrec, "expected 5669687")
    z <- suppressWarnings(fread(tmp, sep = "$", quote = "\"", header = FALSE, select = c(1, 
        2, 3), fill = TRUE, showProgress = TRUE, encoding = "unknown", colClasses = "character"))
    if (nrow(z) != nrec) 
        stop0("Report_Drug reconstructed parse row count mismatch")
    setnames(z, c("report_drug_id", "report_id", "product_id"))
    for (n in names(z)) set(z, j = n, value = clean(z[[n]]))
    z
}
files <- list(reports = onefile("reports.txt"), report_drug = onefile("report_drug.txt"), 
    indication = onefile("report_drug_indication.txt"), reactions = onefile("reactions.txt"), 
    ingredients = onefile("drug_product_ingredients.txt"), priorities = file.path(LOCK_DIR, 
        "LOCKED_9_PRIORITIES.csv"), class_map = CLASS_MAP)
if (!all(file.exists(unlist(files)))) stop0("Missing required input")
pri <- fread(files$priorities)
need <- c("priority_id", "pair_key", "drug_combination", "event")
if (nrow(pri) != 9L || !all(need %chin% names(pri)) || anyDuplicated(pri$priority_id) || 
    anyDuplicated(pri$pair_key)) stop0("Locked 9-priority check failed")
fwrite(pri[, ..need], file.path(PRE, "04_LOCKED_PRIORITY_check.csv"))
logmsg("readiness start: direction blinded.")
reports_dt <- readpos(files$reports, c(1, 2, 3, 5), c("report_id", "report_no", "version_no", 
    "initial_received"))
reports_dt[, `:=`(date, date_parse(initial_received))]
if (any(!is.na(reports_dt$initial_received) & is.na(reports_dt$date))) stop0("Unparsed nonblank report date")
reports_dt <- unique(reports_dt[!is.na(date) & date >= as.Date("2014-01-01") & date <= 
    as.Date("2026-03-31")], by = "report_id")
ids <- reports_dt$report_id
if (!length(ids)) stop0("No study-window reports")
ind <- readpos(files$indication, c(2, 5), c("report_id", "indication"))
pdinc <- c("PARKINSON DISEASE", "PARKINSONS DISEASE", "PARKINSON'S DISEASE", "IDIOPATHIC PARKINSON'S DISEASE", 
    "PARKINSONISM", "PARKINSONIAN")
ind[, `:=`(u, toupper(indication))]
strict <- unique(ind[report_id %chin% ids & u %chin% pdinc, report_id])
if (!length(strict)) stop0("Strict-PD mapping yielded zero reports")
base <- merge(data.table(report_id = strict), reports_dt[, .(report_id, calendar_period = factor(qperiod(date)))], 
    by = "report_id", all.x = TRUE)
ing <- readpos(files$ingredients, c(2, 5), c("product_id", "ingredient"))
ing[, `:=`(z, node_norm(ingredient))]
cm <- fread(files$class_map)
cm[, `:=`(base_ingredient, node_norm(base_ingredient))]
nodefun <- function(x) {
    o <- node_norm(x)
    o[grepl("levodopa|foslevodopa", o)] <- "levodopa_containing_therapy"
    for (n in c("apomorphine", "pramipexole", "rotigotine", "entacapone", "amantadine")) o[startsWith(o, 
        n)] <- n
    for (j in seq_len(nrow(cm))) o[startsWith(o, cm$base_ingredient[j])] <- tolower(cm$class_node[j])
    o
}
ing[, `:=`(node, nodefun(z))]
pn <- unique(ing[!is.na(product_id) & !is.na(node) & nzchar(node), .(product_id, 
    node)])
if (!nrow(pn)) stop0("Ingredient-node map empty")
rd <- read_report_drug_logical(files$report_drug)
if (anyDuplicated(rd$report_drug_id) || nrow(rd) != 5669687L) stop0("DA validated Report_Drug logical-record invariant failed")
rd <- rd[grepl("^[0-9]+$", report_drug_id) & !is.na(report_id) & grepl("^[0-9]+$", 
    report_id) & !is.na(product_id) & grepl("^[0-9]+$", product_id)]
rd <- merge(rd[report_id %chin% strict, .(report_id, product_id)], pn, by = "product_id", 
    all.x = TRUE, allow.cartesian = TRUE)
rd <- rd[!is.na(node) & nzchar(node)]
nc <- rd[, .(n_nodes = uniqueN(node)), by = report_id]
base <- merge(base, nc, by = "report_id", all.x = TRUE)
base[is.na(n_nodes), `:=`(n_nodes, 0L)]
base[, `:=`(n_nodes_group, bgroup(n_nodes))]
rx <- readpos(files$reactions, c(2, 6), c("report_id", "pt"))
rx <- unique(rx[report_id %chin% strict & !is.na(pt), .(report_id, pt = tolower(pt))])
if (!nrow(rx)) stop0("No strict-PD reaction data")
ev <- function(p) unique(rx[grepl(p, pt, perl = TRUE), report_id])
evs <- list(ICD_compulsive = ev("impulse control|pathological gambling|gambling|compulsive|hypersexual|binge eating|compulsive eating|punding|kleptomania|shopping addiction"), 
    Sleepiness = ev("somnolence|sleepiness|excessive daytime sleepiness|sleep attack|hypersomnia|sudden onset of sleep"), 
    Hallucination_psychosis = ev("hallucination|visual hallucination|auditory hallucination|psychosis|psychotic disorder|delusion|paranoia"), 
    Confusion_delirium = ev("^confusional arousal$|^confusional state$|^delirium$|^delirium febrile$|^dementia of the alzheimer's type, with delirium$|^disorientation$|^intensive care unit delirium$|^postoperative delirium$"))
ac <- rx[, .(total_ae_terms = uniqueN(pt)), by = report_id]
base <- merge(base, ac, by = "report_id", all.x = TRUE)
base[is.na(total_ae_terms), `:=`(total_ae_terms, 0L)]
fwrite(data.table(metric = c("Canada_total_study_reports", "reports_with_indication_information", 
    "strict_PD_reports", "reports_with_drug_records", "reports_with_reaction_records"), 
    value = c(length(ids), uniqueN(ind[report_id %chin% ids, report_id]), length(strict), 
        uniqueN(rd$report_id), uniqueN(rx$report_id))), file.path(PRE, "06_STRICT_PD_COUNTS.csv"))
parsekey <- function(x) {
    z <- strsplit(x, " \\| ")[[1]]
    if (length(z) != 3L) 
        stop0("Invalid pair_key")
    list(A = node_norm(z[1]), B = node_norm(z[2]), D = z[3])
}
maskfit <- function(d, sparse) {
    cv <- c("n_nodes_group", "non_target_ae_terms_group", "calendar_period")
    cv <- cv[vapply(cv, function(x) length(unique(d[[x]][!is.na(d[[x]])])) > 1L, 
        logical(1))]
    f <- as.formula(paste("E ~", paste(c("A", "B", "AB", cv), collapse = " + ")))
    mm <- tryCatch(model.matrix(f, d), error = function(e) NULL)
    if (is.null(mm)) 
        return(list(rank = FALSE, est = NA_character_, done = FALSE, fe = FALSE, 
            fs = FALSE, why = "model_matrix_error"))
    if (qr(mm)$rank != ncol(mm)) 
        return(list(rank = FALSE, est = NA_character_, done = FALSE, fe = FALSE, 
            fs = FALSE, why = "rank_deficient"))
    en <- if (sparse) 
        "brglm2_penalized_logistic_wald_CI"
    else "standard_glm"
    fit <- tryCatch(suppressWarnings(if (sparse) 
        glm(f, d, family = binomial("logit"), method = brglm2::brglmFit, control = brglm2::brglmControl(type = "AS_mean"))
    else glm(f, d, family = binomial("logit"))), error = function(e) NULL)
    if (is.null(fit)) 
        return(list(rank = TRUE, est = en, done = FALSE, fe = FALSE, fs = FALSE, 
            why = "fit_error"))
    co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (is.null(co) || !("AB" %in% rownames(co))) 
        return(list(rank = TRUE, est = en, done = FALSE, fe = FALSE, fs = FALSE, 
            why = "AB_coefficient_unavailable"))
    fe <- is.finite(as.numeric(co["AB", 1]))
    fs <- is.finite(as.numeric(co["AB", 2])) && as.numeric(co["AB", 2]) > 0
    list(rank = TRUE, est = en, done = TRUE, fe = fe, fs = fs, why = if (fe && fs) "" else "nonfinite_interaction_or_se")
}
cov <- vector("list", 9)
ready <- vector("list", 9)
for (i in seq_len(nrow(pri))) {
    p <- parsekey(pri$pair_key[i])
    if (!(p$D %chin% names(evs))) 
        stop0("Locked event unavailable")
    d <- copy(base)
    d[, `:=`(A, as.integer(report_id %chin% unique(rd[node == p$A, report_id])))]
    d[, `:=`(B, as.integer(report_id %chin% unique(rd[node == p$B, report_id])))]
    d[, `:=`(AB, as.integer(A == 1L & B == 1L))]
    d[, `:=`(E, as.integer(report_id %chin% evs[[p$D]]))]
    d[, `:=`(non_target_ae_terms, pmax(0, total_ae_terms - E))]
    d[, `:=`(non_target_ae_terms_group, bgroup(non_target_ae_terms))]
    n <- nrow(d)
    nA <- sum(d$A)
    nB <- sum(d$B)
    nAB <- sum(d$AB)
    nE <- sum(d$E)
    nABE <- sum(d$AB == 1L & d$E == 1L)
    cc <- c(n00_event0 = sum(d$A == 0 & d$B == 0 & d$E == 0), n00_event1 = sum(d$A == 
        0 & d$B == 0 & d$E == 1), n10_event0 = sum(d$A == 1 & d$B == 0 & d$E == 0), 
        n10_event1 = sum(d$A == 1 & d$B == 0 & d$E == 1), n01_event0 = sum(d$A == 
            0 & d$B == 1 & d$E == 0), n01_event1 = sum(d$A == 0 & d$B == 1 & d$E == 
            1), n11_event0 = sum(d$A == 1 & d$B == 1 & d$E == 0), n11_event1 = sum(d$A == 
            1 & d$B == 1 & d$E == 1))
    if (sum(cc) != n) 
        stop0("2x2 cells fail sum check for", pri$priority_id[i])
    fr <- maskfit(d, nABE < 30L || min(cc) < 5L)
    ok <- fr$rank && fr$done && fr$fe && fr$fs
    cov[[i]] <- data.table(priority_id = pri$priority_id[i], pair_key = pri$pair_key[i], 
        n_strict_PD = n, nA = nA, nB = nB, nAB = nAB, nE = nE, nAB_E = nABE, as.list(cc))
    ready[[i]] <- data.table(priority_id = pri$priority_id[i], pair_key = pri$pair_key[i], 
        exact_specification_recovered = TRUE, required_covariates_available = TRUE, 
        design_matrix_full_rank = fr$rank, estimator_used = fr$est, fit_completed = fr$done, 
        finite_interaction_estimate = fr$fe, finite_interaction_se = fr$fs, estimable = ok, 
        nonestimable_reason = fr$why)
    logmsg("masked readiness complete:", pri$priority_id[i])
}
cov <- rbindlist(cov)
ready <- rbindlist(ready)
fwrite(cov, file.path(readiness, "canada_coverage_9.csv"))
fwrite(ready, file.path(readiness, "canada_model_readiness_9.csv"))
e <- sum(ready$estimable)
g <- if (e >= 6L) "GO" else if (e == 5L) "BORDERLINE_STOP" else "STOP"
if (g != "GO") stop("model blocked because readiness criterion was not GO")
model <- file.path(WORK_DIR, "canada_model_results")
dir.create(model, recursive = TRUE, showWarnings = FALSE)
faers <- fread(FAERS_REFERENCE)
if (nrow(faers) != 9L || !all(c("priority_id", "pair_key", "FAERS_interaction_OR") %chin% 
    names(faers)) || !setequal(faers$priority_id, pri$priority_id)) stop("Invalid model FAERS reference")
fit_model <- function(d, sparse) {
    cv <- c("n_nodes_group", "non_target_ae_terms_group", "calendar_period")
    cv <- cv[vapply(cv, function(x) length(unique(d[[x]][!is.na(d[[x]])])) > 1L, 
        logical(1))]
    f <- as.formula(paste("E ~", paste(c("A", "B", "AB", cv), collapse = " + ")))
    fit <- if (sparse) 
        glm(f, d, family = binomial("logit"), method = brglm2::brglmFit, control = brglm2::brglmControl(type = "AS_mean"))
    else glm(f, d, family = binomial("logit"))
    co <- summary(fit)$coefficients
    if (!("AB" %in% rownames(co))) 
        stop("AB coefficient unavailable in model")
    b <- as.numeric(co["AB", 1])
    s <- as.numeric(co["AB", 2])
    if (!is.finite(b) || !is.finite(s) || s <= 0) 
        stop("Nonfinite model interaction")
    list(estimator = if (sparse) "brglm2_penalized_logistic_wald_CI" else "standard_glm", 
        OR = exp(b), low = exp(b - 1.96 * s), high = exp(b + 1.96 * s))
}
out <- vector("list", 9)
for (i in seq_len(nrow(pri))) {
    p <- parsekey(pri$pair_key[i])
    d <- copy(base)
    d[, `:=`(A, as.integer(report_id %chin% unique(rd[node == p$A, report_id])))]
    d[, `:=`(B, as.integer(report_id %chin% unique(rd[node == p$B, report_id])))]
    d[, `:=`(AB, as.integer(A == 1L & B == 1L))]
    d[, `:=`(E, as.integer(report_id %chin% evs[[p$D]]))]
    d[, `:=`(non_target_ae_terms, pmax(0, total_ae_terms - E))]
    d[, `:=`(non_target_ae_terms_group, bgroup(non_target_ae_terms))]
    nA <- sum(d$A)
    nB <- sum(d$B)
    nAB <- sum(d$AB)
    nE <- sum(d$E)
    nABE <- sum(d$AB == 1L & d$E == 1L)
    cc <- c(sum(d$A == 0 & d$B == 0 & d$E == 0), sum(d$A == 0 & d$B == 0 & d$E == 
        1), sum(d$A == 1 & d$B == 0 & d$E == 0), sum(d$A == 1 & d$B == 0 & d$E == 
        1), sum(d$A == 0 & d$B == 1 & d$E == 0), sum(d$A == 0 & d$B == 1 & d$E == 
        1), sum(d$A == 1 & d$B == 1 & d$E == 0), sum(d$A == 1 & d$B == 1 & d$E == 
        1))
    z <- fit_model(d, nABE < 30L || min(cc) < 5L)
    fr <- faers[priority_id == pri$priority_id[i]]
    out[[i]] <- data.table(priority_id = pri$priority_id[i], pair_key = pri$pair_key[i], 
        drug_combination = pri$drug_combination[i], event = pri$event[i], nA = nA, 
        nB = nB, nAB = nAB, nE = nE, nAB_E = nABE, estimator = z$estimator, interaction_OR = z$OR, 
        interaction_low95 = z$low, interaction_high95 = z$high, same_direction_as_FAERS = ((z$OR > 
            1) == (fr$FAERS_interaction_OR > 1)), canada_lower_CI_gt_1 = (z$low > 
            1))
}
out <- rbindlist(out)
fwrite(out, file.path(model, "canada_fixed_9_results.csv"))
con <- merge(faers[, .(priority_id, FAERS_interaction_OR)], out[, .(priority_id, 
    Canada_interaction_OR = interaction_OR, same_direction_as_FAERS)], by = "priority_id")
con[, `:=`(log_FAERS_OR = log(FAERS_interaction_OR), log_Canada_OR = log(Canada_interaction_OR), 
    same_direction = same_direction_as_FAERS)]
con[, `:=`(same_direction_as_FAERS, NULL)]
fwrite(con, file.path(model, "faers_canada_concordance.csv"))
rho <- cor(con$log_FAERS_OR, con$log_Canada_OR, method = "spearman")
