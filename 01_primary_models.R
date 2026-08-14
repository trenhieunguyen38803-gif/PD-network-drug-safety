options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) stop("Usage: Rscript 01_primary_models.R EVENT_CSV DRUG_CSV NODE_MAP_CSV COUNT_REFERENCE_CSV TIME_MAP_CSV OUTPUT_DIR")
input_paths <- normalizePath(args[1:5], winslash = "/", mustWork = TRUE)
output_root <- normalizePath(args[6], winslash = "/", mustWork = FALSE)
out <- file.path(output_root, "primary_models")
oa <- file.path(out, "01_TIME_ASSIGNMENT")
or <- file.path(out, "02_SPECIFICATION_RESULTS")
dir.create(oa, recursive = TRUE, showWarnings = FALSE)
dir.create(or, recursive = TRUE, showWarnings = FALSE)
wr <- function(x, p) fwrite(as.data.table(x), p, na = "")
inps <- data.table(item = c("event", "strict_drug", "node_map", "count_reference", 
    "time_map"), path = input_paths)
norm_node <- function(x) {
    x <- as.character(x)
    x <- trimws(tolower(x))
    x <- gsub("[[:space:]]+", " ", x)
    x
}
bin_group <- function(x, breaks = NULL, labels = NULL) {
    x <- suppressWarnings(as.numeric(x))
    if (is.null(breaks)) {
        qs <- unique(as.numeric(stats::quantile(x, probs = c(0, 0.5, 0.8, 0.95, 1), 
            na.rm = TRUE, type = 7)))
        if (length(qs) < 3) 
            return(factor(ifelse(is.na(x), "missing", "one_group")))
        breaks <- qs
    }
    out <- cut(x, breaks = breaks, include.lowest = TRUE, ordered_result = TRUE)
    out <- as.character(out)
    out[is.na(out)] <- "missing"
    factor(out)
}
ev <- fread(inps[item == "event", path], integer64 = "character")
sd <- fread(inps[item == "strict_drug", path], integer64 = "character")
mp <- fread(inps[item == "node_map", path], integer64 = "character")
cr <- fread(inps[item == "count_reference", path], integer64 = "character")
tm <- fread(inps[item == "time_map", path], integer64 = "character")
for (x in list(ev, sd, mp, cr, tm)) if ("primaryid" %in% names(x)) x[, `:=`(primaryid, 
    as.character(primaryid))]
tm[, `:=`(primaryid, as.character(primaryid))]
tm[, `:=`(source_quarter, as.character(source_quarter))]
tm[, `:=`(calendar_period, as.character(calendar_period))]
mult <- tm[, .(time_rows = .N, unique_periods = uniqueN(calendar_period), unique_source_quarters = uniqueN(source_quarter)), 
    by = primaryid]
wr(mult, file.path(oa, "time_map_primaryid_multiplicity.csv"))
dups <- mult[time_rows > 1]
wr(dups, file.path(oa, "duplicated_primaryids.csv"))
dd <- tm[primaryid %chin% dups$primaryid][order(primaryid, source_quarter)]
wr(dd, file.path(oa, "duplicated_rows_detail.csv"))
tm[, `:=`(year, as.integer(substr(source_quarter, 1, 4)))]
tm[, `:=`(quarter, as.integer(sub(".*Q", "", source_quarter)))]
setorder(tm, primaryid, year, quarter, source_quarter)
tmu <- tm[, .SD[1], by = primaryid]
tmu[, `:=`(c("year", "quarter"), NULL)]
if (anyDuplicated(tmu$primaryid)) stop("time assignment not unique")
wr(tmu, file.path(oa, "deterministic_time_assignment.csv"))
reports <- unique(as.character(ev$primaryid))
base0 <- data.table(primaryid = reports)
ev[, `:=`(primaryid, as.character(primaryid))]
base0 <- merge(base0, ev[, c("primaryid", "total_ae_terms_rebuilt", intersect(c("ICD_compulsive", 
    "Sleepiness", "Hallucination_psychosis", "Confusion_delirium"), names(ev))), 
    with = FALSE], by = "primaryid", all.x = TRUE)
base0[is.na(total_ae_terms_rebuilt), `:=`(total_ae_terms_rebuilt, 0)]
card <- data.table(merge_name = c("event_to_base", "node_count_to_base", "time_map_original_to_base", 
    "time_map_deduplicated_to_base"), left_rows = c(nrow(data.table(primaryid = reports)), 
    nrow(base0), nrow(base0), nrow(base0)), right_rows = c(nrow(ev), NA_integer_, 
    nrow(tm), nrow(tmu)), output_rows = c(nrow(base0), NA_integer_, nrow(merge(copy(base0), 
    unique(tm[, .(primaryid, calendar_period)]), by = "primaryid", all.x = TRUE)), 
    nrow(merge(copy(base0), tmu[, .(primaryid, calendar_period)], by = "primaryid", 
        all.x = TRUE))))
wr(card, file.path(oa, "merge_cardinality.csv"))
sd[, `:=`(primaryid, as.character(primaryid))]
mp[, `:=`(drug_std, norm_node(drug_std))]
mp[, `:=`(node_norm, norm_node(node))]
sd[, `:=`(drug_std, norm_node(drug_std))]
mp <- unique(mp[!is.na(drug_std) & drug_std != "" & !is.na(node_norm) & node_norm != 
    "", .(drug_std, node_norm)])
sn <- merge(unique(sd[!is.na(drug_std) & drug_std != "", .(primaryid, drug_std)]), 
    mp, by = "drug_std", all.x = TRUE, allow.cartesian = TRUE)
sn[is.na(node_norm) | node_norm == "", `:=`(node_norm, drug_std)]
nd <- unique(sn[, .(primaryid, node_norm)])
nn <- nd[, .(n_nodes = uniqueN(node_norm)), by = primaryid]
base <- merge(base0, nn, by = "primaryid", all.x = TRUE)
base[is.na(n_nodes), `:=`(n_nodes, 0)]
base[, `:=`(n_nodes_group, bin_group(n_nodes))]
base <- merge(base, tmu[, .(primaryid, calendar_period)], by = "primaryid", all.x = TRUE)
base[is.na(calendar_period) | calendar_period == "", `:=`(calendar_period, "unknown")]
base[, `:=`(calendar_period, factor(calendar_period))]
base[, `:=`(non_target_ae_terms, pmax(0, as.numeric(total_ae_terms_rebuilt)))]
base[, `:=`(non_target_ae_terms_group, bin_group(non_target_ae_terms))]
if (anyDuplicated(base$primaryid)) stop("FINAL ANALYTIC BASE HAS DUPLICATED PRIMARYID")
wr(base, file.path(oa, "final_unique_report_base.csv"))
beforeafter <- data.table(working = c("event_report_base", "after_original_time_merge", 
    "after_deterministic_time_merge"), rows = c(nrow(base0), card[merge_name == "time_map_original_to_base", 
    output_rows], nrow(base)), unique_primaryid = c(uniqueN(base0$primaryid), uniqueN(merge(copy(base0), 
    unique(tm[, .(primaryid, calendar_period)]), by = "primaryid", all.x = TRUE)$primaryid), 
    uniqueN(base$primaryid)))
wr(beforeafter, file.path(oa, "before_after_report_counts.csv"))
node_reports <- split(nd$primaryid, nd$node_norm)
dom <- intersect(c("ICD_compulsive", "Sleepiness", "Hallucination_psychosis", "Confusion_delirium"), 
    names(ev))
ca <- unique(cr[, .(node_a = norm_node(node_a), node_b = norm_node(node_b), domain_col = as.character(domain_col))])
ca <- ca[domain_col %in% dom]
ca[, `:=`(pair_key, paste(node_a, node_b, domain_col, sep = " | "))]
if (nrow(ca) != 60) stop("candidate family not 60")
`%chin%` <- data.table::`%chin%`
fit1 <- function(x, sp) {
    A <- node_reports[[x$node_a]]
    B <- node_reports[[x$node_b]]
    if (is.null(A)) 
        A <- character()
    if (is.null(B)) 
        B <- character()
    d <- copy(base)
    d[, `:=`(A = as.integer(primaryid %chin% A), B = as.integer(primaryid %chin% 
        B))]
    d[, `:=`(E, as.integer(get(x$domain_col)))]
    d[is.na(E), `:=`(E, 0L)]
    d[, `:=`(AB, as.integer(A == 1 & B == 1))]
    d[, `:=`(non_target_ae_terms, pmax(0, as.numeric(total_ae_terms_rebuilt) - E))]
    d[, `:=`(non_target_ae_terms_group, bin_group(non_target_ae_terms))]
    maplev <- function(v) {
        lv <- levels(v)
        setNames(ifelse(seq_along(lv) <= 2, "low", "high"), lv)[as.character(v)]
    }
    if (sp == "C") {
        d[, `:=`(n_nodes_group, factor(maplev(n_nodes_group), levels = c("low", "high")))]
        d[, `:=`(non_target_ae_terms_group, factor(maplev(non_target_ae_terms_group), 
            levels = c("low", "high")))]
    }
    f <- if (sp == "A") 
        E ~ A + B + AB + n_nodes_group + non_target_ae_terms_group + calendar_period
    else if (sp == "B") 
        E ~ A + B + AB + calendar_period
    else E ~ A + B + AB + n_nodes_group + non_target_ae_terms_group + calendar_period
    fs <- paste(deparse(f), collapse = " ")
    nabe <- sum(d$AB == 1 & d$E == 1)
    sparse <- nabe < 30
    fit <- tryCatch(if (sparse) 
        glm(f, data = d, family = binomial("logit"), method = brglm2::brglmFit, control = brglm2::brglmControl(type = "AS_mean"))
    else glm(f, data = d, family = binomial("logit")), error = function(e) e)
    if (inherits(fit, "error")) 
        return(data.table(pair_key = x$pair_key, node_a = x$node_a, node_b = x$node_b, 
            domain_col = x$domain_col, n_A = sum(d$A), n_B = sum(d$B), n_AB = sum(d$AB), 
            n_AB_E = nabe, estimator = if (sparse) "brglm2_AS_mean" else "glm", model_formula = fs, 
            coefficient = NA_real_, OR = NA_real_, low95 = NA_real_, high95 = NA_real_, 
            p = NA_real_, A_OR = NA_real_, B_OR = NA_real_, convergence = FALSE, 
            penalized_trigger = sparse, warning = conditionMessage(fit), model_status = "FAIL"))
    co <- summary(fit)$coefficients
    abn <- if ("AB" %in% rownames(co)) 
        "AB"
    else "A:B"
    data.table(pair_key = x$pair_key, node_a = x$node_a, node_b = x$node_b, domain_col = x$domain_col, 
        n_A = sum(d$A), n_B = sum(d$B), n_AB = sum(d$AB), n_AB_E = nabe, estimator = if (sparse) 
            "brglm2_AS_mean"
        else "glm", model_formula = fs, coefficient = co[abn, 1], OR = exp(co[abn, 
            1]), low95 = exp(co[abn, 1] - 1.96 * co[abn, 2]), high95 = exp(co[abn, 
            1] + 1.96 * co[abn, 2]), p = 2 * pnorm(-abs(co[abn, 1]/co[abn, 2])), 
        A_OR = exp(co["A", 1]), B_OR = exp(co["B", 1]), convergence = isTRUE(fit$converged), 
        penalized_trigger = sparse, warning = "", model_status = "OK")
}
run <- function(s) {
    z <- rbindlist(lapply(seq_len(nrow(ca)), function(i) fit1(ca[i], s)), fill = TRUE)
    if (sum(z$model_status == "OK") != 60) 
        stop(paste("fits failed", s))
    z[, `:=`(q, p.adjust(p, "BH"))]
    z[, `:=`(nominal_flag, p < 0.05)]
    z[, `:=`(positive_nominal_flag, OR > 1 & p < 0.05)]
    z[, `:=`(fdr_flag, q < 0.05)]
    z[, `:=`(positive_fdr_flag, OR > 1 & q < 0.05)]
    z[, `:=`(retained, OR > 1 & low95 > 1 & q < 0.05)]
    z[, `:=`(downgraded, !retained)]
    z[, `:=`(downgrade_reason, ifelse(retained, "none", ifelse(OR <= 1, "nonpositive_direction", 
        ifelse(low95 <= 1, "CI_not_above_1", ifelse(q >= 0.05, "BH_not_significant", 
            "other")))))]
    z
}
for (s in c("A", "B", "C")) {
    z <- run(s)
    z[, `:=`(specification, s)]
    wr(z, file.path(or, paste0("specification_", s, "_60_rows.csv")))
    assign(s, z)
}
cat("PRIMARY_MODELS_COMPLETE\n")
