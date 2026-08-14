options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
    library(data.table)
    library(brglm2)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop("Usage: Rscript 02_role_calendar_sensitivity.R SPECIFICATION_A_CSV BASE_CSV DRUG_CSV NODE_MAP_CSV OUTPUT_DIR")
spec_path <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
base_path <- normalizePath(args[2], winslash = "/", mustWork = TRUE)
drug_path <- normalizePath(args[3], winslash = "/", mustWork = TRUE)
map_path <- normalizePath(args[4], winslash = "/", mustWork = TRUE)
out <- normalizePath(args[5], winslash = "/", mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
norm_node <- function(x) trimws(tolower(gsub("[[:space:]]+", " ", as.character(x))))
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
fit_one <- function(dt, sparse) {
    covars <- c("n_nodes_group", "non_target_ae_terms_group", "calendar_period")
    covars <- covars[covars %in% names(dt)]
    covars <- covars[vapply(covars, function(x) uniqueN(dt[[x]][!is.na(dt[[x]])]) > 
        1, logical(1))]
    form <- as.formula(paste("E ~", paste(c("A", "B", "AB", covars), collapse = " + ")))
    fit <- tryCatch(if (sparse) 
        glm(form, data = dt, family = binomial("logit"), method = brglmFit, control = brglmControl(type = "AS_mean"))
    else glm(form, data = dt, family = binomial("logit")), error = function(e) e)
    if (inherits(fit, "error")) 
        return(list(status = "FAIL", estimator = if (sparse) "brglm2" else "glm", 
            warning = conditionMessage(fit), OR = NA_real_, low95 = NA_real_, high95 = NA_real_, 
            p = NA_real_))
    co <- coef(summary(fit))
    if (!"AB" %in% rownames(co)) 
        return(list(status = "FAIL", estimator = if (sparse) "brglm2" else "glm", 
            warning = "AB coefficient missing", OR = NA_real_, low95 = NA_real_, 
            high95 = NA_real_, p = NA_real_))
    est <- co["AB", "Estimate"]
    se <- co["AB", "Std. Error"]
    pcol <- grep("Pr", colnames(co), value = TRUE)[1]
    list(status = if (isTRUE(fit$converged) || sparse) "OK" else "CHECK_CONVERGENCE", 
        estimator = if (sparse) "brglm2_penalized_logistic_wald_CI" else "glm", warning = "", 
        OR = exp(est), low95 = exp(est - 1.96 * se), high95 = exp(est + 1.96 * se), 
        p = co["AB", pcol])
}
spec <- fread(spec_path)
stopifnot(nrow(spec) == 60)
base <- fread(base_path)
base[, `:=`(primaryid, as.character(primaryid))]
stopifnot(uniqueN(base$primaryid) == nrow(base))
drug <- fread(drug_path)
drug[, `:=`(primaryid = as.character(primaryid), drug_std = norm_node(drug_std), 
    role_cod = as.character(role_cod))]
map <- fread(map_path)
map[, `:=`(drug_std = norm_node(drug_std), node_norm = norm_node(node))]
map <- unique(map[drug_std != "" & node_norm != "", .(drug_std, node_norm)])
node_role <- merge(unique(drug[primaryid != "" & drug_std != "", .(primaryid, drug_std, 
    role_cod)]), map, by = "drug_std", all.x = TRUE, allow.cartesian = TRUE)
node_role[is.na(node_norm) | node_norm == "", `:=`(node_norm, drug_std)]
node_role <- unique(node_role[, .(primaryid, node_norm, role_cod)])
node_any <- split(node_role$primaryid, node_role$node_norm)
node_role_ok <- split(node_role[role_cod %chin% c("PS", "SS", "I"), primaryid], node_role[role_cod %chin% 
    c("PS", "SS", "I"), node_norm])
make_reports <- function(node, allowed = FALSE) {
    z <- if (allowed) 
        node_role_ok[[node]]
    else node_any[[node]]
    if (is.null(z)) 
        character()
    else z
}
make_dt <- function(row, variant, period = NULL) {
    ev <- row$domain_col
    Aany <- make_reports(row$node_a)
    Bany <- make_reports(row$node_b)
    Aok <- make_reports(row$node_a, TRUE)
    Bok <- make_reports(row$node_b, TRUE)
    dt <- base[if (is.null(period)) 
        TRUE
    else calendar_period == period, .(primaryid, total_ae_terms_rebuilt, n_nodes_group, 
        calendar_period, E = as.integer(get(ev)))]
    Aa <- dt$primaryid %chin% Aany
    Ba <- dt$primaryid %chin% Bany
    Ar <- dt$primaryid %chin% Aok
    Br <- dt$primaryid %chin% Bok
    if (variant == "at_least_one_PS_SS_I") {
        dt[, `:=`(A = as.integer(Aa), B = as.integer(Ba), AB = as.integer(Aa & Ba & 
            (Ar | Br)))]
    }
    else if (variant == "both_PS_SS_I") {
        dt[, `:=`(A = as.integer(Ar), B = as.integer(Br), AB = as.integer(Ar & Br))]
    }
    else {
        dt[, `:=`(A = as.integer(Aa), B = as.integer(Ba), AB = as.integer(Aa & Ba))]
    }
    dt[, `:=`(non_target_ae_terms_group, bin_group(pmax(0, total_ae_terms_rebuilt - 
        E)))]
    dt[, `:=`(n_nodes_group = factor(n_nodes_group), non_target_ae_terms_group = factor(non_target_ae_terms_group), 
        calendar_period = factor(calendar_period))]
    dt
}
run_row <- function(row, variant, period = NULL) {
    dt <- make_dt(row, variant, period)
    nAB <- sum(dt$AB == 1)
    nABE <- sum(dt$AB == 1 & dt$E == 1)
    nE <- sum(dt$E == 1)
    a <- nABE
    b <- nAB - nABE
    c <- nE - nABE
    d <- nrow(dt) - a - b - c
    sparse <- nABE < 30 || min(a, b, c, d) < 5
    res <- fit_one(dt, sparse)
    data.table(pair_key = row$pair_key, node_a = row$node_a, node_b = row$node_b, 
        event_domain = row$domain_col, variant = variant, calendar_period = if (is.null(period)) 
            "ALL"
        else period, N = nrow(dt), n_AB = nAB, n_AB_E = nABE, sparse_trigger = sparse, 
        status = res$status, estimator = res$estimator, OR = res$OR, low95 = res$low95, 
        high95 = res$high95, p = res$p, warning = res$warning)
}
role_out <- rbindlist(lapply(seq_len(nrow(spec)), function(i) rbind(run_row(spec[i], 
    "at_least_one_PS_SS_I"), run_row(spec[i], "both_PS_SS_I"))))
fwrite(role_out, file.path(out, "role_sensitivity.csv"))
periods <- c("2014_2017", "2018_2021", "2022_2026")
calendar_out <- rbindlist(lapply(seq_len(nrow(spec)), function(i) rbindlist(lapply(periods, 
    function(p) run_row(spec[i], "any_role", p)))))
fwrite(calendar_out, file.path(out, "calendar_period_stability.csv"))
