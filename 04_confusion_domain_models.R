options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({library(data.table); library(brglm2)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 13L) stop("Usage: Rscript 04_confusion_domain_models.R RAW_REACTION_CSV EVENT_FLAGS_CSV STRICT_DRUG_CSV NODE_MAP_CSV TIME_MAP_CSV SPEC_A_CSV SPEC_B_CSV SPEC_C_CSV ROLE_CSV CALENDAR_CSV MAPPING_AFTER_CSV EXPERT_CSV OUTPUT_DIR")
path_arg <- function(x) normalizePath(x, winslash = "/", mustWork = TRUE)
paths <- list(
  raw_reac = path_arg(args[1]), old_event = path_arg(args[2]), strict_drug = path_arg(args[3]),
  node_map = path_arg(args[4]), time_map = path_arg(args[5]), old_A = path_arg(args[6]),
  old_B = path_arg(args[7]), old_C = path_arg(args[8]), old_role = path_arg(args[9]),
  old_calendar = path_arg(args[10]), mapping_after = path_arg(args[11]), expert = path_arg(args[12])
)
out <- normalizePath(args[13], winslash = "/", mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(all(file.exists(unlist(paths))))

norm_node <- function(x) {
  x <- as.character(x)
  x <- trimws(tolower(x))
  x <- gsub("[[:space:]]+", " ", x)
  x
}

bin_group <- function(x, breaks = NULL, labels = NULL) {
  x <- suppressWarnings(as.numeric(x))
  if (is.null(breaks)) {
    qs <- unique(as.numeric(stats::quantile(x, probs = c(0, 0.5, 0.8, 0.95, 1), na.rm = TRUE, type = 7)))
    if (length(qs) < 3) return(factor(ifelse(is.na(x), "missing", "one_group")))
    breaks <- qs
  }
  out <- cut(x, breaks = breaks, include.lowest = TRUE, ordered_result = TRUE)
  out <- as.character(out)
  out[is.na(out)] <- "missing"
  factor(out)
}

old_A <- fread(paths$old_A)
old_B <- fread(paths$old_B)
old_C <- fread(paths$old_C)
stopifnot(
  nrow(old_A) == 60L,
  nrow(old_B) == 60L,
  nrow(old_C) == 60L,
  identical(old_A$pair_key, old_B$pair_key),
  identical(old_A$pair_key, old_C$pair_key),
  uniqueN(old_A$pair_key) == 60L
)
cands <- old_A[, .(pair_key, node_a, node_b, domain_col)]
conf_keys <- cands[domain_col == "Confusion_delirium"]
stopifnot(nrow(conf_keys) == 15L)

old_event <- fread(paths$old_event, integer64 = "character")
old_event[, primaryid := as.character(primaryid)]
stopifnot(
  nrow(old_event) == 126045L,
  uniqueN(old_event$primaryid) == 126045L
)
reac <- fread(
  paths$raw_reac,
  colClasses = c(
    primaryid = "character",
    pt = "character",
    source_file = "character"
  ),
  showProgress = FALSE
)
reac[, `:=`(
  primaryid = as.character(primaryid),
  pt = tolower(trimws(as.character(pt)))
)]
reac <- unique(
  reac[
    primaryid != "" & !is.na(primaryid) &
      pt != "" & !is.na(pt)
  ]
)
stopifnot(nrow(reac) == 499031L, uniqueN(reac$primaryid) == 126045L)
after_map <- fread(paths$mapping_after)
retained_pt <- after_map$preferred_term
stopifnot(length(retained_pt) == 10L, uniqueN(retained_pt) == 10L)
new_conf_ids <- unique(reac[pt %chin% retained_pt, primaryid])

event_new <- copy(old_event)
event_new[, Confusion_delirium := primaryid %chin% new_conf_ids]
old_true <- old_event[Confusion_delirium %in% TRUE, primaryid]
new_true <- event_new[Confusion_delirium %in% TRUE, primaryid]
stopifnot(
  all(new_true %chin% old_true),
  length(setdiff(new_true, old_true)) == 0L
)
removed_report_n <- length(setdiff(old_true, new_true))

strict <- fread(paths$strict_drug, integer64 = "character")
map <- fread(paths$node_map, integer64 = "character")
tm <- fread(paths$time_map, integer64 = "character")
strict[, `:=`(
  primaryid = as.character(primaryid),
  drug_std = norm_node(drug_std),
  role_cod = as.character(role_cod)
)]
map[, `:=`(
  drug_std = norm_node(drug_std),
  node_norm = norm_node(node)
)]
map <- unique(
  map[
    !is.na(drug_std) & drug_std != "" &
      !is.na(node_norm) & node_norm != "",
    .(drug_std, node_norm)
  ]
)
node_role <- merge(
  unique(
    strict[
      !is.na(drug_std) & drug_std != "",
      .(primaryid, drug_std, role_cod)
    ]
  ),
  map,
  by = "drug_std",
  all.x = TRUE,
  allow.cartesian = TRUE
)
node_role[is.na(node_norm) | node_norm == "", node_norm := drug_std]
node_role <- unique(node_role[, .(primaryid, node_norm, role_cod)])
nd <- unique(node_role[, .(primaryid, node_norm)])
node_reports <- split(nd$primaryid, nd$node_norm)

tm[, `:=`(
  primaryid = as.character(primaryid),
  source_quarter = as.character(source_quarter),
  calendar_period = as.character(calendar_period)
)]
tm[, year := as.integer(substr(source_quarter, 1L, 4L))]
tm[, quarter := as.integer(sub(".*Q", "", source_quarter))]
setorder(tm, primaryid, year, quarter, source_quarter)
tmu <- tm[, .SD[1L], by = primaryid]
tmu[, c("year", "quarter") := NULL]
stopifnot(uniqueN(tmu$primaryid) == nrow(tmu))

build_base <- function(event_dt) {
  domains <- c(
    "ICD_compulsive", "Sleepiness",
    "Hallucination_psychosis", "Confusion_delirium"
  )
  b <- data.table(primaryid = unique(as.character(event_dt$primaryid)))
  b <- merge(
    b,
    event_dt[
      ,
      c("primaryid", "total_ae_terms_rebuilt", domains),
      with = FALSE
    ],
    by = "primaryid",
    all.x = TRUE
  )
  b[is.na(total_ae_terms_rebuilt), total_ae_terms_rebuilt := 0]
  nn <- nd[, .(n_nodes = uniqueN(node_norm)), by = primaryid]
  b <- merge(b, nn, by = "primaryid", all.x = TRUE)
  b[is.na(n_nodes), n_nodes := 0]
  b[, n_nodes_group := bin_group(n_nodes)]
  b <- merge(
    b,
    tmu[, .(primaryid, calendar_period)],
    by = "primaryid",
    all.x = TRUE
  )
  b[
    is.na(calendar_period) | calendar_period == "",
    calendar_period := "unknown"
  ]
  b[, calendar_period := factor(calendar_period)]
  b[, non_target_ae_terms := pmax(
    0, as.numeric(total_ae_terms_rebuilt)
  )]
  b[, non_target_ae_terms_group := bin_group(
    non_target_ae_terms
  )]
  stopifnot(nrow(b) == 126045L, uniqueN(b$primaryid) == 126045L)
  b
}
base_old <- build_base(old_event)
base_new <- build_base(event_new)

fit_primary <- function(ca, specification, base_dt) {
  A_ids <- node_reports[[ca$node_a]]
  B_ids <- node_reports[[ca$node_b]]
  if (is.null(A_ids)) A_ids <- character()
  if (is.null(B_ids)) B_ids <- character()
  d <- copy(base_dt)
  d[, `:=`(
    A = as.integer(primaryid %chin% A_ids),
    B = as.integer(primaryid %chin% B_ids)
  )]
  d[, E := as.integer(get(ca$domain_col))]
  d[is.na(E), E := 0L]
  d[, AB := as.integer(A == 1L & B == 1L)]
  d[, non_target_ae_terms := pmax(
    0, as.numeric(total_ae_terms_rebuilt) - E
  )]
  d[, non_target_ae_terms_group := bin_group(
    non_target_ae_terms
  )]
  if (specification == "C") {
    map_level <- function(v) {
      lv <- levels(v)
      setNames(
        ifelse(seq_along(lv) <= 2L, "low", "high"),
        lv
      )[as.character(v)]
    }
    d[, n_nodes_group := factor(
      map_level(n_nodes_group),
      levels = c("low", "high")
    )]
    d[, non_target_ae_terms_group := factor(
      map_level(non_target_ae_terms_group),
      levels = c("low", "high")
    )]
  }
  form <- if (specification == "A") {
    E ~ A + B + AB + n_nodes_group +
      non_target_ae_terms_group + calendar_period
  } else if (specification == "B") {
    E ~ A + B + AB + calendar_period
  } else {
    E ~ A + B + AB + n_nodes_group +
      non_target_ae_terms_group + calendar_period
  }
  formula_string <- paste(deparse(form), collapse = " ")
  n_AB_E <- sum(d$AB == 1L & d$E == 1L)
  sparse <- n_AB_E < 30L
  fit <- tryCatch(
    if (sparse) {
      glm(
        form,
        data = d,
        family = binomial("logit"),
        method = brglm2::brglmFit,
        control = brglm2::brglmControl(type = "AS_mean")
      )
    } else {
      glm(form, data = d, family = binomial("logit"))
    },
    error = function(e) e
  )
  base_result <- data.table(
    pair_key = ca$pair_key,
    node_a = ca$node_a,
    node_b = ca$node_b,
    domain_col = ca$domain_col,
    n_A = sum(d$A),
    n_B = sum(d$B),
    n_AB = sum(d$AB),
    n_AB_E = n_AB_E,
    estimator = if (sparse) "brglm2_AS_mean" else "glm",
    model_formula = formula_string,
    penalized_trigger = sparse
  )
  if (inherits(fit, "error")) {
    return(cbind(
      base_result,
      data.table(
        coefficient = NA_real_,
        OR = NA_real_,
        low95 = NA_real_,
        high95 = NA_real_,
        p = NA_real_,
        A_OR = NA_real_,
        B_OR = NA_real_,
        convergence = FALSE,
        warning = conditionMessage(fit),
        model_status = "FAIL"
      )
    ))
  }
  co <- summary(fit)$coefficients
  ab_name <- if ("AB" %in% rownames(co)) "AB" else "A:B"
  p_value <- 2 * pnorm(-abs(co[ab_name, 1L] / co[ab_name, 2L]))
  cbind(
    base_result,
    data.table(
      coefficient = co[ab_name, 1L],
      OR = exp(co[ab_name, 1L]),
      low95 = exp(co[ab_name, 1L] - 1.96 * co[ab_name, 2L]),
      high95 = exp(co[ab_name, 1L] + 1.96 * co[ab_name, 2L]),
      p = p_value,
      A_OR = exp(co["A", 1L]),
      B_OR = exp(co["B", 1L]),
      convergence = isTRUE(fit$converged),
      warning = "",
      model_status = "OK"
    )
  )
}

finalize_spec <- function(new_conf, old_full, specification) {
  stopifnot(
    nrow(new_conf) == 15L,
    all(new_conf$model_status == "OK")
  )
  z <- rbindlist(
    list(
      copy(old_full[domain_col != "Confusion_delirium"]),
      copy(new_conf)
    ),
    fill = TRUE,
    use.names = TRUE
  )
  z <- z[match(old_full$pair_key, pair_key)]
  stopifnot(
    nrow(z) == 60L,
    identical(z$pair_key, old_full$pair_key)
  )
  z[, q := p.adjust(p, method = "BH")]
  z[, nominal_flag := p < 0.05]
  z[, positive_nominal_flag := OR > 1 & p < 0.05]
  z[, fdr_flag := q < 0.05]
  z[, positive_fdr_flag := OR > 1 & q < 0.05]
  z[, retained := OR > 1 & low95 > 1 & q < 0.05]
  z[, downgraded := !retained]
  z[, downgrade_reason := fifelse(
    retained,
    "none",
    fifelse(
      OR <= 1,
      "nonpositive_direction",
      fifelse(
        low95 <= 1,
        "CI_not_above_1",
        fifelse(q >= 0.05, "BH_not_significant", "other")
      )
    )
  )]
  z[, specification := specification]
  z[, update_source := fifelse(
    domain_col == "Confusion_delirium",
    "REFITTED_updateED_EXACT_PT_FLAG",
    "REUSED_fixed_NON_CONFUSION_ESTIMATE"
  )]
  z
}

run_conf_spec <- function(specification) {
  rbindlist(
    lapply(
      seq_len(nrow(conf_keys)),
      function(i) {
        fit_primary(conf_keys[i], specification, base_new)
      }
    ),
    fill = TRUE
  )
}
new_A_conf <- run_conf_spec("A")
new_B_conf <- run_conf_spec("B")
new_C_conf <- run_conf_spec("C")
A <- finalize_spec(new_A_conf, old_A, "A")
B <- finalize_spec(new_B_conf, old_B, "B")
C <- finalize_spec(new_C_conf, old_C, "C")

fwrite(
  A,
  file.path(out, "SPEC_A_60_ROWS.csv")
)
fwrite(
  B,
  file.path(out, "SPEC_B_60_ROWS.csv")
)
fwrite(
  C,
  file.path(out, "SPEC_C_60_ROWS.csv")
)

compare_A <- merge(
  old_A,
  A,
  by = "pair_key",
  suffixes = c("_old", "_new"),
  sort = FALSE
)
compare_A[, `:=`(
  event_domain = domain_col_new,
  refit_required = domain_col_new == "Confusion_delirium",
  n_AB_changed = n_AB_old != n_AB_new,
  n_AB_E_changed = n_AB_E_old != n_AB_E_new,
  OR_absolute_change = abs(OR_new - OR_old),
  p_absolute_change = abs(p_new - p_old),
  q_absolute_change = abs(q_new - q_old),
  retained_changed = retained_new != retained_old
)]
compare_A <- compare_A[match(old_A$pair_key, pair_key)]
fwrite(
  compare_A,
  file.path(out, "OLD_VS_NEW_PRIMARY_RESULTS.csv")
)

concordance <- Reduce(
  function(x, y) merge(x, y, by = "pair_key", all = TRUE, sort = FALSE),
  list(
    A[, .(
      pair_key, node_a, node_b, event_domain = domain_col,
      spec_A_OR = OR, spec_A_low95 = low95,
      spec_A_high95 = high95, spec_A_p = p, spec_A_q = q,
      spec_A_retained = retained
    )],
    B[, .(
      pair_key,
      spec_B_OR = OR, spec_B_low95 = low95,
      spec_B_high95 = high95, spec_B_p = p, spec_B_q = q,
      spec_B_retained = retained
    )],
    C[, .(
      pair_key,
      spec_C_OR = OR, spec_C_low95 = low95,
      spec_C_high95 = high95, spec_C_p = p, spec_C_q = q,
      spec_C_retained = retained
    )]
  )
)
concordance <- concordance[
  match(old_A$pair_key, pair_key)
]
concordance[, retained_specification_count :=
  as.integer(spec_A_retained) +
    as.integer(spec_B_retained) +
    as.integer(spec_C_retained)]
concordance[, retained_in_at_least_2_of_3 :=
  retained_specification_count >= 2L]
concordance[, retained_in_all_3 :=
  retained_specification_count == 3L]
fwrite(
  concordance,
  file.path(out, "FINAL_SPECIFICATION_CONCORDANCE.csv")
)

expert <- fread(paths$expert)
node_label_map <- c(
  "加巴喷丁类" = "gabapentinoid",
  "含左旋多巴治疗" = "levodopa_containing_therapy",
  "普拉克索" = "pramipexole",
  "补充剂开放类别" = "supplement",
  "恩他卡朋" = "entacapone",
  "苯二氮䓬类" = "benzodiazepine",
  "罗替戈汀" = "rotigotine",
  "金刚烷胺" = "amantadine",
  "非阿片镇痛药类" = "analgesic_nonopioid",
  "阿扑吗啡" = "apomorphine"
)
event_label_map <- c(
  "嗜睡/日间困倦" = "Sleepiness",
  "幻觉/精神病样症状" = "Hallucination_psychosis",
  "意识混乱/谵妄" = "Confusion_delirium",
  "冲动控制/强迫行为" = "ICD_compulsive"
)
expert[, c("drug_1_cn", "drug_2_cn") :=
  tstrsplit(`药物组合`, " \\+ ", fixed = FALSE)]
expert[, pair_key := paste(
  unname(node_label_map[drug_1_cn]),
  unname(node_label_map[drug_2_cn]),
  unname(event_label_map[`事件`]),
  sep = " | "
)]
stopifnot(
  !anyNA(expert$pair_key),
  uniqueN(expert$pair_key) == 22L,
  setequal(expert$pair_key, old_A[retained %in% TRUE, pair_key])
)
expert_key <- expert[, .(
  pair_key,
  fixed_expert_canonical_id = canonical_id,
  fixed_expert_class = `最终分类`,
  fixed_expert_source = basename(paths$expert)
)]

primary <- merge(
  A[retained %in% TRUE],
  concordance[, .(
    pair_key,
    retained_specification_count,
    retained_in_at_least_2_of_3,
    retained_in_all_3
  )],
  by = "pair_key",
  all.x = TRUE,
  sort = FALSE
)
primary <- merge(
  primary,
  expert_key,
  by = "pair_key",
  all.x = TRUE,
  sort = FALSE
)
primary[, previous_primary := pair_key %chin%
  old_A[retained %in% TRUE, pair_key]]
primary[, review_status := fifelse(
  previous_primary & !is.na(fixed_expert_class),
  "fixed_EXPERT_REVIEW_REUSED",
  "NEW_PRIMARY_NO_fixed_EXPERT_REVIEW"
)]
primary[, human_scientific_review_required :=
  review_status == "NEW_PRIMARY_NO_fixed_EXPERT_REVIEW"]
primary[, statistical_clinical_concordance :=
  retained_in_at_least_2_of_3 &
    fixed_expert_class == "B" &
    !human_scientific_review_required]
primary <- primary[match(A[retained %in% TRUE, pair_key], pair_key)]
fwrite(
  primary,
  file.path(out, "FINAL_PRIMARY_SET.csv")
)
clinical <- primary[statistical_clinical_concordance %in% TRUE]
fwrite(
  clinical,
  file.path(out, "FINAL_CLINICAL_PRIORITY_SET.csv")
)

decomp <- primary[, .(
  pair_key,
  node_a,
  node_b,
  event_domain = domain_col,
  n_A,
  n_B,
  n_AB,
  n_AB_E,
  conditional_A_OR = A_OR,
  conditional_B_OR = B_OR,
  interaction_OR = OR,
  interaction_low95 = low95,
  interaction_high95 = high95,
  interaction_p = p,
  interaction_q = q,
  estimator,
  penalized_trigger,
  retained_specification_count,
  fixed_expert_class,
  statistical_clinical_concordance
)]
decomp[, joint_vs_double_negative_OR_point :=
  conditional_A_OR * conditional_B_OR * interaction_OR]
decomp[, joint_point_contrast_below_1 :=
  joint_vs_double_negative_OR_point < 1]
decomp[, joint_point_contrast_reference :=
  "A=0,B=0 reporting-odds reference under the fitted model"]
decomp[, interpretation_boundary := paste(
  "Deterministic point contrast only; no CI, P, or q.",
  "The interaction OR is a multiplicative reporting-scale departure",
  "and is not a clinical-risk contrast."
)]
fwrite(
  decomp,
  file.path(
    out,
    "CONDITIONAL_DECOMPOSITION_WITH_JOINT_POINT_CONTRAST.csv"
  )
)

node_any <- split(node_role$primaryid, node_role$node_norm)
node_role_ok <- split(
  node_role[role_cod %chin% c("PS", "SS", "I"), primaryid],
  node_role[role_cod %chin% c("PS", "SS", "I"), node_norm]
)
get_reports <- function(node, role_restricted = FALSE) {
  z <- if (role_restricted) {
    node_role_ok[[node]]
  } else {
    node_any[[node]]
  }
  if (is.null(z)) character() else z
}

fit_sensitivity <- function(ca, variant, period = NULL) {
  A_any <- get_reports(ca$node_a)
  B_any <- get_reports(ca$node_b)
  A_ok <- get_reports(ca$node_a, TRUE)
  B_ok <- get_reports(ca$node_b, TRUE)
  d <- if (is.null(period)) {
    copy(base_new)
  } else {
    copy(base_new[as.character(calendar_period) == period])
  }
  d[, E := as.integer(get(ca$domain_col))]
  any_A <- d$primaryid %chin% A_any
  any_B <- d$primaryid %chin% B_any
  role_A <- d$primaryid %chin% A_ok
  role_B <- d$primaryid %chin% B_ok
  if (variant == "at_least_one_PS_SS_I") {
    d[, `:=`(
      A = as.integer(any_A),
      B = as.integer(any_B),
      AB = as.integer(any_A & any_B & (role_A | role_B))
    )]
  } else if (variant == "both_PS_SS_I") {
    d[, `:=`(
      A = as.integer(role_A),
      B = as.integer(role_B),
      AB = as.integer(role_A & role_B)
    )]
  } else {
    d[, `:=`(
      A = as.integer(any_A),
      B = as.integer(any_B),
      AB = as.integer(any_A & any_B)
    )]
  }
  d[, non_target_ae_terms_group := bin_group(
    pmax(0, as.numeric(total_ae_terms_rebuilt) - E)
  )]
  d[, `:=`(
    n_nodes_group = factor(n_nodes_group),
    non_target_ae_terms_group = factor(non_target_ae_terms_group),
    calendar_period = factor(calendar_period)
  )]
  covariates <- c(
    "n_nodes_group",
    "non_target_ae_terms_group",
    "calendar_period"
  )
  covariates <- covariates[
    vapply(
      covariates,
      function(v) uniqueN(d[[v]][!is.na(d[[v]])]) > 1L,
      logical(1L)
    )
  ]
  form <- as.formula(paste(
    "E ~",
    paste(c("A", "B", "AB", covariates), collapse = " + ")
  ))
  n_AB <- sum(d$AB == 1L)
  n_AB_E <- sum(d$AB == 1L & d$E == 1L)
  n_E <- sum(d$E == 1L)
  cells <- c(
    n_AB_E,
    n_AB - n_AB_E,
    n_E - n_AB_E,
    nrow(d) - n_AB - n_E + n_AB_E
  )
  sparse <- n_AB_E < 30L || min(cells) < 5L
  fit <- tryCatch(
    if (sparse) {
      glm(
        form,
        data = d,
        family = binomial("logit"),
        method = brglm2::brglmFit,
        control = brglm2::brglmControl(type = "AS_mean")
      )
    } else {
      glm(form, data = d, family = binomial("logit"))
    },
    error = function(e) e
  )
  common <- data.table(
    pair_key = ca$pair_key,
    node_a = ca$node_a,
    node_b = ca$node_b,
    event_domain = ca$domain_col,
    variant = variant,
    calendar_period = if (is.null(period)) "ALL" else period,
    N = nrow(d),
    n_AB = n_AB,
    n_AB_E = n_AB_E,
    sparse_trigger = sparse,
    estimator = if (sparse) {
      "brglm2_penalized_logistic_wald_CI"
    } else {
      "glm"
    }
  )
  if (inherits(fit, "error")) {
    return(cbind(
      common,
      data.table(
        status = "FAIL",
        OR = NA_real_,
        low95 = NA_real_,
        high95 = NA_real_,
        p = NA_real_,
        warning = conditionMessage(fit)
      )
    ))
  }
  co <- summary(fit)$coefficients
  if (!"AB" %in% rownames(co)) {
    return(cbind(
      common,
      data.table(
        status = "FAIL",
        OR = NA_real_,
        low95 = NA_real_,
        high95 = NA_real_,
        p = NA_real_,
        warning = "AB coefficient missing"
      )
    ))
  }
  estimate <- co["AB", "Estimate"]
  se <- co["AB", "Std. Error"]
  p_col <- grep("Pr", colnames(co), value = TRUE)[1L]
  cbind(
    common,
    data.table(
      status = if (isTRUE(fit$converged) || sparse) {
        "OK"
      } else {
        "CHECK_CONVERGENCE"
      },
      OR = exp(estimate),
      low95 = exp(estimate - 1.96 * se),
      high95 = exp(estimate + 1.96 * se),
      p = co["AB", p_col],
      warning = ""
    )
  )
}

role_variants <- c("at_least_one_PS_SS_I", "both_PS_SS_I")
role_conf <- rbindlist(
  lapply(
    seq_len(nrow(conf_keys)),
    function(i) {
      rbindlist(lapply(
        role_variants,
        function(v) fit_sensitivity(conf_keys[i], v)
      ))
    }
  ),
  fill = TRUE
)
old_role <- fread(paths$old_role)
role_new <- rbindlist(
  list(
    old_role[event_domain != "Confusion_delirium"],
    role_conf
  ),
  fill = TRUE,
  use.names = TRUE
)
role_order <- old_role[, .(pair_key, variant)]
setkeyv(role_new, c("pair_key", "variant"))
role_new <- role_new[role_order, on = .(pair_key, variant)]

periods <- c("2014_2017", "2018_2021", "2022_2026")
calendar_conf <- rbindlist(
  lapply(
    seq_len(nrow(conf_keys)),
    function(i) {
      rbindlist(lapply(
        periods,
        function(p) fit_sensitivity(conf_keys[i], "any_role", p)
      ))
    }
  ),
  fill = TRUE
)
old_calendar <- fread(paths$old_calendar)
calendar_new <- rbindlist(
  list(
    old_calendar[event_domain != "Confusion_delirium"],
    calendar_conf
  ),
  fill = TRUE,
  use.names = TRUE
)
calendar_order <- old_calendar[, .(pair_key, calendar_period)]
setkeyv(calendar_new, c("pair_key", "calendar_period"))
calendar_new <- calendar_new[
  calendar_order,
  on = .(pair_key, calendar_period)
]

primary_reference <- A[, .(
  pair_key,
  primary_n_AB = n_AB,
  primary_n_AB_E = n_AB_E,
  primary_OR = OR,
  primary_low95 = low95,
  primary_high95 = high95,
  primary_p = p,
  primary_q = q,
  primary_retained = retained
)]
clinical_keys <- clinical$pair_key
annotate_sensitivity <- function(x, analysis_type) {
  drop_cols <- intersect(
    c(
      "primary_n_AB", "primary_n_AB_E", "primary_OR",
      "primary_low95", "primary_high95", "primary_p",
      "primary_q", "primary_retained", "analysis_type",
      "estimable", "direction_agrees", "positive_ci",
      "primary_positive_ci", "n_AB_change", "n_AB_E_change",
      "final_clinical_priority"
    ),
    names(x)
  )
  if (length(drop_cols)) x[, (drop_cols) := NULL]
  y <- merge(x, primary_reference, by = "pair_key", all.x = TRUE)
  stopifnot(nrow(y) == nrow(x), !anyNA(y$primary_OR))
  y[, `:=`(
    analysis_type = analysis_type,
    estimable = status == "OK" & is.finite(OR),
    direction_agrees = fifelse(
      status == "OK" & is.finite(OR),
      (OR > 1) == (primary_OR > 1),
      NA
    ),
    positive_ci = fifelse(
      status == "OK" & is.finite(low95),
      low95 > 1,
      NA
    ),
    primary_positive_ci = primary_low95 > 1,
    n_AB_change = n_AB - primary_n_AB,
    n_AB_E_change = n_AB_E - primary_n_AB_E,
    final_clinical_priority = pair_key %chin% clinical_keys
  )]
  y
}
role_new <- annotate_sensitivity(role_new, "role")
calendar_new <- annotate_sensitivity(calendar_new, "calendar_period")
stopifnot(
  nrow(role_new) == 120L,
  nrow(calendar_new) == 180L,
  uniqueN(role_new$pair_key) == 60L,
  uniqueN(calendar_new$pair_key) == 60L
)
fwrite(
  role_new,
  file.path(out, "FINAL_ROLE_SENSITIVITY_AFTER_update.csv")
)
fwrite(
  calendar_new,
  file.path(out, "FINAL_CALENDAR_SENSITIVITY_AFTER_update.csv")
)
