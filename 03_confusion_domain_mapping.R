suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript 03_confusion_domain_mapping.R EVENT_DICTIONARY_CSV OUTPUT_DIR")
dict_path <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
out_dir <- normalizePath(args[2], winslash = "/", mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

excluded <- c(
  "device use confusion", "product appearance confusion", "product confusion",
  "product dosage form confusion", "product dose confusion", "product label confusion",
  "product name confusion", "product packaging confusion", "product regimen confusion"
)

dict <- fread(dict_path)
before <- dict[event_domain == "Confusion_delirium"]
stopifnot(nrow(before) == 19L)
stopifnot(setequal(before$preferred_term, c(
  excluded, "altered state of consciousness", "confusion postoperative",
  "confusional arousal", "confusional state", "delirium", "delirium febrile",
  "dementia of the alzheimer's type, with delirium", "disorientation",
  "intensive care unit delirium", "postoperative delirium"
)))

before[, membership := "BEFORE_INCLUDED"]
after <- copy(before[!preferred_term %chin% excluded])
after[, `:=`(
  matching_rule = "exact normalized PT membership after tolower(trimws(pt))",
  regex_from_formal_script = NA_character_,
  formal_status = "UPDATED_EXACT_MEMBERSHIP",
  membership = "AFTER_INCLUDED",
  matched_report_count = NA_integer_,
  matched_reac_rows = NA_integer_
)]
stopifnot(nrow(after) == 10L)

diff <- before[, .(
  event_domain, preferred_term, before_included = TRUE,
  after_included = !preferred_term %chin% excluded,
  action = fifelse(preferred_term %chin% excluded,
    "EXCLUDE_NONCLINICAL_PRODUCT_OR_DEVICE_CONFUSION",
    "RETAIN_CLINICAL_CONFUSION_OR_DELIRIUM"),
  decision_basis = fifelse(preferred_term %chin% excluded,
    "product/device or medication-use confusion is outside the target clinical event",
    "existing clinical term retained; no new PT added")
)]
stopifnot(sum(!diff$after_included) == 9L, sum(diff$after_included) == 10L)
 setorder(before, preferred_term)
 setorder(after, preferred_term)
 setorder(diff, preferred_term)
fwrite(before, file.path(out_dir, "CONFUSION_PT_MAPPING_BEFORE.csv"))
fwrite(after, file.path(out_dir, "CONFUSION_PT_MAPPING_AFTER.csv"))
fwrite(diff, file.path(out_dir, "CONFUSION_PT_MAPPING_DIFF.csv"))
