# =============================================================================
# test-pcov_gates.R -- sanity checks for classify_pcov_gate()
#
# Not a full testthat suite (adjust to your package's actual test harness),
# but exercises every branch this function has, since this exact class of
# logic is what the whole pcov_pass saga was about getting wrong quietly.
# Run interactively after source()-ing pcov_gates.R, or wire into
# testthat::test_that() blocks as-is.
# =============================================================================

stopifnot_all_equal <- function(actual, expected, label) {
  ok <- isTRUE(all.equal(actual, expected)) || identical(actual, expected)
  if (!ok) stop(sprintf("FAIL [%s]: got %s, expected %s",
                        label, paste(actual, collapse=","), paste(expected, collapse=",")))
  cat(sprintf("PASS [%s]\n", label))
}

# ---- 1. Both LOQs defined: pure positional classification, pcov irrelevant -
r <- classify_pcov_gate(
  x_log10 = c(0.5, 1.0, 1.5, 2.0, 2.5),
  pcov    = c(999, 999, 999, 999, 999),   # deliberately garbage/noisy pcov
  pcov_threshold = 20,
  lloq_log10 = 1.0, uloq_log10 = 2.0
)
stopifnot_all_equal(as.character(r$pcov_gate_class),
                    c("belowLLOQ","inDynamicRange","inDynamicRange","inDynamicRange","aboveULOQ"),
                    "both-LOQ: positional, ignores garbage pcov")
stopifnot_all_equal(r$pcov_pass, c(FALSE, TRUE, TRUE, TRUE, FALSE),
                    "both-LOQ: pcov_pass matches gate_class")

# ---- 2. Boundary inclusivity: exactly at LLOQ/ULOQ counts as inDynamicRange
r <- classify_pcov_gate(x_log10 = c(1.0, 2.0), pcov = c(5, 5),
                        pcov_threshold = 20, lloq_log10 = 1.0, uloq_log10 = 2.0)
stopifnot_all_equal(as.character(r$pcov_gate_class), c("inDynamicRange","inDynamicRange"),
                    "both-LOQ: boundary points included")

# ---- 3. Neither LOQ defined: fallback to pcov vs threshold, classed by inflection
r <- classify_pcov_gate(
  x_log10 = c(0.0, 1.0, 2.0, 3.0),
  pcov    = c(50,  10,  10,  50),          # fails at the two ends
  pcov_threshold = 20,
  lloq_log10 = NA, uloq_log10 = NA,
  inflect_x  = 1.5
)
stopifnot_all_equal(as.character(r$pcov_gate_class),
                    c("belowLLOQ","inDynamicRange","inDynamicRange","aboveULOQ"),
                    "fallback: threshold pass/fail, side by inflection")
stopifnot_all_equal(r$pcov_pass, c(FALSE, TRUE, TRUE, FALSE),
                    "fallback: pcov_pass matches gate_class")

# ---- 4. One LOQ missing (e.g. ULOQ undefined): still falls back fully -----
r <- classify_pcov_gate(
  x_log10 = c(0.5, 1.5), pcov = c(50, 5), pcov_threshold = 20,
  lloq_log10 = 1.0, uloq_log10 = NA, inflect_x = 1.0
)
stopifnot_all_equal(as.character(r$pcov_gate_class), c("belowLLOQ","inDynamicRange"),
                    "one-LOQ-missing: fallback branch used, not partial LOQ logic")

# ---- 5. Fallback with NO inflection point: failing points stay NA class,
#         but pcov_pass is still definitively FALSE for them --------------
r <- classify_pcov_gate(
  x_log10 = c(0.5, 1.5), pcov = c(50, 5), pcov_threshold = 20,
  lloq_log10 = NA, uloq_log10 = NA, inflect_x = NA
)
stopifnot_all_equal(is.na(r$pcov_gate_class), c(TRUE, FALSE),
                    "fallback, no inflection: failing point's class is NA")
stopifnot_all_equal(r$pcov_pass, c(FALSE, TRUE),
                    "fallback, no inflection: pcov_pass still correct")

# ---- 6. NA x_log10 propagates to NA everywhere, regardless of branch -----
r <- classify_pcov_gate(
  x_log10 = c(NA, 1.5), pcov = c(5, 5), pcov_threshold = 20,
  lloq_log10 = 1.0, uloq_log10 = 2.0
)
stopifnot_all_equal(is.na(r$pcov_gate_class), c(TRUE, FALSE), "NA x_log10: propagates (both-LOQ)")
stopifnot_all_equal(is.na(r$pcov_pass), c(TRUE, FALSE), "NA x_log10: propagates to pcov_pass")

# ---- 7. NA pcov in fallback branch propagates too -------------------------
r <- classify_pcov_gate(
  x_log10 = c(1.5, 1.5), pcov = c(NA, 5), pcov_threshold = 20,
  lloq_log10 = NA, uloq_log10 = NA, inflect_x = 1.0
)
stopifnot_all_equal(is.na(r$pcov_gate_class), c(TRUE, FALSE), "NA pcov (fallback): propagates")

cat("\nAll classify_pcov_gate() sanity checks passed.\n")
