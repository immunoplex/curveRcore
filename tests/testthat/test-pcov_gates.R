# Every branch of classify_pcov_gate() -- this class of logic is what the
# pcov_pass saga was about getting wrong quietly.

test_that("both LOQs defined: positional classification, pcov ignored", {
  r <- classify_pcov_gate(x_log10 = c(0.5, 1.0, 1.5, 2.0, 2.5),
                          pcov = rep(999, 5), pcov_threshold = 20,
                          lloq_log10 = 1.0, uloq_log10 = 2.0)
  expect_equal(as.character(r$pcov_gate_class),
               c("belowLLOQ", "inDynamicRange", "inDynamicRange", "inDynamicRange", "aboveULOQ"))
  expect_equal(r$pcov_pass, c(FALSE, TRUE, TRUE, TRUE, FALSE))
})

test_that("points exactly at LLOQ/ULOQ are in the dynamic range", {
  r <- classify_pcov_gate(x_log10 = c(1.0, 2.0), pcov = c(5, 5),
                          pcov_threshold = 20, lloq_log10 = 1.0, uloq_log10 = 2.0)
  expect_equal(as.character(r$pcov_gate_class), c("inDynamicRange", "inDynamicRange"))
})

test_that("neither LOQ defined: pcov threshold, side from the inflection point", {
  r <- classify_pcov_gate(x_log10 = c(0, 1, 2, 3), pcov = c(50, 10, 10, 50),
                          pcov_threshold = 20, lloq_log10 = NA, uloq_log10 = NA,
                          inflect_x = 1.5)
  expect_equal(as.character(r$pcov_gate_class),
               c("belowLLOQ", "inDynamicRange", "inDynamicRange", "aboveULOQ"))
  expect_equal(r$pcov_pass, c(FALSE, TRUE, TRUE, FALSE))
})

test_that("one LOQ missing uses the full fallback, not partial LOQ logic", {
  r <- classify_pcov_gate(x_log10 = c(0.5, 1.5), pcov = c(50, 5), pcov_threshold = 20,
                          lloq_log10 = 1.0, uloq_log10 = NA, inflect_x = 1.0)
  expect_equal(as.character(r$pcov_gate_class), c("belowLLOQ", "inDynamicRange"))
})

test_that("fallback with no inflection: failing class is NA but pcov_pass is FALSE", {
  r <- classify_pcov_gate(x_log10 = c(0.5, 1.5), pcov = c(50, 5), pcov_threshold = 20,
                          lloq_log10 = NA, uloq_log10 = NA, inflect_x = NA)
  expect_equal(is.na(r$pcov_gate_class), c(TRUE, FALSE))
  expect_equal(r$pcov_pass, c(FALSE, TRUE))
})

test_that("NA x_log10 propagates to both columns", {
  r <- classify_pcov_gate(x_log10 = c(NA, 1.5), pcov = c(5, 5), pcov_threshold = 20,
                          lloq_log10 = 1.0, uloq_log10 = 2.0)
  expect_equal(is.na(r$pcov_gate_class), c(TRUE, FALSE))
  expect_equal(is.na(r$pcov_pass), c(TRUE, FALSE))
})

test_that("NA pcov in the fallback branch propagates", {
  r <- classify_pcov_gate(x_log10 = c(1.5, 1.5), pcov = c(NA, 5), pcov_threshold = 20,
                          lloq_log10 = NA, uloq_log10 = NA, inflect_x = 1.0)
  expect_equal(is.na(r$pcov_gate_class), c(TRUE, FALSE))
})
