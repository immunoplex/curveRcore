test_that("weights_contract has the expected shape and keys", {
  ct <- weights_contract()
  expect_setequal(names(ct), c("table", "column", "type", "key", "nullable",
                               "added_in", "models", "description"))
  s <- weights_contract("sample"); f <- weights_contract("fit")
  # primary keys are unchanged by 0.6.0: the sample hierarchy equals calib_samples'
  expect_equal(s$column[s$key],
               c("curve_id", "method", "sampleid", "patientid", "timeperiod", "dilution"))
  expect_equal(f$column[f$key], c("multiplate_group_id", "method"))
  # 14 / 18 pre-0.6.0 columns (verbatim from the live DB) + the 0.6.0 additions
  expect_equal(sum(s$added_in == "0.5"), 14L); expect_equal(sum(s$added_in == "0.6.0"), 9L)
  expect_equal(sum(f$added_in == "0.5"), 18L); expect_equal(sum(f$added_in == "0.6.0"), 14L)
  expect_false(anyDuplicated(paste(ct$table, ct$column)) > 0)
  expect_true(all(nzchar(ct$description)))
  # every 0.6.0 column is nullable except `model`, so pre-0.6.0 rows stay valid
  new <- ct[ct$added_in == "0.6.0", ]
  expect_equal(new$column[!new$nullable], c("model", "model"))
})

test_that("weights_contract can be restricted to a model", {
  add  <- weights_contract("fit", model = "additive")
  mult <- weights_contract("fit", model = "multiplicative")
  expect_true(all(c("tau", "kappa", "reliability") %in% add$column))
  expect_false("phi" %in% add$column)
  expect_true("phi" %in% mult$column)
  expect_false("tau" %in% mult$column)
  expect_true(all(c("model", "beta1", "n_eff") %in% add$column))
  expect_error(weights_contract(model = "bogus"))
})

# ---- DDL <-> contract sync ---------------------------------------------------

ddl_path <- function() {
  p <- testthat::test_path("..", "..", "R", "calib_schema.sql")
  if (file.exists(p)) p else NA_character_
}

ddl_columns <- function(sql, table) {
  start <- grep(paste0("CREATE TABLE IF NOT EXISTS madi_results.", table, " \\("), sql)
  stopifnot(length(start) == 1L)
  end <- start + which(grepl("^\\);", sql[(start + 1):length(sql)]))[1]
  body <- sql[(start + 1):(end - 1)]
  body <- body[grepl("^    [a-z_0-9]+ ", body) & !grepl("^    CONSTRAINT", body)]
  sub("^    ([a-z_0-9]+) .*$", "\\1", body)
}

ddl_added <- function(sql, table) {
  start <- grep(paste0("^ALTER TABLE madi_results.", table, "$"), sql)
  stopifnot(length(start) == 1L)
  rest <- sql[(start + 1):length(sql)]
  end  <- which(grepl(";$", rest))[1]
  blk  <- rest[seq_len(end)]
  sub("^ *ADD COLUMN IF NOT EXISTS ([a-z_0-9]+) .*$", "\\1", blk)
}

test_that("calib_schema.sql matches weights_contract, column for column and in order", {
  p <- ddl_path()
  testthat::skip_if(is.na(p), "calib_schema.sql not available (installed package)")
  sql <- readLines(p, warn = FALSE)
  for (tb in c("sample", "fit")) {
    sqlname <- if (tb == "sample") "calib_weights" else "calib_weights_fit"
    ct <- weights_contract(tb)
    expect_equal(ddl_columns(sql, sqlname), ct$column[ct$added_in == "0.5"],
                 info = paste(sqlname, "CREATE TABLE"))
    expect_equal(ddl_added(sql, sqlname), ct$column[ct$added_in == "0.6.0"],
                 info = paste(sqlname, "ALTER TABLE additions"))
  }
  # keys in the DDL
  expect_true(any(grepl("calib_weights_pkey PRIMARY KEY \\(curve_id, method, sampleid, patientid, timeperiod, dilution\\)", sql)))
  expect_true(any(grepl("calib_weights_fit_pkey PRIMARY KEY \\(multiplate_group_id, method\\)", sql)))
})

# ---- validate_weights_table --------------------------------------------------

sample_rows <- function(model = "multiplicative") {
  x <- data.frame(curve_id = 1:2, method = "frequentist", sampleid = c("a", "b"),
                  patientid = "p", timeperiod = "t", dilution = "d",
                  se = c(0.1, 0.2), sigma = c(0.3, 0.4), w = c(10, 6), w_norm = c(1.2, 0.8),
                  stringsAsFactors = FALSE)
  if (model == "additive") {
    x$v_meas <- c(0.01, 0.04); x$tau2 <- 0.09; x$reliability <- c(0.9, 0.7)
    x$censor_type <- "none"
  }
  x
}

test_that("validate_weights_table accepts conforming rows for each model", {
  expect_silent(validate_weights_table(sample_rows(), "sample", "multiplicative"))
  expect_silent(validate_weights_table(sample_rows("additive"), "sample", "additive"))
  x <- sample_rows("additive"); x$model <- "additive"
  expect_silent(validate_weights_table(x, "sample", "additive"))
})

test_that("validate_weights_table rejects bad rows with a useful message", {
  expect_error(validate_weights_table(transform(sample_rows(), oops = 1), "sample",
                                      "multiplicative"), "unknown column")
  expect_error(validate_weights_table(sample_rows()[, -1], "sample", "multiplicative"),
               "missing required column.*curve_id")
  x <- sample_rows(); x$sampleid[1] <- NA
  expect_error(validate_weights_table(x, "sample", "multiplicative"), "sampleid.*missing")
  x <- sample_rows("additive"); x$model <- "multiplicative"
  expect_error(validate_weights_table(x, "sample", "additive"), "model column")
  # additive-only columns filled in multiplicative rows
  expect_error(validate_weights_table(sample_rows("additive"), "sample", "multiplicative"),
               "belongs to the additive model")
})

test_that("fit rows: phi is multiplicative-only, tau/kappa additive-only", {
  fit_mult <- data.frame(multiplate_group_id = "g", method = "bayesian",
                         phi = 1.1, beta1 = 0.9, n_eff = 100, stringsAsFactors = FALSE)
  expect_silent(validate_weights_table(fit_mult, "fit", "multiplicative"))
  expect_error(validate_weights_table(transform(fit_mult, tau = 0.4), "fit",
                                      "multiplicative"), "additive model")
  fit_add <- data.frame(multiplate_group_id = "g", method = "bayesian", tau = 0.4,
                        kappa = 1, reliability = 0.8, beta1 = 1, fix_beta1 = TRUE,
                        stringsAsFactors = FALSE)
  expect_silent(validate_weights_table(fit_add, "fit", "additive"))
  expect_error(validate_weights_table(transform(fit_add, phi = 1), "fit", "additive"),
               "multiplicative model")
})
