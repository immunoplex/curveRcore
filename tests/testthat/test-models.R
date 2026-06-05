# test-models.R — Tests for forward models, inverses, derivatives, gradients,
#                  formulas, grid, transforms, settings, and result_class

# =============================================================================
# Shared test parameters
# =============================================================================
test_params <- list(a = 1.0, b = 0.5, c = 2.0, d = 4.5, g = 1.5)
x_log <- seq(-1, 4, length.out = 20)   # log10(conc) scale
x_pos <- 10^x_log                       # raw conc scale (positive)


# =============================================================================
# Forward models
# =============================================================================
test_that("all five forward models return correct length and are monotonic", {
  p <- test_params

  models <- list(
    logistic4    = function(x) logistic4(x, p$a, p$b, p$c, p$d),
    logistic5    = function(x) logistic5(x, p$a, p$b, p$c, p$d, p$g),
    loglogistic4 = function(x) loglogistic4(x_pos, p$a, p$b, p$c, p$d),
    loglogistic5 = function(x) loglogistic5(x, p$a, p$b, p$c, p$d, p$g),
    gompertz4    = function(x) gompertz4(x, p$a, p$b, p$c, p$d)
  )

  for (nm in names(models)) {
    y <- if (nm == "loglogistic4") models[[nm]](x_pos) else models[[nm]](x_log)
    expect_length(y, 20)
    expect_true(all(is.finite(y)), info = paste(nm, "has non-finite values"))
    expect_true(all(diff(y) > 0), info = paste(nm, "is not monotonically increasing"))
  }
})

test_that("logistic5 with g=1 equals logistic4", {
  p <- test_params
  y4 <- logistic4(x_log, p$a, p$b, p$c, p$d)
  y5 <- logistic5(x_log, p$a, p$b, p$c, p$d, g = 1)
  expect_equal(y4, y5, tolerance = 1e-12)
})

test_that("all models approach asymptotes", {
  p <- test_params
  expect_equal(logistic4(-100, p$a, p$b, p$c, p$d), p$a, tolerance = 1e-6)
  expect_equal(logistic4(100, p$a, p$b, p$c, p$d), p$d, tolerance = 1e-6)
  expect_equal(gompertz4(-100, p$a, p$b, p$c, p$d), p$a, tolerance = 1e-6)
  expect_equal(gompertz4(100, p$a, p$b, p$c, p$d), p$d, tolerance = 1e-6)
})


# =============================================================================
# Inverses: forward → inverse round-trip
# =============================================================================
test_that("inverse functions recover x from y for all five models", {
  p <- test_params
  tol <- 1e-8

  # logistic4
  y4 <- logistic4(x_log, p$a, p$b, p$c, p$d)
  x_back4 <- inv_logistic4(y4, p$a, p$b, p$c, p$d)
  ok4 <- !is.na(x_back4)
  expect_true(sum(ok4) > 10)
  expect_equal(x_back4[ok4], x_log[ok4], tolerance = tol)

  # logistic5
  y5 <- logistic5(x_log, p$a, p$b, p$c, p$d, p$g)
  x_back5 <- inv_logistic5(y5, p$a, p$b, p$c, p$d, p$g)
  ok5 <- !is.na(x_back5)
  expect_true(sum(ok5) > 10)
  expect_equal(x_back5[ok5], x_log[ok5], tolerance = tol)

  # loglogistic4 (positive x only)
  y_ll4 <- loglogistic4(x_pos, p$a, p$b, p$c, p$d)
  x_back_ll4 <- inv_loglogistic4(y_ll4, p$a, p$b, p$c, p$d)
  expect_equal(x_back_ll4, x_pos, tolerance = tol)

  # loglogistic5
  y_ll5 <- loglogistic5(x_log, p$a, p$b, p$c, p$d, p$g)
  x_back_ll5 <- inv_loglogistic5(y_ll5, p$a, p$b, p$c, p$d, p$g)
  expect_equal(x_back_ll5, x_log, tolerance = tol)

  # gompertz4
  y_gom <- gompertz4(x_log, p$a, p$b, p$c, p$d)
  x_back_gom <- inv_gompertz4(y_gom, p$a, p$b, p$c, p$d)
  expect_equal(x_back_gom, x_log, tolerance = tol)
})


# =============================================================================
# Derivatives: dy/dx positive for all models
# =============================================================================
test_that("first derivatives are positive for all five models", {
  p <- test_params

  expect_true(all(dydx_logistic4(x_log, p$a, p$b, p$c, p$d) > 0))
  expect_true(all(dydx_logistic5(x_log, p$a, p$b, p$c, p$d, p$g) > 0))
  expect_true(all(dydx_loglogistic4(x_pos, p$a, p$b, p$c, p$d) > 0))
  expect_true(all(dydx_loglogistic5(x_log, p$a, p$b, p$c, p$d, p$g) > 0))
  expect_true(all(dydx_gompertz4(x_log, p$a, p$b, p$c, p$d) > 0))
})

test_that("analytical derivatives match numerical for logistic4", {
  p <- test_params
  h <- 1e-7
  x0 <- 2.0
  analytical <- dydx_logistic4(x0, p$a, p$b, p$c, p$d)
  numerical  <- (logistic4(x0 + h, p$a, p$b, p$c, p$d) -
                   logistic4(x0 - h, p$a, p$b, p$c, p$d)) / (2 * h)
  expect_equal(analytical, numerical, tolerance = 1e-5)
})


# =============================================================================
# Gradients: make_inv_and_grad_fixed dispatch
# =============================================================================
test_that("make_inv_and_grad_fixed works for all five models (a free)", {
  p <- test_params
  y_mid <- (p$a + p$d) / 2

  for (model in available_models()) {
    # Factory takes only model and fixed_a — no y argument
    fns <- make_inv_and_grad_fixed(model, fixed_a = NULL)
    expect_type(fns, "list")
    expect_named(fns, c("inv", "grad", "grad_y"))

    pv <- if (model %in% c("logistic5", "loglogistic5")) {
      c(a = p$a, b = p$b, c = p$c, d = p$d, g = p$g)
    } else {
      c(a = p$a, b = p$b, c = p$c, d = p$d)
    }

    # Closures now take (y, p) — y_mid passed as first argument
    x_val <- fns$inv(y_mid, pv)
    expect_true(is.finite(x_val), info = paste(model, "inv failed"))

    grad_theta <- fns$grad(y_mid, pv)
    expect_true(all(is.finite(grad_theta)), info = paste(model, "grad failed"))

    grad_y_val <- fns$grad_y(y_mid, pv)
    expect_true(is.finite(grad_y_val), info = paste(model, "grad_y failed"))
  }
})

test_that("make_inv_and_grad_fixed works for all five models (a fixed)", {
  p <- test_params
  y_mid <- (p$a + p$d) / 2

  for (model in available_models()) {
    # Factory takes fixed_a — no y argument
    fns <- make_inv_and_grad_fixed(model, fixed_a = p$a)
    expect_type(fns, "list")

    # When a is fixed, p must not contain a
    pv <- if (model %in% c("logistic5", "loglogistic5")) {
      c(b = p$b, c = p$c, d = p$d, g = p$g)
    } else {
      c(b = p$b, c = p$c, d = p$d)
    }

    # Closures take (y, p)
    x_val <- fns$inv(y_mid, pv)
    expect_true(is.finite(x_val), info = paste(model, "inv failed"))
  }
})

# =============================================================================
# Formulas: build_nls_formulas
# =============================================================================
test_that("build_nls_formulas returns formulas for all five models", {
  formulas <- build_nls_formulas(
    model_names = available_models(),
    response_variable = "y",
    independent_variable = "x"
  )
  expect_length(formulas, 5)
  expect_named(formulas, available_models())
  for (f in formulas) {
    expect_s3_class(f, "formula")
    expect_true("a" %in% all.vars(f))  # a is free
  }
})

test_that("build_nls_formulas with fixed_a removes a from free params", {
  formulas <- build_nls_formulas(
    model_names = "logistic4",
    response_variable = "y",
    fixed_a = 1.5,
    is_log_response = FALSE
  )
  vars <- all.vars(formulas$logistic4)
  expect_false("a" %in% vars)
  expect_true("b" %in% vars)
})


# =============================================================================
# Grid generation
# =============================================================================
test_that("generate_prediction_grid works with flat args", {
  grid <- generate_prediction_grid(
    std_curve_conc = 10000,
    n_grid = 50,
    is_log_independent = TRUE
  )
  expect_equal(nrow(grid), 50)
  expect_named(grid, c("log10_concentration", "concentration", "x_fit"))
  expect_true(all(diff(grid$x_fit) > 0))
  expect_equal(grid$x_fit, grid$log10_concentration)
})

test_that("generate_prediction_grid works with S3 backwards compat", {
  ac <- new_antigen_constraints("test", std_curve_conc = 10000)
  fo <- new_fit_options(n_grid = 100)
  grid <- generate_prediction_grid(ac, fo, is_log_independent = TRUE)
  expect_equal(nrow(grid), 100)
})


# =============================================================================
# Preprocessing
# =============================================================================
test_that("preprocess_standards computes concentration and log-transforms", {
  df <- data.frame(dilution = c(1, 10, 100, 1000), mfi = c(20000, 5000, 500, 50))
  result <- preprocess_standards(
    data = df,
    antigen_settings = list(standard_curve_concentration = 10000),
    response_variable = "mfi",
    independent_variable = "concentration",
    is_log_response = TRUE,
    is_log_independent = TRUE,
    apply_prozone = FALSE
  )
  prepped <- result$data
  expect_true("concentration" %in% names(prepped))
  # log10(10000/1) = 4, log10(10000/10) = 3, etc.
  expect_equal(prepped$concentration[1], 4, tolerance = 1e-10)
  # Response should be log10-transformed
  expect_equal(prepped$mfi[1], log10(20000), tolerance = 1e-10)
})


# =============================================================================
# Settings: resolve_fixed_lower_asymptote
# =============================================================================
test_that("resolve_fixed_lower_asymptote returns NULL for default method", {
  ac <- new_antigen_constraints("test", l_asy_method = "default")
  expect_null(resolve_fixed_lower_asymptote(ac))
})

test_that("resolve_fixed_lower_asymptote returns value for user_defined", {
  ac <- new_antigen_constraints("test", l_asy_min = 50, l_asy_max = 50,
                                l_asy_method = "user_defined")
  expect_equal(resolve_fixed_lower_asymptote(ac), 50)
})

test_that("resolve_fixed_lower_asymptote returns value for blank methods", {
  ac <- new_antigen_constraints("test", l_asy_min = 120,
                                l_asy_method = "geometric_mean_of_blanks")
  expect_equal(resolve_fixed_lower_asymptote(ac), 120)
})


# =============================================================================
# Result class
# =============================================================================
test_that("new_calibration_result validates required meta fields", {
  expect_error(
    new_calibration_result(meta = list(method = "frequentist")),
    "missing required fields"
  )
})

test_that("new_calibration_result creates valid object with minimal meta", {
  meta <- list(
    method = "frequentist", package = "test", curve_id = "1",
    response_var = "y", independent_var = "x",
    is_log_response = TRUE, is_log_independent = TRUE
  )
  cr <- new_calibration_result(meta = meta)
  expect_s3_class(cr, "calibration_result")
  expect_equal(cr$meta$curve_id, "1")
})

test_that("agreement_metrics computes correctly", {
  x <- 1:10
  y <- 1:10 + 0.5
  m <- agreement_metrics(x, y)
  expect_equal(m$bias, 0.5)
  expect_equal(m$cor, 1.0)
  expect_true(m$ccc > 0.9)
})


# =============================================================================
# Utils
# =============================================================================
test_that("available_models returns all five", {
  expect_equal(
    available_models(),
    c("logistic4", "loglogistic4", "gompertz4", "logistic5", "loglogistic5")
  )
})

test_that("model_params returns correct params for each model", {
  expect_equal(model_params("logistic4"), c("a", "b", "c", "d"))
  expect_equal(model_params("logistic5"), c("a", "b", "c", "d", "g"))
  expect_equal(model_params("gompertz4"), c("a", "b", "c", "d"))
  expect_error(model_params("bad_model"))
})
