test_that("release_manifest records the versions that determine a result", {
  lock <- tempfile(fileext = ".lock"); writeLines("{}", lock)
  m <- release_manifest("i-spi-compute-worker", "1.1.0", "abc123", cmdstan = FALSE,
                        lockfile = lock, extra = list(image = list(tag = "1.1.0")))
  expect_equal(m$calib_schema, calib_schema_version())
  expect_equal(m$packages$curveRcore, as.character(utils::packageVersion("curveRcore")))
  expect_true(is.na(m$packages[["curveRweights"]]) || is.character(m$packages[["curveRweights"]]))
  expect_equal(m$renv_lock_md5, unname(tools::md5sum(lock)))
  expect_true(is.na(m$cmdstan))
  expect_equal(m$image$tag, "1.1.0")
  expect_match(m$created_at, "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}Z$")
  expect_true(all(c("manifest_version", "component", "version", "commit", "r_version",
                    "platform", "packages", "cmdstan", "calib_schema",
                    "renv_lock_md5") %in% names(m)))
  expect_error(release_manifest(extra = list(1)), "must be named")
})

test_that("manifest_json writes valid, stable JSON", {
  j <- manifest_json(list(a = "x\"y", b = NA, c = list(d = 1.5, e = TRUE), f = c(1, 2),
                          g = NULL, h = Inf, i = list()))
  expect_equal(j, '{"a":"x\\"y","b":null,"c":{"d":1.5,"e":true},"f":[1,2],"g":null,"h":null,"i":[]}')
  m <- release_manifest(cmdstan = FALSE)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    back <- jsonlite::fromJSON(manifest_json(m, pretty = TRUE), simplifyVector = FALSE)
    expect_equal(back$calib_schema, calib_schema_version())
    expect_equal(names(back), names(m))
  }
})

test_that("calib_schema.sql records calib_schema_version() and adds calib_run.manifest", {
  p <- testthat::test_path("..", "..", "R", "calib_schema.sql")
  testthat::skip_if(!file.exists(p), "calib_schema.sql not available (installed package)")
  sql <- paste(readLines(p, warn = FALSE), collapse = "\n")
  expect_match(sql, sprintf("VALUES \\('calib_schema', '%s'", gsub(".", "\\.", calib_schema_version(), fixed = TRUE)))
  expect_match(sql, "ALTER TABLE madi_results.calib_run\n    ADD COLUMN IF NOT EXISTS manifest jsonb;")
  expect_match(sql, "CREATE TABLE IF NOT EXISTS madi_results.schema_version")
})
