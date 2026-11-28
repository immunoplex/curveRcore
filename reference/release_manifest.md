# Composite release manifest

Collects the versions that determine a result: R, the curveR packages,
CmdStan (when cmdstanr and an installation are available), the calib
schema version, the hash of an renv lockfile, and whatever the caller
adds (the component's own name, version and commit). The worker stores
it with every job in `calib_run.manifest`; apps show it and attach it to
downloads.

## Usage

``` r
release_manifest(
  component = NA_character_,
  version = NA_character_,
  commit = NA_character_,
  packages = c("curveRcore", "curveRfreq", "curveRbayes", "curveRweights"),
  cmdstan = TRUE,
  lockfile = NULL,
  extra = list()
)
```

## Arguments

- component:

  Optional name of the component building the manifest, e.g.
  `"i-spi-compute-worker"` or `"i-spi"`.

- version:

  Optional version of that component (e.g. its image tag).

- commit:

  Optional commit (git SHA) of that component.

- packages:

  Packages whose installed versions are recorded (`NA` when not
  installed).

- cmdstan:

  Logical: try to record the CmdStan version through cmdstanr.

- lockfile:

  Optional path to an `renv.lock`; its MD5 is recorded.

- extra:

  Optional named list of further entries (scalars or named lists).

## Value

A named list with entries `manifest_version`, `created_at` (UTC, ISO
8601), `component`, `version`, `commit`, `r_version`, `platform`,
`packages` (named list), `cmdstan`, `calib_schema`, `renv_lock_md5`, and
the entries of `extra`. Serialise with
[`manifest_json()`](https://immunoplex.github.io/curveRcore/reference/manifest_json.md).

## See also

[`manifest_json()`](https://immunoplex.github.io/curveRcore/reference/manifest_json.md),
[`calib_schema_version()`](https://immunoplex.github.io/curveRcore/reference/calib_schema_version.md)

## Examples

``` r
m <- release_manifest("example", "1.0.0", cmdstan = FALSE)
m$calib_schema
#> [1] "0.7.1"
cat(manifest_json(m))
#> {"manifest_version":1,"created_at":"2026-11-28T11:20:08Z","component":"example","version":"1.0.0","commit":null,"r_version":"4.6.2","platform":"x86_64-pc-linux-gnu","packages":{"curveRcore":"0.7.1.9000","curveRfreq":null,"curveRbayes":null,"curveRweights":null},"cmdstan":null,"calib_schema":"0.7.1","renv_lock_md5":null}
```
