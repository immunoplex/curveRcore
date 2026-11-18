# Serialise a release manifest to JSON

Minimal, dependency-free JSON writer for manifests: named lists become
objects, unnamed vectors of length \> 1 arrays, `NA`/`NULL` `null`.
Suitable for `calib_run.manifest` (jsonb) and JSON sidecars.

## Usage

``` r
manifest_json(x, pretty = FALSE)
```

## Arguments

- x:

  A manifest from
  [`release_manifest()`](https://immunoplex.github.io/curveRcore/reference/release_manifest.md)
  (or any list of scalars and named lists).

- pretty:

  Logical: indent with two spaces.

## Value

Character scalar.
