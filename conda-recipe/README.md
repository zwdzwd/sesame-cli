# Conda recipe for sesame-cli

Builds the `sesame` binary and publishes it to the **zhou-lab** channel.

## Install (end users)

```sh
conda install -c zhou-lab -c conda-forge sesame-cli
```

This installs the `sesame` command. (The conda package is named `sesame-cli` to
avoid confusion with the Bioconductor R package `sesame`; the binary is `sesame`.)

## Build & upload (maintainers)

Requires `conda-build` and `anaconda-client`:

```sh
conda install -n base conda-build anaconda-client

# Build (Linux + macOS; Windows is skipped). conda-forge supplies zlib + libcurl.
conda build -c conda-forge conda-recipe/

# Upload the built package to the zhou-lab channel
anaconda login
anaconda upload -u zhou-lab $(conda build -c conda-forge conda-recipe/ --output)
```

To upload automatically on every successful build:

```sh
conda config --set anaconda_upload yes
```

## Automated builds (CI)

A `.github/workflows/conda-build.yml` (mirroring cinderplot's) can run
`conda build` on every push/PR across **linux-64**, **osx-64**, and
**osx-arm64**, uploading each package as a build artifact, and publish to the
`zhou-lab` channel on a `vX.Y.Z` tag (guarded so the tag must match the recipe
version). Publishing needs an `ANACONDA_TOKEN` repository secret with upload
rights to the org. Cutting a release:

```sh
# bump `version` in conda-recipe/meta.yaml and SESAME_VERSION in include/sesame.h,
# commit, then:
git tag -a v0.1.0 -m "sesame-cli 0.1.0" && git push origin v0.1.0
```

## Notes

- **In-repo build.** The recipe builds the checked-out tree (`source: path: ..`),
  including the vendored `YAME/` + `htslib` submodule, so a release build is just:
  check out the `vX.Y.Z` tag, then `conda build`. Keep `version` in `meta.yaml`
  and `SESAME_VERSION` in `include/sesame.h` in sync with the tag.
- **Submodule.** `YAME/` (and its nested `htslib/`) must be checked out before
  building (`git submodule update --init --recursive`); `path:` copies whatever
  is on disk. `htslib`'s `config.h` and `version.h` are committed, so no autoconf
  / git step is needed at build time.
- **Rebuilds the vendored static libs.** `YAME/libyame.a` and
  `YAME/htslib/libhts.a` are committed for the author's platform; `build.sh`
  deletes them (and the stale `.o`) so they rebuild for the target arch.
- **Flag handling.** The vendored Makefiles have rigid `CFLAGS` and ignore
  `CPPFLAGS`, so `build.sh` bakes conda's compile flags into `CC` (which every
  sub-make uses verbatim). `LDFLAGS` reaches the final link via the environment.
- **Dependencies.** Only zlib and libcurl (both from conda-forge, with
  `run_exports` that pin the runtime automatically). Make sure conda-forge is on
  your channel list when building.
```
