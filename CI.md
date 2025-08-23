CI layout

This project uses a professional CI/CD layout with the following workflows in `.github/workflows`:

- `ci.yml` — primary continuous integration: runs tests across platforms (Linux, macOS, Windows) and Julia versions, caches packages, captures logs and coverage artifacts.
- `docs.yml` — builds documentation if `docs/make.jl` exists.
- `release.yml` — runs tests when a tag `v*.*.*` is pushed and creates a GitHub release with a source artifact.
- `Unit_tests.yml` — removed; use `ci.yml` as single source of truth.

Local reproduction

Run tests locally in the project environment:

```bash
julia --project=@. -e 'import Pkg; Pkg.instantiate(); Pkg.test()'
```

If you want to run the same steps as CI (instantiate and precompile):

```bash
julia --project=@. -e 'import Pkg; Pkg.instantiate(); Pkg.precompile()'
```

Notes

- 32-bit i386 emulation is available via the `Unit_tests.yml`-derived emulation job (QEMU + Docker), but it's experimental and non-blocking.
- For production-grade 32-bit testing, provisioning a 32-bit self-hosted runner is recommended.

Optional integrations

- Code coverage upload to Codecov is available when you set the `CODECOV_TOKEN` secret in the repository. The CI step is conditional and will be skipped when the token is missing.
- Documentation build uses Documenter; a minimal `docs/make.jl` and `docs/Project.toml` are included in this repo. The docs job will run only if `docs/make.jl` exists.

Registry readiness

- To prepare this package for registry submission, ensure `Project.toml` contains `name`, `uuid`, `version`, `license`, `authors`, and `compat` for Julia.
