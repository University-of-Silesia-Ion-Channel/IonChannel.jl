# IonChannel

IonChannel.jl — a Julia package for detecting and analyzing ion channel events in time-series data.

This repository contains:

- src/: package source code
- test/: unit and integration tests
- scripts/: helper scripts for CI and validation

Quick start

1. Activate project and instantiate:

   julia --project=@. -e 'using Pkg; Pkg.instantiate()'

2. Run tests:

   julia --project=@. -e 'using Pkg; Pkg.test()'

Contributing

Please open issues or PRs. See `CI.md` for how CI is configured.
