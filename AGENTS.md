# Repository Guidelines

## Project Structure & Module Organization
Cloudy’s C++17 sources live in `source/`, with compiler-specific wrappers under `source/sys_*`. Building drops artifacts (`cloudy.exe`, `runtests.exe`) in place, so keep working inside that directory. Runtime atomic and molecular datasets are in `data/`; respect the existing naming so lookups in `atmdat_*` stay valid. Long-form docs are in `docs/`, while `scripts/` contains operational helpers (Perl runners, plotting). Integration inputs and expected outputs sit in `tsuite/` (`tsuite/auto` mirrors nightly coverage). Vendor dependencies (UnitTest++) are staged via `library/` and copied into `source/lib` by the Makefile.

## Build, Test, and Development Commands
- `make -C source -j4` builds `cloudy.exe` using the default GCC/Clang toolchain.
- `make -C source test` compiles the UnitTest++ harness and runs `runtests.exe`.
- `perl tsuite/auto/run_parallel.pl -j4` executes the regression suite; rerun failures with `perl tsuite/auto/rerun_parallel.pl`.
- `./source/cloudy.exe < tsuite/auto/<case>.in` runs a single scenario; capture output with shell redirection.
Use `make -C source clean` before pushing large changes; `make -C source valgrind-test` is available when chasing leaks.

## Coding Style & Naming Conventions
Follow the existing tab-based indentation and brace-on-new-line layout in `source/*.cpp`. Functions and types generally use PascalCase (`InitDefaultsPreparse`), while module-level globals stay lowercase (`nzone`). Prefer descriptive namespaces over new globals, and keep diagnostic hooks (`DEBUG_ENTRY`, `ASSERT`) intact. Do not reformat generated headers such as `cloudyconfig.h`; they are produced by `configure.sh`.

## Testing Guidelines
Unit tests must remain green via `make -C source test`. For physics regressions, add or update monitors in `tsuite/auto` scripts and document intent after the blank line separator. Nightly automation expects no `botched monitor` or `PROBLEM` markers—review `serious.txt`/`minor.txt` if `checkall.pl` emits them. When optimization flags differ, run the suite once without defining `NDEBUG` to exercise assertions.

## Commit & Pull Request Guidelines
Recent history favors short, imperative commit subjects (`Fix typo`, `Update release date`). Group related changes, reference issue IDs when applicable, and note any data/table updates. Pull requests should include: 1) a concise problem statement and solution summary, 2) test evidence (`make test`, `run_parallel.pl` output), and 3) verification notes for any science-facing results (plots, changed monitors). Screenshots or tables help reviewers triage substantial output diffs.
