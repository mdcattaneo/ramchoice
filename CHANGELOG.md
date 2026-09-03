# Changelog

Notable project changes are listed from newest to oldest.

## 3.0.0 (2026-09-03)

- Added subject-clustered inference to `aomTest()`, `hlaoTest()`, and
  `hlaoNoPITest()`, together with a tidy cluster-capable `ramTest()` interface.
  The implementation estimates the joint covariance of all
  menu--outcome frequencies from cluster influence vectors, preserves
  cross-menu dependence, uses cluster multiplier calibration, and supplies a
  cluster-Hoeffding fallback for sparse cells.
- Fixed two-step AOM moment selection when exactly one effective inequality is
  retained; the selected simulated statistics now preserve their matrix
  dimension.
- Corrected the documentation of the legacy GMS tuning default. A null
  `MNRatioGMS` uses `1/log(N)`, so the Gaussian recentering multiplies the
  estimated moments by `1/sqrt(log(N))`; this clarification does not change
  numerical behavior.
- Added `aomIdentify()` with an exact mixed-integer characterization of the
  sharp homogeneous-AOM compatible-preference set. The routine reports model
  feasibility and pairwise revealed preference without factorial ranking
  enumeration, while retaining the legacy candidate-ranking interfaces.
- Added `hlaoEvent()` and exact restricted-master/MILP-pricing column
  generation for structured benchmark H-LAO events. `hlaoModel()` now switches
  automatically from direct enumeration when the ranking count exceeds
  `max_rankings` and reports active-column and pricing diagnostics.
- Added full-attention agreement to `hlaoModel()`. Complete positive-reach
  designs use the direct full-attention formula; incomplete independent,
  preference--stopping-dependent, and no-SPI designs use sharp linear-program
  bounds.
- Corrected the benchmark H-LAO population LP to impose preference-probability
  adding up explicitly. The equation is redundant under positive reach but is
  necessary for sharp event bounds when every observed menu has zero reach.
- Added studentized Bonferroni inversion of undivided H-LAO pairwise moments,
  including exact quadratic confidence-set components and an uninformative
  convention at degenerate zero reach.
- Added the sharp path-independence-robust H-LAO population polytope through
  `hlaoModel(..., dependence = "noPI")` and `dependence = "all"`.
- Added `hlaoNoPITest()` for exact linear projection of simultaneous primitive
  probability bands without Sequential Path Independence or suffix closure.
- Added simultaneous H-LAO specification diagnostics for recovered-attention
  overload and Block--Marschak violations, including conservative recursive
  bands for the recovered full-attention rule and an opt-in direct
  delta-Gaussian diagnostic for regular complete-menu designs.
- Added covariance-aware correlated-Gaussian H-LAO probability bands with an
  exact-binomial safeguard for sparse and degenerate menu--outcome cells and a
  common error-budget split when both components are active.
- Added `hlaoModel()`, `hlaoTest()`, and `hlaoRankings()` for population
  H-LAO identification, weak-reach pairwise inference, and finite-sample
  dependence-robust projection bounds.
- Added an internal one-entry cache for `genMat()` so repeated simulations with
  unchanged menus and preferences reuse the same restriction matrix.
- Added `aomTest()` and `aomModel()` as model-specific homogeneous-AOM
  interfaces with tidy inference output and population violation diagnostics,
  while retaining the complete legacy calculation for auditing.
- Started the version 3.0 development line.
- Added automated tests for the version 2.2 public API, result contracts,
  seeded numerical fixtures, and the `rAtte()` compatibility alias.
- Recorded the exact 2020 JPE replication checkpoint and legacy AOM numerical
  baselines in `COMPATIBILITY.md`.
- Established the GitHub development repository with an R-only package layout,
  continuous integration, and standard project-maintenance templates.
- Imported `ramchoice` 2.2 as the backward-compatibility baseline for future
  RAM and AOM development.

## 2.2

- Baseline CRAN release predating this development repository.
