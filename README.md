# Random Limited Attention Methods

The `ramchoice` package implements revealed-preference, attention, estimation,
inference, and specification procedures for random limited-attention models,
including the Random Attention Model (RAM) and the Attention Overload Model
(AOM).

- `aomTest`: conducts sample inference for candidate preferences under
  homogeneous AOM and returns tidy hypothesis-level results. An optional
  cluster identifier enables cluster-robust covariance estimation and
  multiplier critical values.
- `aomModel`: evaluates population AOM compatibility and reports the number
  and magnitude of preference-specific inequality violations.
- `aomIdentify`: represents the sharp homogeneous-AOM preference identified set
  through a mixed-integer formulation. It tests model feasibility and every
  pairwise comparison without enumerating all strict rankings, while
  distinguishing certified infeasibility from other solver failures.
- `hlaoModel`: recovers list-based attention and computes sharp independent,
  dependence-robust, or path-independence-robust bounds for heterogeneous-
  preference events and full-attention agreement. The no-SPI mode permits
  arbitrary observed-menu domains. For benchmark independence, structured
  events use status-checked column generation when ranking enumeration exceeds
  the requested cap. Returned diagnostics include solver statuses, reduced
  costs, primal and dual residuals, an optimality-gap bound, tolerance, and
  certificate checks.
- `hlaoEvent`: describes pairwise, joint-above, and top-choice preference
  events in a form that can be optimized by H-LAO column generation.
- `hlaoTest`: provides weak-reach pairwise intervals and dependence-robust
  projection inference for H-LAO using finite-sample Hoeffding or
  covariance-aware correlated-Gaussian probability bands, with the common
  error budget split between Gaussian and exact-binomial cells when both are
  used. It also reports a studentized undivided-moment inversion for pairwise
  preference shares, together with simultaneous attention-overload and
  Block--Marschak specification diagnostics. Conservative outer-region
  diagnostics remain the default; an opt-in direct delta-Gaussian method is
  available on complete menu domains with positive terminal reach. Supplying
  a cluster identifier preserves cross-menu dependence within clusters and
  replaces Gaussian calibration with cluster multiplier draws and a
  cluster-Hoeffding fallback.
- `hlaoNoPITest`: provides exact finite-sample LP projection inference for the
  path-independence-robust H-LAO model, with the same optional clustered
  probability region.
- `hlaoRankings`: enumerates strict rankings for direct H-LAO calculations on
  small and moderate universes.
- `revealPref`: tests candidate preference orderings and conducts
  revealed-preference inference under RAM and AOM.
- `ramTest`: returns tidy RAM candidate-ranking inference with optional
  cluster-robust covariance and multiplier critical values.
- `revealAtte`: computes bounds on attention frequencies under AOM.
- `revealPrefModel`: checks whether population choice probabilities are
  compatible with RAM and/or AOM.
- `sumData`: constructs menu frequencies and empirical choice probabilities
  from individual choice data.
- `genMat`: generates the constraint matrices used in preference analysis.
- `logitAtte`: computes choice probabilities and attention frequencies under
  the logit attention rule.
- `logitSimu`: simulates choice data under the logit attention rule.
- `rAtte`: legacy interface retained for backward compatibility; new code
  should use `revealPref`.

The package also includes the simulated dataset `ramdata` for illustration.

## R Implementation

To install or update the released R package, type:

```r
install.packages("ramchoice")
```

- Help: [R manual](https://cran.r-project.org/web/packages/ramchoice/ramchoice.pdf),
  [CRAN repository](https://cran.r-project.org/package=ramchoice).

To install the development version from GitHub, type:

```r
remotes::install_github("mdcattaneo/ramchoice", subdir = "R/ramchoice")
```

## References

### Random Attention Model

- Cattaneo, Ma, Masatlioglu, and Suleymanov (2020):
  [A Random Attention Model](https://doi.org/10.1086/706861).<br>
  *Journal of Political Economy* 128(7): 2796-2836.<br>
  [Supplemental Appendix](https://mdcattaneo.github.io/papers/Cattaneo-Ma-Masatlioglu-Suleymanov_2020_JPE--Supplement.pdf)

### Attention Overload Model

- Cattaneo, Cheung, Ma, and Masatlioglu (2026):
  [Attention Overload](https://arxiv.org/abs/2110.10650).<br>
  Working paper.
