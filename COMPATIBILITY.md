# Backward-Compatibility Baseline

This document records the version 2.2 behavior that the 3.0 development line
must preserve unless a change is explicitly documented. The baseline was
audited on July 13, 2026, using R 4.6.1 and MASS 7.3-65.

## Version 3.0 Release Verification

On September 3, 2026, the unchanged published replication script
`replication-CMMS_2020_JPE/Cattaneo-Ma-Masatlioglu-Suleymanov_2020_JPE.R`
was rerun against the installed `ramchoice` 3.0.0 source build using R 4.5.2
on Ubuntu 26.04.1. The five rejection frequencies exactly matched the
published checkpoints recorded below. The 3.0.0 source package also passed
all 272 automated tests, examples, and the PDF-manual check.

## Source Baseline

The tracked package source was compared with the archived `v_2.2/ramchoice`
source. The only differences are removal of trailing whitespace and a final
blank line. Package code, documentation, data, and generated namespace content
otherwise match the archived source.

## Public API

| Function | Version 2.2 role | Compatibility requirement |
|---|---|---|
| `sumData()` | Aggregate menu counts, choice probabilities, and covariance estimates | Preserve arguments and the six named result components |
| `genMat()` | Construct RAM/AOM moment-inequality matrices | Preserve arguments and `R`/`ConstN` outputs |
| `logitAtte()` | Compute logit-attention choice probabilities | Preserve arguments and `choiceProb`/`atteFreq` outputs |
| `logitSimu()` | Simulate logit-attention choice data | Preserve arguments, RNG behavior, and `menu`/`choice` outputs |
| `revealPref()` | Sample preference inference under RAM/AOM | Preserve arguments, class, and nested result fields |
| `rAtte()` | Legacy alias for `revealPref()` | Keep as an exact compatibility alias |
| `revealAtte()` | Attention-frequency bounds | Preserve arguments, class, and bound/critical-value fields |
| `revealPrefModel()` | Population model compatibility | Preserve arguments, class, and inequality fields |

The registered S3 classes are `ramchoiceRevealPref`,
`ramchoiceRevealAtte`, and `ramchoiceRevealPrefModel`, each with `print()` and
`summary()` methods.

## Replication Bindings

- The 2020 JPE replication calls `logitSimu()` and `revealPref()`. It reads
  `logitSimu()$menu`, `logitSimu()$choice`, `revealPref()$Tstat`, and
  `revealPref()$critVal$GMS`.
- The archived homogeneous-AOM replication calls `rAtte()` and reads the same
  test-statistic and GMS critical-value fields.
- The archived heterogeneous-preference replication calls only `sumData()`;
  its lower-bound routine is replication-local and is not a version 2.2
  package interface.

## Numerical Checkpoints

The unchanged 2020 JPE simulation was run from
`replication-CMMS_2020_JPE/Cattaneo-Ma-Masatlioglu-Suleymanov_2020_JPE.R`.
With seed 42, 1,000 Monte Carlo repetitions, sample size 200 per menu, and
2,000 critical-value draws, version 2.2 returned:

```text
H1     H2     H3     H4     H5
0.000  0.004  0.144  0.275  0.351
```

These values exactly equal the targets printed in the published replication
script. Elapsed time was 350.77 seconds on the audit machine.

The archived homogeneous-AOM files `AOM/Result-1-.txt` through
`AOM/Result-7-.txt` each contain 2,000 rows and 12 rejection indicators. They
are the panel-level numerical baselines for the new homogeneous-AOM pipeline.
For panel 1, their rejection frequencies are:

```text
n=50:   0.0135  0.1590  0.0385  0.0795
n=100:  0.0060  0.2595  0.0645  0.1390
n=200:  0.0040  0.4680  0.1355  0.3440
```

Automated tests additionally freeze exported names, function signatures,
seeded choice-simulation output, constraint dimensions, result classes and
field names, deterministic test statistics, and the `rAtte()` alias. Simulated
critical values are checked for their contracts, ordering, admissible ranges,
and same-platform seeded reproducibility rather than exact cross-platform
draws. The eigenvector-based covariance square root used by `MASS::mvrnorm()`
can legitimately change signs across BLAS/LAPACK implementations while
preserving the simulated distribution.

## Known Version 2.2 Issues

These issues are not compatibility promises. Fix them only in focused changes
with tests and changelog entries:

- `revealPref()` tests `is.character("method")` instead of the supplied
  `method` object.
- `revealPrefModel()` has its probability-row-sum validation commented out.
- The one-menu early return from `genMat()` uses `constN`, whereas regular
  returns use `ConstN`.

Printed whitespace and internal helper organization may change, provided the
documented objects and numerical results remain compatible.
