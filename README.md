# AdaMantCV: Adaptive Mantel Test with Cross-Validation

Code for

> Gauran, I. I., Xue, G., Chen, C., Ombao, H., and Yu, Z. (2022). Ridge
> Penalization in High Dimensional Testing with Applications to Imaging
> Genetics. *Frontiers in Neuroscience* 16:836100.
> [doi.org/10.3389/fnins.2022.836100](https://doi.org/10.3389/fnins.2022.836100)

## The idea

Is a response associated with a large set of predictors taken together, such
as a brain measure with hundreds of genetic variants? A ridge-penalized
Mantel test answers this with one statistic, but its power depends on the
ridge penalty, and the best penalty is unknown in advance. AdaMantCV chooses
the penalty by generalized or leave-one-out cross-validation, repeats that
choice on every permutation of the response, and ranks the observed result
among the permutations. The paper studies this test in simulations, on liver
toxicity data, and on EEG coherence and SNP data. The method is set out in
the paper.

## What this repository does

It rebuilds every table and figure of the paper that comes from code
(Tables 1 to 4; Figures 1 to 3 and 5 to 8; Figure 4 is a diagram) from the
archived runs of the original analysis, and checks each one against the
published paper.

The code exists twice, in R (`R/`, `scripts/`) and in Python (`python/`).
The Python version is an independent second implementation: it gives the same
numbers as R, cell by cell for the tables and mark by mark for the figures, so
the results do not depend on one language or one set of libraries. A Nix
flake pins R, Python and every package, so the checks run the same way on any
machine.

## Running

With [Nix](https://nixos.org) (flakes enabled), from the repository root:

```sh
nix run .#table4        # one output: table1 to table4, fig1 to fig3, fig5 to fig8
nix run .#all           # every output, in R
nix run .#py-table4     # the same output from the Python code (py-all for every one)
nix flake check         # build everything, compare with the paper and R with Python, run the tests
nix run .#help          # every command
```

Tables go to `results/` (csv) and figures to `figures/` (pdf), named after
the paper. `checks/RESULTS.md` records the comparison with the paper and
`checks/PYTHON.md` the comparison of R with Python. Every function in R and
Python has unit tests, and GitHub runs `nix flake check` on every push.

## Data

- **Liver toxicity data** (`data/liver/`): the public `liver.toxicity` data
  of the Bioconductor package mixOmics (Bushel et al. 2007; Heinloth et al.
  2004). Used by Figures 1 and 2.
- **EEG coherence and SNPs** of 350 subjects: subject-level data, not in the
  repository. They are needed only to recompute Table 4 and Figures 5 to 8
  from scratch (`--rerun`); the archived runs in `archive/` rebuild them
  without the data.

## More detail

[`docs/DETAILS.md`](docs/DETAILS.md) has the file layout, every output with
its script, input and check, run times, the tests, the data pipeline and the
pinned software versions.

## Citation

Please cite the paper above when you use this code.
