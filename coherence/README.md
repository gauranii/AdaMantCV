# How the EEG coherence was computed

These scripts record how the EEG input of the paper was made. They are kept as
they were written (2018 and 2021) and are not run by the flake: they need the
raw data in `Data/01_raw/` and `Data/02_coherence/`, which are private and not
in this repository.

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `coherence_functions.R` | | Functions: periodogram matrix of the 64 channels of one trial by FFT (`get_pgram`), averages within the five bands with edges 0, 4, 8, 16, 32, 45 Hz (`get_spectmat`), coherence \|S_jk\|^2 / (\|S_jj\| \|S_kk\|) (`coherence`) |
| 2 | `compute_band_coherence.R` | EEG trials of one subject | The 5 x 64 x 64 band coherence of that subject; over all 350 subjects this gives `Data/02_coherence/vwm_coh_clean.RData` |
| 3 | `build_channel_band_files.R` | `Data/02_coherence/vwm_coh_clean.RData`, `Data/01_raw/eeg_dat_snp_all.RData`, `Data/01_raw/cups_top_500_SNPs.txt` | One file per channel and band: the 64 coherences of that channel and the 497 candidate SNPs (3 of the 500 dropped) |
| 4 | | the files of step 3, `Data/03_triples_channels/triple.csv` | `Data/04_coherence_SNPs/dat_ch<k>.csv`, one per channel-channel-band triple k (coherence, then the 497 SNPs), the input of Table 4 and Figures 5 to 8 |

The script of step 4 was not found; `triple.csv` lists the 250 triples it kept.
The paths inside the scripts (`setwd()`, `load()`) are those of the machine
they were written on.

## The same steps as functions

`R/coherence.R` and `python/adamant/coherence.py` (the same function names)
do what these scripts do, with the settings as arguments. The scripts stay
here unchanged as the record.

| Step | Functions | Agreement with the original |
|---|---|---|
| 1 | `get_pgram()`; `pgram_frequencies()`, `band_of_frequency()`, `band_rows()` (the rows of the periodogram matrix averaged in each band); `coherence_matrix()` | R: identical to `get_pgram()`, `get_spectmat()` and `coherence()` on made-up signals (`tests/testthat/test-coherence.R`). Python: the same band rows; the coherence agrees with R to about 3e-16 (the two Fourier transforms round differently) |
| 2 | `band_spectrum()`, `band_coherence()` | R: identical to the loop of `compute_band_coherence.R` (`freq_cutoff = 0.4`) and to `get_spectmat()` (`freq_cutoff = 50 / time points`) on made-up signals. Not checked on the data: the EEG series are not in the data folder, so `vwm_coh_clean.RData` cannot be rebuilt |
| 3 | `read_snp_names()`, `candidate_snps()`, `channel_band_table()`, `write_channel_band_files()` | R: the files of the loop of `build_channel_band_files.R`, byte for byte. Python: the bytes of R |
| 4 | `triple_table()`, `write_triple_files()` | Reconstructed from the files. Rebuilt from the data, all 250 files have the numbers of `Data/04_coherence_SNPs/` (largest difference 5.3e-15, the rounding of the files to 15 significant digits), and 249 have their bytes when the numbers are written correctly rounded; one field of `dat_ch195.csv` differs in its last digit |

Three things the functions keep from the scripts, because the data were made
this way:

- The band edges are divided by the number of time points and compared with
  the frequencies of `spec.pgram()` (1/N, 2/N, ..., N the padded length),
  while row i of the periodogram matrix is the Fourier frequency (i - 1)/n.
  Every band therefore lies one frequency below its edges: delta is 0 to 3
  cycles per trial, theta 4 to 7, alpha 8 to 15, beta 16 to 31, gamma 32 to
  44.
- The first line of `cups_top_500_SNPs.txt` is read as a header, and rows 164
  and 357 of the remaining 499 are skipped (these two SNPs are not among the
  genotyped ones): 497 SNPs.
- Triple k = (chan1, chan2, band) is column chan1 of the file of step 3 of
  channel chan2: `dat_array[, chan1, chan2, band]`. The coherence array is
  symmetric only up to rounding, and this orientation is the one whose bytes
  match (the other matches 117 of the 250 files).

## Checking steps 3 and 4 against the data

```sh
ADAMANT_DATA="<data folder>" nix develop -c bash checks/check_coherence.sh
```

rebuilds the files of steps 3 and 4 with R and with Python in a temporary
folder, compares the 250 files per triple with `Data/04_coherence_SNPs/` and
the files of R with those of Python, prints counts and the largest
difference (never data), and deletes the folder. It takes about 4 minutes
and 2.7 GB of memory (the genotypes). It is not part of `nix flake check`, because the data
are private.

R writes a number with 15 significant digits by rounding in the precision of
the machine. On x86-64, where the files of the paper were written, it has
extended precision and the digits are almost always the correctly rounded
ones; on Apple silicon it has not, and about 1 % of the numbers are written
with fewer digits (rounded once more). The check therefore reports
the bytes twice: of the files as R and Python write them on this machine,
and of the rebuilt numbers written correctly rounded (`%.15g`).

The unit tests (`tests/testthat/test-coherence.R`,
`python/tests/test_coherence.py`) use made-up signals and run with the other
tests.
