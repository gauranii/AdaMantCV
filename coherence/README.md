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
