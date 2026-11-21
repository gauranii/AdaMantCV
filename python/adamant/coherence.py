# -----------------------------------------------------------------------------
# File:    python/adamant/coherence.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Input of Table 4 and Figures 5 to 8.
# Purpose: How the EEG input of the paper was made, as functions: the
#          periodogram matrix of the channels of one trial, its averages
#          within the five frequency bands and over the trials, the band
#          coherence, the files per channel and band (coherences and the 497
#          candidate SNPs) and the files per channel-channel-band triple.
#          Channels, bands and rows are numbered from 1, as in R.
# Inputs:  None (the data are private; see coherence/README.md).
# Outputs: None.
# Usage:   from adamant.coherence import band_coherence
#          band_coherence(eeg, freq_cutoff=0.4)
# Source:  R/coherence.R.
# -----------------------------------------------------------------------------
import functools
import os

import numpy as np

from .setup import format_number, read_csv

# Edges of the bands delta, theta, alpha, beta, gamma, in cycles per trial
# (Hz for trials of one second).
COHERENCE_BAND_EDGES = np.array([0.0, 4.0, 8.0, 16.0, 32.0, 45.0])

# Rows of the table of candidate SNPs that are used: its 499 rows without the
# two SNPs (rows 164 and 357) that are not among the genotyped SNPs.
SNP_ROWS = [*range(1, 164), *range(165, 357), *range(358, 500)]


# ----- Step 1: periodogram, band averages, coherence -----

def get_pgram(eeg_trial):
    """Periodogram matrix of the channels of one trial.

    R's `mvfft()` is the plain discrete Fourier transform, sum over t of
    x_t exp(-2 pi i k t / n), with no taper, detrending, padding or scaling;
    `numpy.fft.fft()` is the same transform. The two use different
    algorithms, so the values agree to rounding (about 1e-13 relative), not
    bit for bit.

    Parameters: eeg_trial (array, channels x time points).
    Returns: complex array, time points x channels x channels: at Fourier
    frequency k / time points (k from 0), the outer product d_k d_k^H of the
    transforms of the channels.
    """
    d = np.fft.fft(np.asarray(eeg_trial, dtype=np.float64), axis=1)
    return d.T[:, :, None] * np.conj(d.T)[:, None, :]


def nextn(n):
    """`nextn(n)`: the smallest integer >= n that is a product of powers of
    2, 3 and 5. Returns: int."""
    while True:
        m = n
        for factor in (2, 3, 5):
            while m % factor == 0:
                m //= factor
        if m == 1:
            return n
        n += 1


def pgram_frequencies(n_time):
    """Frequencies `spec.pgram()` reports for a series.

    The original code took its frequencies from
    `spec.pgram(x, plot = FALSE)$freq`. Only the frequencies were used, so
    the taper and the detrending of `spec.pgram()` play no part; its padding
    does: the series is padded to `nextn(n_time)` points, and the
    frequencies are `seq.int(from = 1/N, by = 1/N, length.out = N %/% 2)`,
    computed here as R computes them.

    Parameters: n_time (number of time points).
    Returns: float array, 1/N, 2/N, ..., (N // 2)/N in cycles per time
    point, with N = nextn(n_time).
    """
    N = nextn(n_time)
    return 1 / N + np.arange(N // 2) * (1 / N)


def band_of_frequency(freq, n_time):
    """Band of each frequency, as `cut(freq, breaks = edges / n_time)`.

    Parameters: freq (frequencies in cycles per time point), n_time (number
    of time points of a trial).
    Returns: int array: band b (1 to 5) when
    edge[b] / n_time < freq <= edge[b + 1] / n_time, otherwise 0.
    """
    band = np.searchsorted(COHERENCE_BAND_EDGES / n_time, freq, side="left")
    return np.where(band < len(COHERENCE_BAND_EDGES), band, 0)


def band_rows(n_time, freq_cutoff):
    """Rows of the periodogram matrix that are averaged in each band.

    As in the original code: the rows of `get_pgram()` are first restricted
    to those whose frequency is below `freq_cutoff`, the frequencies being
    recycled along the rows, and the rows of band b are then those at the
    positions of the frequencies of band b. Row i of `get_pgram()` is the
    Fourier frequency (i - 1) / n_time while `pgram_frequencies()` starts at
    1 / N, so every band lies one row below its edges, and the delta band
    includes frequency zero.

    Parameters: n_time (number of time points of a trial), freq_cutoff
    (frequencies at or above it are dropped, in cycles per time point: 0.4
    in compute_band_coherence.R, 50 / n_time in `get_spectmat()`).
    Returns: list of five int arrays of row numbers, counted from 1.
    """
    freq = pgram_frequencies(n_time)
    kept = np.flatnonzero(np.resize(freq < freq_cutoff, n_time)) + 1
    band = band_of_frequency(freq, n_time)
    return [kept[np.flatnonzero(band == b)] for b in range(1, len(COHERENCE_BAND_EDGES))]


# ----- Step 2: band coherence of one subject -----

def band_spectrum(eeg, freq_cutoff):
    """Spectral matrix of one subject in the five bands.

    The means add their terms one after the other, as R's `mean()` does
    (where R adds in extended precision, on x86-64, the last bit can
    differ).

    Parameters: eeg (array, channels x time points x trials), freq_cutoff
    (see `band_rows()`).
    Returns: complex array, 5 x channels x channels: the periodogram matrix
    averaged over the rows of each band, then over the trials.
    """
    eeg = np.asarray(eeg, dtype=np.float64)
    n_chans, n_time, n_trials = eeg.shape
    rows = band_rows(n_time, freq_cutoff)
    spectmat = np.zeros((len(rows), n_chans, n_chans), dtype=np.complex128)
    for trial in range(n_trials):
        pgram = get_pgram(eeg[:, :, trial])
        for band, band_rows_ in enumerate(rows):
            total = np.zeros((n_chans, n_chans), dtype=np.complex128)
            for row in band_rows_:
                total += pgram[row - 1]
            spectmat[band] += total / len(band_rows_)
    return spectmat / n_trials


def coherence_matrix(spectmat):
    """Coherence from a spectral matrix.

    Parameters: spectmat (complex array, bands x channels x channels).
    Returns: float array of the same size: |S_jk|^2 / (|S_jj| |S_kk|) in
    each band.
    """
    spectmat = np.asarray(spectmat)
    power = np.abs(np.diagonal(spectmat, axis1=1, axis2=2))
    return np.abs(spectmat) ** 2 / (power[:, :, None] * power[:, None, :])


def band_coherence(eeg, freq_cutoff):
    """Band coherence of one subject.

    Parameters: eeg (array, channels x time points x trials), freq_cutoff
    (see `band_rows()`).
    Returns: float array, 5 x channels x channels.
    """
    return coherence_matrix(band_spectrum(eeg, freq_cutoff))


# ----- Step 3: one file per channel and band -----

def read_snp_names(file):
    """Names of the candidate SNPs that are used.

    Parameters: file (path of cups_top_500_SNPs.txt, one SNP per line; as in
    the original code its first line is taken as a header, which leaves
    499).
    Returns: list of the 497 names in the rows `SNP_ROWS`.
    """
    with open(file) as handle:
        names = handle.read().split()[1:]
    return [names[row - 1] for row in SNP_ROWS]


def candidate_snps(genotypes, snp_names):
    """Genotypes of the candidate SNPs.

    Parameters: genotypes (DataFrame or dict of columns, subjects x SNPs,
    by SNP name), snp_names (names of the SNPs to take).
    Returns: float array, subjects x SNPs in the order of `snp_names`.
    """
    return np.column_stack([np.asarray(genotypes[name], dtype=np.float64) for name in snp_names])


def channel_band_table(coherences, snps, channel, band):
    """Coherences of one channel in one band, and the SNPs.

    Parameters: coherences (array, subjects x channels x channels x bands),
    snps (array of `candidate_snps()`), channel, band (numbers from 1).
    Returns: float array: the coherence of every channel with `channel`,
    `coherences[:, :, channel - 1, band - 1]`, then the SNPs.
    """
    return np.column_stack([coherences[:, :, channel - 1, band - 1], snps])


def write_table(table, file, formatted=format_number):
    """`write.csv(as.data.frame(table), file)`: columns "V1", "V2", ..., row
    names "1", "2", ..., numbers with 15 significant digits as R writes them.

    Parameters: table (2-d float array), file (path), formatted (function
    giving the text of a number).
    Returns: None.
    """
    with open(file, "w") as handle:
        handle.write(",".join(['""'] + [f'"V{j}"' for j in range(1, table.shape[1] + 1)]) + "\n")
        for i, row in enumerate(table.tolist(), start=1):
            handle.write(",".join([f'"{i}"'] + [formatted(value) for value in row]) + "\n")


def write_channel_band_files(coherences, snps, band, out_dir):
    """Write the files dat_ch<channel>b<band>.csv of one band.

    Parameters: coherences, snps, band (as for `channel_band_table()`),
    out_dir (folder the files are written to).
    Returns: None.
    """
    formatted = functools.lru_cache(maxsize=None)(format_number)  # the SNPs repeat in every file
    for channel in range(1, coherences.shape[2] + 1):
        write_table(channel_band_table(coherences, snps, channel, band),
                    os.path.join(out_dir, f"dat_ch{channel}b{band}.csv"), formatted)


# ----- Step 4: one file per triple -----

def triple_table(channel_band, channel, n_chans):
    """Coherence of one channel pair and the SNPs, from a table of step 3.

    Parameters: channel_band (array of `channel_band_table()`, or the file
    written from it without its row names, of one channel and band), channel
    (number of the other channel, from 1), n_chans (number of channels: the
    columns before the SNPs).
    Returns: float array: the coherence in column `channel` of the table,
    then its SNPs.
    """
    return np.column_stack([channel_band[:, channel - 1], channel_band[:, n_chans:]])


def write_triple_files(triples, n_chans, channel_band_dir, out_dir):
    """Write the files dat_ch<k>.csv of a set of triples.

    Triple k (chan1, chan2, band) is cut from the file of step 3 of channel
    chan2 and that band, column chan1: `coherences[:, chan1, chan2, band]`
    counted from 1. (The coherence array is symmetric only up to rounding;
    this is the orientation of the files of the paper.)

    Parameters: triples (DataFrame with columns `chan1`, `chan2`, `band`,
    indexed by the triple numbers k: `triple.csv` read with `index_col=0`),
    n_chans (number of channels), channel_band_dir (folder with the files
    of step 3), out_dir (folder the files are written to).
    Returns: None.
    """
    formatted = functools.lru_cache(maxsize=None)(format_number)  # the SNPs repeat in every file
    for k, triple in triples.iterrows():
        channel_band = read_csv(os.path.join(channel_band_dir,
                                             f"dat_ch{triple['chan2']}b{triple['band']}.csv"),
                                index_col=0).to_numpy(dtype=np.float64)
        write_table(triple_table(channel_band, int(triple["chan1"]), n_chans),
                    os.path.join(out_dir, f"dat_ch{k}.csv"), formatted)
