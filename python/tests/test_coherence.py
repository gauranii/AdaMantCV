# -----------------------------------------------------------------------------
# File:    python/tests/test_coherence.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Unit tests of python/adamant/coherence.py: the bands of the
#          frequencies, the properties of the coherence, the files of steps 3
#          and 4 (the cases of tests/testthat/test-coherence.R), and
#          agreement with R/coherence.R on made-up input: the same band rows,
#          the same coherence up to the rounding of the two Fourier
#          transforms, the same files byte for byte.
# Inputs:  R/coherence.R (run with Rscript by the tests marked against_r);
#          made-up signals and files in temporary folders (the EEG and SNP
#          data are not used).
# Outputs: None.
# Usage:   python -m pytest python/tests/test_coherence.py   # inside `nix develop`
# -----------------------------------------------------------------------------
import os
import subprocess

import numpy as np
import pandas as pd
import pytest

from adamant.coherence import (band_coherence, band_of_frequency, band_rows, band_spectrum,
                               candidate_snps, channel_band_table, coherence_matrix, get_pgram,
                               nextn, pgram_frequencies, read_snp_names, triple_table,
                               write_channel_band_files, write_table, write_triple_files)
from adamant.setup import REPO_ROOT
from inputs import close

R_SCRIPT = """
a <- commandArgs(trailingOnly = TRUE)
source(file.path(a[1], "R", "coherence.R"))
work <- a[2]
read_array <- function(name, dims) array(readBin(file.path(work, name), "double", prod(dims)), dims)

# band rows and coherence of the recordings
sizes <- matrix(scan(file.path(work, "sizes.txt"), quiet = TRUE), ncol = 4, byrow = TRUE)
out <- file(file.path(work, "coherence.bin"), "wb")
for (i in seq_len(nrow(sizes))) {
  eeg <- read_array(sprintf("eeg%d.bin", i), sizes[i, 1:3])
  cutoff <- sizes[i, 4] / sizes[i, 2]
  rows <- band_rows(sizes[i, 2], cutoff)
  writeBin(c(lengths(rows), unlist(rows)), out, size = 4L)
  writeBin(as.vector(band_coherence(eeg, cutoff)), out)
}
close(out)

# steps 3 and 4
coherences <- read_array("coherences.bin", c(6, 4, 4, 5))
genotypes <- as.data.frame(read_array("genotypes.bin", c(6, 498)))
colnames(genotypes) <- readLines(file.path(work, "genotyped.txt"))
snps <- candidate_snps(genotypes, read_snp_names(file.path(work, "snps.txt")))
dir.create(file.path(work, "r3")); dir.create(file.path(work, "r4"))
for (band in 1:5) write_channel_band_files(coherences, snps, band, file.path(work, "r3"))
triples <- read.csv(file.path(work, "triple.csv"), row.names = 1)
write_triple_files(triples, 4, file.path(work, "r3"), file.path(work, "r4"))
"""


def made_up_eeg(n_chans, n_time, n_trials, cycles, seed):
    """A made-up recording: white noise in every channel, and one sinusoid
    shared by channels 1 and 2 with a phase that changes from trial to trial.

    Parameters: n_chans, n_time, n_trials (size of the recording), cycles
    (cycles per trial of the shared sinusoid), seed.
    Returns: float array, channels x time points x trials.
    """
    rng = np.random.default_rng(seed)
    eeg = rng.standard_normal((n_chans, n_time, n_trials))
    t = np.arange(1, n_time + 1)
    for trial in range(n_trials):
        wave = 30 * np.sin(2 * np.pi * cycles * t / n_time + rng.uniform(0, 2 * np.pi))
        eeg[:2, :, trial] += wave
    return eeg


def made_up_inputs(folder):
    """A made-up coherence array, SNP table and SNP list for steps 3 and 4.

    Parameters: folder (the list of SNP names is written to it as snps.txt).
    Returns: dict with `coherences` (6 subjects x 4 x 4 channels x 5 bands),
    `genotypes` (DataFrame of 498 SNPs: the 500 listed without two), `file`
    (the list of 500 names), `names`.
    """
    rng = np.random.default_rng(6)
    names = [f"rs{1000 + j}" for j in range(1, 501)]
    file = os.path.join(folder, "snps.txt")
    with open(file, "w") as handle:
        handle.write("\n".join(names) + "\n")
    genotyped = [name for j, name in enumerate(names, start=1) if j not in (165, 358)]
    genotypes = pd.DataFrame(rng.standard_normal((6, 498)), columns=rng.permutation(genotyped))
    return {"coherences": rng.uniform(size=(6, 4, 4, 5)), "genotypes": genotypes,
            "file": file, "names": names}


def made_up_triples():
    """Two triples, numbered 7 and 9, of the channels of `made_up_inputs()`."""
    return pd.DataFrame({"chan1": [1, 2], "chan2": [4, 3], "band": [2, 5]}, index=[7, 9])


def test_get_pgram_is_the_outer_product_of_the_fourier_transforms():
    trial = made_up_eeg(4, 60, 1, cycles=5, seed=1)[:, :, 0]
    pgram = get_pgram(trial)
    d = np.fft.fft(trial, axis=1)
    assert pgram.shape == (60, 4, 4)
    assert np.isclose(pgram[6, 1, 2], d[1, 6] * np.conj(d[2, 6]), rtol=1e-14)
    assert np.allclose(pgram[6], np.conj(pgram[6].T), rtol=1e-14)


def test_pgram_frequencies_are_those_of_a_series_padded_to_nextn():
    assert [nextn(n) for n in (60, 101, 127, 200, 256)] == [60, 108, 128, 200, 256]
    freq = pgram_frequencies(127)
    assert freq.shape == (64,)
    assert freq[0] == 1 / 128 and freq[-1] == 0.5


def test_frequencies_are_assigned_to_the_five_bands_by_their_right_closed_edges():
    n_time = 256
    cycles = np.array([0, 0.5, 4, 4.5, 8, 9, 16, 17, 32, 33, 45, 46])
    assert band_of_frequency(cycles / n_time, n_time).tolist() == [0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0]
    band = band_of_frequency(pgram_frequencies(n_time), n_time)
    assert np.bincount(band)[1:].tolist() == [4, 4, 8, 16, 13]


def test_band_rows_lie_one_row_below_the_band_edges():
    rows = [r.tolist() for r in band_rows(256, 0.4)]
    assert rows == [list(range(a, b)) for a, b in ((1, 5), (5, 9), (9, 17), (17, 33), (33, 46))]
    assert [r.tolist() for r in band_rows(256, 50 / 256)] == rows
    assert band_rows(101, 0.4)[4].tolist() == [*range(35, 44), *range(55, 60)]
    assert band_rows(256, 40.5 / 256)[4].tolist() == [*range(33, 41), *range(129, 134)]


def test_channels_sharing_a_sinusoid_are_coherent_in_its_band_only():
    # 10 cycles per trial is row 11 of the periodogram: the alpha band (rows 9 to 16)
    coh = band_coherence(made_up_eeg(5, 256, 40, cycles=10, seed=3), 0.4)
    assert coh[2, 0, 1] > 0.99
    assert np.delete(coh[:, 0, 1], 2).max() < 0.2
    assert coh[2, :2, 2:].max() < 0.2
    assert coh[2, 2:, 2:][np.triu_indices(3, 1)].max() < 0.2


def test_the_coherence_is_symmetric_between_0_and_1_and_1_on_the_diagonal():
    coh = band_coherence(made_up_eeg(6, 200, 5, cycles=20, seed=4), 0.4)
    assert coh.shape == (5, 6, 6)
    assert np.all((coh >= 0) & (coh <= 1 + 1e-12))
    close(coh, coh.transpose(0, 2, 1))
    close(np.diagonal(coh, axis1=1, axis2=2), np.ones((5, 6)))


def test_one_trial_has_coherence_1_where_a_band_is_a_single_frequency():
    # a single outer product d d^H has |S_jk|^2 = S_jj S_kk
    spectmat = get_pgram(made_up_eeg(4, 64, 1, cycles=3, seed=5)[:, :, 0])[2:4]
    close(coherence_matrix(spectmat), np.ones((2, 4, 4)))


def test_band_spectrum_averages_the_band_rows_and_the_trials():
    eeg = made_up_eeg(3, 256, 4, cycles=10, seed=8)
    rows = band_rows(256, 0.4)[2]
    expected = np.mean([get_pgram(eeg[:, :, trial])[rows - 1].mean(axis=0) for trial in range(4)],
                       axis=0)
    close(band_spectrum(eeg, 0.4)[2].real, expected.real)
    close(band_spectrum(eeg, 0.4)[2].imag, expected.imag)


def test_the_497_candidate_snps_skip_the_header_line_and_two_rows(tmp_path):
    inputs = made_up_inputs(str(tmp_path))
    snp_names = read_snp_names(inputs["file"])
    assert snp_names == [name for j, name in enumerate(inputs["names"], start=1)
                         if j not in (1, 165, 358)]
    snps = candidate_snps(inputs["genotypes"], snp_names)
    assert snps.shape == (6, 497)
    assert np.array_equal(snps[:, 163], inputs["genotypes"][inputs["names"][165]])
    with pytest.raises(KeyError):
        candidate_snps(inputs["genotypes"], [inputs["names"][164]])


def test_the_file_of_a_triple_is_a_row_id_its_coherence_then_the_snps(tmp_path):
    step3, step4 = tmp_path / "step3", tmp_path / "step4"
    step3.mkdir()
    step4.mkdir()
    inputs = made_up_inputs(str(tmp_path))
    snps = candidate_snps(inputs["genotypes"], read_snp_names(inputs["file"]))
    for band in range(1, 6):
        write_channel_band_files(inputs["coherences"], snps, band, str(step3))
    assert len(os.listdir(step3)) == 20
    table = pd.read_csv(step3 / "dat_ch3b5.csv", index_col=0)
    assert table.shape == (6, 4 + 497)
    close(table["V2"], inputs["coherences"][:, 1, 2, 4])

    write_triple_files(made_up_triples(), 4, str(step3), str(step4))
    assert sorted(os.listdir(step4)) == ["dat_ch7.csv", "dat_ch9.csv"]
    with open(step4 / "dat_ch9.csv") as handle:
        text = handle.read()
    assert text.split("\n")[0] == ",".join(['""'] + [f'"V{j}"' for j in range(1, 499)])
    table = pd.read_csv(step4 / "dat_ch9.csv")
    assert table.shape == (6, 499)
    assert table.iloc[:, 0].tolist() == [1, 2, 3, 4, 5, 6]
    close(table.iloc[:, 1], inputs["coherences"][:, 1, 2, 4])
    close(table.iloc[:, 2:], snps)
    # the same file written from the table in memory
    direct = tmp_path / "direct.csv"
    write_table(triple_table(channel_band_table(inputs["coherences"], snps, 3, 5), 2, 4), direct)
    assert direct.read_text() == text


@pytest.mark.against_r
def test_python_agrees_with_r(tmp_path):
    work = str(tmp_path)
    # recordings: channels, time points, trials, cutoff in cycles per trial
    sizes = [(5, 256, 6, 102.4), (3, 200, 4, 50), (3, 127, 3, 50), (4, 101, 3, 40.4),
             (3, 256, 2, 40.5)]
    recordings = [made_up_eeg(*size[:3], cycles=10, seed=i) for i, size in enumerate(sizes)]
    for i, eeg in enumerate(recordings, start=1):
        eeg.ravel(order="F").tofile(os.path.join(work, f"eeg{i}.bin"))
    np.savetxt(os.path.join(work, "sizes.txt"), np.array(sizes), fmt="%.17g")
    inputs = made_up_inputs(work)
    inputs["coherences"].ravel(order="F").tofile(os.path.join(work, "coherences.bin"))
    inputs["genotypes"].to_numpy().ravel(order="F").tofile(os.path.join(work, "genotypes.bin"))
    with open(os.path.join(work, "genotyped.txt"), "w") as handle:
        handle.write("\n".join(inputs["genotypes"].columns) + "\n")
    made_up_triples().to_csv(os.path.join(work, "triple.csv"))
    script = os.path.join(work, "coherence.R")
    with open(script, "w") as handle:
        handle.write(R_SCRIPT)
    subprocess.run(["Rscript", "--vanilla", script, REPO_ROOT, work], check=True)

    with open(os.path.join(work, "coherence.bin"), "rb") as handle:
        for size, eeg in zip(sizes, recordings):
            n_chans, n_time = size[:2]
            cutoff = size[3] / n_time
            lengths = np.fromfile(handle, dtype=np.int32, count=5)
            r_rows = np.fromfile(handle, dtype=np.int32, count=lengths.sum())
            rows = band_rows(n_time, cutoff)
            assert [len(r) for r in rows] == lengths.tolist()
            assert np.concatenate(rows).tolist() == r_rows.tolist()
            r_coh = np.fromfile(handle, dtype=np.float64, count=5 * n_chans * n_chans)
            np.testing.assert_allclose(band_coherence(eeg, cutoff).ravel(order="F"), r_coh,
                                       rtol=1e-10, atol=1e-12)

    snps = candidate_snps(inputs["genotypes"], read_snp_names(inputs["file"]))
    for step in ("py3", "py4"):
        os.mkdir(os.path.join(work, step))
    for band in range(1, 6):
        write_channel_band_files(inputs["coherences"], snps, band, os.path.join(work, "py3"))
    write_triple_files(made_up_triples(), 4, os.path.join(work, "py3"), os.path.join(work, "py4"))
    for ours, theirs, count in (("py3", "r3", 20), ("py4", "r4", 2)):
        files = sorted(os.listdir(os.path.join(work, theirs)))
        assert len(files) == count and files == sorted(os.listdir(os.path.join(work, ours)))
        for name in files:
            with open(os.path.join(work, ours, name), "rb") as a, \
                    open(os.path.join(work, theirs, name), "rb") as b:
                assert a.read() == b.read(), name
