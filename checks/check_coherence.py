# -----------------------------------------------------------------------------
# File:    checks/check_coherence.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Rebuild steps 3 and 4 of the coherence pipeline with
#          python/adamant/coherence.py from the private data and compare the
#          files of the triples with those the paper used.
# Inputs:  <data>/01_raw/cups_top_500_SNPs.txt,
#          <data>/03_triples_channels/triple.csv,
#          <data>/04_coherence_SNPs/dat_ch<k>.csv; from the work folder
#          genotypes.bin, genotyped.txt and coherences.bin, written by
#          checks/check_coherence.R from the two .RData files.
# Outputs: In the work folder: python/step4/dat_ch<k>.csv, python/step3.md5
#          (checksums of the files of step 3, which are deleted band by
#          band) and python/summary.csv. Prints counts and the largest
#          difference, never data.
# Usage:   python checks/check_coherence.py <data folder> <work folder>
#          (called by checks/check_coherence.sh from the repository root)
# -----------------------------------------------------------------------------
import hashlib
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))

import numpy as np
import pandas as pd

from adamant.coherence import (candidate_snps, channel_band_table, read_snp_names, triple_table,
                               write_channel_band_files, write_triple_files)
from adamant.setup import read_csv

N_SUBJECTS, N_CHANS, N_BANDS = 350, 64, 5


def compare_triple(k, triple, coherences, snps, data_dir, step4):
    """Compare the rebuilt file of one triple with the file of the paper.

    Parameters: k (triple number), triple (its row of triple.csv),
    coherences, snps (the inputs of step 3), data_dir (data folder), step4
    (folder of the rebuilt files).
    Returns: dict: whether the two files have the same bytes, whether the
    file of the paper has the bytes of the rebuilt numbers written correctly
    rounded to 15 significant digits, the number of fields that then
    differ, and the largest absolute difference of the numbers.
    """
    name = f"dat_ch{k}.csv"
    with open(os.path.join(data_dir, "04_coherence_SNPs", name)) as handle:
        paper = handle.read()
    with open(os.path.join(step4, name)) as handle:
        written = handle.read()
    rebuilt = triple_table(channel_band_table(coherences, snps, triple["chan2"], triple["band"]),
                           triple["chan1"], N_CHANS)
    lines = paper.splitlines()
    fields = np.array([line.split(",") for line in lines[1:]])
    rounded = np.array([["%.15g" % value for value in row] for row in rebuilt.tolist()])
    assert fields[:, 1:].shape == rounded.shape and lines[0] == written.splitlines()[0]
    assert fields[:, 0].tolist() == [f'"{i}"' for i in range(1, N_SUBJECTS + 1)]
    return {"triple": k, "same_bytes": int(paper == written),
            "same_bytes_rounded": int(np.array_equal(fields[:, 1:], rounded)),
            "fields_differing": int((fields[:, 1:] != rounded).sum()),
            "max_difference": np.abs(fields[:, 1:].astype(np.float64) - rebuilt).max()}


def main(data_dir, work):
    """Rebuild the files, compare them, and print the summary line.

    Parameters: data_dir (data folder), work (work folder).
    Returns: None.
    """
    with open(os.path.join(work, "genotyped.txt")) as handle:
        genotyped = handle.read().split()
    genotypes = np.fromfile(os.path.join(work, "genotypes.bin")).reshape(
        (N_SUBJECTS, len(genotyped)), order="F")
    snps = candidate_snps(dict(zip(genotyped, genotypes.T)),
                          read_snp_names(os.path.join(data_dir, "01_raw", "cups_top_500_SNPs.txt")))
    coherences = np.fromfile(os.path.join(work, "coherences.bin")).reshape(
        (N_SUBJECTS, N_CHANS, N_CHANS, N_BANDS), order="F")
    triples = read_csv(os.path.join(data_dir, "03_triples_channels", "triple.csv"), index_col=0)

    step3 = os.path.join(work, "python", "step3")
    step4 = os.path.join(work, "python", "step4")
    os.makedirs(step3)
    os.makedirs(step4)
    checksums = []
    for band in range(1, N_BANDS + 1):
        write_channel_band_files(coherences, snps, band, step3)
        for name in sorted(os.listdir(step3)):
            with open(os.path.join(step3, name), "rb") as handle:
                checksums.append(f"{hashlib.md5(handle.read()).hexdigest()} {name}")
        write_triple_files(triples[triples["band"] == band], N_CHANS, step3, step4)
        for name in os.listdir(step3):
            os.remove(os.path.join(step3, name))
    with open(os.path.join(work, "python", "step3.md5"), "w") as handle:
        handle.write("\n".join(checksums) + "\n")

    summary = pd.DataFrame([compare_triple(k, triple, coherences, snps, data_dir, step4)
                            for k, triple in triples.iterrows()])
    summary.to_csv(os.path.join(work, "python", "summary.csv"), index=False)
    print(f"Python: {len(checksums)} files of step 3, {len(triples)} of step 4; files with the "
          f"bytes of the paper: {summary['same_bytes'].sum()} as written, "
          f"{summary['same_bytes_rounded'].sum()} written correctly rounded "
          f"({summary['fields_differing'].sum()} fields differ); "
          f"largest difference {summary['max_difference'].max():.3g}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
