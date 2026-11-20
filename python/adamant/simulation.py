# -----------------------------------------------------------------------------
# File:    python/adamant/simulation.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Ridge Penalization in High
#          Dimensional Testing with Applications to Imaging Genetics.
#          Front. Neurosci. 16:836100. Tables 1 to 3, Figure 3.
# Purpose: Simulation study of AdaMantCV: simulated designs, one replication,
#          the driver that runs replications in parallel and stores their
#          p-values, and the rejection-rate table of a set of p-values.
# Inputs:  archive/run0204b/table<N>_pvalues.csv, fig3_lambdas.csv (unless
#          rerun).
# Outputs: results/<name>.csv (the table); in a rerun also
#          results/<name>_<settings>/scenarioSS_runRRRR.csv (one per
#          replication).
# Usage:   from adamant.simulation import simulation_study
#          simulation_study("table1_type1_error", "zero", opts, PAPER, PATHS)
# Source:  R/simulation.R.
# -----------------------------------------------------------------------------
import glob
import os
import re

import numpy as np
import pandas as pd

from . import rbase, rng
from .adamant_cv import adamant_cv
from .setup import (cached_csv, log_message, output_file, parallel_map, parse_index_set,
                    read_csv, settings_tag, write_csv)

# Simulation design of Tables 1 to 3 and Figure 3, passed as `design`.
SIMULATION_DESIGN = {
    "n": 350,                                               # sample size
    "rho": 0.025,                                           # correlation, compound symmetric
    "beta_mean": 1,                                         # value or mean of the coefficients
    "quantiles": {"none": 0, "Q1": 0.25, "Q2": 0.5},        # truncation tau
    "covariances": {"CS": "CSYM", "Heteroskedastic": "NCDIAG"},
    "dimensions": [500, 1000],                              # numbers of predictors p
}


def make_covariance(type, p, rho):
    """Covariance matrix of the simulated predictors.

    Parameters: type ("CSYM", compound symmetric, or "NCDIAG", diagonal with
    variances log(2), ..., log(p + 1)), p (dimension), rho (correlation of
    the compound symmetric matrix).
    Returns: p x p matrix.
    """
    if type == "CSYM":
        Sigma = np.full((p, p), float(rho), order="F")
        np.fill_diagonal(Sigma, 1.0)
        return Sigma
    if type == "NCDIAG":
        return np.asfortranarray(np.diag(np.log(np.arange(2, p + 2))))
    raise ValueError(f"unknown covariance type {type}")


def simulate_design(n, p, covariance, rho, beta_type, mu, sigma, seed):
    """Simulated standardized design and coefficients.

    Parameters: n, p (numbers of observations and predictors), covariance,
    rho (as for `make_covariance()`), beta_type ("zero" for the null,
    "fixed" for all equal to `mu`, or "random" for N(mu, sigma), drawn right
    after X from the same random stream), mu, sigma (mean and sd of the
    coefficients), seed (seed of the design).
    Returns: dict with `X` (n x p, standardized) and `beta` (length p).
    """
    rng.set_seed(seed)
    X = rbase.scale(rbase.mvrnorm(n, make_covariance(covariance, p, rho)))
    if beta_type == "zero":
        beta = np.zeros(p)
    elif beta_type == "fixed":
        beta = np.full(p, float(mu))
    elif beta_type == "random":
        beta = rng.rnorm(p, mu, sigma)
    else:
        raise ValueError(f"unknown beta type {beta_type}")
    return {"X": X, "beta": beta}


def simulation_scenarios(beta_fixed, design):
    """The 12 settings of a table in the row order of the paper.

    Parameters: beta_fixed (True for fixed coefficients, scenario numbers 20
    to 36 of the original runs; False for zero or random ones, 2 to 18),
    design (simulation design: quantiles, covariances, dimensions).
    Returns: DataFrame with columns scenario, covariance, p, tau.
    """
    taus, covariances = list(design["quantiles"]), list(design["covariances"])
    rows = [{"scenario": 18 * beta_fixed + 9 * (p == 1000) + 3 * taus.index(tau)
             + covariances.index(covariance) + 2,
             "covariance": covariance, "p": p, "tau": tau}
            for covariance in covariances for p in design["dimensions"] for tau in taus]
    return pd.DataFrame(rows)


def simulation_run(simulated, run, q, sigma, n_negative, n_positive, n_perms, seed,
                   upper=1000):
    """One replication: a new response from the design, then one AdaMantCV
    test.

    The permutation seed is the same in every replication, as in the
    original.

    Parameters: simulated (output of `simulate_design()`), run (replication
    number; the response is drawn with seed `10 * run + seed`), q,
    n_negative, n_positive, n_perms, seed, upper (as for `adamant_cv()`),
    sigma (sd of the noise).
    Returns: output of `adamant_cv()`.
    """
    rng.set_seed(10 * run + seed)
    X = simulated["X"]
    Y = rbase.matprod(X, simulated["beta"])[:, 0] + rng.rnorm(X.shape[0], 0, sigma)
    return adamant_cv(X, Y, q, n_negative, n_positive, n_perms, seed, upper)


def simulate_designs(scenarios, beta_type, settings, design):
    """Simulated designs of a set of scenarios, one per covariance and
    dimension.

    Parameters: scenarios (rows of `simulation_scenarios()`), beta_type (as
    for `simulate_design()`), settings (settings of the paper: seed, sigma),
    design (simulation design).
    Returns: dict of `simulate_design()` outputs keyed "<covariance> <p>".
    """
    keys = scenarios[["covariance", "p"]].drop_duplicates()
    return {f"{covariance} {p}": simulate_design(
        design["n"], p, design["covariances"][covariance], design["rho"], beta_type,
        design["beta_mean"], settings["sigma"], settings["seed"])
        for covariance, p in zip(keys["covariance"], keys["p"])}


def run_replications(run_dir, beta_type, scenarios, opts, settings, design):
    """Run the replications that are not yet stored and return all stored
    ones.

    Each replication seeds itself, so the result does not depend on the
    number of parallel processes or on the order of the replications.

    Parameters: run_dir (folder with one csv file per replication),
    beta_type (as for `simulate_design()`), scenarios (rows of
    `simulation_scenarios()` to run), opts (options: runs, perms, negative,
    positive, cores), settings (settings of the paper: seed, sigma, upper),
    design (simulation design).
    Returns: DataFrame with columns scenario, run and the four p-values
    p_gcv, p_loo, p_ga_gcv, p_ga_loo of every stored replication.
    """
    os.makedirs(run_dir, exist_ok=True)

    def run_file(scenario, run):
        return os.path.join(run_dir, "scenario%02d_run%04d.csv" % (scenario, run))

    simulated = simulate_designs(scenarios, beta_type, settings, design)
    settings_rows = scenarios.to_dict("records")
    jobs = [(run, setting) for setting in settings_rows
            for run in parse_index_set(opts["runs"], count=True)
            if not os.path.exists(run_file(setting["scenario"], run))]
    log_message("%s: %d replications to run", os.path.basename(run_dir), len(jobs))

    def replicate(job):
        run, setting = job
        ada = simulation_run(simulated[f"{setting['covariance']} {setting['p']}"], run,
                             design["quantiles"][setting["tau"]], settings["sigma"],
                             opts["negative"], opts["positive"], opts["perms"],
                             settings["seed"], settings["upper"])
        p = ada["p_values"]
        write_csv(pd.DataFrame({"scenario": [setting["scenario"]], "run": [run],
                                "p_gcv": [p["gcv"]], "p_loo": [p["loo"]],
                                "p_ga_gcv": [p["ga_gcv"]], "p_ga_loo": [p["ga_loo"]]}),
                  run_file(setting["scenario"], run))

    parallel_map(replicate, jobs, opts["cores"])
    return pd.concat([read_csv(file) for file in sorted(glob.glob(os.path.join(run_dir, "*")))],
                     ignore_index=True)


def rejection_rates(scenarios, pvalues, alpha):
    """Rejection rates of the four tests in every setting.

    A test rejects when its stored p-value is below `alpha`.

    Parameters: scenarios (rows of `simulation_scenarios()`), pvalues
    (DataFrame as returned by `run_replications()`), alpha (significance
    level).
    Returns: `scenarios` with the columns runs (number of replications) and
    the rates gcv, loocv, ga_gcv, ga_loocv.
    """
    columns = {"gcv": "p_gcv", "loocv": "p_loo", "ga_gcv": "p_ga_gcv", "ga_loocv": "p_ga_loo"}
    complete = pvalues.dropna(subset=list(columns.values()))
    rows = []
    for setting in scenarios.to_dict("records"):
        chosen = complete[complete["scenario"] == setting["scenario"]]
        if len(chosen) == 0:
            continue
        rates = {rate: int((chosen[column] < alpha).sum()) / len(chosen)
                 for rate, column in columns.items()}
        rows.append({**setting, "runs": len(chosen), **rates})
    return pd.DataFrame(rows)


def simulation_study(name, beta_type, opts, settings, paths, design=SIMULATION_DESIGN):
    """Build one of Tables 1 to 3.

    Rejection rates from the p-values of the archived run
    (archive/run0204b/table<N>_pvalues.csv) or, with `opts["rerun"]`, of
    replications run here. The table is written to results/<name>.csv.

    Parameters: name ("Table1_type1-error", "Table2_power_fixed-effects" or
    "Table3_power_random-effects"), beta_type (as for `simulate_design()`), opts
    (options: rerun, scenarios, "all" or an index text, and those of
    `run_replications()`), settings (settings of the paper), paths (project
    paths: archive, results), design (simulation design).
    Returns: the table.
    """
    scenarios = simulation_scenarios(beta_type == "fixed", design)
    if opts["scenarios"] != "all":
        scenarios = scenarios[scenarios["scenario"].isin(parse_index_set(opts["scenarios"]))]
    if opts["rerun"]:
        run_dir = os.path.join(paths["results"], f"{name}_{settings_tag(opts)}")
        pvalues = run_replications(run_dir, beta_type, scenarios, opts, settings, design)
    else:
        pvalues = read_csv(os.path.join(paths["archive"], "run0204b",
                                        re.sub("_.*", "_pvalues.csv", name).lower()))
    table = rejection_rates(scenarios, pvalues, settings["alpha"])
    write_csv(table, output_file(paths["results"], f"{name}.csv"))
    print(table.to_string(index=False))
    return table


def chosen_penalties(panels, replication, p, covariance, opts, settings, paths,
                     design=SIMULATION_DESIGN):
    """Penalties chosen by GCV and LOOCV in one replication (Figure 3).

    From the archived run (archive/run0204b/fig3_lambdas.csv) or, with
    `opts["rerun"]`, simulated and chosen here (stored in results/).

    Parameters: panels (DataFrame with columns panel, beta_type, tau),
    replication (replication number), p (number of predictors), covariance
    (covariance type, as for `make_covariance()`), opts (options: rerun,
    perms, negative, positive, cores), settings (settings of the paper),
    paths (project paths: archive, results), design (simulation design).
    Returns: DataFrame with columns panel, b (1 = observed data, then the
    permutations), GCV, LOOCV.
    """
    if not opts["rerun"]:
        return read_csv(os.path.join(paths["archive"], "run0204b", "fig3_lambdas.csv"))

    def compute():
        simulated = {beta_type: simulate_design(design["n"], p, covariance, design["rho"],
                                                beta_type, design["beta_mean"],
                                                settings["sigma"], settings["seed"])
                     for beta_type in ("zero", "fixed")}

        def choose(panel):
            ada = simulation_run(simulated[panel["beta_type"]], replication,
                                 design["quantiles"][panel["tau"]], settings["sigma"],
                                 opts["negative"], opts["positive"], opts["perms"],
                                 settings["seed"], settings["upper"])
            return pd.DataFrame({"panel": panel["panel"], "beta_type": panel["beta_type"],
                                 "tau": panel["tau"],
                                 "b": np.arange(1, len(ada["lambda_gcv"]) + 1),
                                 "GCV": ada["lambda_gcv"], "LOOCV": ada["lambda_loo"]})

        return pd.concat(parallel_map(choose, panels.to_dict("records"), opts["cores"]),
                         ignore_index=True)

    return cached_csv(output_file(paths["results"], f"fig3_lambdas_{settings_tag(opts)}.csv"),
                      compute)
