"""Audit the plate-wise normalization used for the cyst-swelling figures.

The reviewer correctly noted that applying different normalizations to separate
plates can change a Mann–Whitney test if values from those plates are pooled.
The manuscript analyses do not pool plates: each reported comparison is made
within one physical plate. This script makes that design auditable by:

1. reproducing the median-control normalization for both follow-up experiments;
2. exporting each well's raw and normalized cyst-area value with the plate
   control medians used for normalization; and
3. verifying for every within-plate pair of experimental conditions that the
   two-sided Mann–Whitney p-value is identical on raw and normalized values.

Run from the repository root:
    python 03_Target_Validation/normalization_audit.py

Outputs:
    data/normalization_audit/per_well_raw_and_normalized.csv
    data/normalization_audit/within_plate_mann_whitney_audit.csv
"""

from __future__ import annotations

from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = DATA / "normalization_audit"
READOUT = "obj.Mean(area).um2.meas"
NORMALIZED = "normalized_cyst_swelling_pct"
EXPECTED_WELL_ROWS = 1_265
EXPECTED_PAIRWISE_TESTS = 4_485


def load_target_validation() -> pd.DataFrame:
    path = DATA / "target_validation" / "ADPKD-TargetValidationScreen_Batch3791_and_Batch3753.csv"
    df = pd.read_csv(path)
    # Match the analysis notebook's quality-control filtering.
    df = df.assign(QC_problems=df["QC_problems"].fillna("-"))
    df = df.query("QC_problems == '-'").dropna(subset=[READOUT]).copy()
    df["experiment"] = "target_validation"
    df["batch"] = "Batch3791_and_Batch3753"
    return df


def load_compound_exploration() -> pd.DataFrame:
    path = DATA / "compound_exploration" / "ADPKD-CpdExplorationScreen_Batch4042_and_Batch4064.csv"
    df = pd.read_csv(path)
    plate_mapping = {"4042_1": 1, "4042_2": 2, "4064_1": 3}
    strip_prefix = {
        column: column.replace("plate.layout.info.", "")
        for column in df.columns
        if column.startswith("plate.layout.info.")
    }
    df = df.rename(columns=strip_prefix).rename(columns={"Type": "treatment_type"})
    df["plate_id"] = df["ID"].map(plate_mapping)
    if df["plate_id"].isna().any():
        unknown = sorted(df.loc[df["plate_id"].isna(), "ID"].dropna().unique())
        raise ValueError(f"Unmapped exploration plate IDs: {unknown}")
    df["plate_id"] = df["plate_id"].astype(int)
    df["experiment"] = "compound_exploration"
    df["batch"] = df["batch.code"].astype(str)
    return df.dropna(subset=[READOUT]).copy()


def normalize_within_plate(df: pd.DataFrame) -> pd.DataFrame:
    """Apply the normalization used in the analysis notebooks.

    DMSO median maps to 0% and the 2.5 µM FSK median maps to 100%:
        100 * (raw - DMSO median) / (FSK median - DMSO median)
    """
    pieces = []
    for (experiment, plate_id), plate in df.groupby(["experiment", "plate_id"], sort=True):
        dmso = plate.loc[
            (plate["treatment_type"] == "solvent_ctrl") & (plate["Stimulant dose"] == 0.1), READOUT
        ].median()
        fsk = plate.loc[
            (plate["treatment_type"] == "stim_only") & (plate["Stimulant dose"] == 2.5), READOUT
        ].median()
        if not np.isfinite(dmso) or not np.isfinite(fsk) or fsk <= dmso:
            raise ValueError(
                f"Invalid control window for {experiment}, plate {plate_id}: DMSO={dmso}, FSK={fsk}"
            )
        plate = plate.copy()
        plate["dmso_median_raw"] = dmso
        plate["fsk_2_5uM_median_raw"] = fsk
        plate[NORMALIZED] = 100.0 * (plate[READOUT] - dmso) / (fsk - dmso)
        pieces.append(plate)
    return pd.concat(pieces, ignore_index=True)


def condition_label(row: pd.Series) -> str:
    """Create an exact condition label for grouping replicate wells."""
    fields = [
        ("type", row.get("treatment_type")),
        ("stimulant", row.get("Stimulant")),
        ("stimulant_dose", row.get("Stimulant dose")),
        ("treatment", row.get("Treatment")),
        ("treatment_dose", row.get("Treatment concentration")),
        ("treatment_2", row.get("Treatment 2")),
        ("treatment_2_dose", row.get("Treatment 2 concentration")),
    ]
    return ";".join(f"{name}={value}" for name, value in fields if pd.notna(value))


def audit_pairwise_tests(df: pd.DataFrame) -> pd.DataFrame:
    """Compare raw and normalized p-values for all condition pairs per plate."""
    records = []
    testable = df[df["treatment_type"].isin(["solvent_ctrl", "stim_only", "treatment"])].copy()
    testable["condition"] = testable.apply(condition_label, axis=1)

    for (experiment, plate_id), plate in testable.groupby(["experiment", "plate_id"], sort=True):
        grouped = {name: group for name, group in plate.groupby("condition", sort=True) if len(group) >= 2}
        for (condition_1, group_1), (condition_2, group_2) in combinations(grouped.items(), 2):
            p_raw = mannwhitneyu(group_1[READOUT], group_2[READOUT], alternative="two-sided").pvalue
            p_normalized = mannwhitneyu(
                group_1[NORMALIZED], group_2[NORMALIZED], alternative="two-sided"
            ).pvalue
            records.append(
                {
                    "experiment": experiment,
                    "plate_id": plate_id,
                    "condition_1": condition_1,
                    "n_1": len(group_1),
                    "condition_2": condition_2,
                    "n_2": len(group_2),
                    "p_raw": p_raw,
                    "p_normalized": p_normalized,
                    "absolute_difference": abs(p_raw - p_normalized),
                    "identical_within_tolerance": bool(np.isclose(p_raw, p_normalized, rtol=0, atol=1e-15)),
                }
            )

    audit = pd.DataFrame(records)
    if audit.empty or not audit["identical_within_tolerance"].all():
        failures = audit.loc[~audit["identical_within_tolerance"]] if not audit.empty else audit
        raise AssertionError(f"Raw/normalized Mann–Whitney audit failed:\n{failures}")
    return audit


def export_columns(df: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "experiment",
        "batch",
        "plate_id",
        "well.name",
        "treatment_type",
        "Stimulant",
        "Stimulant dose",
        "Stimulant unit",
        "Treatment",
        "Treatment concentration",
        "Treatment unit",
        "Treatment 2",
        "Treatment 2 concentration",
        "Treatment 2 unit",
        READOUT,
        "dmso_median_raw",
        "fsk_2_5uM_median_raw",
        NORMALIZED,
    ]
    return df.reindex(columns=columns).sort_values(["experiment", "plate_id", "well.name"])


def main() -> None:
    combined = pd.concat([load_target_validation(), load_compound_exploration()], ignore_index=True)
    normalized = normalize_within_plate(combined)
    audit = audit_pairwise_tests(normalized)

    if len(normalized) != EXPECTED_WELL_ROWS:
        raise AssertionError(
            f"Expected {EXPECTED_WELL_ROWS:,} QC-filtered wells, found {len(normalized):,}"
        )
    if len(audit) != EXPECTED_PAIRWISE_TESTS:
        raise AssertionError(
            f"Expected {EXPECTED_PAIRWISE_TESTS:,} within-plate comparisons, found {len(audit):,}"
        )

    OUT.mkdir(parents=True, exist_ok=True)
    values_path = OUT / "per_well_raw_and_normalized.csv"
    audit_path = OUT / "within_plate_mann_whitney_audit.csv"
    export_columns(normalized).to_csv(values_path, index=False)
    audit.to_csv(audit_path, index=False)

    print(f"Exported {len(normalized):,} well-level rows to {values_path.relative_to(ROOT)}")
    print(f"Audited {len(audit):,} within-plate condition pairs")
    print(f"Maximum |p_raw - p_normalized|: {audit['absolute_difference'].max():.3g}")
    print(f"Exported audit to {audit_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
