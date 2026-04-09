from pathlib import Path

import pandas as pd

root_path = Path(__file__).parent.absolute()

df = pd.read_csv(root_path / "David_screen_exp1_Batch3791_exp2_aggr_Fullset.csv")
df["plate.name"] = df["plate.name"].astype(int) + 1

df = pd.concat(
    [pd.read_csv(root_path / "David_screen_exp1_Batch3753_exp1_aggr_Fullset.csv"), df], ignore_index=True
)

nuclei_cols = [c for c in df.columns if "nc." in c]

df = df.drop(columns=nuclei_cols)

std_cols = [c for c in df.columns if "Standard deviation" in c]
dapi_cols = [c for c in df.columns if "DAPI." in c]
path_cols = [c for c in df.columns if "fpath" in c]
tritc_cols = [c for c in df.columns if "TRITC." in c]

dont_need = [
    "obj.Mean(perimeter).um.meas",
    "obj.Mean(area_BoundingBox).um2.meas",
    "obj.Mean(ratio_Area_BoundingBox_Area).um2.meas",
    "obj.Mean(major_Axis).um.meas",
    "obj.Mean(minor_Axis).um.meas",
    "obj.Standard deviation(area).um2.meas",
    "obj.Standard deviation(perimeter).um.meas",
    "obj.Standard deviation(area_BoundingBox).um2.meas",
    "obj.Standard deviation(major_Axis).um.meas",
    "obj.Standard deviation(minor_Axis).um.meas",
    "obj.Sum(area).um2.meas",
    "obj.Mean(area).meas",
    "obj.Mean(perimeter).meas",
    "obj.Mean(area_BoundingBox).meas",
    "obj.Mean(ratio_Area_BoundingBox_Area).meas",
    "obj.Mean(major_Axis).meas",
    "obj.Mean(minor_Axis).meas",
    "obj.Mean(axis_Ratio_Minor_Major).meas",
    "obj.Mean(angle_major_to_x_axis).meas",
    "obj.Mean(feret).meas",
    "obj.Mean(minFeret).meas",
    "obj.Mean(feretAngle).meas",
    "obj.Mean(roundness).meas",
    "obj.Mean(solidity).meas",
    "obj.Mean(circularity).meas",
    "obj.Mean(equivdiameter).meas",
    "obj.Mean(eccentricity).meas",
    "obj.Mean(number_of_end_points).meas",
    "obj.Mean(number_of_junction_points).meas",
    "obj.Mean(number_of_single_junction_points).meas",
    "obj.Mean(number_of_triple_points).meas",
    "obj.Mean(number_of_quadruple_points).meas",
    "obj.Mean(number_of_branches).meas",
    "obj.Mean(maximum_length_of_branches).meas",
    "obj.Mean(average_length_of_branches).meas",
    "obj.Mean(accumulated_length_of_branches).meas",
    "obj.Mean(accumulated_intensity).meas",
    "obj.Mean(mean_intensity).meas",
    "obj.Mean(std_intensity).meas",
    "obj.Mean(maximum_intensity).meas",
    "obj.Mean(minimum_intensity).meas",
]

internal_qc_cols = [
    "QC",
    "Fibres",
    "2D growth",
    "Bubbles",
    "Pipetting effects",
    "Replicate inconsistency",
    "Autofluorescence of compound",
    "Crowded clusters",
    "Collapsed gel",
    "Staining of background",
    "IA_problem",
    "Other",
]

df = df.drop(
    columns=[
        "proj.code",
        "proj.name",
        "exp.code",
        "exp.name",
        "batch.code",
        "batch.name",
        "owner",
        "root.name",
        "site.name",
        "time.point",
        "plate.layout.info.Project",
        "image.fname.DAPI",
        "image.fname.TRITC",
        "plate.layout.info.Pipeting location",
        "plate.layout.info.Model",
        "plate.layout.info.Tuning set",
        "plate.layout.info.Step-size",
        "plate.layout.info.Exposure duration",
        "plate.layout.info.Culture duration",
        *std_cols,
        *dapi_cols,
        *dont_need,
        *path_cols,
        *tritc_cols,
        *internal_qc_cols,
    ]
)

strip_columns = {c: c.replace("plate.layout.info.", "") for c in df.columns if "plate.layout.info." in c}

df = (
    df.rename(columns=strip_columns)
    .rename(columns={"Type": "treatment_type", "plate.name": "plate_id"})
    .drop(columns=["Stimulation solvent", "Treatment solvent", "Treatment 2 solvent"])  # it's all DMSO
    .query("plate_id != 4")  # remove plate 4, not our compounds!
)

df.to_csv("ADPKD-TargetValidationScreen_Batch3791_and_Batch3753.csv", index=False)
