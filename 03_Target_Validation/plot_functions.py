import re
from itertools import product
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statannotations.Annotator import Annotator

ibm_colors = [
    "#648FFF",
    "#785EF0",
    "#DC267F",
    "#FE6100",
    "#FFB000",
]

rgb_colors = [mcolors.to_rgba(c) for c in ibm_colors]
plot_colors = [rgb_colors[0], rgb_colors[2], rgb_colors[3], rgb_colors[4]]


def add_treatment_2(row):
    if pd.notna(row["Treatment_2"]):
        if row["treatment_type"] == "treatment" and row["Treatment"] == "DMSO":
            # In this case, we have like DMSO 0.1% + CPA 1µM
            return row["Treatment_2"] + " " + str(row["Treatment_2_concentration"]) + row["Treatment_2_unit"]
        elif row["Treatment_2"] == "DMSO":
            return (
                row["treat_string"]
                + " + "
                + row["Stimulant"]
                + " "
                + str(row["Stimulant dose"])
                + row["Stimulant unit"]
            )
        else:
            return (
                row["treat_string"]
                + " + "
                + row["Treatment_2"]
                + " "
                + str(row["Treatment_2_concentration"])
                + row["Treatment_2_unit"]
            )
    else:
        return row["treat_string"]


def change_control_row(row):
    if row["treatment_type"] == "solvent_ctrl":
        return "DMSO 0.2%"  # two doses of DMSO 0.1%
    elif row["treatment_type"] == "stim_only":
        return f"FSK {row['Stimulant dose']}µM"
    else:
        return row["treat_string"]


def clean_treatment_string(text, single_treatment: bool = False):
    if single_treatment:
        return re.sub(r"\s\d\.\d+µM", "", text)
    else:
        return re.sub(r"([A-Za-z0-9]+)\s\d\.\d+µM(?=.+)", r"\1", text)


def plot_with_curves(
    df,
    plate_id,
    toplot_compound,
    stim_dose=2.5,
    single_treatment: bool = False,
    plot_normalized=True,
    target=None,
    savefig=True,
    figs_root=None,
):
    # Some variables...
    spheresize_col = f"obj.Mean(area).um2.meas{'_norm' if plot_normalized else ''}"
    rm_fsk_pattern = re.compile(r"\s\+\sFSK \d\.\d+µM")
    target_str = f"-{target}" if target is not None else ""

    # prepare the data by cleaning the treatment strings as we want  them for the plot
    df = df.query(f"plate_id == {plate_id}")
    df = (
        df.replace({"uM": "µM", "pct": "%"})
        .rename(
            columns={
                "Treatment 2": "Treatment_2",
                "Treatment 2 concentration": "Treatment_2_concentration",
                "Treatment 2 unit": "Treatment_2_unit",
            }
        )
        .assign(
            treat_string=lambda x: x["Treatment"]
            + " "
            + x["Treatment concentration"].astype(str)
            + x["Treatment unit"]
        )
    )
    df["treat_string"] = df.apply(add_treatment_2, axis=1)
    df["treat_string"] = df.apply(change_control_row, axis=1)

    # Colors for the treatment groups
    colors = ["#648FFF", "#DC267F", "#FFB000"]
    rgb_colors = [mcolors.to_rgba(c) for c in colors]

    # Define subsets for the control groups
    unstimulated = df.query("treat_string == 'DMSO 0.2%'").assign(whatever=1)
    stimulated = df.query("treat_string == 'FSK 2.5µM'").assign(whatever=1)
    half_stimulated = df.query("treat_string == 'FSK 0.79µM'").assign(whatever=1)

    controls = [unstimulated, half_stimulated, stimulated]
    control_colors = ["silver", "rosybrown", "maroon"]

    fig, axs = plt.subplots(ncols=2, figsize=(5, 3), sharey=True, gridspec_kw={"width_ratios": [1, 5]})
    control_ax = axs[0]
    for control, color in zip(controls, control_colors):
        sns.boxplot(
            data=control,
            x="whatever",
            y=spheresize_col,
            hue="treat_string",
            ax=control_ax,
            palette=[mcolors.to_rgb(color)],
        )

    ax = axs[1]  # where the treatment groups will be plotted
    seaborn_graphic = sns.pointplot(
        data=(
            df.rename(columns={"Stimulant dose": "Stimulant_dose"})
            .query("treat_string.str.contains(@toplot_compound)")
            .query(f"Stimulant_dose == {stim_dose}")
            .sort_values(["Treatment", "Treatment concentration"])
            .assign(logC=lambda x: np.log10(x["Treatment concentration"]))
            .assign(
                hue_string=lambda x: x["treat_string"].apply(
                    clean_treatment_string, single_treatment=single_treatment
                )
            )
            .assign(
                hue_string=lambda x: x["hue_string"]
                .apply(lambda s: rm_fsk_pattern.sub("", s))
                .str.replace("Aldosterone", "ALD")
            )
        ),
        x="logC",
        y=spheresize_col,
        hue="hue_string",
        palette=rgb_colors,
        native_scale=True,
        ax=ax,
        zorder=10,
    )
    seaborn_graphic.legend().set_visible(False)

    ax = plt.gca()

    # make a dotted line for the median of each control group
    median_unstimulated = unstimulated[spheresize_col].median()
    median_stimulated = stimulated[spheresize_col].median()
    median_half_stimulated = half_stimulated[spheresize_col].median()
    ax.axhline(median_unstimulated, color="silver", linestyle="--")  # label="median(DMSO 0.2%)")
    ax.axhline(median_half_stimulated, color="rosybrown", linestyle="--")  # label="median(FSK 0.79µM)")
    ax.axhline(median_stimulated, color="maroon", linestyle="--")  # label="median(FSK 2.5µM)")

    # make the legends with all the labels
    handles1, labels1 = control_ax.get_legend_handles_labels()
    handles2, labels2 = ax.get_legend_handles_labels()  # This will now include the horizontal lines
    all_handles = handles2 + handles1
    all_labels = labels2 + labels1
    if control_ax.get_legend() is not None:
        control_ax.get_legend().remove()
    if ax.get_legend() is not None:
        ax.get_legend().remove()

    ax.legend(
        all_handles, all_labels, bbox_to_anchor=(1.04, 0), loc="lower left", borderaxespad=0, frameon=False
    )

    ax.set_xlabel(r"log$_{10}$" + f"(Concentration) {toplot_compound} [µM]")
    ax.set_ylabel(r"Mean-Aggregated Cyst Size ($\mu$m$^2$)")
    plt.xticks([-3.0, -1.0, 0.0])
    ax.set_title(f"Plate {plate_id} - {toplot_compound} + FSK {stim_dose}µM")
    ax.spines[["right", "top"]].set_visible(False)

    control_ax.legend().set_visible(False)
    control_ax.set_ylabel(r"Mean-Aggregated Cyst Size ($\mu$m$^2$)")
    control_ax.xaxis.set_ticklabels([])
    control_ax.set_xlabel("Control Groups")
    control_ax.spines[["right", "top"]].set_visible(False)

    if savefig:
        if figs_root is None:
            figs_root = Path(__file__).parents[1] / "figures"
        fig.savefig(
            figs_root
            / f"lineplot-plate{plate_id}{target_str}_{toplot_compound}_FSK{str(stim_dose).replace('.', '-')}.png",
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
        )
    return seaborn_graphic


def box_mann_whitney_u(
    df,
    compounds,
    plate_number,
    double_treat_cpd,
    stimulant_dose,
    savefig=True,
    plot_normalized=True,
    target=None,
    figs_root=None,
):
    """Method for creating the boxplots used in the paper. If you're actually reading this code,
    I'm so sorry for the mess 🥲... The data transformation is a bit weird to get the desired plot.

    But you can do it!! I believe in you!"""
    target_str = f"-{target}" if target is not None else ""
    readout = "obj.Mean(area).um2.meas_norm" if plot_normalized else "obj.Mean(area).um2.meas"
    subset = df.query(f"plate_id == {plate_number}")
    other_dose_idxs = subset[
        (subset["Stimulant dose"] != stimulant_dose)
        & (subset["treatment_type"].isin(["stim_only", "treatment"]))
    ].index.tolist()
    subset = subset.drop(index=other_dose_idxs)

    idxs = subset.loc[
        (subset["treatment_type"].isin(["solvent_ctrl", "stim_only"])) & (subset["Stimulant"] == "FSK")
    ].index.tolist()
    subset.loc[idxs, "Treatment"] = subset.loc[idxs, "Stimulant"]
    subset.loc[idxs, "Treatment concentration"] = stimulant_dose
    subset.loc[idxs, "Treatment unit"] = "uM"

    subset = subset.assign(
        treat_string=lambda x: x["Treatment"]
        + x["Treatment concentration"].astype(str)
        + " "
        + x["Treatment 2"]
        + x["Treatment 2 concentration"].astype(str)
    ).assign(
        treat_string=lambda x: (
            x["treat_string"]
            .str.replace("DMSO0.1 DMSO0.1", "DMSO 0.2%")  # two doses of DMSO 0.1%
            .str.replace(f"FSK{stimulant_dose}", f"FSK {stimulant_dose}µM")
            .str.replace("DMSO0.1", "")
            .str.replace("0.001", " 0.001µM ")
            .str.replace("1.0", " 1µM ")
            .str.replace("0.1 ", " 0.1µM ")
            .str.replace(f"{double_treat_cpd}0.1", f"{double_treat_cpd} 0.1µM")
            .str.strip()
        ).apply(
            lambda s: " ".join(s.split())  # remove double spaces
        )
    )

    single_treatment = subset[
        (subset["Treatment 2"] == f"{double_treat_cpd}") & (subset["Treatment"] == "DMSO")
    ].copy()
    single_treatment["Treatment"] = f"{double_treat_cpd}"

    for cpd in compounds:
        subset_controls = pd.concat(
            [
                subset.query("treatment_type.isin(['solvent_ctrl'])"),
                subset.query("treatment_type.isin(['stim_only'])"),
                single_treatment.sort_values(["Treatment concentration", "Treatment 2 concentration"]),
            ],
            ignore_index=True,
        )
        subset_controls["Treatment 2"] = "Control"

        subset_cpd = subset.query(f"Treatment == '{cpd}'")
        subset_cpd["Treatment 2 Sort Key"] = subset_cpd["Treatment 2"].apply(
            lambda x: 0 if x == "DMSO" else 1
        )
        subset_cpd = subset_cpd.sort_values(
            ["Treatment 2 Sort Key", "Treatment 2 concentration", "Treatment concentration"]
        ).drop(columns=["Treatment 2 Sort Key"])

        fig, ax = plt.subplots(figsize=(6, 4))
        plot_data = pd.concat(
            [
                subset_controls,
                subset_cpd,
            ]
        )
        plot_data["Treatment 2"] = plot_data["Treatment 2"].apply(
            lambda x: f"+ FSK {stimulant_dose}µM" if x == "DMSO" else x
        )
        plot_data["Treatment 2"] = plot_data["Treatment 2"].apply(
            lambda x: (
                f"+ {x} + FSK {stimulant_dose}µM" if x not in ["Control", f"+ FSK {stimulant_dose}µM"] else x
            )
        )

        def update_treatment_2(row, double_treat_cpd):
            """Add the concentration of the double treatment!! This is for the legend only"""
            return row["Treatment 2"].replace(
                double_treat_cpd, f"{double_treat_cpd} {row['Treatment 2 concentration']}µM"
            )

        doubletreat_search_pattern = re.compile(f"(?=.{cpd})|(?=.{double_treat_cpd})")
        mask = plot_data.treat_string.str.contains(doubletreat_search_pattern, regex=True)
        plot_data.loc[mask, "Treatment 2"] = plot_data.loc[mask].apply(
            update_treatment_2, double_treat_cpd=double_treat_cpd, axis=1
        )
        sns.boxplot(
            data=plot_data,
            x="treat_string",
            y=readout,
        )

        labels = plot_data[["treat_string", "Treatment 2"]].drop_duplicates()["Treatment 2"]
        uniq_labels = (
            labels.drop_duplicates()
            .str.replace("Aldosterone", "ALD")
            .str.replace(" + FSK ", "\n+ FSK ")
            .tolist()
        )
        uniq_colors = plot_colors
        # 2 controls: DMSO, FSK
        # 6 single treatments: 3 concentrations of cpd, 3 concentrations of double_treat_cpd
        # 3 double treatments: 3 concentrations of cpd + double_treat_cpd 0.1µM
        # 3 double treatments: 3 concentrations of cpd + double_treat_cpd 1µM
        colors = [uniq_colors[0]] * 2 + [uniq_colors[1]] * 6 + [uniq_colors[2]] * 3 + [uniq_colors[3]] * 3

        legend_handles = [  # handcrafted legend. Horrible, I know
            mpatches.Patch(facecolor=color, label=label, edgecolor="black", linewidth=0.5)
            for label, color in zip(uniq_labels, uniq_colors)
        ]

        bars = ax.patches
        for bar, color in zip(bars, colors):
            bar.set_facecolor(color)

        ax.legend(
            handles=legend_handles,
            bbox_to_anchor=(1.04, 0),
            loc="lower left",
            borderaxespad=0,
            frameon=False,
        )

        # Prepare pairs for statistical annotation
        control_condition = f"FSK {stimulant_dose}µM"
        concentrations = ["0.001µM", "0.1µM", "1µM"]
        treatment_conditions = [f"{a} {b}" for a, b in product([cpd, double_treat_cpd], concentrations)]
        double_treat_comparisons = (
            # DTC (double treatment compound) v.s. PC (plate compound) + DTC
            [(f"{double_treat_cpd} 1µM", f"{cpd} {c} {double_treat_cpd} 1µM") for c in concentrations]
            + [(f"{double_treat_cpd} 0.1µM", f"{cpd} {c} {double_treat_cpd} 0.1µM") for c in concentrations]
            # compound alone v.s. double treatment
            + [(f"{cpd} 1µM", f"{cpd} 1µM {double_treat_cpd} {c}") for c in concentrations[1:]]
            + [(f"{cpd} 0.1µM", f"{cpd} 0.1µM {double_treat_cpd} {c}") for c in concentrations[1:]]
            + [(f"{cpd} 0.001µM", f"{cpd} 0.1µM {double_treat_cpd} {c}") for c in concentrations[1:]]
            # double treatment v.s. double treatment
            + [
                (f"{cpd} 1µM {double_treat_cpd} 1µM", f"{cpd} 1µM {double_treat_cpd} 0.1µM"),
                (f"{cpd} 0.1µM {double_treat_cpd} 1µM", f"{cpd} 0.1µM {double_treat_cpd} 0.1µM"),
                (f"{cpd} 0.001µM {double_treat_cpd} 1µM", f"{cpd} 0.1µM {double_treat_cpd} 0.1µM"),
            ]
        )
        pairs = [(control_condition, cond) for cond in treatment_conditions] + double_treat_comparisons
        # Add statistical annotations
        annotator = Annotator(ax, pairs, data=plot_data, x="treat_string", y=readout)
        annotator.configure(
            test="Mann-Whitney",
            text_format="star",
            loc="inside",
            hide_non_significant=True,
        )
        annotator.apply_and_annotate()

        labels = [tick.get_text() for tick in ax.get_xticklabels()]
        # get rid of the + doublee_treat_cpd <concentration> part of the string (already in the legend)
        new_labels = [
            re.sub(double_treat_cpd + r"\s\d\.?\d*\s?µM", "", x) if x not in labels[:6] else x for x in labels
        ]
        ax.set_xticklabels(new_labels, rotation=45, ha="right", rotation_mode="anchor", fontsize=10)

        ax.set_xlabel("")
        ax.set_ylabel(r"Mean-Aggregated Cyst Size ($\mu$m$^2$)")
        ax.set_title(cpd)
        # ax.set_title(f"Plate {plate_number} - {cpd}")

        ax.grid(axis="y", alpha=0.5)
        ax.spines[["right", "top"]].set_visible(False)
        ax.set_axisbelow(True)

        if savefig:
            if figs_root is None:
                figs_root = Path(__file__).parents[1] / "figures"
            fig.savefig(
                figs_root
                / (
                    f"boxplot-MannWhitneyU-plate{plate_number}{target_str}_{cpd}"
                    f"_FSK{str(stimulant_dose).replace('.', '-')}_{double_treat_cpd}.png"
                ),
                dpi=300,
                bbox_inches="tight",
                facecolor="white",
            )
            fig.savefig(
                figs_root
                / (
                    f"boxplot-MannWhitneyU-plate{plate_number}{target_str}_{cpd}"
                    f"_FSK{str(stimulant_dose).replace('.', '-')}_{double_treat_cpd}.svg"
                ),
                bbox_inches="tight",
                facecolor="none",
            )

        plt.show()


def box_mann_whitney_u_no_dt(
    df,
    compounds,
    plate_number,
    stimulant_dose,
    savefig=True,
    plot_normalized=True,
    target=None,
    figs_root=None,
):
    """Similar to box_mann_whitney_u but without double treatment comparisons"""
    target_str = f"-{target}" if target is not None else ""
    readout = "obj.Mean(area).um2.meas_norm" if plot_normalized else "obj.Mean(area).um2.meas"
    subset = df.query(f"plate_id == {plate_number}")
    other_dose_idxs = subset[
        (subset["Stimulant dose"] != stimulant_dose)
        & (subset["treatment_type"].isin(["stim_only", "treatment"]))
    ].index.tolist()
    subset = subset.drop(index=other_dose_idxs)

    idxs = subset.loc[
        (subset["treatment_type"].isin(["solvent_ctrl", "stim_only"])) & (subset["Stimulant"] == "FSK")
    ].index.tolist()
    subset.loc[idxs, "Treatment"] = subset.loc[idxs, "Stimulant"]
    subset.loc[idxs, "Treatment concentration"] = stimulant_dose
    subset.loc[idxs, "Treatment unit"] = "uM"

    subset = subset.assign(
        treat_string=lambda x: (
            x["Treatment"]
            + " "
            + x["Treatment concentration"].astype(str)
            + x["Treatment unit"].str.replace("uM", "µM").str.replace("pct", "%")
        ).str.replace("DMSO 0.1%", "DMSO 0.2%")
    )

    for cpd in compounds:
        subset_controls = pd.concat(
            [
                subset.query("treatment_type.isin(['solvent_ctrl'])"),
                subset.query("treatment_type.isin(['stim_only'])"),
            ],
            ignore_index=True,
        )
        subset_controls["Treatment 2"] = "Control"

        subset_cpd = (
            subset.query("treatment_type.isin(['treatment'])")
            .sort_values(["Treatment concentration"])
            .query(f"Treatment == '{cpd}'")
        )

        fig, ax = plt.subplots(figsize=(3, 4))
        plot_data = pd.concat(
            [
                subset_controls,
                subset_cpd,
            ]
        )

        def remove_single_treatment_C(row, rm_pattern):
            if row["treatment_type"] == "treatment":
                return rm_pattern.sub("", row["treat_string"])
            else:
                return "Control"

        remove_pattern = re.compile(r"\s\d\.\d+µM$")
        plot_data["Treatment 2"] = plot_data.apply(
            remove_single_treatment_C, rm_pattern=remove_pattern, axis=1
        )

        sns.boxplot(
            data=plot_data,
            x="treat_string",
            y=readout,
        )

        labels = (
            plot_data[["treat_string", "Treatment 2"]]
            .drop_duplicates()["Treatment 2"]
            .str.replace(cpd, f"+ FSK {stimulant_dose}µM")
            .str.replace("Aldosterone", "ALD")
        )
        uniq_labels = labels.drop_duplicates().tolist()
        uniq_colors = plot_colors
        colors = [uniq_colors[0]] * 2 + [uniq_colors[1]] * 6 + [uniq_colors[2]] * 6

        legend_handles = [
            mpatches.Patch(facecolor=color, label=label, edgecolor="black", linewidth=0.5)
            for label, color in zip(uniq_labels, uniq_colors)
        ]

        bars = ax.patches
        for bar, color in zip(bars, colors):
            bar.set_facecolor(color)

        ax.legend(
            handles=legend_handles, bbox_to_anchor=(1.04, 0), loc="lower left", borderaxespad=0, frameon=False
        )

        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=10)
        # Prepare pairs for statistical annotation
        control_condition = f"FSK {stimulant_dose}µM"
        concentrations = ["0.001µM", "0.1µM", "1.0µM"]
        treatment_conditions = [f"{a} {b}" for a, b in product([cpd], concentrations)]
        pairs = [(control_condition, cond) for cond in treatment_conditions]
        ### Add statistical annotations
        annotator = Annotator(ax, pairs, data=plot_data, x="treat_string", y=readout)
        annotator.configure(
            test="Mann-Whitney",
            text_format="star",
            loc="inside",
            hide_non_significant=True,
        )
        annotator.apply_and_annotate()

        ax.set_xlabel("")
        ax.set_ylabel(r"Mean-Aggregated Cyst Size ($\mu$m$^2$)")
        ax.set_title(cpd)
        # ax.set_title(f"Plate {plate_number} - {cpd}")

        ax.grid(axis="y", alpha=0.5)
        ax.spines[["right", "top"]].set_visible(False)
        ax.set_axisbelow(True)

        if savefig:
            if figs_root is None:
                figs_root = Path(__file__).parents[1] / "figures"
            fig.savefig(
                figs_root
                / (
                    f"boxplot-MannWhitneyU-plate{plate_number}{target_str}_{cpd}"
                    f"_FSK{str(stimulant_dose).replace('.', '-')}.png"
                ),
                dpi=300,
                bbox_inches="tight",
                facecolor="white",
            )
            fig.savefig(
                figs_root
                / (
                    f"boxplot-MannWhitneyU-plate{plate_number}{target_str}_{cpd}"
                    f"_FSK{str(stimulant_dose).replace('.', '-')}.svg"
                ),
                bbox_inches="tight",
                facecolor="none",
            )

        plt.show()
