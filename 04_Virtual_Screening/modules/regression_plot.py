import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score


def plot_errors(ax, linemin, linemax):
    ax.plot(
        [linemin, linemax],
        [linemin, linemax],
        color="black",
        linewidth=2,
        linestyle="--",
        alpha=1,
        zorder=0,
        label="Perfect fit",
    )
    toadd_elements = [
        [1, "--", "1 log unit error", 0.7],
        [2, ":", "2 log unit error", 0.2],
    ]
    for elements in toadd_elements:
        factor, linestyle, label, alpha = elements
        ax.plot(
            [linemin + factor, linemax + factor],
            [linemin, linemax],
            color="grey",
            linewidth=2,
            linestyle=linestyle,
            alpha=alpha,
            zorder=0,
        )
        ax.plot(
            [linemin, linemax],
            [linemin + factor, linemax + factor],
            color="grey",
            linewidth=2,
            linestyle=linestyle,
            alpha=alpha,
            zorder=0,
            label=label,
        )
    return ax


def plot_regression_performance(
    df,
    xaxis: str,
    yaxis: str,
    hue: str = None,
    cmap="plasma",
    alpha=0.5,
    s=75,
    vmax=1,
    vmin=0.4,
):
    """Plots the performance of the regression model including the fitted line, the
    perfect fit line and the 1 and 2 log unit error lines.

    Args:
        df: dataframe containing the data to plot.
        xaxis: name of the column to plot on the x axis.
        yaxis: name of the column to plot on the y axis.
        hue: hue to use for the color of the points. Usually similarity metric.
        cmap: Colormap to be used in the plot. Defaults to "plasma".
        alpha: alpha value for transparency of the points. Defaults to 0.5.
        s: size of the points in the scatterplot. Defaults to 75.
        vmax: maximum value for the `hue` metric (usually similarity). Defaults to 1.
        vmin: minimum value for the `hue` metric (usually similarity). Defaults to 0.4.

    Returns:
        Matplotlib figure and axis.
    """
    # Make the plot
    fig, ax = plt.subplots(figsize=(6, 6))
    sc = plt.scatter(
        df[xaxis].values,
        df[yaxis].values,
        c=(df[hue].values if hue is not None else "lightblue"),
        cmap=cmap,
        s=s,
        edgecolors="black",
        alpha=alpha,
        vmax=vmax,
        vmin=vmin,
    )
    if hue is not None:
        cax = plt.axes([0.95, 0.3, 0.05, 0.55])  # [left, bottom, width, height]
        cbar = plt.colorbar(
            sc,
            cax=cax,
            label="Max similarity to training set",
        )
        cbar.set_ticks([0.4, 0.6, 0.8, 1])
    # Set the scale and limits for the plot
    ax.set(xscale="linear", yscale="linear")
    _max = max(df[xaxis].max(), df[yaxis].max())
    _min = min(df[xaxis].min(), df[yaxis].min())
    ax.set_xlim(_min - 0.5, _max + 0.5)
    ax.set_ylim(_min - 0.5, _max + 0.5)

    # Define how the ticks are formatted
    ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

    lr = LinearRegression()
    lr.fit(df[[xaxis]].values, df[yaxis].values)
    xmin, xmax = ax.get_xlim()
    linspace = np.linspace(xmin, xmax, 100).ravel()
    predictions = lr.predict(linspace.reshape(-1, 1))
    r2 = r2_score(df[xaxis], df[yaxis])
    ax.plot(
        linspace,
        predictions,
        color="tab:red",
        label=f"Fitted Line (R² = {r2:.2f})",
        alpha=0.8,
    )
    linemin = 0
    linemax = max(df[xaxis].max(), df[yaxis].max()) + 0.5
    # Add lines at plus and minus 5-fold error in log scale
    plot_errors(ax, linemin, linemax)
    ax.legend(bbox_to_anchor=(1.04, 0), loc="lower left", borderaxespad=0)
    # Define the axis labels
    ax.set_xlabel(xaxis)
    ax.set_ylabel(yaxis)
    # Add some text
    ax.grid(axis="both", alpha=0.5, linestyle="--")
    ax.set_axisbelow(True)

    # Add a title
    ax.set_title("Model Performance")
    return fig, ax
