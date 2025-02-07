import re
import warnings
from functools import partial
from itertools import combinations, compress
from multiprocessing import Pool, cpu_count

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import xgboost as xgb
from matplotlib.cm import ScalarMappable
from rdkit import Chem
from rdkit.Chem import Descriptors
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.covariance import EllipticEnvelope
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import IsolationForest, RandomForestClassifier
from sklearn.feature_selection import RFECV
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.neighbors import KernelDensity, LocalOutlierFactor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import (
    LabelEncoder,
    MaxAbsScaler,
    MinMaxScaler,
    Normalizer,
    PowerTransformer,
    QuantileTransformer,
    RobustScaler,
    StandardScaler,
)
from sklearn.svm import SVC, LinearSVC
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
from venn import venn


def get_spearman_corrs(feat1: str, feat2: str, subset_df: pd.DataFrame = None) -> dict:
    """
    Function to run the non-redundant spearman correlations in parallel.

    Args:
        col_names: columns (features) to be correlated.
        df: dataframe subset for calculating the correlation values.
        Defaults to None.

    Returns:
        dictionary with the correlation values held by the key [f'{feat1}~{feat2}'].
    """
    values1 = subset_df.loc[:, feat1]
    values2 = subset_df.loc[:, feat2]

    result_dict = dict()
    result_dict[f"{feat1}~{feat2}"] = spearmanr(values1, values2)
    return result_dict


def discrete_transform(row):
    """
    Function to convert row valeus into discrete.
    1 -> statistically significant (p-value < 0.05)
    0 -> not statistically significant (p-value > 0.05)
    """

    def transform_function(value):
        t_function = lambda x: 1 if x < 0.05 else 0
        if np.isnan(value):
            return value
        else:
            return t_function(value)

    new_row = row.apply(transform_function)
    return new_row


def scale_features(data, scaling_method="standard"):
    """Scales (standardizes) the data in the dataframe and returns
    it with the same column names. Numpy arrays also suported.

    Args:
        data: Data to be standardized (pandas dataframe or numpy array)
        scaling_method: Standardization method. . Defaults to "standard".

    Raises:
        TypeError: when the input data is not a pandas dataframe or a numpy array.
        ValueError: when the scaling method of choice is not available.

    Returns:
        X_stand (pandas dataframe or numpy array): standardized data
    """
    if not any([isinstance(data, pd.DataFrame), isinstance(data, np.ndarray)]):
        raise TypeError("Data should be a pandas dataframe or a numpy array.")

    scaling_method = scaling_method.lower()
    available = [
        "standard",
        "minmax",
        "maxabs",
        "robust",
        "quantile",
        "power",
        "normalizer",
        "none",
    ]
    if isinstance(data, pd.DataFrame):
        assert (
            all(np.issubdtype(dtype, np.number) for dtype in data.dtypes),
            "Input data: non-numeric data found.",
        )
        columns = data.columns
        X = data.to_numpy()
    else:
        X = data

    if scaling_method == "none":
        return X
    elif scaling_method == "standard":
        X_stand = StandardScaler().fit_transform(X)
    elif scaling_method == "minmax":
        X_stand = MinMaxScaler().fit_transform(X)
    elif scaling_method == "maxabs":
        X_stand = MaxAbsScaler().fit_transform(X)
    elif scaling_method == "robust":
        X_stand = RobustScaler().fit_transform(X)
    elif scaling_method == "quantile":
        X_stand = QuantileTransformer().fit_transform(X)
    elif scaling_method == "power":
        X_stand = PowerTransformer().fit_transform(X)
    elif scaling_method == "normalizer":
        X_stand = Normalizer().fit_transform(X)
    else:
        raise ValueError(f"Type not recognized. Try from {available}")

    if isinstance(data, pd.DataFrame):
        return pd.DataFrame(columns=columns, data=X_stand)
    else:
        return X_stand


def GetMolFromSmi(smi):
    """Simple wrapper function to make it pickleable"""
    return Chem.MolFromSmiles(smi)


def get_physchem_descriptors(
    smiles: list, njobs: int = 1, scaling=None, return_names: bool = False
) -> pd.DataFrame:
    """
    Get physicochemical descriptors from a list of SMILES strings.
    """
    # get the descriptors
    descriptors = [Descriptors._descList[i] for i in range(len(Descriptors._descList))]
    descriptors = np.vstack(descriptors)
    descriptors = dict(zip(descriptors[:, 0], descriptors[:, 1]))
    if return_names:
        return descriptors
    # Removing some of the descriptors that are not my interest:
    desired_descriptors = {}
    dont_want_pattern = re.compile(
        r"PEOE_|Kappa\d|Chi\d|BCUT2D_|Morgan\d|_VSA|TPSA|fr_"
    )
    for desc in descriptors.keys():
        if dont_want_pattern.findall(desc):
            continue
        else:
            if callable(descriptors[desc]):
                desired_descriptors[desc] = descriptors[desc]
            else:
                continue

    with Pool(njobs) as pool:
        mols = pool.map(GetMolFromSmi, smiles)

    # return desired_descriptors
    properties = {key: None for key in desired_descriptors.keys()}
    df = pd.DataFrame(columns=["SMILES"] + list(desired_descriptors.keys()))

    for key in properties.keys():
        properties[key] = [desired_descriptors[key](mol) for mol in mols]

    properties_df = pd.DataFrame.from_dict(properties)
    columns = properties_df.columns
    if scaling is not None:
        df[columns] = scale_features(properties_df, scaling_method=scaling)
    else:
        df[columns] = pd.DataFrame.from_dict(properties)
    return df


def lda_explainability(
    df,
    features,
    y_column,
    max_components: int = 10,
    standardization=None,
    model_kwargs=None,
    ax=None,
):
    """Plots the step explained variance ration for the LDA project

    Args:
        data: whole dataframe with [features, y_column] in it.
        features: the numerical features used to train the LDA.
        y_column: column with the class labels. Encoded by LabelEncoder().
        max_components: max components for the plot. Defaults to 10.
        model_kwargs: kwargs for the PCA. Defaults to None.

    Returns:
        ax, lda, X_train_lda -> matplotlib axis, lda object,
        lda projection values.
    """
    label_encoder = LabelEncoder()

    if standardization is not None:
        X = scale_features(df[features], standardization).to_numpy()
    else:
        X = df[features].to_numpy()

    if model_kwargs is None:
        lda = LinearDiscriminantAnalysis(n_components=max_components)
    else:
        lda = LinearDiscriminantAnalysis(n_components=max_components, **model_kwargs)

    y = label_encoder.fit_transform(df[y_column].to_numpy().ravel())
    X_train_lda = lda.fit_transform(X, y)

    if ax is None:
        ax = plt.gca()
    ax.bar(
        range(1, max_components + 1),
        lda.explained_variance_ratio_,
        align="center",
    )
    ax.step(
        range(1, max_components + 1),
        np.cumsum(lda.explained_variance_ratio_),
        where="mid",
    )
    ax.set_ylabel("Explained variance ratio")
    ax.set_xlabel("Principal components")

    return (
        ax,
        lda,
        X_train_lda,
    )


def pca_explainability(
    data,
    max_components: int = 10,
    scaling_type: str = "standard",
    ax=None,
    model_kwargs=None,
):
    """
    Plots the step explained variance ration for the PCA projection.

    Args:
        data: numpy array of your model data. It should be standardized.
        max_components: max components for the plot. Defaults to 10.
        model_kwargs: kwargs for the PCA. Defaults to None.

    Returns:
        ax, X_fitted_pa -> matplotlib axis, PCA values
    """
    if model_kwargs is None:
        model_kwargs = {"n_components": max_components}
    elif all([model_kwargs is not None, "n_components" not in model_kwargs]):
        model_kwargs["n_components"] = max_components

    pca = PCA(**model_kwargs)
    scaled_data = scale_features(data, scaling_method=scaling_type)
    X_fitted_pca = pca.fit_transform(scaled_data)

    if ax is None:
        ax = plt.gca()
    ax.bar(
        range(1, max_components + 1),
        pca.explained_variance_ratio_,
        align="center",
    )
    ax.step(
        range(1, max_components + 1),
        np.cumsum(pca.explained_variance_ratio_),
        where="mid",
    )
    ax.set_ylabel("Explained variance ratio")
    ax.set_xlabel("Principal components")

    return ax, X_fitted_pca


def fit_pca(
    df: pd.DataFrame,
    features: list,
    n_components: int = 2,
    scaling_method: str = "standard",
):
    """Fits a PCA to the dataset and returns the
    pca object and the transformed dataset.
 
    Args:
        df: dataframe with the data to be projected.
        features: Features used for the dimensionalty reduction
        n_components: Final number of dimensions. Defaults to 2.
        scaling_method: Method used by function `scale_features`.
            Defaults to 'standard'.

    Returns:
        pca (object), projected_x (np.array)
    """
    X = df[features].to_numpy()
    X_stand = scale_features(X, scaling_method)

    pca = PCA(n_components=n_components)
    X_fitted_pca = pca.fit_transform(X_stand)
    return pca, X_fitted_pca


def fit_lda(
    df: pd.DataFrame,
    features: list,
    y_encoding: str,
    n_components: int = 2,
    scaling_method: str = "standard",
    lda_kwargs=None,
):
    """
    Fits a LDA model to the dataset and returns the
    lda object and the transformed dataset.

    Args:
        df: dataframe with the data to be projected.
        features: Features used for the dimensionalty reduction
        y_encoding: Column name with strings to be encoded or with numerical values.
        n_components: Final number of dimensions. Defaults to 2.
        scaling_method: Method used by function `scale_features`. Defaults to 'standard'.

    Raises:
        ValueError: When the y_encoding is not a column in the dataframe.

    Returns:
        pca (object), projected_x (np.array), encoding_dict (dict|None if numerical)
    """
    if all([isinstance(y_encoding, str), y_encoding in df.columns]):
        if df[y_encoding].dtype == "object":
            encoder = LabelEncoder().fit(df[y_encoding])
            y = encoder.transform(df[y_encoding])
            encoding_dict = dict(
                zip(encoder.classes_, encoder.transform(encoder.classes_))
            )
        else:
            y = df[y_encoding].to_numpy()
            encoding_dict = None
    else:
        raise ValueError(
            "The `y_encoding` paramenter should be a column within the dataframe."
        )

    X = df[features].to_numpy()
    X_stand = scale_features(X, scaling_method)

    lda_params = {
        "n_components": n_components,
    }
    if lda_kwargs is not None:
        lda_params.update(lda_kwargs)
    lda = LinearDiscriminantAnalysis(**lda_params)

    X_fitted_lda = lda.fit_transform(X_stand, y)
    return lda, X_fitted_lda, encoding_dict


def fit_kde(df, X_pca, pca, percent=1, sign="<", show_idx=False, **kde_kwargs):
    """
    Function to fit a kde to the data and display a selection
    within the distribution based on the `percent` and the `sign`
    parameters.

    Args:
        X_pca: X-trained vector on the PCA.
        pca: pca object returned from sklearn.decomposition.PCA
        percent: Percentage of the data to be selected from density.
            Defaults to 1.
        sign: Either "<" or ">" to select the data from the density.
            Defaults to "<".

    Returns:
        fig, ax, selected_idxs
    """
    if kde_kwargs is None:
        kde_kwargs = {"kernel": "gaussian", "bandwidth": 0.5}
    kde = KernelDensity(**kde_kwargs).fit(X_pca)
    kde.score_samples(X_pca)

    reset_idx_df = df.reset_index(drop=True)

    # Highlight the points that are in the highest or lowest density region
    # < 1 for lowest 1% or > 99 for highest 1%, for example
    log_dens = kde.score_samples(X_pca)
    threshold = np.percentile(log_dens, percent)
    if sign == "<":
        selection = X_pca[np.exp(log_dens) < np.exp(threshold)]
        selected_idxs = np.where(np.exp(log_dens) < np.exp(threshold))[0]
    elif sign == ">":
        selection = X_pca[np.exp(log_dens) > np.exp(threshold)]
        selected_idxs = np.where(np.exp(log_dens) > np.exp(threshold))[0]

    fig, axs = plt.subplots(ncols=2, figsize=(10, 4))

    # Plot the selection in red
    kde_plot = sns.kdeplot(x=X_pca[:, 0], y=X_pca[:, 1], cmap="Reds", ax=axs[0])

    axs[0].scatter(
        selection[:, 0],
        selection[:, 1],
        color="red",
        label="selection",
        alpha=0.4,
    )

    axs[0].set_xlabel(
        f"PC1, explained variance: {pca.explained_variance_ratio_[0]*100:.2f}%"
    )
    axs[0].set_ylabel(
        f"PC2, explained variance: {pca.explained_variance_ratio_[1]*100:.2f}%"
    )
    axs[0].set_title(
        "Outlier detection through Kernel Density Estimation:\n"
        f"Highlighted: points with exp(log_dens) {sign} {np.exp(threshold):.4f}"
    )
    if show_idx == True:
        for idx in selected_idxs:
            axs[0].text(
                X_pca[idx, 0],
                X_pca[idx, 1],
                str(idx),
            )
    elif isinstance(show_idx, str):
        labels = reset_idx_df.iloc[selected_idxs][show_idx]
        for idx, label in zip(selected_idxs, labels):
            axs[0].text(
                X_pca[idx, 0],
                X_pca[idx, 1],
                label,
            )
    axs[0].legend()  # (bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)
    axs[1].hist(np.exp(log_dens), bins=50)
    axs[1].set_ylabel("Frequency")
    axs[1].set_xlabel("exp(log density)")
    axs[1].set_title("Histogram of exp(log density)")
    axs[1].vlines(
        x=np.exp(threshold),
        ymin=axs[1].get_ylim()[0],
        ymax=axs[1].get_ylim()[1],
        colors=["red"],
        linestyles="dashed",
        label="Selection\nthreshold",
    )
    for tick in axs[1].get_xticklabels():
        tick.set_rotation(45)
    axs[1].legend()  # (bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)
    return fig, axs, selected_idxs, log_dens


def fit_kmeans(df, features, n_clusters: int, **kmeans_kwargs):
    """
    Function to fit a k-means clustering to the data and return
    the km object and the labels (y_km).

    Args:
        df: Dataframe with the data to be clustered.
        features: List of features to be used in the clustering.
        n_clusters: Number of clusters to be used in the clustering.
        kmeans_kwargs: Keyword arguments to be passed to the kmeans object.

    Returns:
        km (object), y_km (labels array)
    """
    if kmeans_kwargs is None:
        kmeans_kwargs = {
            "init": "random",
            "n_init": 10,
            "max_iter": 300,
            "tol": 1e-04,
            "random_state": 0,
        }

    X = df[features].to_numpy()
    X_stand = scale_features(X, scaling_method="standard")

    km = KMeans(n_clusters=n_clusters, **kmeans_kwargs)
    y_km = km.fit_predict(X_stand)

    return km, y_km


def pca_scatter(
    df,
    id_col: str,
    X_fitted_pca: np.ndarray,
    pca: PCA,
    write_label=None,
    cmap="viridis",
    ax=None,
    scatter_kwargs=None,
):
    """Make a PCA scatter plot with the data from the `df` dataframe.
    The `id_col` will be used to color the points based on the `toplot_classes`.

    Args:
        df: pd.DataFrame containg id_col and toplot_classes
        id_col: column with with the identification of the datapoints.
            Will plot all unique values in the column.
        X_fitted_pca: X fitted vactor after PCA
        pca: pca object from sklearn.decomposition.PCA
        cmap: colormap to be used in the plot. Defaults to 'viridis'.

    Returns:
        [ax, class_idx_dict] (dictionary with classes and the indexes
        of the points)
    """
    toplot_classes = df[id_col].unique()

    df = df.copy()
    df_idxs = df.index.to_numpy()

    class_pca_idx_dict = dict()
    class_df_idx_dict = dict()
    # Creating a dictionary with format {class:idxs}
    for class_ in toplot_classes:
        pca_idx = np.where(df[id_col] == class_)[0]
        df_idx = np.take(df_idxs, pca_idx)
        class_pca_idx_dict.update({class_: pca_idx.tolist()})
        class_df_idx_dict.update({class_: df_idx.tolist()})
    # Sort the classes in the dictionary by alphabetical order:
    # class_pca_idx_dict = dict(sorted(class_pca_idx_dict.items()))
    # class_df_idx_dict = dict(sorted(class_df_idx_dict.items()))
    # Setting colors for scatter plot
    scalar_mappable = ScalarMappable(cmap=cmap)
    colors = scalar_mappable.to_rgba(range(len(toplot_classes)), alpha=0.8).tolist()

    if ax is None:
        ax = plt.gca()

    if scatter_kwargs is None:
        scatter_kwargs = dict()

    for color, key in zip(colors, class_pca_idx_dict.keys()):
        ax.scatter(
            X_fitted_pca[:, 0][class_pca_idx_dict[key]],  # x axis
            X_fitted_pca[:, 1][class_pca_idx_dict[key]],  # y axis
            color=color,
            label=key,
            **scatter_kwargs
            # alpha=0.4
        )
    if write_label is not None:
        if write_label.lower() == "index":
            for idx, label in enumerate(df.index):
                ax.text(
                    X_fitted_pca[idx, 0],
                    X_fitted_pca[idx, 1],
                    label,
                )
    ax.set_xlabel(
        f"PC1, explained variance: {pca.explained_variance_ratio_[0]*100:.2f}%"
    )
    ax.set_ylabel(
        f"PC2, explained variance: {pca.explained_variance_ratio_[1]*100:.2f}%"
    )
    ax.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)
    return ax, class_df_idx_dict


def kmeans_inertia(X, n_clusters=range(1, 11), ax=None):
    """Function to plot the inertia of the kmeans clustering
    for the range of clusters given by `n_clusters`.

    Args:
        X: X vector for performing the clustering.
        n_clusters: range of clusters to be tested. Defaults to range(1, 10).
        ax:  matplotlib axis object. Defaults to None.

    Returns:
         ax -> matplotlib axis
    """
    inertia = []
    for n in n_clusters:
        kmeans = KMeans(
            n_clusters=n,
            init="k-means++",
            n_init=10,
            max_iter=300,
            tol=1e-04,
            random_state=0,
        )
        kmeans.fit(X)
        inertia.append(kmeans.inertia_)
    if ax is None:
        ax = plt.gca()
    ax.plot(n_clusters, inertia, marker="o")
    ax.set_xlabel("Number of clusters")
    ax.set_ylabel("Inertia")
    return ax


def plot_multidim_kmeans_clustering(km, y_km, X_fitted_pca, pca, ax=None):
    """Plots the kmeans clustering. Function made for 2D PCA.

    Args:
        km: kmeans object from sklearn.cluster.KMeans
        y_km: y_kmean returned from km.predict(X)
        X_fitted_pca: X fitted vactor after PCA
        pca: pca object from sklearn.decomposition.
        ax: matplotlib axis object. Defaults to None.

    Returns:
        ax -> matplotlib axis
    """
    cluster_idxs = list(range(km.n_clusters))

    scalar_mappable = ScalarMappable(cmap="plasma")
    colors = scalar_mappable.to_rgba(cluster_idxs, alpha=0.8).tolist()
    shapes = ["s", "o", "v", "^", "P", "h", "X"]
    if ax is None:
        ax = plt.gca()

    for clu_idxs, color, shape in zip(cluster_idxs, colors, shapes):
        ax.scatter(
            X_fitted_pca[y_km == clu_idxs, 0],  # x axis
            X_fitted_pca[y_km == clu_idxs, 1],  # y axis
            s=50,
            color=color,
            marker=shape,
            edgecolor="black",
            label=f"Cluster {clu_idxs}",
            alpha=0.7,
        )

    ax.scatter(
        km.cluster_centers_[:, 0],
        km.cluster_centers_[:, 1],
        s=250,
        marker="*",
        color="red",
        edgecolor="black",
        label="Centroids",
    )

    ax.set_title("K-means clustering of the treatment conditions")
    ax.set_xlabel(
        f"PC1, explained variance: {pca.explained_variance_ratio_[0]*100:.2f}%"
    )
    ax.set_ylabel(
        f"PC2, explained variance: {pca.explained_variance_ratio_[1]*100:.2f}%"
    )
    ax.legend()
    return ax


def simple_spearman_corr(
    df,
    features: list,
    discrete: bool = True,
    njobs=5,
):
    """
    Simpler function equevalemt to the `calculate_spearman_corr` function.

    Features should contain the numerical features to perform the validation.

    Usage example:
    >>> corr_heatmap_df, pv_heatmap_df = simple_spearman_corr(df, features, discrete=True)
    >>> fig, ax = selchem.plot_spearman_corr_heatmap(corr_heatmap_df,pv_heatmap_df)
    """
    df = df[features].copy()

    # Calculate spearman correlation for both non-redundant feature combinations
    combi = list(combinations(features, 2))
    with Pool(njobs) as pool:
        results = pool.starmap(partial(get_spearman_corrs, subset_df=df), combi)

    # Nested dictionaries to hold the values for the dataframe
    pvals_dict = {val: {v: np.NaN for v in features} for val in features}
    corrs_dict = {val: {v: np.NaN for v in features} for val in features}

    pvals_to_correct = list()
    spear_corrs = list()
    key_names = list()

    # unpacking results from multiprocessing
    # values are held by a key in the form of (feature1~feature2)
    for d in results:
        key = list(d.keys())[0]
        col1, col2 = key.split("~")
        corr, pval = [i for i in d[key]]
        pvals_dict[col1][col2] = pval
        corrs_dict[col1][col2] = corr
        pvals_to_correct.append(pval)
        spear_corrs.append(corr)
        key_names.append(key)

    # Correcting pvalues for multiple testing and updating original dictionary
    multitest_result = multipletests(pvals=pvals_to_correct, method="fdr_bh")
    print(
        "Number of null hypothesis that were rejected:",
        np.count_nonzero(multitest_result[0] == True),
    )
    print(
        "Number of null hypothesis that were accepted:",
        np.count_nonzero(multitest_result[0] == False),
    )

    # Updating the dictionary with the corrected pvalues
    corrected_pvals = multitest_result[1]
    for pv, k in zip(corrected_pvals, key_names):
        col1, col2 = k.split("~")
        pvals_dict[col1][col2] = pv

    pv_heatmap_df = pd.DataFrame.from_dict(pvals_dict, orient="index")
    corr_heatmap_df = pd.DataFrame.from_dict(corrs_dict, orient="index")
    signif_heatmap_df = pv_heatmap_df.apply(discrete_transform, axis=0)

    if discrete:
        return corr_heatmap_df, signif_heatmap_df
    else:
        return corr_heatmap_df, pv_heatmap_df


def drop_corr_features(
    df: pd.DataFrame,
    feature_cols: list,
    corr_threshold: float = 0.9,
    corr_method: str = "pearson",
):
    """
    Function to drop highly correlated features from a dataframe.

    Args:
        df: pd.DataFrame with the features to be evaluated.
        feature_cols: subset of features to be evaluated.
        corr_threshold: correlation threshold for feature selection.
            positive and negative values will be evaluated. Defaults to 0.9.
        corr_method: supported: `pearson`, `kendall`, `spearman`.
            Defaults to "pearson".
        plot_corr: Whether to plot the data on current axis. Defaults to False.

    Raises:
        AttributeError: if `corr_method` is not supported.

    Returns:
        df[selected features] (pd.DataFrame), correlation matrix (pd.DataFrame)
    """
    t = corr_threshold
    if corr_method.lower() not in ["pearson", "kendall", "spearman"]:
        raise AttributeError(
            "Invalid corr_method. Supported methods:  'pearson', 'kendall', 'spearman'"
        )

    contains_cols = all(np.isin(np.array(feature_cols), df.columns))
    assert contains_cols, "features_cols not in df.columns"

    data = df[feature_cols].copy()
    corr_matrix = data.corr(method=corr_method)
    # Takes only the upper triangle of the similarity matrix
    upper_bool = np.triu(np.ones(corr_matrix.shape), k=1).astype(bool)
    upper_tri = corr_matrix.where(upper_bool)

    to_drop = [
        c for c in upper_tri.columns if any((upper_tri[c] > t) | (upper_tri[c] < -t))
    ]
    new_df = df.copy().drop(columns=to_drop)
    return new_df, corr_matrix


def outlier_flagging(
    df,
    features: list,
    toflag_treatments: list = ["neg_control", "pos_control"],
    contamination: float = "auto",
    n_jobs: int = "auto",
    verbose: bool = False,
):
    """
    Function to perform outlier flagging on the dataframe using a Isolation Forest
    algorithm. The function will produce an additional column named "Outlier" with
    the boolean values of the outliers. This function will use `self.desired_cols`
    to perform the flagging.

    Args:
        toflag_treatments: Values from 'TreatmentType' column to be flaggged. If
        set to 'all', will apply to all treatment conditions.
            Defaults to ["neg_control", "pos_control"].
        contamination: parameter for the sklearn function. Defaults to "auto".
        n_jobs: number of jobs to run in parallel. Defaults to "auto".
        verbose: Show print statements, and progress bar. Defaults to False.

    Returns:
        toflag_list: List with the indexes of the outliers.
    """
    if n_jobs == "auto":
        n_jobs = cpu_count() - 2
    if verbose:
        disable = False
    else:
        disable = True
    subset_df = df.query(f"TreatmentType.isin(@toflag_treatments)").query("QC == 'OK'")
    iso_forest = IsolationForest(
        n_jobs=n_jobs, random_state=0, contamination=contamination
    )
    toflag_list = []
    for t in tqdm(toflag_treatments, disable=disable):
        grouped_plates = subset_df.query(f"TreatmentType == @t").groupby("PlateID")
        for pnumber, plate_df in grouped_plates:
            X = plate_df[features].to_numpy()
            df_index = plate_df.index.to_numpy()
            y_pred = iso_forest.fit(X).predict(X)

            outlier_idx = np.where(y_pred == -1)[0]
            outlier_df_idx = df_index[outlier_idx].tolist()
            toflag_list = toflag_list + outlier_df_idx

    df = df.assign(Outlier=df.index.isin(toflag_list))
    if verbose:
        print(f"Flagged {len(toflag_list)} out of {len(subset_df)} wells as outliers.")
    return df, toflag_list


def multivar_kde_outlier(
    df: pd.DataFrame,
    feature_names: list,
    contamination: float,
    bandwidth: float,
    col4scatter="Identifier",
):
    """
    Function to visualize the outlier samples through multivariate
    kernel density estimation. The way I'm scoring and plotting the
    outlier_idx came from the following source:
    https://towardsdatascience.com/kernel-density-estimation-for-anomaly-detection-in-python-part-1-452c5d4c32ec


    Args:
        df: pd.DataFrame
        feature_names: list of column names with desired features.
        contamination: fraction of the data will be flagged.
        bandwidth: bandwidth for the kernel density estimation.
        col4scatter: column for the scatterplot point labels.

    Raises:
        ValueError: whe the contamination is set to higher than 0.5.

    Returns:
        fig, axs, df.iloc[outliers_idx] -> matplotlib figure, axes, and
        the dataframe with the outlier_idx.
    """
    if contamination >= 0.5:
        raise ValueError("Contamination must be less than 0.5")
    if col4scatter == "Identifier":
        df = df.assign(
            Identifier=lambda x: x["Screening"]
            + x["PlateID"].astype(str)
            + "_"
            + x["PlateRow"]
            + x["PlateColumn"].astype(str)
        ).copy()  # Assigning an identifier for the scatterplot

    X = df[feature_names].to_numpy()
    labels = df["Compound"].values
    X_stand = scale_features(X, scaling_method="standard")

    pca = PCA(n_components=2)
    X_fitted_pca = pca.fit_transform(X_stand)

    kde = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(X_stand)

    """
    For performing the outlier detection, we'll normalize the log density by dividing it
    by its norm.
    
    This is achieved by using `numpy.linalg.norm()` which, if used like this, returns
    the Frobenius norm (also known as the euclidian or the L2 norm). This norm, when 
    calculated for a $m\times n$ matrix $A$, is defined by the square root of the sum 
    of the absolute squares of its elements.
    
    $||A||_F = \sqrt{\sum_{i=1}^{m}\sum_{j=1}^{n}|a_{ij}|^2}$
    
    (Text from [MathWorld](https://mathworld.wolfram.com/FrobeniusNorm.html))
    """
    log_dens = kde.score_samples(X_stand)
    scores = log_dens / np.linalg.norm(log_dens)

    threshold = np.percentile(scores, contamination * 100)
    # Highlight the points that are in the lowest density region
    # (i.e. the points that are most likely to be outlier_idx)
    outlier_idx = X_fitted_pca[scores < threshold]
    outliers_idx = np.where(scores < threshold)[0]
    other_samples = X_fitted_pca[scores >= threshold]

    fig, axs = plt.subplots(ncols=2, figsize=(12, 4))

    axs[0].scatter(
        outlier_idx[:, 0],
        outlier_idx[:, 1],
        color="darkorange",
        alpha=1,
        label="outlier_idx",
    )
    axs[0].scatter(
        other_samples[:, 0],
        other_samples[:, 1],
        color="lime",
        alpha=0.2,
        label="other samples",
    )

    axs[0].set_xlabel(
        f"PC1, explained variance: {pca.explained_variance_ratio_[0]*100:.2f}%"
    )
    axs[0].set_ylabel(
        f"PC2, explained variance: {pca.explained_variance_ratio_[1]*100:.2f}%"
    )
    axs[0].legend()

    if col4scatter is not None:
        labels = df.iloc[outliers_idx][col4scatter].values
        # rotation = 9
        for idx, label in zip(outliers_idx, labels):
            axs[0].text(
                X_fitted_pca[idx, 0], X_fitted_pca[idx, 1], label, rotation=0, size=8
            )
    axs[1].hist(scores, bins=50)
    axs[1].set_ylabel("Frequency")
    axs[1].set_xlabel("$\\frac{log(density)}{norm(density)}$", fontsize=14)
    axs[1].set_title("Histogram of $log(density)/norm(density)$")
    axs[1].vlines(
        x=threshold,
        ymin=axs[1].get_ylim()[0],
        ymax=axs[1].get_ylim()[1],
        colors=["red"],
        linestyles="dashed",
        label="Selection\nthreshold",
    )
    axs[1].set_ylim(0, 100)
    for tick in axs[1].get_xticklabels():
        tick.set_rotation(45)
    # axs[0].legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.0)

    print(f"Total number of samples: {(len(df))}")
    print(f"{len(outliers_idx)} Outliers detected")
    print(df.iloc[outliers_idx]["Identifier"].values)
    return fig, axs, df.iloc[outliers_idx]


def multivar_sklearn_outlier(
    df: pd.DataFrame,
    feature_names: list,
    method: str,
    contamination="auto",
    col4scatter=None,
    suppress_warnings: bool = True,
    ax=None,
    n_jobs: int = "auto",
):
    """
    Function to apply and visualize results from a outlier (or novelty) detection
    method from scikit-learn.

    Args:
        df: dataframe with the data.
        feature_names: list of columns -> features to be used for outlier detection.
        method: name of the sklearn algortihm to be used.
        contamination: fraction of contamination [0 to 1]. Defaults to "auto".
        col4scatter: not implemented yet. Defaults to "Identifier".
        suppress_warnings: supresses sklearn warnings. Defaults to True.
        ax: matplotlib axis object. Defaults to None.
        n_jobs: number of jobs. Defaults to "auto".

    Raises:
        ValueError: when the method is not avaiable/not implemented.

    Returns:
        ax, class_pca_idx_dict [dict], algorithm: matplotlib ax, dictionary with
        the indexes belonging to each class, trained outlier prediction agorithm.
    """
    # To add documentation
    if suppress_warnings:
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        warnings.filterwarnings("ignore", category=UserWarning)

    if n_jobs == "auto":
        n_jobs = cpu_count() - 2
    anomaly_algorithms = {
        "robust covariance scaled": make_pipeline(
            StandardScaler(),
            EllipticEnvelope(
                contamination=contamination,
                assume_centered=True,
                random_state=0,
            ),
        ),
        "robust covariance": EllipticEnvelope(contamination=contamination),
        "isolation forest": IsolationForest(
            contamination=contamination, random_state=0, n_jobs=n_jobs
        ),
        "local outlier factor scaled": make_pipeline(
            StandardScaler(),
            LocalOutlierFactor(
                n_neighbors=35, contamination=contamination, n_jobs=n_jobs
            ),
        ),
        "local outlier factor": LocalOutlierFactor(
            n_neighbors=35,
            contamination=contamination,
            n_jobs=n_jobs,
        ),
    }
    X = df[feature_names].to_numpy()
    df_index = df.index.to_numpy()

    if method.lower() not in anomaly_algorithms.keys():
        raise ValueError(
            f"Method {method} not recognized. Use one of the following: {list(anomaly_algorithms.keys())}"
        )
    else:
        algorithm = anomaly_algorithms[method.lower()]

    algorithm.fit(X)

    # These method don't have the .predict() method
    if method.lower() in ["local outlier factor", "local outlier factor scaled"]:
        y_pred = algorithm.fit_predict(X)
    else:
        y_pred = algorithm.fit(X).predict(X)
    outlier_idx = np.where(y_pred == -1)[0]
    outlier_df_idx = df_index[outlier_idx]

    df = df.assign(
        Selection=["Outlier" if idx in outlier_df_idx else "Normal" for idx in df.index]
    )

    pca, X_fitted_PCA = fit_pca(df, feature_names)

    ax, class_df_idx_dict = pca_scatter(
        df, "Selection", X_fitted_PCA, pca, cmap="coolwarm", ax=ax
    )

    if col4scatter is not None:
        labels = df.iloc[outlier_df_idx][col4scatter].values
        # rotation = 9
        for idx, label in zip(outlier_df_idx, labels):
            ax.text(
                X_fitted_PCA[idx, 0], X_fitted_PCA[idx, 1], label, rotation=0, size=8
            )

    ax.set_title(f"{method.capitalize()} - {len(outlier_idx)} outlier_idx")

    return ax, class_df_idx_dict, algorithm


class feature_selector:
    """
    Class for feature selection. Currently only supports the recursive
    feature elimination with cross validation, but the idea is that it
    could support other methods as well. To initialize the an object of
    this class, you should provide a dataframe and the numeric features
    that will be used for training.

    Code example:
    >>> encoder =  {'neg_control': 1, 'pos_control': 0}
    >>> selector = feature_selector(raw_df, feature_cols)
    >>> selector.estimator_from_known('svc')
    >>> task_rfecv = selector.rfecv_feature_selector(
    >>>     'pos_x_neg', 'TreatmentType',
    >>>     query_strings=['pos_control', 'neg_control'],
    >>>     encoder=encoder)
    >>> fig, axs = selector.plot_rfecv_results(task='pos_x_neg')
    """

    def __init__(self, df, features) -> None:
        self.df = df.copy()
        self.all_features = features
        self.id_cols = [c for c in self.df.columns if c not in features]
        self.estimator = None
        self.rfecv = dict()
        self.rfecv_folds = dict()
        self.rfecv_scoring = dict()
        self._selected_features = dict()

    def initialize_estimator(self, estimator, estimator_kwargs):
        """
        Load a estimator from the sklearn.feature_selection module
        and initialize it with the given parameters. This estimator
        will be used by other methodsof this class.

        Args:
            estimator: sklearn.feature_selection estimator
            **estimator_kwargs: keyword arguments for the estimator.
        """
        self.estimator = estimator(**estimator_kwargs)
        return self.estimator

    def estimator_from_known(self, estimator_name: str, estimator_kwargs: dict = None):
        """
        Same as initialize_estimator but with estimator already implemented
        within this module.

        Args:
            estimator_name: the name of the estimator already implemented.

        Returns:
            The estimator object with the initialized parameters.
        """
        cores = cpu_count() - 2

        if estimator_name.lower() == "svc":
            if estimator_kwargs is None:
                print(f"Using internal {estimator_name} parameters.")
                estimator_kwargs = {
                    "kernel": "linear",
                    "C": 1,
                    "random_state": 0,
                    "class_weight": "balanced",
                }
            self.estimator = SVC(**estimator_kwargs)

        elif estimator_name.lower() == "linearsvc_l1":
            if estimator_kwargs is None:
                print(f"Using internal {estimator_name} parameters.")
                estimator_kwargs = {
                    "penalty": "l1",
                    "C": 1,
                    "random_state": 0,
                    "class_weight": "balanced",
                    "dual": False,
                    "max_iter": 3000,
                }
            self.estimator = LinearSVC(**estimator_kwargs)
        elif estimator_name.lower() == "linearsvc_l2":
            if estimator_kwargs is None:
                print(f"Using internal {estimator_name} parameters.")
                estimator_kwargs = {
                    "penalty": "l2",
                    "C": 1,
                    "random_state": 0,
                    "class_weight": "balanced",
                    "dual": False,
                    "max_iter": 3000,
                }
            self.estimator = LinearSVC(**estimator_kwargs)

        elif estimator_name.lower() == "xgboost":
            if estimator_kwargs is None:
                print(f"Using internal {estimator_name} parameters.")
                estimator_kwargs = {
                    "eval_metric": "logloss",
                    "n_jobs": cores,
                    "random_state": 0,
                    "use_label_encoder": False,
                }
            self.estimator = xgb.XGBClassifier(**estimator_kwargs)

        elif estimator_name.lower() == "logistic_regression":
            if estimator_kwargs is None:
                print(f"Using internal {estimator_name} parameters.")
                estimator_kwargs = {
                    "penalty": "l2",
                    "solver": "lbfgs",
                    "random_state": 0,
                    "n_jobs": cores,
                }
            self.estimator = LogisticRegression(**estimator_kwargs)

        elif estimator_name.lower() == "random_forest":
            if estimator_kwargs is None:
                print(f"Using internal {estimator_name} parameters.")
                estimator_kwargs = {
                    "n_estimators": 100,
                    "random_state": 0,
                    "n_jobs": cores,
                }
            self.estimator = RandomForestClassifier(**estimator_kwargs)

        else:
            raise ValueError(f"Unknown estimator: {estimator_name}")
        return self.estimator

    def rfecv_feature_selector(
        self,
        task: list,
        column: str,
        encoder: dict,
        query_strings: list = None,
        cv_folds: int = 8,
        scoring: str = "roc_auc",
        standardize: bool = True,
        rfecv_kwargs: dict = None,
    ):
        """
        Should be used for a binary classification task. This method
        will select features that best discriminate between the two
        classes by performing a recursive feature elimination with
        8x cross validation, using a support vector classifier as the
        estimator.

        Args:
            task: task name to be saved under `self._selected_features`
            column: dataframe column to be used as the target variable
            query_strings: list of classes within `column` to use pd.DataFrame.query()
            encoder: dictionary with the mapping of the classes
            cv_folds: Folds for the cross validation. Defaults to 8.
            scoring: sklearn performance metric. Defaults to "roc_auc".
            standardize: To standardized or not. Defaults to True.

        Returns:
            task_rfecv -> sklearn.feature_selection.RFECV object
        """
        if query_strings is not None:
            df = self.df.query(f"{column}.isin({query_strings})").copy()
            print(
                f"Number of {query_strings[0]} samples: ",
                df.query(f'{column} == "{query_strings[0]}"').shape[0],
            )
            print(
                f"Number of {query_strings[1]} samples: ",
                df.query(f'{column} == "{query_strings[1]}"').shape[0],
            )
        else:
            df = self.df
        X = df[self.all_features]
        y = df[column].map(encoder).values
        X_stand = scale_features(X, scaling_method="standard")

        # Loading the estimator for the rfecv
        if self.estimator is None:
            raise ValueError(
                "Estimator not initialized. Use initialize_estimator()"
                "or estimator_from_known()"
            )

        rfecv_parameters = {
            "estimator": self.estimator,
            "step": 1,
            "cv": StratifiedKFold(cv_folds),
            "scoring": scoring,
        }
        if rfecv_kwargs is not None:
            print("Seeting default parameters for RFECV.")
            rfecv_parameters.update(rfecv_kwargs)
        # For available scoring functions, check
        # : https://scikit-learn.org/stable/modules/model_evaluation.html
        task_rfecv = RFECV(**rfecv_parameters)

        if standardize:
            task_rfecv.fit(X_stand, y.ravel())
        else:
            task_rfecv.fit(X, y.ravel())

        # Updating dictionaries with the results
        self.rfecv.update({task: task_rfecv})
        self.rfecv_scoring.update({task: scoring})
        self.rfecv_folds.update({task: cv_folds})
        self._selected_features.update(
            {task: list(compress(self.all_features, task_rfecv.support_))}
        )

        print("---Done---")
        print("Optimal number of features : %d" % task_rfecv.n_features_)
        print(f"Best features : {self._selected_features[task]}")
        return task_rfecv

    def plot_rfecv_results(self, task: str):
        """Plots the results of the RFECV for a given task.

        Args:
            task: name of the task from rfecv_feature_selector

        Returns:
            fig, axs: matplotlib figure and axes
        """

        fig, axs = plt.subplots(ncols=2, figsize=(10, 4))

        split_pattern = re.compile("split(\d+)_test")

        task_rfecv = self.rfecv[task]
        cv_folds = self.rfecv_folds[task]
        scoring = self.rfecv_scoring

        for key in task_rfecv.cv_results_.keys():
            if split_pattern.findall(key):
                axs[0].plot(task_rfecv.cv_results_[key], label=key.split("_")[0])
        axs[0].set_xlabel("Number of features")
        axs[0].set_ylabel(f"Cross validation score ({scoring[task]})")
        axs[0].set_title(f"RFE {cv_folds}x CV - {scoring[task]} by number of features")
        axs[0].legend(ncol=2)

        # Plot features VS. cross-validation scores using task_rfecv.cv_results_
        axs[1].plot(
            range(0, len(task_rfecv.cv_results_["mean_test_score"]) + 0),
            task_rfecv.cv_results_["mean_test_score"],
        )
        axs[1].set_xlabel("Number of features")
        axs[1].set_ylabel(f"Cross validation score ({scoring[task]})")
        axs[1].set_title(f"RFE {cv_folds}x CV - Mean {scoring[task]} values")
        fig.tight_layout()
        return fig, axs

    def venn_best_features(self, tasks: list, ax=None):
        if ax is None:
            ax = plt.gca()
        plot_dict = {task: set(self._selected_features[task]) for task in tasks}
        venn(plot_dict, ax=ax)
        print("Features in common: ", list(set.intersection(*plot_dict.values())))
        return ax
