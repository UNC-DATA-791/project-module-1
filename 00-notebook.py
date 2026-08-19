# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo",
#     "pooch",
#     "polars",
#     "altair>=5.4",
#     "vegafusion",
#     "vl-convert-python",
#     "pyarrow",
#     "scikit-learn==1.9.0",
# ]
# ///

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")

with app.setup:
    import altair as alt
    import marimo as mo
    import numpy as np
    import polars as pl
    import pooch
    from sklearn.decomposition import PCA

    # This allows us to make Altair charts with lots of data points
    alt.data_transformers.enable("vegafusion")


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Week 1 — get the data and look at the counts

    In this week's notebook you will access gene expression data from [Bottomly et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0017820) and explore the number of reads we have in each sample and the sparsity of our gene x sample matrix.

    ### Experimental details

    Bottomly and colleagues studied 21 mice from two inbred strains. Mice from the same strain are identical genetically. Ten mice were of the **C57BL/6J** strain, and eleven mice were from strain **DBA/2J**.

    From each mouse the team sampled a part of the brain called the striatum. The team measured the expression of every gene in the striatum samples.
    """)
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Let's get the data!

    We will get the read counts for all samples in the study from [ReCount](https://rna.recount.bio/) which is a public archive of count tables for RNA-seq studies.

    Let's use this [pooch](https://github.com/fatiando/pooch) library to download all the data. Pooch 1) downloads the data, and 2) checks the integrity of the file via the `sha256` hash.
    """)
    return


@app.cell
def _():
    RECOUNT_URL = "https://bowtie-bio.sourceforge.net/recount"

    URLS = {
        "bottomly_count_table.txt": f"{RECOUNT_URL}/countTables/bottomly_count_table.txt",
        "bottomly_phenodata.txt": f"{RECOUNT_URL}/phenotypeTables/bottomly_phenodata.txt",
    }

    count_df = pl.read_csv(URLS['bottomly_count_table.txt'], separator='\t')
    pheno_df = pl.read_csv(URLS['bottomly_phenodata.txt'], separator=' ').rename({'sample.id': 'sample_id'})
    count_df
    return count_df, pheno_df


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    Let's do some quick sanity checks since we know the study design beforehand.
    """)
    return


@app.cell
def _(count_df):
    assert count_df.shape[1] - 1  == 21, "Expected 21 samples."
    assert count_df.shape[0] == 36536, "Expected 36536 genes."
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    Did everything pass? Why are we subtracting one in the first assertion above? **Answer in a markdown cell below.**
    """)
    return


@app.cell
def _():
    # answer here, delete this comment and make this a markdown cell

    return


@app.cell
def _(pheno_df):
    pheno_df
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    How many samples do we have from each strain? (Hint: use `pheno_df`).
    """)
    return


@app.cell
def _():
    # your code here

    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    Join `pheno_df` with `count_df`. You have to unpivot `count_df` into long format.

    Make a barchart of read depth by sample. Color the bars by `strain`.
    """)
    return


@app.cell
def _():
    # your code here

    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## Sparsity

    For each gene and strain, what fraction of samples have a zero count? Let's look at the distribution.

    Starting from `count_df` and `pheno_df`, build a `sparsity_df`:

    1. Unpivot `count_df` into long format, so each row is one (gene, sample) pair with its count.
    2. Join in `strain` from `pheno_df`, matching on the sample ID.
    3. Group by `gene` and `strain`.
    4. Within each group, compute the fraction of samples with a zero count.
    """)
    return


@app.cell
def _():
    # your code here to make the df

    return


@app.cell
def _():
    # your code here for the sparsity histogram

    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    Suppose a gene has zero counts in every sample of one strain, but nonzero counts in the other strain. **How would you quantify the difference in expression between the two strains for this gene?** Think about what happens if you try a standard metric like log fold-change. 🤔
    """).callout(kind="info")
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    ## PCA

    Do the two strains separate out if we look at overall gene expression? Raw counts aren't comparable across samples because library sizes differ (see the read depth chart above), so we normalize first:

    1. Convert counts to CPM (counts per million), so every sample sums to the same total.
    2. Log-transform (`log2(CPM + 1)`) since expression is heavily right-skewed.
    3. Run PCA on the log-CPM matrix (samples x genes) and keep the first two components.
    """)
    return


@app.cell
def _(count_df, pheno_df):
    sample_ids = count_df.columns[1:]
    X = count_df.select(pl.exclude("gene")).to_numpy().T  # samples x genes

    lib_sizes = X.sum(axis=1, keepdims=True)
    log_cpm = np.log2(X / lib_sizes * 1e6 + 1)

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(log_cpm)

    pca_df = pl.DataFrame(
        {
            "sample_id": sample_ids,
            "PC1": pcs[:, 0].tolist(),
            "PC2": pcs[:, 1].tolist(),
        }
    ).join(pheno_df, on="sample_id")

    pca_df
    return pca, pca_df


@app.cell
def _(pca, pca_df):
    var_explained = pca.explained_variance_ratio_

    alt.Chart(pca_df).mark_point(size=80, filled=True).encode(
        x=alt.X("PC1", title=f"PC1 ({var_explained[0]:.0%} variance)"),
        y=alt.Y("PC2", title=f"PC2 ({var_explained[1]:.0%} variance)"),
        color="strain",
        tooltip=["sample_id", "strain"],
    )
    return


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    What is this PCA chart telling us? Answer in a markdown cell.
    """)
    return


@app.cell
def _():
    # answer here, delete this comment and make this a markdown cell

    return


if __name__ == "__main__":
    app.run()
