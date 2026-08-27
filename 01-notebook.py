# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "altair==6.2.2",
#     "marimo>=0.24.0",
#     "numpy==2.5.2",
#     "polars==1.44.1",
#     "pyarrow==25.0.1",
#     "statsmodels==0.15.0",
# ]
# ///

import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Unit 2 assignment

    In this assignment we are going to calculate log fold change of gene expression. We'll use the experiment and data from [Bottomly et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3063777/).
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Unit 2 assignment

    In this assignment we are going to calculate log fold change of gene expression. We'll use the experiment and data from [Bottomly et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3063777/).
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    This cell is the same task as our first assignment. Load both files into polars tables from the URLs below, and rename the phenodata's `sample.id` column to `sample_id`.

    ### 🚧 Needs your help
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Unit 2 assignment

    In this assignment we are going to calculate log fold change of gene expression. We'll use the experiment and data from [Bottomly et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC3063777/).
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    We also have a table `pheno_df` that holds sample metadata. We'll randomly split the samples into two groups, making sure both mouse strains are approximately equally represented in each group. The idea is that we can see how stable our log fold change estimates are by comparing the values we get for each gene across each sample group.

    **Look at the `split` column. Does it make sense?**
    """)
    return


@app.cell
def _(pheno_df, pl):
    seed = 42

    pheno_df_split = pheno_df.with_columns(
        (pl.int_range(pl.len()).shuffle(seed=seed).over("strain") >= 5)
        .cast(pl.Int8)
        .alias("split")
    )

    pheno_df_split
    return pheno_df_split, seed


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Convenience functions
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    This function takes a gene ID and returns the counts for just that gene. **Run it yourself to see how it works.**
    """)
    return


@app.cell
def _(count_df, pheno_df_split, pl):
    def get_gene_table(gene_id='ENSMUSG00000055010', count_df=count_df, pheno_df=pheno_df_split):
        return(
            count_df
            .unpivot(index=['gene'], variable_name='sample_id', value_name='count')
            .with_columns(
                pl.col('count').sum().over('sample_id').alias('total_count')
            )
            .filter(pl.col('gene')==gene_id)
            .join(pheno_df_split, on=['sample_id'])
        )

    return (get_gene_table,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    `get_gene_lfc` takes a count table for one gene (produced with `get_gene_table`) and computes its log2 fold change between strains, separately for each of the two sample splits.

    We're fitting a Poisson regression where `count` depends on `strain`, using `total_count` as an offset to account for each sample's sequencing depth:

    $$
    \log\big(E[\texttt{count}]\big) = \beta_0 + \beta_1 \cdot \texttt{strain} + \log(\texttt{total\_count})
    $$

    Here $\beta_1$ is the **log fold change (LFC)** between strains, on the natural-log scale. This is the value you convert to log2 fold change in step 3 below.

    You've been given the model formula and the `statsmodels` code to create and fit the model. **You need to finish this function**.

    For each split in the input gene count table, the function should:

    1. Filter for the rows in that split.
    2. Fit a Poisson regression of `count` on `strain`, using the log of `total_count` as an offset (to adjust for sequencing depth).
    3. Convert log fold change from natural log to log2 fold change.

    It should return a dictionary with the gene's ID and its log2 fold change for each split, e.g. `{'gene_id': ..., 'log2fc_0': ..., 'log2fc_1': ...}`.

    ### 🚧 Needs your help
    """)
    return


@app.cell
def _(pl):
    def get_gene_lfc(gene_df):
        formula = "count ~ C(strain, Treatment(reference='C57BL/6J'))"

        rec = {}
        for split in [0, 1]:
        
            split_df = gene_df.filter(pl.col('split')==split).to_pandas()
        
            # This is the model code assuming your split specific dataframe is called "split_df" 
            # res = smf.poisson(formula, data=split_df, offset=np.log(split_df['total_count'])).fit(disp=False)
        
            raise NotImplementedError("Fit the Poisson model, convert the coef. to log2fc, and store it in rec.")

        return rec

    return (get_gene_lfc,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    The code below runs the `get_gene_table` and `get_gene_lfc` functions to produce a table with log2 fold change estimates for each gene (within each "split").
    """)
    return


@app.cell
def _(count_df, get_gene_lfc, get_gene_table, mo, pl, seed):
    with mo.persistent_cache(name="lfc_split_half"):
        selected_genes = (
            count_df
            .with_columns(pl.sum_horizontal(pl.exclude('gene')).alias('total_count'))
            .filter(pl.col('total_count') >= 10)
            .get_column('gene')
            .sample(n=2000, seed=seed)
        )

        lfc_df = pl.DataFrame(
            [get_gene_lfc(get_gene_table(gene_id=gene_id)) for gene_id in selected_genes]
        )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Make a scatter chart of log2 fold change for split1 versus split2.

    🚧 Needs your help
    """)
    return


@app.cell
def _():
    # your chart code here


    return


if __name__ == "__main__":
    app.run()
