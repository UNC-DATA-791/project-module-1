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
    ## Unit 3 assignment notebook


    In this notebook, we will calculate log2-fold change of gene expression across mouse strains from the [Bottomly et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0017820) dataset using **negative binomial regression**.
    """)
    return


@app.cell
def _():
    import altair as alt
    alt.data_transformers.enable("vegafusion")

    import numpy as np
    import polars as pl
    import statsmodels.formula.api as smf
    import warnings

    return alt, np, pl, smf


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Load the data
    """)
    return


@app.cell
def _(pl):
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
def _(mo):
    mo.md(r"""
    Below is a histogram of mean count across all samples for all the genes in the dataset. For the log2-fold change analysis, we'll select genes for which we have moderate counts. The orange ribbon shows the range of log10 mean count of genes we will select for analysis.
    """)
    return


@app.cell(hide_code=True)
def _(alt, count_df, pl):
    source = (
        count_df
        .with_columns(
            pl.mean_horizontal([c for c in count_df.columns if c != 'gene']).alias('mean_count'),
        )
        .filter(pl.col('mean_count')>10)
        .with_columns(
            pl.col('mean_count').log10().alias('log_mean_count')
        )
    )

    upper_threshold = 2.25
    lower_threshold = 2
    moderately_expressed = source.filter(pl.col('log_mean_count')>lower_threshold, pl.col('log_mean_count')<upper_threshold).get_column('gene')

    _ribbon = alt.Chart(
        pl.DataFrame({'lower': [lower_threshold], 'upper': [upper_threshold]})
    ).mark_rect(opacity=0.5, color='orange').encode(
        x=alt.X('lower').scale(zero=False),
        x2=alt.X2('upper'),
    )

    _hist = alt.Chart(source, width=600, height=200).mark_bar().encode(
        x=alt.X('log_mean_count').bin(maxbins=50),
        y=alt.Y('count()')
    )

    _bars = _ribbon + _hist

    _bars
    return (moderately_expressed,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    This `Series` holds the IDs of all the genes we will consider.
    """)
    return


@app.cell
def _(moderately_expressed):
    moderately_expressed
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Convenience functions.

    These functions get a gene-specific count table and calculate log2-fold change for a gene count table using negative binomial regression.
    """)
    return


@app.cell
def _(count_df, np, pheno_df, pl, smf):
    def get_gene_table(gene_id='ENSMUSG00000055010', count_df=count_df, pheno_df=pheno_df):
        return(
            count_df
            .unpivot(index=['gene'], variable_name='sample_id', value_name='count')
            .with_columns(
                pl.col('count').sum().over('sample_id').alias('total_count')
            )
            .with_columns(
                (pl.col('count')/pl.col('total_count')).alias('count_norm')
            )
            .filter(pl.col('gene')==gene_id)
            .join(pheno_df, on=['sample_id'])
        )


    def get_gene_lfc(gene_df):
        gene_df = gene_df.to_pandas()
        gene_id = gene_df['gene'].iloc[0]
    
        formula = "count ~ C(strain, Treatment(reference='C57BL/6J'))"
        fit_result = smf.negativebinomial(formula, data=gene_df, offset=np.log(gene_df['total_count'])).fit(disp=False)
    
        coef = "C(strain, Treatment(reference='C57BL/6J'))[T.DBA/2J]"
        log2fc = fit_result.params[coef]
        p_value = fit_result.pvalues[coef]
    
        return {'gene': gene_id, 'log2fc': log2fc, 'p_value': p_value}

    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    For each gene in the `moderately_expressed` list,

    1. Find log2 fold-change in expression in strain *DBA/2J* vs. *C57BL/6J* with the `get_gene_table` and `get_gene_lfc` functions.
    2. Collect the results into a `pl.DataFrame` called `lfc_df`, and
    3. Show the dataframe as a sortable, single-selection `mo.ui.table` called `lfc_table` (select the first row by default).

    Some genes won't fit cleanly and you can silence statsmodels' `ConvergenceWarning` and `HessianInversionWarning` while fitting (see solution for how to do this).

    <details>
    <summary>👀 Take a peek at the solution</summary>

    ```python
    from statsmodels.tools.sm_exceptions import ConvergenceWarning, HessianInversionWarning

    with warnings.catch_warnings():

        warnings.simplefilter('ignore', category=ConvergenceWarning)  # this silences the warning from statsmodels
        warnings.simplefilter('ignore', category=HessianInversionWarning)

        lfc_df = pl.DataFrame(
            [get_gene_lfc(get_gene_table(gene_id=gene_id)) for gene_id in moderately_expressed]
        )

    lfc_table = mo.ui.table(
        lfc_df.sort(pl.col('log2fc').abs(), descending=True),
        selection='single',
        initial_selection=[0]
    )

    lfc_table
    ```

    </details>
    """)
    return


@app.cell
def _():
    # Your code goes here



    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Using the gene selected in `lfc_table`, look up its ID and pull that gene's per-sample data with `get_gene_table` to plot normalized counts (`count_norm`) against `strain`.

    **You should be able to select a gene in the table above and see the counts in the chart below.**

    <details>
    <summary>👀 Take a peek at the solution</summary>

    ```python
    gene_id = lfc_table.value.get_column('gene')[0]

    alt.Chart(get_gene_table(gene_id)).mark_point().encode(
        x=alt.X('count_norm'),
        y=alt.Y('strain')
    )
    ```

    </details>
    """)
    return


@app.cell
def _():
    # Your code goes here



    return


if __name__ == "__main__":
    app.run()
