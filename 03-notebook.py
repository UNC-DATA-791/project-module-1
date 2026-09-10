import marimo

__generated_with = "0.24.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell
def _():
    import polars as pl
    import altair as alt
    alt.data_transformers.enable("vegafusion")
    return


@app.cell
def _():
    DATA_URLS = {
        'LogFoldChange': 'https://raw.githubusercontent.com/UNC-DATA-791/project-module-1/refs/heads/main/data/lfc_actual_vs_moderated.csv',
        'LFC_experimental_splits': 'https://raw.githubusercontent.com/UNC-DATA-791/project-module-1/refs/heads/main/data/lfc_split1_vs_split2.csv'
    }

    DATA_URLS
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## MLE versus MAP estimates of log2 fold change
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Use the URL from the `DATA_URLS` dictionary to read in the log fold change values for gene expression over the mouse strains in the Bottomly data set. The `key` is `'LogFoldChange'`.

    <details><summary>Hint</summary>
    ```python
    lfc = pl.read_csv(DATA_URLS['LogFoldChange'])
    lfc
    ```
    </details>
    """)
    return


@app.cell
def _():
    # Code goes here


    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    You are going to plot `baseMean` versus the log fold change estimate from the MLE method versus the empirical Bayes method. The MLE estimate is in the `lfc_mle` column and the "shrunk" estimate is in the `lfc_shrunk` column.

    Use the dropdown below to select the log fold change value you are visualizing in the chart.

    <details>
    <summary>Hint</summary>
    ```python
    # this is a chart for the MLE
    alt.Chart(lfc, width=600, height=300).mark_point(size=15, fill='steelblue', fillOpacity=0.2, strokeOpacity=0.5).encode(
        x=alt.X('baseMean').scale(type='log'),
        y=alt.Y('lfc_mle')  # change the argument here to refer to the dropdown value
    )
    ```
    </details>
    """)
    return


@app.cell
def _(mo):
    y_axis_col = mo.ui.dropdown(options=['lfc_mle', 'lfc_shrunk'], value='lfc_mle')
    y_axis_col
    return (y_axis_col,)


@app.cell
def _(y_axis_col):
    y_axis_col.value
    return


@app.cell
def _():
    # Code goes here


    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Comparing log2 fold change when the full data set is split in half.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Use the URL from the `DATA_URLS` dictionary to read in the log fold change values for gene expression calculated from two halves of the data set. The `key` is `'LFC_experimental_splits'`.

    <details>
    <summary>Hint</summary>
    ```python
    lfc_splits = pl.read_csv(DATA_URLS['LFC_experimental_splits'])
    lfc_splits
    ```
    </details>
    """)
    return


@app.cell
def _():
    # Code goes here


    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    Plot `split1` versus `split2`. These columns hold the log2 fold change estimates. The `method` column indicates if the values are estimated via the **MLE** or **MAP** method. I want you to use two scatter charts to display the comparison (`split1` values versus `split2`). You have several options:

    1. You can use `alt.hconcat` and put two individual charts as arguments.
    2. You can create each chart (one for MLE and one for MAP) in individual cells.
    3. You can encode `method` with the `column` channel (that's the method in your hint).

    <details>
    <summary>Hint</summary>
    ```python
    alt.Chart(lfc_splits).mark_point().encode(
        x=alt.X('split1'),
        y=alt.Y('split2'),
        column=alt.Column('method')
    ).resolve_scale(x="independent", y="independent")
    ```
    </details>
    """)
    return


@app.cell
def _():
    # Code goes here


    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    In a markdown cell below. Briefly describe why the MAP estimates are more consistent between experiments.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    *put your answer here*
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""

    """)
    return


if __name__ == "__main__":
    app.run()
