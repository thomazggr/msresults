import pandas as pd

dfresults = pd.read_csv(
    "./mRNA-analysis/GSE66676/GSE66676.top.table.tsv",
    sep="\t",
    usecols=["adj.P.Val", "P.Value", "logFC", "Gene.symbol"],
)
# dfresults.drop(columns=["t", "B", "Description", "ID", "GB_ACC"], inplace=True)
colunas = ["adjpval", "pval", "logfc", "id"]
dfresults.columns = colunas
print(dfresults.shape)

frames = [
    dfresults.query("pval <= 0.01 & logfc > 1"),
    dfresults.query("pval <= 0.01 & logfc < -1"),
]
results_de = pd.concat(frames)
results_de.reset_index(inplace=True, drop=True)
# results_de['id'].to_clipboard(index=False)
print(results_de)