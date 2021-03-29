import pandas as pd
import pyperclip

dfresults = pd.read_csv(
    "./mRNA-analysis/GSE48452/GSE48452.top.table.tsv",
    sep="\t",
    usecols=["adj.P.Val", "P.Value", "logFC", "Gene.symbol"],
)
# dfresults.drop(columns=["t", "B", "Description", "ID", "GB_ACC"], inplace=True)
colunas = ["adjpval", "pval", "logfc", "id"]
dfresults.columns = colunas

frames = [
    dfresults.query("pval <= 0.01 & logfc > 1"),
    dfresults.query("pval <= 0.01 & logfc < -1"),
]
results_de = pd.concat(frames)
results_de.reset_index(inplace=True, drop=True)
results_de.dropna(inplace=True, axis="index")
idlist = results_de["id"]
res = []
for it in idlist:
    res.extend(it.split("///"))
pyperclip.copy("\n".join(res))
# results_de['id'].to_clipboard(index=False)
# print(results_de)
