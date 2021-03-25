import pandas as pd

dfresults = pd.read_csv(
    "./mRNA-analysis/GSE151158/GSE151158-CTvsNASH.top.table.tsv",
    sep="\t",
)
# dfresults.drop(columns=["t", "B", "Description", "ID", "GB_ACC"], inplace=True)
colunas = ["id", "adjpval", "pval", "logfc"]
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

# results = pd.read_csv("./mRNA-analysis/GSE37031/genes_david_conversion.txt", sep="\t")
# results.To.to_clipboard(index=False)