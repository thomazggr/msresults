import pandas as pd

dfresults = pd.read_csv(
    "./mRNA-analysis/GSE83452/GSE83452-diet-NASHBASEvsNASHFOLLOW.top.table.tsv",
    sep="\t",
    usecols=["ID", "adj.P.Val", "P.Value", "logFC", "SPOT_ID"],
)
# dfresults.drop(columns=["t", "B", "Description", "ID", "GB_ACC"], inplace=True)
colunas = ["id", "adjpval", "pval", "logfc", "acc"]
dfresults.columns = colunas
print(dfresults.shape)
frames = [
    dfresults.query("pval <= 0.01 & logfc > 1"),
    dfresults.query("pval <= 0.01 & logfc < -1"),
]
results_de = pd.concat(frames)
results_de.reset_index(inplace=True, drop=True)
results_de["acc"].to_clipboard(index=False)
print(results_de)

# results = pd.read_csv("./mRNA-analysis/GSE37031/genes_david_conversion.txt", sep="\t")
# results.To.to_clipboard(index=False)