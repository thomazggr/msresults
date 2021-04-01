import pandas as pd

gens = [
    "GADD45G",
    "FADS2",
    "P4HA1",
    "SOCS2",
    "FOSB",
    "IGFBP2",
    "TRHDE",
    "ME1",
]

# get data from datasets results
gs1 = pd.read_csv(
    "mRNa-analysis/GSE37031/results_converted.txt",
    sep="\t",
    usecols=["logFC", "SPOT_ID"],
)
gs1.columns = ["logFC", "ID"]
gs1 = gs1[["ID", "logFC"]]

gs2 = pd.read_csv(
    "mRNa-analysis/GSE48452/GSE48452.top.table.tsv",
    sep="\t",
    usecols=["logFC", "Gene.symbol"],
)
gs2.columns = ["logFC", "ID"]
gs2 = gs2[["ID", "logFC"]]

gs3 = pd.read_csv(
    "mRNa-analysis/GSE89632/results.txt",
    sep="\t",
    usecols=["logFC", "Symbol"],
    skiprows=1,
)
gs3.columns = ["logFC", "ID"]
gs3 = gs3[["ID", "logFC"]]

gs1x = gs1[gs1["ID"].isin(gens)].query("logFC >= 1 | logFC <= -1")
gs2x = gs2[gs2["ID"].isin(gens)].query("logFC >= 1 | logFC <= -1")
gs3x = gs3[gs3["ID"].isin(gens)].query("logFC >= 1 | logFC <= -1")

print(gs1x)
print(gs2x)
print(gs3x)

# gs1x.to_clipboard()
# gs2x.to_clipboard()
gs3x.to_clipboard()
