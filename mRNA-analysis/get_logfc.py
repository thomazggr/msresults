import pandas as pd
from scripts.reverse import rvs

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

mir1 = []
logfc1 = []
for i, j in gs1x.iterrows():
    mir1.append(j["ID"])
    logfc1.append(rvs(j["logFC"]))

gs1x = pd.DataFrame(data=[mir1, logfc1])
gs1x = gs1x.T
gs1x.columns = ["ID", "GSE37031"]

mir2 = []
logfc2 = []
for i, j in gs2x.iterrows():
    mir2.append(j["ID"])
    logfc2.append(rvs(j["logFC"]))

gs2x = pd.DataFrame(data=[mir2, logfc2])
gs2x = gs2x.T
gs2x.columns = ["ID", "GSE48452"]

mir3 = []
logfc3 = []
for i, j in gs3x.iterrows():
    mir3.append(j["ID"])
    logfc3.append(rvs(j["logFC"]))

gs3x = pd.DataFrame(data=[mir3, logfc3])
gs3x = gs3x.T
gs3x.columns = ["ID", "GSE89632"]


res = pd.merge(gs1x, gs2x, on="ID", how="inner")

res = pd.merge(res, gs3x, on="ID", how="inner")

print(gs1x, gs2x, gs3x)