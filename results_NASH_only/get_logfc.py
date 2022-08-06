import pandas as pd
from scripts.reverse import rvs

mirs = [
    "hsa-miR-200b",
    "hsa-miR-483",
    "hsa-miR-155",
    "hsa-miR-940",
    "hsa-miR-222",
    "hsa-miR-221",
    "hsa-miR-181d",
    "hsa-miR-200a",
    "hsa-miR-193b",
    "hsa-miR-31",
    "hsa-miR-143",
    "hsa-miR-199b",
    "hsa-miR-378",
    "hsa-miR-150",
    "hsa-miR-182",
    "hsa-miR-10a",
    "hsa-miR-21",
    "hsa-miR-146b",
    "hsa-miR-200c",
    "hsa-miR-1234",
]

# get data from datasets results
gs1 = pd.read_csv(
    "datasets/GSE49012/GSE49012_w_logfc.txt", sep="\t", usecols=["ID", "logFC"]
)
gs2 = pd.read_csv(
    "datasets/GSE59492/GSE59492_w_logfc.txt",
    sep="\t",
    usecols=["logFC", "miRNA_ID_LIST"],
)
gs2.columns = ["logFC", "ID"]
gs2 = gs2[["ID", "logFC"]]

gs1x = gs1[gs1["ID"].isin(mirs)].query("logFC >= 1 | logFC <= -1")
gs2x = gs2[gs2["ID"].isin(mirs)].query("logFC >= 1 | logFC <= -1")

mir = []
logfc = []
for i, j in gs1x.iterrows():
    mir.append(j["ID"])
    logfc.append(rvs(j["logFC"]))

gs1x = pd.DataFrame(data=[mir, logfc])
gs1x = gs1x.T
gs1x.columns = ["ID", "GSE49012"]

mir2 = []
logfc2 = []
for i, j in gs2x.iterrows():
    mir.append(j["ID"])
    logfc.append(rvs(j["logFC"]))

gs2x = pd.DataFrame(data=[mir, logfc])
gs2x = gs2x.T
gs2x.columns = ["ID", "GSE59492"]


res = pd.merge(gs1x, gs2x, on="ID", how="inner")
print(res)