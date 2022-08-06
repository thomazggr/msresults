import pandas as pd

# mirs from NASH datasets
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
    "hsa-miR-141",
    "hsa-miR-378",
    "hsa-miR-150",
    "hsa-miR-182",
    "hsa-miR-10a",
    "hsa-miR-21",
    "hsa-miR-146b",
    "hsa-miR-200c",
    "hsa-miR-551b",
    "hsa-miR-1234",
]

# mirdip results from genes from NASH datasets
genes_mirs = pd.read_csv(
    "results_nashonly/mRNA_nash/mirdip_n.csv",
    skiprows=35,
    usecols=["Gene Symbol", "MicroRNA", "Integrated Score", "Number of Sources"],
)
mirslst = []
for it in genes_mirs.MicroRNA:
    if it.endswith("5p") | it.endswith("3p"):
        it = it[:-3]
        mirslst.append(it)
    else:
        mirslst.append(it)
genes_mirs.MicroRNA = mirslst
filtrado = genes_mirs[genes_mirs["MicroRNA"].isin(mirs)]
filtrado.columns = ["gene", "mir", "score", "sources"]
# filtrado2 = filtrado.query("sources >= 14")
filtrado.reset_index(drop=True, inplace=True)

filtrado["score"].to_clipboard(index=False, header=False)
