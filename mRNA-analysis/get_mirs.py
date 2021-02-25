import pandas as pd

mirdip = pd.read_csv(
    "./mRNA-analysis/mirdip_results.csv",
    skiprows=35,
    usecols=[
        "Gene Symbol",
        "Uniprot",
        "MicroRNA",
        "Integrated Score",
        "Number of Sources",
    ],
)
print(
    mirdip.groupby("MicroRNA")["Integrated Score"]
    .mean()
    .sort_values(ascending=False)
    .filter(
        regex="hsa-miR-221|hsa-miR-155|hsa-miR-222|hsa-miR-422a|hsa-miR-150|hsa-miR-378|hsa-miR-182|hsa-miR-21.*",
        axis="index",
    )
)
print(
    mirdip.MicroRNA.value_counts().filter(
        regex="hsa-miR-221|hsa-miR-155|hsa-miR-222|hsa-miR-422a|hsa-miR-150|hsa-miR-378|hsa-miR-182|hsa-miR-21.*",
        axis="index",
    )
)
