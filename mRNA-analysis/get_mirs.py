import pandas as pd

mirdip = pd.read_csv(
    "results_nashonly/mRNA_nash/mirdip.csv",
    skiprows=20,
    usecols=[
        "Gene Symbol",
        "Uniprot",
        "MicroRNA",
        "Integrated Score",
        "Number of Sources",
    ],
)
"""print(
    mirdip.groupby("MicroRNA")["Integrated Score"]
    .mean()
    .sort_values(ascending=False)
    .filter(
        regex="hsa-miR-200b|hsa-miR-483|hsa-miR-155|hsa-miR-940|hsa-miR-222|hsa-miR-221|\
            hsa-miR-181d|hsa-miR-200a|hsa-miR-193b|hsa-miR-31|hsa-miR-143|hsa-miR-199b|\
                hsa-miR-378|hsa-miR-150|hsa-miR-182|hsa-miR-10a|hsa-miR-21|hsa-miR-146b|\
                    hsa-miR-200c|hsa-miR-1234.*",
        axis="index",
    )
)
print(
    mirdip.MicroRNA.value_counts().filter(
        regex="hsa-miR-200b|hsa-miR-483|hsa-miR-155|hsa-miR-940|hsa-miR-222|hsa-miR-221|\
            hsa-miR-181d|hsa-miR-200a|hsa-miR-193b|hsa-miR-31|hsa-miR-143|hsa-miR-199b|\
                hsa-miR-378|hsa-miR-150|hsa-miR-182|hsa-miR-10a|hsa-miR-21|hsa-miR-146b|\
                    hsa-miR-200c|hsa-miR-1234.*",
        axis="index",
    )
)"""
mirdip = mirdip[
    mirdip.MicroRNA.str.contains(
        "hsa-miR-200b|hsa-miR-483|hsa-miR-155|hsa-miR-940|hsa-miR-222|hsa-miR-221|\
            hsa-miR-181d|hsa-miR-200a|hsa-miR-193b|hsa-miR-31|hsa-miR-143|hsa-miR-199b|\
                hsa-miR-378|hsa-miR-150|hsa-miR-182|hsa-miR-10a|hsa-miR-21|hsa-miR-146b|\
                    hsa-miR-200c|hsa-miR-1234"
    )
]

mirdip[["Gene Symbol", "MicroRNA", "Integrated Score"]].to_clipboard(
    index=False, header=False
)
