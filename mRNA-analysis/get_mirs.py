import pandas as pd

mirdip = pd.read_csv(
    "./mRNA-analysis/mirdip_results_2de3.csv",
    skiprows=156,
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
    .head(10)
)
print(mirdip.MicroRNA.value_counts().head(10))
