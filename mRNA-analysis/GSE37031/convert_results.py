import pandas as pd

init_results = pd.read_csv("mRNa-analysis/GSE37031/results.txt", sep="\t", skiprows=1)

converted_genes = pd.read_csv(
    "mRNa-analysis/GSE37031/converted_new.txt", sep="\t", usecols=["From", "To"]
)

dictgens = ((row["From"], row["To"]) for it, row in converted_genes.iterrows())
dictgens = dict(dictgens)
init_results["SPOT_ID"] = init_results["SPOT_ID"].map(dictgens)
init_results.to_csv(
    "mRNa-analysis/GSE37031/results_converted.txt", sep="\t", index=False
)
