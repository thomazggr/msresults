import pandas as pd

datasets = {
    "GSE33814": "Gene.ID",
    "GSE37031": "SPOT_ID", 
    "GSE48452": "Gene.ID", 
    "GSE63067": "Gene.ID", 
    "GSE89632": "Entrez_Gene_ID"}

FILTER_STATE = "`adj.P.Val` <= 0.05 and (`logFC` < -1 or `logFC` > 1)"

GENES_TO_FILTER = ["51280", "3485", "10537", "10912", "7057", "2353", 
                    "3484", "80303", "8835", "8660", "83716", "55908", "29947"]

for dataset, column in zip(datasets.keys(), datasets.values()):
    initial_results = pd.read_csv(f"{dataset}/{dataset}.tsv", sep="\t")
    initial_results.dropna(subset=[column], inplace=True)
    initial_results[column] = initial_results[column].astype(str)
    initial_results = initial_results[[column, "logFC", "P.Value", "adj.P.Val"]]
    results_filtered = (
        initial_results
        .query(FILTER_STATE)
        # .drop_duplicates(
            # subset=[column],
            # keep="first"
        # )
        .rename(
            columns={
                column:"Gene.ID"
            }
        )
    )
    gene_string = "|".join(GENES_TO_FILTER)
    results_filtered = results_filtered[results_filtered["Gene.ID"].str.contains(gene_string)]
    results_filtered["Gene.ID"] = results_filtered["Gene.ID"].str.replace(".0", "", regex=False)

    print(dataset)
    print(results_filtered.reset_index(drop=True))
    print("")
    print("")

