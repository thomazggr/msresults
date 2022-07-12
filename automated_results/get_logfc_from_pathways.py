import os

import pandas as pd

kegg_results = pd.read_csv("padj005_logfc1/results__kegg_full_table.txt", sep="\t")

pathways_possible = kegg_results[["Description", "geneID"]]
print(pathways_possible)

print("\nNumber of pathway to read:")
index_needed = int(input(">>> "))

try:
    pathway = pathways_possible.iloc[index_needed]
    print("")
    print(pathway["Description"])
    print("")
    print(pathway["geneID"].split("/"))
    list_gene_ids = pathway["geneID"].split("/")
except IndexError as ie:
    print(f"ERROR >>> No pathway found for the index {index_needed}")
    quit()

datasets = {
    "GSE33814": "Gene.ID",
    "GSE37031": "SPOT_ID", 
    "GSE48452": "Gene.ID", 
    "GSE63067": "Gene.ID", 
    "GSE89632": "Entrez_Gene_ID"}

FILTER_STATE = "`adj.P.Val` <= 0.05 and (`logFC` < -1 or `logFC` > 1)"

resulted_dataframe = pd.DataFrame(
    data = {
        "Gene.ID":[],
        "logFC":[],
        "P.Value":[],
        "adj.P.Val":[],
        "Dataset":[]
    }
)

for dataset, column in zip(datasets.keys(), datasets.values()):
    ds_results = pd.read_csv(f"{dataset}/{dataset}.tsv", sep="\t")
    ds_results.dropna(subset=[column], inplace=True)
    ds_results[column] = ds_results[column].astype(str)
    string_ids = "|".join(list_gene_ids)
    
    ds_results_filtered = ds_results[ds_results[column].str.contains(string_ids)]
    ds_results_filtered = ds_results_filtered[[column, "logFC", "P.Value", "adj.P.Val"]]
    ds_final = (
        ds_results_filtered
        .query(FILTER_STATE)
        .drop_duplicates(
            subset=[column],
            keep="first"
        )
        .rename(
            columns={
                column:"Gene.ID"
            }
        )
    )

    ds_final["Gene.ID"] = ds_final["Gene.ID"].str.replace(".0", "", regex=False)
    ds_final["Dataset"] = dataset

    resulted_dataframe = pd.concat([resulted_dataframe, ds_final])

    print("")
    print("")
    print(dataset)
    print("")
    print(ds_final)
    print("")
    print("")

df_dict_stan = {
    "GSE33814" :{},
    "GSE37031" :{},
    "GSE89632" :{}
}

for row in resulted_dataframe.itertuples():
    df_dict_stan[row[5]][row[1]] = row[2]

print("DATASET FINAL")
print("")
# print(resulted_dataframe)
print(pd.DataFrame(df_dict_stan).transpose())
print("")
print("")

to_R_script = ",".join(list_gene_ids)
os.system(f"Rscript get_gene.R --organism hsa --geneids {to_R_script}")