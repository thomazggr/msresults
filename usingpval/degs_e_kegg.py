import pandas as pd

kegg = pd.read_csv("results_nashonly/miRNA_nash/kegg_pathway.txt", sep="\t")

list_genes = "GADD45G|FADS2|P4HA1|SOCS2|FOSB|IGFBP2|TRHDE|ME1"

r = kegg[kegg.Genes.str.contains(list_genes)]

list_genes2 = [
    "GADD45G",
    "FADS2",
    "P4HA1",
    "SOCS2",
    "FOSB",
    "IGFBP2",
    "TRHDE",
    "ME1",
]

# r2 = kegg[~kegg.Genes.isin(list_genes2)]

r2 = r[["Term", "Genes", "PValue", "Fold Enrichment"]]
r2["Term"].replace(to_replace="^.*?:", value="", regex=True, inplace=True)
r2.reset_index(inplace=True, drop=True)
gens = r2["Genes"]
idx = 0

for row in gens:
    a = row.split(", ")
    rs = [x for x in a if any(word in x for word in list_genes2)]
    r2.at[idx, "Genes"] = str(rs)
    idx += 1

print(r2)
