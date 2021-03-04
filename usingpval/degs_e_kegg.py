import pandas as pd

kegg = pd.read_csv("usingpval/keggo/kegg_david.txt", sep="\t")

list_genes = "NR4A2|NFIL3|CEBPD|RCAN1|TP53I3|THBS1|FOSB|ME1|GOLM1|EFHD1|TMEM45B|GADD45G|DNMT3L|IGFBP1|FOS|SOCS2|CRISPLD2"

r = kegg[kegg.Genes.str.contains(list_genes)]

list_genes2 = [
    "NR4A2",
    "NFIL3",
    "CEBPD",
    "RCAN1",
    "TP53I3",
    "THBS1",
    "FOSB",
    "ME1",
    "GOLM1",
    "EFHD1",
    "TMEM45B",
    "GADD45G",
    "DNMT3L",
    "IGFBP1",
    "FOS",
    "SOCS2",
    "CRISPLD2",
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
