import pandas as pd

genes = pd.read_csv("results_nashonly/miRNA_nash/mirdip.csv", sep=",", skiprows=39)
genes.drop(columns=["Pseudogene", "Score Class"], inplace=True)
colunas = ["gene", "uniprot", "mir", "score", "nsources", "sources"]
genes.columns = colunas
genes_sources = genes.query("nsources >= 14")

a = pd.DataFrame(genes_sources["gene"].unique())
a.columns = ["genes"]
a["genes"].to_clipboard(index=False)
