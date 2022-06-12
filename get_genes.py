import pandas as pd
# from operator import is_not
# from functools import partial
# import rpy2.robjects as ro
# from rpy2.robjects.packages import importr
# from rpy2.robjects import pandas2ri
# from rpy2 import robjects

# from rpy2.robjects.conversion import localconverter
# biomart = importr("biomaRt")
# dplyr = importr("dplyr")
# rprint = robjects.globalenv.find("print")
# bmart = robjects.r("""
# (biomaRt::getBM(attributes = c(
# 'external_gene_name',
# 'description',
# 'entrezgene_id'),
# mart = (biomaRt::useMart(biomart = 'ensembl',
# dataset = 'mmusculus_gene_ensembl',
# host = 'https://dec2021.archive.ensembl.org')
# )
# )) %>% mutate_all(as.character)
# """)
#
# with localconverter(ro.default_converter + pandas2ri.converter):
# pd_from_r_df = ro.conversion.rpy2py(bmart)

FILTER_STATE = "pval < 0.05 and (logfc < -1 or logfc > 1)"
print("FILTERING USED FOR THIS RUN")
print(FILTER_STATE)
print()
print()


gse_33814 = pd.read_csv("automated_results/GSE33814/GSE33814.tsv", sep="\t")
# print(gse_33814.columns)
gse_33814_cols = gse_33814[["Gene.symbol", "logFC", "P.Value", "adj.P.Val"]]

gse_33814_cols.columns = ["gene", "logfc", "pval", "adjpval"]
gse_33814_filtered = gse_33814_cols.query(FILTER_STATE)
# print(gse_33814_filtered.tail())
print("GSE33814 GENE LIST")
print([x for x in list(gse_33814_filtered["gene"]) if str(x) != 'nan'])
print()
print()


gse_37031 = pd.read_csv("automated_results/GSE37031/GSE37031.tsv", sep="\t")
# print(gse_37031.columns)
gse_37031_cols = gse_37031[["SPOT_ID", "logFC", "P.Value", "adj.P.Val"]]

gse_37031_cols.columns = ["gene", "logfc", "pval", "adjpval"]

gse_37031_filtered = gse_37031_cols.query(FILTER_STATE)

# gse_37031_filtered["gene"] = gse_37031_filtered["gene"].to_string()
# pd_from_r_df["entrezgene_id"] = pd_from_r_df["entrezgene_id"].to_string()

# df_final = gse_37031_filtered.merge(
#    pd_from_r_df, how="left", left_on="gene", right_on="entrezgene_id")

# print(list(df_final["description"]))
print("GSE37031 GENE LIST")
# print(list(gse_37031_filtered["gene"]))  # ENTREZ_GENE_ID
print([x for x in list(gse_37031_filtered["gene"]) if str(x) != 'nan'])
# ENTREZ_GENE_ID
print()
print()

gse_48452 = pd.read_csv("automated_results/GSE48452/GSE48452.tsv", sep="\t")
# print(gse_48452.columns)
gse_48452_cols = gse_48452[["Gene.symbol", "logFC", "P.Value", "adj.P.Val"]]

gse_48452_cols.columns = ["gene", "logfc", "pval", "adjpval"]

gse_48452_filtered = gse_48452_cols.query(FILTER_STATE)
print("GSE48452 GENE LIST")
# print(list(gse_48452_filtered["gene"]))
print([x for x in list(gse_48452_filtered["gene"]) if str(x) != 'nan'])
print()
print()


gse_63067 = pd.read_csv("automated_results/GSE63067/GSE63067.tsv", sep="\t")
# print(gse_63067.columns)
gse_63067_cols = gse_63067[["Gene.symbol", "logFC", "P.Value", "adj.P.Val"]]

gse_63067_cols.columns = ["gene", "logfc", "pval", "adjpval"]

gse_63067_filtered = gse_63067_cols.query(FILTER_STATE)
print("GSE63067 GENE LIST")
# print(list(gse_63067_filtered["gene"]))
print([x for x in list(gse_63067_filtered["gene"]) if str(x) != 'nan'])
print()
print()

gse_89632 = pd.read_csv("automated_results/GSE89632/GSE89632.tsv", sep="\t")
# print(gse_89632.columns)
gse_89632_cols = gse_89632[["Symbol", "logFC", "P.Value", "adj.P.Val"]]

gse_89632_cols.columns = ["gene", "logfc", "pval", "adjpval"]

gse_89632_filtered = gse_89632_cols.query(FILTER_STATE)
print("GSE89632 GENE LIST")
# print(list(gse_89632_filtered["gene"]))
print([x for x in list(gse_89632_filtered["gene"]) if str(x) != 'nan'])
print()
print()
