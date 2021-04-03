# SCRIPT TO CREATE TABLE WITH DATA FROM SCORE BETWEEN miRNA AND IT'S TARGET GENE THAT IS PRESENT IN A SPECIFIC PATHWAY FROM KEGG
# DATA WILL BE USED IN CIRCOS PLOT WHICH IT'S CONFIG FILE HAS TO BE FOLLOWED TO GET A BETTER VISUALIZATION USING BLACK COLOR FROM
# MULTINDEX THAT HAS BEEN CREATED. ROW WITH COL COLORS HAS TO BE ALWAYS ON.

import pandas as pd

# Same path for mirDIP results, KEGG pathway files and building final table
# This function creates the dataframe from mirdip score with genes
def mirdipresults():
    mirsscoretargets = pd.read_csv(
        "results_nashonly/miRNA_nash/mirdip.csv", sep=",", skiprows=39
    )
    mirsscoretargets.drop(
        columns=["Uniprot", "Score Class", "Pseudogene"], inplace=True
    )
    return mirsscoretargets


mirdip = mirdipresults()

# Gets kegg data and sets a pathway to be searched for
def keggpathways(mirdip):
    keggpathways = pd.read_csv("results_nashonly/miRNA_nash/kegg_pathway.txt", sep="\t")
    keggpathways.drop(
        columns=[
            "Category",
            "Bonferroni",
            "Benjamini",
            "Pop Total",
            "%",
            "PValue",
            "List Total",
        ],
        inplace=True,
    )
    kegg_filtered = keggpathways[
        keggpathways.Term.str.contains("Sphingolipid signaling pathway")
    ]
    kegg_filtered.reset_index(inplace=True)
    genes_filtered = str(kegg_filtered.at[0, "Genes"]).replace(", ", "|")
    mirstargets = mirdip[mirdip["Gene Symbol"].str.contains(genes_filtered)]
    targetscan = mirstargets[mirstargets.Sources.str.contains("TargetScan")]
    targetscan.columns = ["gene", "mirna", "score", "nsources", "sources"]
    targetscan.sort_values(by="nsources", ascending=False, inplace=True)
    targetscan.reset_index(inplace=True, drop=True)
    # targetscan.drop(index=74, inplace=True)
    return targetscan


tablefilter = keggpathways(mirdip)

# Build the table that will be used in Circos Plot with multi index with black color for genes and let miRNAs be colored from circos
def buildtable(tablefiltered):
    pd.set_option("precision", 0)
    dataframe = tablefiltered.pivot(index="mirna", columns="gene", values="score")
    dataframe.fillna(0, inplace=True)
    df = (dataframe * 100).astype(int)
    df.loc["Total", :] = df.sum(axis=0)
    df.sort_values(by=["Total"], axis=1, ascending=False, inplace=True)
    df.drop(index=["Total"], inplace=True)
    # df.loc[(len(df.index.values))] = ['0,0,0 '] * len(df.columns)
    df.columns = pd.MultiIndex.from_product([["0,0,0"], df.columns])
    df2 = df.astype(int)
    # print(df2)
    df2.to_csv(
        "results_nashonly/miRNA_nash/data_circos_sphingo.txt",
        sep="\t",
    )


buildtable(tablefilter)