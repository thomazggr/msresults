from bioinfokit import visuz
import dash_bio as dashbio
import pandas as pd
import dash_core_components as dcc

# df = pd.read_csv("datasets/GSE10694/GSE10694_results.txt", sep="\t", skiprows=[0])
# df = pd.read_csv("datasets/GSE49012/GSE49012_w_logfc.txt", sep="\t")
# df = pd.read_csv("datasets/GSE59492/GSE59492_w_logfc.txt", sep="\t")
df = pd.read_csv("datasets/GSE74618/GSE74618_results.txt", sep="\t")

component = dcc.Graph(
    figure=dashbio.VolcanoPlot(
        dataframe=df,
        gene="miRNA_ID_LIST",
        effect_size="logFC",
        p="P.Value",
        snp="miRNA_ID_LIST",
    )
)

# df = pd.read_csv('', sep='\t', skiprows=[0])

"""visuz.gene_exp.volcano(
    df=df,
    lfc="logFC",
    pv="P.Value",
    show=True,
    plotlegend=True,
    legendpos="upper right",
    color=("#E10600FF", "grey", "#00239CFF"),
    sign_line=False,
)"""
