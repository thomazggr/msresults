import pandas as pd

ans = input("Print any miR to check? Y/N\n").upper()

if ans == "N":
    confirm = input(
        "If there is no print it will continue to write/overwrite the miR file. \n Do you want to continue? Y/N\n"
    ).upper()
    if confirm == "Y":
        pass
    else:
        exit()
elif ans == "Y":
    pass
else:
    exit()

mir = "miR-21"


def gse10694():
    complete_df = pd.read_csv(
        "./datasets/GSE10694/GSE10694_results.txt", sep="\t", skiprows=[0]
    )
    complete_df.drop(columns=["t", "B", "SPOT_ID", "ID"], inplace=True)
    colunas = ["adjpval", "pval", "logfc", "id"]
    complete_df.columns = colunas
    filtra_hsa = complete_df[complete_df.id.str.contains("hsa", na=False)]
    if ans == "Y":
        filtra_p = filtra_hsa[filtra_hsa.id.str.contains(mir)]
        print(filtra_p)
    elif ans == "N":
        filtra_p = filtra_hsa.query("pval <= 0.05 & logfc < -1 | logfc > 1")
        filtra_p["id"].to_csv("./tbu/GSE10694/mirs-pval.txt", index=False)


# gse10694()


def gse59492():
    complete_df = pd.read_csv("./datasets/GSE59492/GSE59492_w_logfc.txt", sep="\t")
    complete_df.drop(columns=["t", "B", "SPOT_ID", "ID"], inplace=True)
    colunas = ["adjpval", "pval", "logfc", "id"]
    complete_df.columns = colunas
    filtra_hsa = complete_df[complete_df.id.str.contains("hsa", na=False)]
    if ans == "Y":
        filtra_p = filtra_hsa[filtra_hsa.id.str.contains(mir)]
        print(filtra_p)
    elif ans == "N":
        filtra_p = filtra_hsa.query("pval <= 0.01 & logfc < -1 | logfc > 1")
        filtra_p["id"].to_csv("./datasets/GSE59492/mirs-pval.txt", index=False)


gse59492()


def gse74618():
    complete_df = pd.read_csv("./datasets/GSE74618/GSE74618_results.txt", sep="\t")
    complete_df.drop(columns=["t", "B", "SPOT_ID", "ID"], inplace=True)
    colunas = ["adjpval", "pval", "logfc", "id"]
    complete_df.columns = colunas
    filtra_hsa = complete_df[complete_df.id.str.contains("hsa", na=False)]
    if ans == "Y":
        filtra_p = filtra_hsa[filtra_hsa.id.str.contains(mir)]
        print(filtra_p)
    elif ans == "N":
        filtra_p = filtra_hsa.query("pval <= 0.05 & logfc < -1 | logfc > 1")
        filtra_p["id"].to_csv("./tbu/GSE74618/mirs-pval.txt", index=False)


# gse74618()


def gse49012():
    complete_df = pd.read_csv("./datasets/GSE49012/GSE49012_w_logfc.txt", sep="\t")
    complete_df.drop(columns=["t", "B", "SPOT_ID"], inplace=True)
    colunas = ["id", "adjpval", "pval", "logfc"]
    complete_df.columns = colunas
    filtra_hsa = complete_df[complete_df.id.str.contains("hsa", na=False)]
    if ans == "Y":
        filtra_p = filtra_hsa[filtra_hsa.id.str.contains(mir)]
        print(filtra_p)
    elif ans == "N":
        filtra_p = filtra_hsa.query("pval <= 0.01 & logfc < -1 | logfc > 1")
        filtra_p["id"].to_csv("./datasets/GSE49012/mirs-pval.txt", index=False)


gse49012()