import pandas as pd

mirdip = pd.read_csv(
    "mirDIP_E_2022_06_13_16 08 41.txt",
    skiprows=117,
    sep="\t"
)

three_intersect_genes = """'FOS', 'IGFBP1', 'IGFBP2', 'THBS1', 'IRS2', 'SOCS2', 'UBD', 'GADD45G', 'DNMT3L', 'GOLM1', 'ANGPTL8', 'EFHD1', 'CRISPLD2'"""

# mirdip = mirdip.query(f"`Gene Symbol` in ({three_intersect_genes})")

summ = mirdip["MicroRNA"].value_counts()#.sort_values("count", ascending=False)

print(summ.head(40))