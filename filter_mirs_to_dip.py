import pandas as pd

def gse10694():
    complete_df = pd.read_csv('./tbu/GSE10694/GSE10694_results.txt', sep='\t', skiprows=[0])
    complete_df.drop(columns=['t', 'B', 'SPOT_ID', 'ID'], inplace=True)
    colunas = ['adjpval', 'pval', 'logfc', 'id']
    complete_df.columns = colunas
    filtra_hsa = complete_df[complete_df.id.str.contains('hsa', na=False)]
    filtra_p = filtra_hsa.query('adjpval <= 0.05')
    filtra_p['id'].to_csv('./tbu/GSE10694/mirs.txt', index=False)

#gse10694()

def gse59492():
    complete_df = pd.read_csv('./tbu/GSE59492/GSE59492_w_logfc.txt', sep='\t')
    complete_df.drop(columns=['t', 'B', 'SPOT_ID', 'ID'], inplace=True)
    colunas = ['adjpval', 'pval', 'logfc', 'id']
    complete_df.columns = colunas
    filtra_hsa = complete_df[complete_df.id.str.contains('hsa', na=False)]
    filtra_p = filtra_hsa.query('adjpval <= 0.05')
    filtra_p['id'].to_csv('./tbu/GSE59492/mirs.txt', index=False)

gse59492()