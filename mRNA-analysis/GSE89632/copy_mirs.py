import pandas as pd

dfresults = pd.read_csv('./mRNA-analysis/GSE89632/results.txt', sep='\t', skiprows=1)
dfresults.drop(columns=['t', 'B', 'ID'], inplace=True)
colunas = ['adjpval', 'pval', 'logfc', 'id']
dfresults.columns = colunas
frames = [dfresults.query("adjpval <= 0.05 & logfc > 1"), dfresults.query("adjpval <= 0.05 & logfc < -1")]
results_de = pd.concat(frames)
results_de.reset_index(inplace=True, drop=True)
results_de['id'].to_clipboard(index=False)
