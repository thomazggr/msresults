import pandas as pd

dfresults = pd.read_csv('./mRNA-analysis/GSE37031/results.txt', sep='\t', skiprows=1)
dfresults.drop(columns=['t', 'B', 'Description', 'ID'], inplace=True)
colunas = ['adjpval', 'pval', 'logfc', 'id']
dfresults.columns = colunas
frames = [dfresults.query("adjpval <= 0.05 & logfc > 1"), dfresults.query("adjpval <= 0.05 & logfc < -1")]
results_de = pd.concat(frames)
results_de.reset_index(inplace=True, drop=True)
#results_de['id'].to_clipboard(index=False)

results = pd.read_csv('./mRNA-analysis/GSE37031/genes_david_conversion.txt', sep='\t')
results.To.to_clipboard(index=False)
