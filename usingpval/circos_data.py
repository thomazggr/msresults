import pandas as pd



#Same path for mirDIP results, KEGG pathway files and building final table

def mirdipresults():
    mirsscoretargets = pd.read_csv('./usingpval/mirdip/mirdip_results.csv', sep=',', skiprows=27)
    mirsscoretargets.drop(columns=['Uniprot', 'Score Class', 'Pseudogene'], inplace=True)
    return mirsscoretargets

mirdip = mirdipresults()

def keggpathways(mirdip):
    keggpathways = pd.read_csv('./usingpval/keggo/kegg_david.txt', sep='\t')
    keggpathways.drop(columns=['Category', 'Bonferroni', 'Benjamini', 'Pop Total', '%', 'PValue', 'List Total'], inplace=True)
    kegg_filtered = keggpathways[keggpathways.Term.str.contains('Sphingolipid signaling pathway')]
    kegg_filtered.reset_index(inplace=True)
    genes_filtered = str(kegg_filtered.at[0, 'Genes']).replace(', ','|')
    mirstargets = mirdip[mirdip['Gene Symbol'].str.contains(genes_filtered)]
    targetscan = mirstargets[mirstargets.Sources.str.contains('TargetScan')]
    targetscan.columns = ['gene', 'mirna', 'score', 'nsources', 'sources']
    targetscan.sort_values(by='nsources', ascending=False, inplace=True)
    targetscan.reset_index(inplace=True, drop=True)
    #targetscan.drop(index=74, inplace=True)
    return targetscan

tablefilter = keggpathways(mirdip)

def buildtable(tablefiltered):
    pd.set_option('precision', 0)
    dataframe = tablefiltered.pivot(index='mirna', columns='gene', values='score')
    dataframe.fillna(0, inplace=True)
    df = (dataframe*100).astype(int)
    df.loc['Total',:] = (df.sum(axis=0))
    df.sort_values(by=['Total'], axis=1, ascending=False, inplace=True)
    df.drop(index=['Total'], inplace=True)
    #df.loc[(len(df.index.values))] = ['0,0,0 '] * len(df.columns)
    df.columns = pd.MultiIndex.from_product([df.columns, ['0,0,0']])
    df2 = df.astype(int)
    df2.to_csv('./usingpval/data_circos.txt', sep='\t')

buildtable(tablefilter)