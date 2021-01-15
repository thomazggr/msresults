import pandas as pd
import numpy as np
from functools import reduce

pd.options.display.float_format = '{:,.4f}'.format

def gse10694():
    results = pd.read_csv('./tbu/GSE10694/GSE10694_results.txt', sep='\t', skiprows=[0], usecols=['logFC', 'miRNA_ID'])
    results.dropna(inplace=True)
    resultsfilter = results[results.miRNA_ID.str.contains('hsa-miR-378|hsa-miR-224|hsa-miR-221|hsa-miR-214|hsa-miR-99a|hsa-miR-148a|hsa-miR-34a|hsa-miR-155|hsa-miR-222|hsa-miR-93|hsa-miR-422a|hsa-miR-424|hsa-miR-106b|hsa-miR-200b|hsa-miR-181d')]
    resultsfilter.columns = ['logfc10', 'mir']
    #resultsfilter.set_index('mir', drop=True, inplace=True)
    return resultsfilter

gse10 = gse10694()

def gse49012():
    results = pd.read_csv('./tbu/GSE49012/GSE49012_w_logfc.txt', sep='\t', usecols=['logFC', 'ID'])
    results.dropna(inplace=True)
    resultsfilter = results[results.ID.str.contains('hsa-miR-378|hsa-miR-224|hsa-miR-221|hsa-miR-214|hsa-miR-99a|hsa-miR-148a|hsa-miR-34a|hsa-miR-155|hsa-miR-222|hsa-miR-93|hsa-miR-422a|hsa-miR-424|hsa-miR-106b|hsa-miR-200b|hsa-miR-181d')]
    resultreplace = resultsfilter['ID'].replace(to_replace=r'[*]', value='', regex=True)
    resultsfilter['ID'] = resultreplace
    resultsfilter.columns = ['mir','logfc49']
    resultsfilter.set_index('mir', drop=True, inplace=True)
    return resultsfilter

gse49 = gse49012()

def gse59492():
    results = pd.read_csv('./tbu/GSE59492/GSE59492_w_logfc.txt', sep='\t', usecols=['logFC', 'miRNA_ID_LIST'])
    results.dropna(inplace=True)
    resultsfilter = results[results.miRNA_ID_LIST.str.contains('hsa-miR-378|hsa-miR-224|hsa-miR-221|hsa-miR-214|hsa-miR-99a|hsa-miR-148a|hsa-miR-34a|hsa-miR-155|hsa-miR-222|hsa-miR-93|hsa-miR-422a|hsa-miR-424|hsa-miR-106b|hsa-miR-200b|hsa-miR-181d')]
    resultsfilter.columns = ['logfc59', 'mir']
    resultsfilter.set_index('mir', drop=True, inplace=True)
    return resultsfilter

gse59=gse59492()

def gse74618():
    results = pd.read_csv('./tbu/GSE74618/GSE74618_results.txt', sep='\t', usecols=['logFC', 'miRNA_ID_LIST'])
    results.dropna(inplace=True)
    resultsfilter = results[results.miRNA_ID_LIST.str.contains('hsa-miR-378|hsa-miR-224|hsa-miR-221|hsa-miR-214|hsa-miR-99a|hsa-miR-148a|hsa-miR-34a|hsa-miR-155|hsa-miR-222|hsa-miR-93|hsa-miR-422a|hsa-miR-424|hsa-miR-106b|hsa-miR-200b|hsa-miR-181d')]
    resultsfilter.columns = ['logfc74', 'mir']
    #resultsfilter.set_index('mir', drop=True, inplace=True)
    return resultsfilter

gse74=gse74618()
#print(gse74, gse10)
result = pd.concat([gse74, gse10, gse59, gse49], axis=0, sort=False)
'''result = gse10.merge(gse74, on=['mir'])
result2 = result.merge(gse59, on=['mir'])
result3 = result2.merge(gse49, on=['mir'])
result3.to_csv('result3.csv')'''

index='hsa-miR-221,hsa-miR-424,hsa-miR-222,hsa-miR-422a,hsa-miR-148a,hsa-miR-34a,hsa-miR-93,hsa-miR-106b,hsa-miR-155,hsa-miR-181d,hsa-miR-224,hsa-miR-214,hsa-miR-99a,hsa-miR-200b'
indexl = index.split(',')

GSE10694val = '2.7035,-1.8175,2.3893,-1.0743,-1.3629,1.4803,1.1582,1.3369,-1.0672,-0.6116,0.9426,0.7873,-0.8358,0.8722'
GSE10694l = GSE10694val.split(',')

GSE59492val = '1.07,0.02,1.02,-0.39299999999999996,-0.35100000000000003,0.948,-0.40399999999999997,0.126,1.33,1.3,1.59,0.878,0.115,1.98'
GSE59492l = GSE59492val.split(',')

GSE49012val = '-2.35817,-1.2412100000000001,-1.52788,-1.3582299999999998,-1.38529,-0.78051,-1.76216,-1.6145200000000002,-1.39993,-1.00824,-0.21122,-2.73584,-0.8829600000000001,-1.82435'
GSE49012l = GSE49012val.split(',')

GSE74618val = '1.87,-0.992,1.61,-2.59,-0.289,2.21,0.613,0.73,1.61,0.089,2.76,-2.76,-1.16,-0.354'
GSE74618l = GSE74618val.split(',')

data = {'GSE10694l':GSE10694l, 'GSE59492l':GSE59492l, 'GSE49012l':GSE49012l, 'GSE74618l':GSE74618l}
df = pd.DataFrame(index=indexl, data=data)
print(df)