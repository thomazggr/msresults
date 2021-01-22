import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import re

def kegg():
    dados_kegg = pd.read_csv('./usingpval/keggo/kegg_david.txt', sep='\t', usecols=['Term', 'PValue'])
    dados_kegg.set_index(['Term'], drop=True, inplace=True)
    dados_kegg = dados_kegg.head(10)
    dados_kegg = np.log10(dados_kegg).abs()
    dados_kegg.reset_index(inplace=True)
    dados_kegg.replace(to_replace='^.*?:', value= '', regex=True, inplace=True)
    dados_kegg.sort_values(by='PValue', axis='index', ascending=True, inplace=True)
    ax = sns.barplot(x='PValue',y='Term', orient='h', data=dados_kegg, color='#771919')
    ax.set_ylabel('Termo KEGG')
    ax.set_xlabel('-log10(Valor de P)')
    plt.show()
#kegg()

def go():
    dados_go = pd.read_csv('./usingpval/keggo/go_david.txt', sep='\t', usecols=['Term', 'PValue'])
    dados_go.set_index(['Term'], drop=True, inplace=True)
    dados_go = dados_go.head(10)
    dados_go = np.log10(dados_go).abs()
    dados_go.reset_index(inplace=True)
    dados_go.replace(to_replace='^.*?~', value= '', regex=True, inplace=True)
    dados_go.sort_values(by='PValue', axis='index', ascending=True, inplace=True)
    ax = sns.barplot(x='PValue',y='Term', orient='h', data=dados_go, color='#007CAB')
    ax.set_ylabel('Termo Gene Ontology')
    ax.set_xlabel('-log10(Valor de P)')
    plt.show()
go()