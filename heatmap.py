import pandas as pd
import GEOparse
import scipy
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import re
from scripts.heatmap import generate_heat_map

def gse10694():
    tabela = GEOparse.get_GEO(geo='GSE59492', destdir='./tbu/', silent=True)
    valores = tabela.pivot_samples('VALUE')
    return valores
tabela = gse10694()

def hm(table):
    # id_and_mir = pd.read_csv('./tbu/GSE10694/GSE10694_results.txt', sep='\t', skiprows=[0])
    # id_and_mir.set_index('ID', inplace=True, drop=True)
    # id_and_mir.sort_index(inplace=True)
    # table['newindex'] = id_and_mir['miRNA_ID']
    # table.set_index('newindex', inplace=True, drop=True)
    filterhsa = table.filter(regex='^hsa', axis=0)
    filtermirs = filterhsa.filter(regex='hsa-miR-378|hsa-miR-224|hsa-miR-221|hsa-miR-214|hsa-miR-99a|hsa-miR-148a|hsa-miR-34a|hsa-miR-155|hsa-miR-222|hsa-miR-93|hsa-miR-422a|hsa-miR-424|hsa-miR-106b|hsa-miR-200b|hsa-miR-181d.*', axis=0)
    samples = "XXXXXXXXXXXXXXXXXXXXXX11111000000"
    nindex = splitamostras(samples=samples)
    filtern = filtermirs.T
    filtern['NIND'] = nindex
    filtern.set_index(filtern['NIND'], inplace=True, drop=True)
    filtern.drop(columns= ['NIND'], index=['X'], inplace=True)
    filtern.sort_index(inplace=True)
    tabt = filtern.T
    tablog = np.log2(tabt)
    tablog.drop(index=['hsa-miR-939_st','hsa-miR-938_st','hsa-miR-937_st','hsa-miR-936_st','hsa-miR-935_st','hsa-miR-934_st','hsa-miR-933_st'], inplace=True)
    generate_heat_map(tablog)

def splitamostras(samples):
       resultado = re.findall('.', samples)
       return resultado

hm(table=tabela)