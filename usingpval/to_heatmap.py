import pandas as pd
import GEOparse
from scripts.heatmap import generate_heat_map
import re

def generate_list(path, skip=False, keepid=False):
    mir = 'hsa-miR-221|hsa-miR-155|hsa-miR-222|hsa-miR-422a|hsa-miR-150|hsa-miR-378|hsa-miR-182|hsa-miR-21'
    if skip == False:
        complete_df = pd.read_csv(path, sep='\t')
    elif skip == True:
        complete_df = pd.read_csv(path, sep='\t', skiprows=[0])

    if keepid == False:
        complete_df.drop(columns=['t', 'B', 'SPOT_ID', 'ID'], inplace=True)
        colunas = ['adjpval', 'pval', 'logfc', 'id']
    elif keepid == True:
        complete_df.drop(columns=['t', 'B', 'SPOT_ID'], inplace=True)
        colunas = ['id','adjpval', 'pval', 'logfc']
    complete_df.columns = colunas
    filtra_hsa = complete_df[complete_df.id.str.contains('hsa', na=False)]
    filtra_p = filtra_hsa[filtra_hsa.id.str.contains(mir)].query('pval <= 0.05 & logfc < -1 | logfc > 1')
    lista = list(filtra_p['id'])
    return lista

lista59492 = generate_list(path='./datasets/GSE59492/GSE59492_w_logfc.txt')
lista74618 = generate_list(path='./datasets/GSE74618/GSE74618_results.txt')
lista10694 = generate_list(path='./datasets/GSE10694/GSE10694_results.txt', skip=True)
lista49012 = generate_list(path='./datasets/GSE49012/GSE49012_w_logfc.txt', keepid=True)


datasets = ['GSE59492', 'GSE74618', 'GSE10694', 'GSE49012']

def run_heatmaps(dataset):
    tabela = GEOparse.get_GEO(geo=dataset, destdir='./junk/', silent=True)
    valores = tabela.pivot_samples('VALUE')
    if dataset == 'GSE59492':
        lista = lista59492
    elif dataset == 'GSE74618':
        lista = lista74618
    elif dataset == 'GSE10694':
        lista = lista10694
    elif dataset == 'GSE49012':
        lista = lista49012
    valores.reset_index(inplace=True)
    print(valores)
    #print(valores[valores['ID_REF'].isin(lista)])
    #print(valores.query('ID_REF in @lista'))
    #resultado = re.findall('.', samples)

for ds in datasets:
    run_heatmaps(ds)
