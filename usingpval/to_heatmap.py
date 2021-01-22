import pandas as pd
import GEOparse
from scripts.heatmap import generate_heat_map
import re

def generate_list(path, skip=False, keepid=False, usemircol=False):
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

    if len(complete_df.columns) == 5:
        colunas = ['id','adjpval', 'pval', 'logfc', 'mir']
    else:
        colunas = ['id','adjpval', 'pval', 'logfc']

    complete_df.columns = colunas
    complete_df.sort_values(by='id', axis='index', inplace=True)
    if usemircol == True:
        filtra_hsa = complete_df[complete_df.mir.str.contains('hsa', na=False)]
        filtra_p = filtra_hsa[filtra_hsa.mir.str.contains(mir)].query('pval <= 0.05 & logfc < -1 | logfc > 1')
        lista = list(filtra_p['id'])
        mirs = list(filtra_p['mir'])
        return lista, mirs
    elif usemircol == False:
        filtra_hsa = complete_df[complete_df.id.str.contains('hsa', na=False)]
        filtra_p = filtra_hsa[filtra_hsa.id.str.contains(mir)].query('pval <= 0.05 & logfc < -1 | logfc > 1')
        lista = list(filtra_p['id'])
        return lista

lista59492 = generate_list(path='./datasets/GSE59492/GSE59492_w_logfc.txt', keepid=True)
lista74618 = generate_list(path='./datasets/GSE74618/GSE74618_results.txt', keepid=True)
lista10694, mirs = generate_list(path='./datasets/GSE10694/GSE10694_results.txt', skip=True, 
    keepid=True, usemircol=True)
lista49012 = generate_list(path='./datasets/GSE49012/GSE49012_w_logfc.txt', keepid=True)


datasets = ['GSE59492', 'GSE74618', 'GSE10694','GSE49012']
order59492 = 'XXXXXXXXXXXXXXXXXXXXXX11111000000'
order74618 = 'X1X1X1X1X1X1X11X1X111X1111111111110000XXXXX11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111000000XXXXXXX'
order10694 = "111111111111111111111111111111111111111111111111111111111111111111111111111111XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX0000000000"
order49012 = 'XXXXXXXXXXXXX111111111000000000000'

def run_heatmaps(dataset):
    path = './junk/' + dataset + '_family.soft.gz'
    tabela = GEOparse.get_GEO(filepath=path, silent=True)
    valores = tabela.pivot_samples('VALUE')
    if dataset == 'GSE59492':
        lista = lista59492
        valores.columns = re.findall('.', order59492)
    elif dataset == 'GSE74618':
        lista = lista74618
        valores.columns = re.findall('.', order74618)
    elif dataset == 'GSE10694':
        lista = lista10694  
        valores.columns = re.findall('.', order10694)      
    elif dataset == 'GSE49012':
        lista = lista49012
        valores.columns = re.findall('.', order49012)
    
    valores.drop(columns=['X'], inplace=True)
    valores.reset_index(inplace=True)

    restable = valores[valores['ID_REF'].isin(lista)]
    if dataset == 'GSE10694':
        restable['ID_REF'] = mirs
    else:
        pass
    print(dataset)
    print(lista)
    pathtocsv = './csv/' + dataset + '_resultado.csv'
    restable.to_csv(pathtocsv, index=False)

for ds in datasets:
    run_heatmaps(ds)
