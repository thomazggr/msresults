import pandas as pd
import GEOparse
import scipy
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import re


def generate_hm(df):
    """
    Function to generate headmaps with dendrogram
    """
    # figure size
    fig = plt.figure(figsize=(10, 8))
    # need to change the line width
    # default line width doesn't look good
    matplotlib.rcParams["lines.linewidth"] = 0.5
    # axis to show where this current axis
    # need to be shown
    ax1 = fig.add_axes([0.09, 0.1, 0.2, 0.6])
    # to generate the array the will be passed
    # to dendrogram, mey
    Y = sch.linkage(df, method="centroid")
    # need to make changes in orientation to
    # make it appear in the left
    # default is top
    Z1 = sch.dendrogram(Y, orientation="left")
    # removing x and y ticks
    ax1.set_xticks([])
    ax1.set_yticks([])
    # to remove the axis
    ax1.axis("off")

    # top side dendogram
    ax2 = fig.add_axes([0.3, 0.71, 0.6, 0.2])
    Y = sch.linkage(df.T, method="ward")
    Z2 = sch.dendrogram(Y)
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.axis("off")

    # main heatmap
    axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.6])
    # getting the index to get our required values
    # that needs to be shown in the heatmap
    idx1 = Z1["leaves"]
    idx2 = Z2["leaves"]
    D = df.values[idx1, :]
    D = df.values[:, idx2]
    #     D = df.values
    # using matshow to display the heatmap
    # colormap: 'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'
    im = axmatrix.matshow(D, aspect="auto", origin="lower", cmap="seismic")
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # xticks to the right (x-axis)
    axmatrix.set_xticks(range(len(df.columns)))
    # this shows the required labels
    # in this case it is the dataframe columns
    axmatrix.set_xticklabels(df.columns, minor=False)
    axmatrix.xaxis.set_label_position("bottom")
    axmatrix.xaxis.tick_bottom()

    # to change the rotation of xticks
    plt.xticks(rotation=90, fontsize=8)

    # xticks to the right (y-axis)
    axmatrix.set_yticks(range(len(df.index)))
    # this shows the required labels
    # in this case it is the dataframe index
    axmatrix.set_yticklabels(df.index, minor=False)
    axmatrix.yaxis.set_label_position("right")
    axmatrix.yaxis.tick_right()

    plt.savefig("teste.png")


def gse():
    tabela = GEOparse.get_GEO(geo="GSE49012", destdir="./junk/", silent=True)
    valores = tabela.pivot_samples("VALUE")
    return valores


tabela = gse()


def hm(table):
    # id_and_mir = pd.read_csv('./tbu/GSE10694/GSE10694_results.txt', sep='\t', skiprows=[0])
    # id_and_mir.set_index('ID', inplace=True, drop=True)
    # id_and_mir.sort_index(inplace=True)
    # table['newindex'] = id_and_mir['miRNA_ID']
    # table.set_index('newindex', inplace=True, drop=True)
    filterhsa = table.filter(regex="^hsa", axis=0)
    filtermirs = filterhsa.filter(
        regex="hsa-miR-200b|hsa-miR-483|hsa-miR-155|hsa-miR-940|hsa-miR-222|hsa-miR-221|\
            hsa-miR-181d|hsa-miR-200a|hsa-miR-193b|hsa-miR-31|hsa-miR-143|hsa-miR-199b|\
                hsa-miR-378|hsa-miR-150|hsa-miR-182|hsa-miR-10a|hsa-miR-21|hsa-miR-146b|\
                    hsa-miR-200c|hsa-miR-1234.*",
        axis=0,
    )
    samples = "XXXXXXXXXXXXX111111111000000000000"
    nindex = splitamostras(samples=samples)
    filtern = filtermirs.T
    filtern["NIND"] = nindex
    filtern.set_index(filtern["NIND"], inplace=True, drop=True)
    filtern.drop(columns=["NIND"], index=["X"], inplace=True)
    filtern.sort_index(inplace=True)
    tabt = filtern.T
    tablog = np.log2(tabt)
    print(tablog)

    generate_hm(tablog)


def splitamostras(samples):
    resultado = re.findall(".", samples)
    return resultado


hm(table=tabela)