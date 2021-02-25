import pandas as pd
from scripts.heatmap import generate_heat_map
import numpy as np

n10694 = pd.read_csv("./csv/new10694.csv")
n49012 = pd.read_csv("./csv/new49012.csv")
n59492 = pd.read_csv("./csv/new59492.csv")
n74618 = pd.read_csv("./csv/new74618.csv")

"""nv2 = n49012.merge(right=n10694, how='left',on='ID_REF', suffixes=('_a', None))
nv3 = nv2.merge(right=n59492, how='left', on='ID_REF')
nv4f = nv3.merge(right=n74618, how='left', on='ID_REF')
nv4f.set_index('ID_REF', drop=True, inplace=True)
nv4flog = np.log2(nv4f)
nv4flog.sort_index(axis='columns', inplace=True)"""

n49012.set_index("ID_REF", drop=True, inplace=True)
n49012x = np.log2(n49012)
n49012x.sort_index(axis="columns", inplace=True)
generate_heat_map(n49012x, nome="hm_n49012.jpg")
