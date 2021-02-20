import pandas as pd
import numpy as np

df = pd.read_csv("csv/new10694.csv", index_col="ID_REF")
df = np.log2(df)
df.to_csv("csv/updated10694.csv")
