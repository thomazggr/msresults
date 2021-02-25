import pandas as pd
import re
import operator
from functools import reduce

o74 = "111111111111111111111111000011111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111000000"
o10 = "1111111111111111111111111111111111111111111111111111111111111111111111111111110000000000"
o49 = "111111111000000000000"
o19 = "11111000000"

df74 = pd.read_csv("./csv/csvn/new74618.csv")
df74.set_index("ID_REF", drop=True, inplace=True)
df74.columns = re.findall(".", o74)

df10 = pd.read_csv("./csv/csvn/updated10694.csv")
df10.set_index("ID_REF", drop=True, inplace=True)
df10.columns = re.findall(".", o10)

df49 = pd.read_csv("./csv/csvn/updated49012.csv")
df49.set_index("ID_REF", drop=True, inplace=True)
df49.columns = re.findall(".", o49)

df59 = pd.read_csv("./csv/csvn/new59492.csv")
df59.set_index("ID_REF", drop=True, inplace=True)
df59.columns = re.findall(".", o19)

dictdf = {"df74": df74, "df10": df10, "df49": df49, "df59": df59}

empdf = pd.DataFrame()

for df in dictdf.values():
    a = df.loc[:, df.columns.isin(["0"])]
    b = a.filter(items=["hsa-miR-422a"], axis="index")
    if b.empty:
        pass
    else:
        c = b.T
        c.to_clipboard(index=False, header=False)
    inp = input("Has been copied. Press N for next\n")
    if inp == "N":
        pass

    empdf.to_clipboard()

# print(mir150_1)
