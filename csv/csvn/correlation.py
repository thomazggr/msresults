import pandas as pd

xp_val = pd.read_csv("results_nashonly/exp-val.csv", sep=",")
xp_val.set_index("ID_REF", inplace=True, drop=True)
xp_val = xp_val.T
# print(xp_val)
xp_val.corr(method="pearson").to_clipboard()
# print(xp_val.corr(method="pearson"))
