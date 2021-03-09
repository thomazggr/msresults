import pandas as pd

xp_val = pd.read_csv("csv/csvn/exp_val.csv", sep=",")

xp_val.corr(method="pearson").to_clipboard()
