import pandas as pd

mirdip_data = pd.read_csv("mirDIP_E_2022_08_07_20 45 44.txt", sep="\t", skiprows=13)

print(mirdip_data.query("`Gene Symbol` == 'UBD'").head())

# print(mirdip_data["MicroRNA"].value_counts())

# print(mirdip_data["Gene Symbol"].unique())