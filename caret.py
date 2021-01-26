from pycaret.datasets import get_data
from pycaret.classification import *
from sklearn.model_selection import train_test_split  # Scikit separar treino/teste
import pandas as pd
import re

data = pd.read_csv('./csv/new10694.csv')
data.set_index('ID_REF', drop=True, inplace=True)

data = data.T
i = "1111111111111111111111111111111111111111111111111111111111111111111111111111110000000000"
data.index = re.findall('.', i)
data.reset_index(inplace=True)
data.drop(columns=['hsa-miR-221', 'hsa-miR-222', 'hsa-miR-155', 'hsa-miR-150'], inplace=True)


'''train, test = train_test_split(data, test_size=0.2,
                                     random_state=42)'''

exp_name = setup(data = data,  target = 'index')
lr = create_model('lr')
plot_model(lr)