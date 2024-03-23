# import requests
# from time import time
#
# baseUrl = 'http://121.40.210.46:8097'
#
# if __name__ == '__main__':
#     api = '/api/admet'
#     url = baseUrl + api
#     param = {
#         # A example SMILES
#         'SMILES': ["O=C(C)Oc1ccccc1C(=O)O"] * 1,
#         'feature': False,
#         'uncertain': True
#     }
#     start = time()
#     # Send a request
#     response = requests.post(url, json=param)
#     if response.status_code == 200:
#         data = response.json()['data']
#     else:
#         print("Request failed!")
#     print(time() - start)
import pandas as pd
from rdkit import Chem

data = pd.read_csv('test.csv')
for item in data.iloc[0]:
    if not isinstance(item, float):
        for i in item.split('|'):
            print(Chem.MolToSmiles(Chem.MolFromSmarts(i.strip())))
