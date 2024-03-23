import pandas as pd
import os

os.chdir("E:\\Papers\\chemfh\\ChemFH")

import static.media.rule.filter_rule as fr

all_rule = ['Aggregators', 'Fluc', 'Blue_fluorescence', 'Green_fluorescence', 'Reactive', 'Other_assay_interference',
            'Promiscuous', 'ALARM_NMR', 'BMS', 'Chelator_Rule', 'GST_FHs_Rule', 'His_FHs_Rule',
            'Luciferase_Inhibitor_Rule', 'NTD', 'PAINS', 'Potential_Electrophilic_Rule']

data = pd.read_csv('service/ChemDiv_pred.csv')
smiles_list = data['smiles'].tolist()
out1, out2, out3, out4 = fr.filter_rule(smiles_list, all_rule)
print(out1)
