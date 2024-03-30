
import os
import re
import chemprop
import pandas as pd

pattern1 = re.compile(r'\[n\+H\]', re.IGNORECASE)
pattern2 = re.compile(r'\[n\+H2\]', re.IGNORECASE)
pattern3 = re.compile(r'\[n\+H3\]', re.IGNORECASE)
pattern4 = re.compile(r'\[n\+@H\]', re.IGNORECASE)


def main(input_path, output_path):
    file_list = os.listdir(
        f'{input_path}')
    for file in file_list:
        if file == 'reactive_lzz.csv' or file == 'promiscuous_lzz.csv':
            df = pd.read_csv(f'{input_path}/{file}')
            smiles_list = df['SMILES']
            labels = df['Label']
            new_dict = {
                'SMILES': [],
                'Label': []
            }
            invalid_smiles_l = []
            for i, label in zip(smiles_list, labels):
                i = str(i)
                if i is not None and pattern1.search(i):
                    i = pattern1.sub('[NH+]', i)
                if i is not None and pattern2.search(i):
                    i = pattern2.sub('[NH2+]', i)
                if i is not None and pattern3.search(i):
                    i = pattern3.sub('[NH3+]', i)
                if i is not None and pattern4.search(i):
                    i = pattern4.sub('[N@H+]', i)
                try:
                    invalid_smiles = chemprop.data.utils.get_invalid_smiles_from_list(
                        [[i]])
                except:
                    invalid_smiles = []
                invalid_smiles_l.append(invalid_smiles)
                if invalid_smiles == []:
                    new_dict['SMILES'].append(i)
                    new_dict['Label'].append(label)

            new_df = pd.DataFrame(new_dict)
            print(f"{file}, len: {len(invalid_smiles_l)}")
            new_df.to_csv(
                f'{output_path}/{file}', index=False)
            print(f"Save {file} successfully!")


if __name__ == "__main__":
    main(input_path="FH_data/data", output_path="FH_data/filter")
