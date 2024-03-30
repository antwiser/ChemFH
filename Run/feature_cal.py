import chemprop
import os
import numpy as np
import pandas as pd
from chemprop.features.features_generators import rdkit_2d_normalized_features_generator, morgan_binary_features_generator
from multiprocessing import Pool
from tqdm import tqdm


def feature_cal_process(wash_smiles):
    feature_list = []
    for s in wash_smiles:
        # feature_list.append(rdkit_2d_normalized_features_generator(s[0]))
        feature_list.append(morgan_binary_features_generator(s[0]))
    return feature_list



def multiple_pool(
        input_path="FH_data/model_data",
        output_path="FH_data/feature_FP"
):
    for i in range(1, 11):
        os.makedirs(f'{output_path}/{i}', exist_ok=True)

        csv_list = os.listdir(f"{input_path}/{i}")
        for csv_name in csv_list:
            if os.path.exists(f'{output_path}/{i}/{csv_name.split(".")[0]}_2d.npy'):
                continue
            df = pd.read_csv(f'{input_path}/{i}/{csv_name}')
            smiles_list = df['SMILES']
            smiles = [[i] for i in smiles_list]
            invalid_smiles = chemprop.data.utils.get_invalid_smiles_from_list(
                smiles)
            print('len', len(invalid_smiles))
            indices = [i for value in invalid_smiles for i,
                       x in enumerate(smiles) if x == value]
            wash_smiles = [x for x in smiles if x not in invalid_smiles]
            num_processes = 100  # 指定进程数，可以根据需要进行调整

            chunk_size = len(wash_smiles) // num_processes

            with Pool(num_processes) as pool:
                sub_smiles = [wash_smiles[i:i+chunk_size]
                              for i in range(0, len(wash_smiles), chunk_size)]
                feature_lists = list(
                    tqdm(pool.map(feature_cal_process, sub_smiles)))

            feature_list = [
                feature for sublist in feature_lists for feature in sublist]
            feature_array = np.array(feature_list)

            np.save(
                f'{output_path}/{i}/{csv_name.split(".")[0]}_2d.npy', feature_array)


def multiple_pool2(
        input_path="FH_data/model_data",
        output_path="FH_data/feature_data"
):

    csv_list = os.listdir(input_path)
    for csv_name in csv_list:
        if os.path.exists(f'{output_path}/{csv_name.split(".")[0]}_2d.npy'):
            continue
        df = pd.read_csv(f'{input_path}/{csv_name}')
        smiles_list = df['SMILES']
        smiles = [[i] for i in smiles_list]
        invalid_smiles = chemprop.data.utils.get_invalid_smiles_from_list(
            smiles)
        print('len', len(invalid_smiles))
        indices = [i for value in invalid_smiles for i,
                   x in enumerate(smiles) if x == value]
        wash_smiles = [x for x in smiles if x not in invalid_smiles]
        num_processes = 100  # 指定进程数，可以根据需要进行调整

        chunk_size = len(wash_smiles) // num_processes

        with Pool(num_processes) as pool:
            sub_smiles = [wash_smiles[i:i+chunk_size]
                          for i in range(0, len(wash_smiles), chunk_size)]
            feature_lists = list(
                tqdm(pool.map(feature_cal_process, sub_smiles)))

        feature_list = [
            feature for sublist in feature_lists for feature in sublist]
        feature_array = np.array(feature_list)

        np.save(
            f'{output_path}/{csv_name.split(".")[0]}_2d.npy', feature_array)


if __name__ == "__main__":
    multiple_pool(
        input_path="/home/fuli/my_code/ChemFH/FH_modeldata/model_data",
        output_path="/home/fuli/my_code/ChemFH/FH_modeldata/FP"
    )
