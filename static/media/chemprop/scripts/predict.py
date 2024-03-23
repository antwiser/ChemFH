import time
import re
import os
import sys
import json
import pandas as pd
import numpy as np

from tqdm import tqdm
from multiprocessing import Pool

import chemprop
from chemprop.features.features_generators import rdkit_2d_normalized_features_generator
from .Process_runner import ProcessRunner

#
# chemfh_colname = [
#     "Aggregators",
#     "FLuc inhibitors",
#     "Blue fluorescence",
#     "Green fluorescence",
#     "Reactive compounds",
#     "Other assay interference",
#     "Promiscuous compounds",
# ]
#
# chemfh_uncolname = [
#     "Aggregators uncertainty",
#     "FLuc inhibitors uncertainty",
#     "Blue fluorescence uncertainty",
#     "Green fluorescence uncertainty",
#     "Reactive compounds uncertainty",
#     "Other assay interference uncertainty",
#     "Promiscuous compounds uncertainty",
# ]

chemfh_colname = [
    "Other assay interference",
    "Blue fluorescence",
    "FLuc inhibitors",
    "Promiscuous compounds",
    "Colloidal aggregators",
    "Reactive compounds",
    "Green fluorescence"
]

chemfh_uncolname = [
    "Other assay interference uncertainty",
    "Blue fluorescence uncertainty",
    "FLuc inhibitors uncertainty",
    "Promiscuous compounds uncertainty",
    "Colloidal aggregators uncertainty",
    "Reactive compounds uncertainty",
    "Green fluorescence uncertainty"
]


class suppress_stdout_stderr(object):
    '''
    A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
       This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).
    '''

    def __init__(self):
        # Open a pair of null files
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = (os.dup(1), os.dup(2))

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close the null files
        os.close(self.null_fds[0])
        os.close(self.null_fds[1])


class ChempropError(Exception):
    pass


# 计算特征
def feature_cal_process(wash_smiles):
    feature_list = []
    for s in wash_smiles:
        feature_list.append(rdkit_2d_normalized_features_generator(s[0]))
    return feature_list


# 并行计算特征
def multiple_pool_feature(smiles_list, num_processes=8):
    smiles = [[i] for i in smiles_list]
    invalid_smiles = chemprop.data.utils.get_invalid_smiles_from_list(smiles)
    indices = [i for value in invalid_smiles for i,
    x in enumerate(smiles) if x == value]
    wash_smiles = [x for x in smiles if x not in invalid_smiles]
    if len(wash_smiles) <= 50:
        feature_list = feature_cal_process(wash_smiles)
    else:
        with Pool(num_processes) as pool:
            chunk_size = len(wash_smiles) // num_processes
            sub_smiles = [wash_smiles[i:i + chunk_size]
                          for i in range(0, len(wash_smiles), chunk_size)]
            feature_lists = list(
                tqdm(pool.map(feature_cal_process, sub_smiles)))
            feature_list = [
                feature for sublist in feature_lists for feature in sublist]
    feature_array = np.array(feature_list)
    indices = list(set(indices))
    return feature_array, indices


def chemfh_pred(
        smiles_file,
        pred_file,
        feature_cmd,
        MODEL_PATH,
        smi_df
):
    arguments = [
        '--test_path', smiles_file,
        '--preds_path', pred_file,
        "--checkpoint_dir", f"{MODEL_PATH}/fold_0/model_0/",
        "--num_workers", "0",
        "--no_cuda",
        "--uncertainty_method", "dropout"
    ]
    arguments += feature_cmd
    # print(arguments)

    args = chemprop.args.PredictArgs().parse_args(arguments)
    preds, uncertainty = chemprop.train.make_predictions(
        args=args,
        return_uncertainty=True,
    )

    sub_df = pd.DataFrame(
        preds,
        columns=chemfh_colname,
        index=list(smi_df['smiles'])
    )

    sub_df_un = pd.DataFrame(
        uncertainty,
        columns=chemfh_uncolname,
        index=smi_df['smiles']

    )
    cols = ['smiles'] + [col for cols in zip(
        chemfh_colname, chemfh_uncolname) for col in cols]
    # 合并 DataFrame
    sub_df = pd.concat([sub_df, sub_df_un], axis=1)
    sub_df = sub_df.reindex(cols, axis=1)
    sub_df['smiles'] = sub_df.index
    return sub_df


def main(smiles_list, model_path, num_processes=20):
    '''
    smiles_list: []
    num_processes: 计算特征时并行进程数， default=1
    feature: 是否计算特征， default=True
    uncertain: 是否计算不确定性， default=False
    group: 需要计算的任务组， default='all'计算所有的性质共 77 个，
            可选 'pcp', 'absorption', 'distribution',
                'metabolism', 'excretion', 'toxicity', 'tox21'
    '''
    with ProcessRunner() as P:
        feature_file = 'feature.npy'
        smiles_file = 'input.csv'
        pred_file = 'output.csv'
        df = pd.DataFrame()

        indices = []

        feature_cmd = [
            "--features_path", f'{os.getcwd()}/{feature_file}',
            "--no_features_scaling"
        ]
        feature_array, indices = multiple_pool_feature(
            smiles_list, num_processes=num_processes)
        with open(feature_file, 'wb') as file:
            np.save(file, feature_array)

        smi_df = pd.DataFrame({'smiles': smiles_list})
        original_index = smi_df.index.tolist()
        removed_rows = pd.DataFrame(
            data={'smiles': smi_df['smiles'][indices]})
        smi_df = smi_df.drop(indices)
        after_index = smi_df.index.tolist()
        smi_df.to_csv(smiles_file, index=False)

        sub_df = chemfh_pred(
            smiles_file,
            pred_file,
            feature_cmd,
            model_path,
            smi_df
        )
        if indices == []:
            sub_df.set_index('smiles', inplace=True)
            df = pd.concat([df, sub_df], axis=1)
        else:
            # 插入无效的 smiles 行
            sub_df = sub_df.merge(
                removed_rows, on='smiles', how='outer').fillna('Invalid SMILES')
            index = after_index + indices
            sub_df['index'] = index
            sub_df.set_index('index', inplace=True)
            sub_df = sub_df.reindex(original_index)
            sub_df.set_index('smiles', inplace=True)
            df = pd.concat([df, sub_df], axis=1)
        return df


# def main(smiles_list, MODEL_PATH, num_processes=1):
#     '''
#     smiles_list: []
#     num_processes: 计算特征时并行进程数， default=1
#     feature: 是否计算特征， default=True
#     uncertain: 是否计算不确定性， default=False
#     group: 需要计算的任务组， default='all'计算所有的性质共 77 个，
#             可选 'pcp', 'absorption', 'distribution',
#                 'metabolism', 'excretion', 'toxicity', 'tox21'
#     '''
#     with ProcessRunner() as P:
#         feature_file = 'feature.npy'
#         smiles_file = 'input.csv'
#         pred_file = 'output.csv'
#         df = pd.DataFrame()
#
#         indices = []
#         feature_cmd = [
#             "--features_path", f'{os.getcwd()}/{feature_file}',
#             "--no_features_scaling"
#         ]
#         feature_array, indices = multiple_pool_feature(
#             smiles_list, num_processes=num_processes)
#         with open(feature_file, 'wb') as file:
#             np.save(file, feature_array)
#
#         smi_df = pd.DataFrame({'smiles': smiles_list})
#         original_index = smi_df.index.tolist()
#         removed_rows = pd.DataFrame(
#             data={'smiles': smi_df['smiles'][indices]})
#         smi_df = smi_df.drop(indices)
#         after_index = smi_df.index.tolist()
#         smi_df.to_csv(smiles_file, index=False)
#
#         sub_df = chemfh_pred(
#             smiles_file,
#             pred_file,
#             feature_cmd,
#             MODEL_PATH,
#             smi_df
#         )
#         print(sub_df.to_dict())
#         if not indices:
#             sub_df.set_index('smiles', inplace=True)
#             df = pd.concat([df, sub_df], axis=1)
#         else:
#             # 插入无效的 smiles 行
#             sub_df = sub_df.merge(
#                 removed_rows, on='smiles', how='outer').fillna('Invalid SMILES')
#             index = after_index + indices
#             sub_df['index'] = index
#             sub_df.set_index('index', inplace=True)
#             sub_df = sub_df.reindex(original_index)
#             sub_df.set_index('smiles', inplace=True)
#             df = pd.concat([df, sub_df], axis=1)
#         print('df: ', df)
#         df.to_csv('df.csv')
#         return df


from django.conf import settings


def predict(smiles_list):
    startTime = time.time()
    path = os.path.join(settings.SITE_ROOT,
                        'static') + '/media/chemprop/checkpoint'
    with suppress_stdout_stderr():  # suppress_stdout_stderr() 用于屏蔽 chemprop 的输出
        result = main(smiles_list,
                      model_path=path,
                      num_processes=20)
    result = result.reset_index(drop=True)
    result, unResult = result[chemfh_colname], result[chemfh_uncolname]
    print(f"time:{time.time() - startTime}")
    return result, unResult


if __name__ == "__main__":
    ss = ['C1=CC=C(C=C1)C(=O)NC2=C(N=C(S2)N)C3=CC=CS3'] * 1
    smiles_list = [
        i for i in ss
    ]
    print(len(smiles_list))

    # 1 个分子总耗时 0.9s, 500个分子耗时15.8s，1000个分子耗时27s
    start = time.time()
    with suppress_stdout_stderr():  # suppress_stdout_stderr() 用于屏蔽 chemprop 的输出
        df = main(smiles_list,
                  model_path='/static/media/chemprop/checkpoint',
                  num_processes=20)
    df.to_csv('result.csv')
    end = time.time()
    print(f"time:{end - start}")
