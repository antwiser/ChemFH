from ninja import Router, Schema
from typing import List, Union
from service.views import wash_input_mol
import pandas as pd
import time
import tempfile
import os
from .mapping import explanation

SITE_ROOT = os.path.abspath(os.path.dirname(__file__))

router = Router()


class SMILESSchema(Schema):
    SMILES: Union[str, List[str]] = 'CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N'
    feature: bool = False


class SuccessDict:
    def __init__(self, message, data):
        self.status = "success"
        self.code = 200
        self.data = data


class ErrorDict:
    def __init__(self, message):
        self.status = "error"
        self.code = 200
        self.msg = message
        self.data = None


all_rule = ['Aggregators', 'Fluc', 'Blue_fluorescence', 'Green_fluorescence', 'Reactive', 'Other_assay_interference',
            'Promiscuous', 'ALARM_NMR', 'BMS', 'Chelator_Rule', 'GST_FHs_Rule', 'His_FHs_Rule',
            'Luciferase_Inhibitor_Rule', 'NTD', 'PAINS', 'Potential_Electrophilic_Rule', 'Lilly']


def returnHandle(taskid, result, noTotal, noValid):
    return_result = {'taskid': taskid, 'number of all molecules': noTotal, 'number of valid molecules': noValid,
                     'data': result.to_dict(orient="records"), 'explanation': {key: value.replace(
            '<br/>', '') for key, value in explanation.items()}}
    return return_result


@router.post('/fh', description="")
def chemfhCal(request, data: SMILESSchema):
    smiles_list, invalidIdx = wash_input_mol(
        data.SMILES, issmiles=True, invalidStr='invalid', returnInvalidIdx=True)
    if True not in invalidIdx:
        return 200, ErrorDict('invalid molecule!').__dict__
    if len(smiles_list) > 5000:
        return 200, ErrorDict('Excessive number of requested molecules!').__dict__
    # 计算模型预测结果
    from static.media.chemprop.scripts.predict import predict
    result, unResult = predict(smiles_list)
    result = result.round(3)
    unResult = unResult.round(6)
    # 计算规则匹配结果
    import static.media.rule.filter_rule as fr
    out1, out2, out3, out4 = fr.filter_rule(smiles_list, all_rule)
    result = pd.concat([result, unResult, out4], axis=1)
    # 合并结果
    smiles_df = pd.DataFrame(smiles_list, columns=['smiles'])
    result = pd.concat([smiles_df, result], axis=1)
    # 整合无效的结果
    for _, flag in enumerate(invalidIdx):
        if not flag:
            new_row = pd.DataFrame([{col: 'Invalid Molecule' if col !=
                                                                'smiles' else data.SMILES[_] for col in
                                     result.columns}])
            df1 = result.iloc[:_]  # 第一部分
            df2 = result.iloc[_:]  # 第二部分
            result = pd.concat(
                [df1, new_row, df2]).reset_index(drop=True)
    # 保存文件
    tmpf = tempfile.NamedTemporaryFile()
    file_name = tmpf.name.split('/')[-1]
    file_name = file_name.split('\\')[-1]
    t = int(time.time())
    file_name = file_name + str(t)
    file_path = os.path.join(SITE_ROOT, '../static/results') + \
                '/csv/' + file_name + '.csv'
    # result.to_csv(file_path, index=False)
    result = returnHandle(file_name, result, len(
        smiles_list), len(smiles_list) - invalidIdx.count(False))
    return 200, SuccessDict('operation success!', result).__dict__
