from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render, get_object_or_404, redirect
from django.urls import reverse
from django.views import generic
from django.utils import timezone
import tempfile
import time
import os
import math
from django.conf import settings
import pandas as pd
import numpy as np
from django.http import FileResponse
from rdkit.Chem import PandasTools as pt
import json
from random import randrange
from django.http import HttpResponse

from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor

import static.media.rule.filter_rule as fr
from rdkit.Chem import Draw
import base64
from io import BytesIO

from service.utils.handle_pdf import gen_pdf


def HighlightAtoms(mol, highlightAtoms, figsize=[400, 200], kekulize=True):
    """This function is used for showing which part of fragment matched the SMARTS by the id of atoms.

    :param mol: The molecule to be visualized
    :type mol: rdkit.Chem.rdchem.Mol
    :param highlightAtoms: The atoms to be highlighted
    :type highlightAtoms: tuple
    :param figsize: The resolution ratio of figure
    :type figsize: list
    :return: a figure with highlighted molecule
    :rtype: IPython.core.display.SVG

    """

    def _revised(svg_words):
        """
        """
        svg_words = svg_words.replace(
            'stroke-width:2px', 'stroke-width:1.5px').replace(
            'fonts-size:17px', 'fonts-size:15px').replace(
            'stroke-linecap:butt', 'stroke-linecap:square').replace(
            'fill:#FFFFFF', 'fill:none').replace(
            'svg:', '')
        return svg_words

    mc = Chem.Mol(mol.ToBinary())

    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(*figsize)
    drawer.DrawMolecule(mc, highlightAtoms=highlightAtoms)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # It seems that the svg renderer used doesn't quite hit the spec.
    # Here are some fixes to make it work in the notebook, although I think
    # the underlying issue needs to be resolved at the generation step
    return _revised(svg)


un_threhold = [
    0.0010918042703688,
    0.004291940348845,
    0.0009872382265281,
    0.0023240092127062,
    0.0016870113071982,
    0.0061948686233965,
    9.781904685012926e-05,
]

all_rule = ['Aggregators', 'Fluc', 'Blue_fluorescence', 'Green_fluorescence', 'Reactive', 'Other_assay_interference',
            'Promiscuous', 'ALARM_NMR', 'BMS', 'Chelator_Rule', 'GST_FHs_Rule', 'His_FHs_Rule',
            'Luciferase_Inhibitor_Rule', 'NTD', 'PAINS', 'Potential_Electrophilic_Rule', 'Lilly']
from rdkit import Chem
from rdkit import RDConfig
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Descriptors, Lipinski, QED
from rdkit.Chem.Scaffolds import MurckoScaffold
from math import pi, log10
from collections.abc import Iterable
from itertools import combinations


def evaluationIndex(request):
    return render(request, "services/evaluation/index.html", {})


def filterIndex(request):
    return render(request, "services/screening/index.html")


from service.utils.utils import getMD5
from .models import chemfh

name2dbName = {
    'id': 'id', 'taskid': 'taskId', 'smiles': 'smiles', 'aggregators': 'Colloidal aggregators',
    'fluc': 'FLuc inhibitors',
    'blue': 'Blue fluorescence',
    'green': 'Green fluorescence', 'reactive': 'Reactive compounds', 'other': 'Other assay interference',
    'promiscuous': 'Promiscuous compounds', 'aggregators_un': 'Colloidal aggregators uncertainty',
    'fluc_un': 'FLuc inhibitors uncertainty', 'blue_un': 'Blue fluorescence uncertainty',
    'green_un': 'Green fluorescence uncertainty', 'reactive_un': 'Reactive compounds uncertainty',
    'other_un': 'Other assay interference uncertainty', 'promiscuous_un': 'Promiscuous compounds uncertainty',
    'aggregators_idx': 'Aggregators_index',
    'fluc_idx': 'Fluc_index', 'blue_idx': 'Blue_fluorescence_index', 'green_idx': 'Green_fluorescence_index',
    'reactive_idx': 'Reactive_index',
    'other_idx': 'Other_assay_interference_index', 'promiscuous_idx': 'Promiscuous_index',
    'alarm': 'ALARM_NMR_index', 'bms': 'BMS_index', 'chelator': 'Chelator_Rule_index',
    'gst': 'GST_FHs_Rule_index', 'his': 'His_FHs_Rule_index',
    'luciferase': 'Luciferase_Inhibitor_Rule_index', 'ntd': 'NTD_index', 'pains': 'PAINS_index',
    'potential': 'Potential_Electrophilic_Rule_index', 'lilly': 'Lilly_index'
}

dbName2name = {v: k for k, v in name2dbName.items()}

other_rule_explanation = [
    'Thiol reactive compounds;<br/> 75 substructures <br/>J Am Chem Soc 2005;127:217–24',
    'Undesirable, reactive compounds;<br/> 176 substructures <br/>J Chem Inf Model 2006;46:1060–8',
    'Chelators;<br/> 55 substructures <br/>ChemMedChem 2010;5:195–9',
    'GST/GSH FHs;<br/> 34 substructures <br/>J Biomol Screen 2016;21:596–607',
    'Ni2+ chelators;<br/> 19 substructures <br/>J Biomol Screen 2014;19:715–26',
    'FLuc inhibitors;<br/> 3 substructures <br/>J Chem Inf Model 2018;58:933–42',
    'Unwanted groups, reactive groups and possible HTS interferences;<br/> 105 substructures <br/>ChemMedChem '
    '2008;3:435–44',
    'Frequent hitters, Alpha-screen artifacts and reactive compound;<br/> 480 substructures <br/>J Med Chem '
    '2010;53:2719–40',
    'Reactive compounds;<br/> 119 substructures <br/>J Chem Inf Model 2012;52:2310–6',
    'Interfere with biological assays; <br/> 273 substructures <br/>J Med Chem 2012;55:9763-9772'
]

mechanism_explanation = [
    'Category 0: Non-aggregators; Category 1: aggregators. <br/>The output value is the probability of being '
    'aggregators, within the range of 0 to 1.',
    'Firefly luciferase (FLuc) inhibitors. Category 0: Non-FLuc inhibitors; Category 1: FLuc inhibitors. '
    '<br/>The output value is the probability of being FLuc inhibitors, within the range of 0 to 1.',
    'Category 0: Non-fluorescence; Category 1: fluorescence. <br/>The output value is the probability of being '
    'fluorescence, within the range of 0 to 1.',
    'Category 0: Non-fluorescence; Category 1: fluorescence. <br/>The output value is the probability of being '
    'fluorescence, within the range of 0 to 1.',
    'Category 0: Non-reactive; Category 1: reactive.<br/> The output value is the probability of being '
    'reactive compounds, within the range of 0 to 1.',
    'Category 0: Non-promiscuous; Category 1: promiscuous. <br/>The output value is the probability of being '
    'promiscuous, within the range of 0 to 1.',
    'Category 0: Non-assay interferences; Category 1: assay interferences. <br/>The output value is the probability '
    'of being assay interferences, within the range of 0 to 1.'
]

other_rule_display = ['ALARM NMR Rule', 'BMS Rule', 'Chelator Rule', 'GST FHs Rule',
                      'His FHs Rule', 'Luciferase Inhibitor Rule', 'NTD Rule', 'PAINS Rule',
                      'Potential Electrophilic Rule', 'Lilly Medchem Rules']


def evaluationCal(request):
    score_level = 'High Risk'
    if request.method == 'POST':
        message = ''
        method = request.POST.get('method')
        if method == '1':
            smiles = request.POST.get('smiles')
        else:
            smiles = str(request.POST.get('mol'))
        if not smiles:
            message = 'The SMILES is invalid! please check!'
            return render(request, 'services/evaluation/index.html', locals())
        structure = ''
        mol = Chem.MolFromSmiles(smiles)
        if not mol:  # mol为空，说明输入的SMILES字符串是错误的
            message = 'The SMILES is invalid! please check!'
            tab_choice = 'input'
            return render(request, 'services/evaluation/index.html', locals())

        token = getMD5(smiles)
        task = chemfh.objects.filter(taskid=token).exists()
        t = token
        if task:
            # 表示该smiles存在数据库
            result = chemfh.objects.get(taskid=token)
            data = {field: value for field, value in result.__dict__.items() if not field.startswith('_')}
            df = pd.DataFrame([data])
            result = df.rename(columns=name2dbName)
        else:
            from static.media.chemprop.scripts.predict import predict
            result, unResult = predict([smiles])
            result = result.round(3)
            out1, out2, out3, out4 = fr.filter_rule([smiles], all_rule)
            smiles_df = pd.DataFrame([smiles], columns=['smiles'])
            result = pd.concat([smiles_df, result, unResult, out4], axis=1)
            # 保存进数据库
            task_id = pd.DataFrame([t], columns=['taskId'])
            result = pd.concat([task_id, result], axis=1)
            save_result = result.rename(columns=dbName2name)
            chemfh.objects.create(
                **(save_result.to_dict(orient="records")[0]))
        structure = HighlightAtoms(mol, highlightAtoms=(), figsize=[400, 400])
        mechanisms_name = ['Colloidal aggregators', 'FLuc inhibitors',
                           'Blue fluorescence', 'Green fluorescence', 'Reactive compounds', 'Promiscuous compounds',
                           'Other assay interference', ]
        uncertainty_name = [item + ' uncertainty' for item in mechanisms_name]
        score_list = result[mechanisms_name].iloc[0]
        score_result = score_list.transpose().reset_index()
        score_result.columns = ['mechanism', 'score']
        uncertainty_list = result[uncertainty_name].iloc[0].values.tolist()
        score_result['uncertainty'] = [item < un_threhold[_] for _, item in enumerate(uncertainty_list)]
        return_result = score_result

        def func_decision(row):
            if row['score'] > 0.5:
                return False
            else:
                return True

        return_result['decision'] = return_result.apply(func_decision, axis=1)
        return_result['attention'] = ['-'] * len(return_result)
        rule_name = [item + '_index' for item in
                     ['Aggregators', 'Fluc', 'Blue_fluorescence', 'Green_fluorescence', 'Reactive', 'Promiscuous',
                      'Other_assay_interference', ]]
        return_result['atom_list'] = result[rule_name].iloc[0].values.tolist()
        return_result['cnt'] = return_result['atom_list'].map(lambda x: len(eval(x) if isinstance(x, str) else x))
        return_result['svg'] = return_result['atom_list'].map(
            lambda x: ''.join([HighlightAtoms(mol, highlightAtoms=item, figsize=[200, 200]) for item in (eval(x) if
                                                                                                         isinstance(x,
                                                                                                                    str) else x)])
            if (eval(x) if isinstance(x, str) else x) != [] else '')
        return_result['explanation'] = mechanism_explanation
        other_rule_name = [item + '_index' for item in
                           ['ALARM_NMR', 'BMS', 'Chelator_Rule', 'GST_FHs_Rule', 'His_FHs_Rule',
                            'Luciferase_Inhibitor_Rule', 'NTD', 'PAINS',
                            'Potential_Electrophilic_Rule', 'Lilly']]
        other_rule_result = result[other_rule_name].iloc[0].transpose().reset_index()
        other_rule_result.columns = ['mechanism', 'atom_list']
        other_rule_result['cnt_match'] = other_rule_result['atom_list'].map(lambda x: len(eval(x) if isinstance(x,
                                                                                                                str)
                                                                                          else x))
        other_rule_result['img_match'] = other_rule_result['atom_list'].map(
            lambda x: ''.join([HighlightAtoms(mol, highlightAtoms=item, figsize=[200, 200]) for item in (eval(x) if
                                                                                                         isinstance(x,
                                                                                                                    str) else x)])
            if (eval(x) if isinstance(x, str) else x) != [] else '')
        other_rule_result['mechanism'] = other_rule_display
        other_rule_result['explanation'] = other_rule_explanation

        # 渲染雷达图
        # radar_view = getPlot()

        return render(request, 'services/evaluation/result_index.html', {
            'smiles': smiles,
            'structure': structure,
            'mechanisms': mechanisms_name,
            'score': return_result['score'].tolist,
            'result': return_result,
            'other_rules': other_rule_result,
            'download': True,
            'taskid': t,
        })
    else:
        return render(request, 'services/evaluation/index.html', locals())


def handle_uploaded_file(f, path):
    with open(path, 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)


def save_to_file(file_name, contents):
    fh = open(file_name, 'wb')
    fh.write(contents)
    fh.close()


SITE_ROOT = os.path.abspath(os.path.dirname(__file__))
import csv
from service.utils.handle_mol import wash_input_mol


def screeningCal(request):
    if request.method == 'POST':
        message = ''
        method = request.POST.get('method')
        file_obj = request.FILES.get('uploadfile', None)
        is_example = request.POST.get('is_example', '1')
        invalid_datas = []  # 表示无效的分子列表
        molecule_cnt = 0
        tmpf = tempfile.NamedTemporaryFile()
        file_name = tmpf.name.split('/')[-1]
        file_name = file_name.split('\\')[-1]
        t = int(time.time())
        file_name = file_name + str(t)
        smiles_list = []
        if method == '1':
            smiles_list = []
            if is_example == '1':  # 使用示例文件
                csv_file = pd.read_csv(
                    os.path.join(SITE_ROOT, '../static') + '/screening/files/example.csv')
                smiles_list = csv_file['SMILES'].tolist()
            else:  # 用户自己上传的文件
                if not file_obj:  # 上传文件无效
                    message = 'The file is invalid! please check!'
                    return render(request, 'services/screening/index.html', locals())
                file_size = file_obj.size
                if file_size == 0:  # 上传文件为空
                    message = 'The file is empty! please check!'
                    return render(request, 'services/screening/index.html', locals())
                if file_size >= 3000000:  # 上传文件太大
                    message = 'Your file size: ' + \
                              str(file_size / 1024) + \
                              'KB. The max file size is 3000KB!'
                    return render(request, 'services/screening/index.html', locals())
                filename, filetype = os.path.splitext(file_obj.name)
                if str.lower(filetype) == '.sdf':
                    file_path = os.path.join(
                        SITE_ROOT, '../static/results') + '/sdf/' + file_name + '.sdf'
                    handle_uploaded_file(file_obj, file_path)
                    suppl = Chem.SDMolSupplier(file_path)
                    smiles_list = [Chem.MolToSmiles(
                        mol) if mol else 'Invalid Molecule' for _, mol in enumerate(suppl)]
                elif str.lower(filetype) == '.txt':
                    for line in file_obj.readlines():
                        line = str(line, encoding='utf-8').strip()
                        if str.lower(line) == 'smiles':
                            continue
                        smiles_list.append(line)
                elif str.lower(filetype) == '.csv':
                    file_data = file_obj.read().decode('utf-8').splitlines()
                    reader = csv.reader(file_data)
                    for row in reader:
                        print(row)
                        if 'smiles' in row[0].lower():
                            continue
                        smiles_list.append(row[0])
                else:
                    message = 'Invalid file type!'
                    return render(request, 'services/screening/index.html', locals())
        elif method == '2':
            smiles_list = request.POST.get('smiles-list')
            smiles_list = smiles_list.strip().split('\r\n')
        print(smiles_list)

        if len(smiles_list) > 5000:
            message = 'The number of molecules exceeds a predetermined number!'
            return render(request, 'services/screening/index.html', locals())

        washed_smiles_list, invalidIdx = wash_input_mol(
            smiles_list, issmiles=True, invalidStr='invalid', returnInvalidIdx=True)

        if True not in invalidIdx:
            message = 'Invalid SMILES!'
            return render(request, 'services/screening/index.html', locals())

        from static.media.chemprop.scripts.predict import predict
        result, unResult = predict(smiles_list)

        # print(result)
        result = result.round(3)
        out1, out2, out3, out4 = fr.filter_rule(smiles_list, all_rule)
        out1 = out1.drop(columns=['Molecule', 'SMILES'])
        out2 = out2.drop(columns=['Molecule', 'SMILES'])
        out3 = out3.map(lambda x: x.replace(' ', '').split('|') if x != '' else '')
        result = pd.concat([result, unResult, out4], axis=1)
        smiles_df = pd.DataFrame(smiles_list, columns=['smiles'])
        result = pd.concat([smiles_df, result], axis=1)
        for _, flag in enumerate(invalidIdx):
            if not flag:
                new_row = pd.DataFrame([{col: 'Invalid Molecule' if col !=
                                                                    'smiles' else smiles_list[_] for col in
                                         result.columns}])
                df1 = result.iloc[:_]  # 第一部分
                df2 = result.iloc[_:]  # 第二部分
                result = pd.concat(
                    [df1, new_row, df2]).reset_index(drop=True)
        # 保存文件
        file_path = os.path.join(SITE_ROOT, '../static/results') + \
                    '/csv/' + file_name + '.csv'
        result.to_csv(file_path, index=False)
        return HttpResponseRedirect(
            reverse("service:result", kwargs={
                "filename": str(file_name)})
        )


def result_file(request, filename):
    result = pd.read_csv(os.path.join(SITE_ROOT, '../static/results') +
                         '/csv/' + filename + '.csv')
    data = result.to_dict(orient='records')
    smiles_list = result['smiles'].tolist()
    washed_smiles_list, invalidIdx = wash_input_mol(
        smiles_list, issmiles=True, invalidStr='invalid', returnInvalidIdx=True)
    return_result = {'data': data, 'taskid': filename, 'success': len(washed_smiles_list), 'successrate': len(
        washed_smiles_list) / len(smiles_list) * 100, 'invalid': len(smiles_list) - len(washed_smiles_list)}
    return_result['invalidrate'] = return_result['invalid'] / \
                                   len(smiles_list) * 100
    return_result['filepath'] = '/static/results/csv/' + filename + '.csv'
    return render(request, "services/screening/result_index.html", {
        'result': return_result,
        'filename': filename,
    })


def screening_detail(request, filename, index):
    index = int(index)
    result = pd.read_csv(os.path.join(SITE_ROOT, '../static/results') +
                         '/csv/' + filename + '.csv')
    smiles = result.iloc[index]['smiles']
    mol = Chem.MolFromSmiles(smiles)
    structure = HighlightAtoms(Chem.MolFromSmiles(smiles), highlightAtoms=(), figsize=[300, 300])
    mechanisms_name = ['Colloidal aggregators', 'FLuc inhibitors',
                       'Blue fluorescence', 'Green fluorescence', 'Reactive compounds', 'Promiscuous compounds',
                       'Other assay interference', ]
    uncertainty_name = [item + ' uncertainty' for item in mechanisms_name]
    score_list = result[mechanisms_name].iloc[index]
    score_result = score_list.transpose().reset_index()
    score_result.columns = ['mechanism', 'score']
    uncertainty_list = result[uncertainty_name].iloc[index].values.tolist()
    score_result['uncertainty'] = uncertainty_list
    return_result = score_result
    return_result['uncertainty'] = return_result['uncertainty'].round(6)

    def func_decision(row):
        if row['score'] > 0.5:
            return False
        else:
            return True

    return_result['decision'] = return_result.apply(func_decision, axis=1)
    return_result['attention'] = ['-'] * len(return_result)

    rule_name = [item + '_index' for item in
                 ['Aggregators', 'Fluc', 'Blue_fluorescence', 'Green_fluorescence', 'Reactive', 'Promiscuous',
                  'Other_assay_interference', ]]
    return_result['atom_list'] = result[rule_name].iloc[index].values.tolist()
    return_result['cnt'] = return_result['atom_list'].map(lambda x: len(eval(x)))
    return_result['svg'] = return_result['atom_list'].map(
        lambda x: ''.join([HighlightAtoms(mol, highlightAtoms=item) for item in eval(x)])
        if eval(x) != [] else '')
    return_result['explanation'] = mechanism_explanation
    # basic = Basic([Chem.MolFromSmiles(mol) for mol in [smiles]])
    # basic_prop = basic.CalculateBasicProperties().round(3)

    other_rule_name = [item + '_index' for item in ['ALARM_NMR', 'BMS', 'Chelator_Rule', 'GST_FHs_Rule', 'His_FHs_Rule',
                                                    'Luciferase_Inhibitor_Rule', 'NTD', 'PAINS',
                                                    'Potential_Electrophilic_Rule', 'Lilly']]
    other_rule_result = result[other_rule_name].iloc[index].transpose().reset_index()
    other_rule_result.columns = ['mechanism', 'atom_list']
    other_rule_result['cnt_match'] = other_rule_result['atom_list'].map(lambda x: len(eval(x)))
    other_rule_result['img_match'] = other_rule_result['atom_list'].map(
        lambda x: ''.join([HighlightAtoms(mol, highlightAtoms=item) for item in eval(x)])
        if eval(x) != [] else '')
    other_rule_result['explanation'] = other_rule_explanation
    other_rule_result['mechanism'] = other_rule_display
    # other_rule_result.to_csv('test.csv', index=False)

    return render(request, 'services/evaluation/result_index.html', {
        'smiles': smiles,
        'structure': structure,
        'mechanisms': mechanisms_name,
        'score': return_result['score'].tolist,
        'result': return_result,
        'other_rules': other_rule_result,
        'download': False,
    })


def datasource(request):
    if request.method == 'POST':
        draw = int(request.POST.get('draw'))  # 記錄操作次數
        start = int(request.POST.get('start'))  # 起始位置
        length = int(request.POST.get('length'))  # 每頁長度
        filename = request.POST.get('filename')
        result_path = os.path.join(settings.SITE_ROOT,
                                   'static') + '/uploadfile/molecules/' + filename + '.csv'
        mol_datas = pd.read_csv(result_path)
        results = mol_datas[start: start + length]
        counts = len(mol_datas) - 1
        datas = []
        for idx, row in results.iterrows():
            res = dict()
            mol = Chem.MolFromSmiles(row['SMILES'])
            res['smiles'] = row['SMILES']
            res['structure'] = HighlightAtoms(mol, highlightAtoms=(), figsize=[150, 150])
            res['status'] = row['status']
            datas.append(res)
        response = dict()
        response['draw'] = draw
        response['recordsTotal'] = counts
        response['recordsFiltered'] = counts
        response['data'] = datas
        return HttpResponse(json.dumps(response), content_type='application/json')
    else:
        return render(request, "services/screening/index.html", locals())


def screening_result(request):
    if request.method == 'POST':
        draw = int(request.POST.get('draw'))  # 記錄操作次數
        start = int(request.POST.get('start'))  # 起始位置
        length = int(request.POST.get('length'))  # 每頁長度
        filename = request.POST.get('filename')
        result_path = pd.read_csv(os.path.join(SITE_ROOT, '../static/results') +
                                  '/csv/' + filename + '.csv')
        results = result_path[start: start + length]
        counts = len(result_path)
        datas = []
        for idx, row in results.iterrows():
            res = dict()
            mol = Chem.MolFromSmiles(row['smiles'])
            res['smiles'] = row['smiles']
            res['structure'] = HighlightAtoms(mol, highlightAtoms=(), figsize=[150, 150])
            res['index'] = idx
            # res['status'] = row['status']
            datas.append(res)
        response = dict()
        response['draw'] = draw
        response['recordsTotal'] = counts
        response['recordsFiltered'] = counts
        response['data'] = datas
        return HttpResponse(json.dumps(response), content_type='application/json')
    else:
        return render(request, "services/screening/index.html", locals())


def handleDownload(result):
    result = result.rename(columns=name2dbName)
    rule_name = [item + '_index' for item in
                 ['Aggregators', 'Fluc', 'Blue_fluorescence', 'Green_fluorescence', 'Reactive', 'Promiscuous',
                  'Other_assay_interference',
                  'ALARM_NMR', 'BMS', 'Chelator_Rule', 'GST_FHs_Rule', 'His_FHs_Rule',
                  'Luciferase_Inhibitor_Rule', 'NTD', 'PAINS',
                  'Potential_Electrophilic_Rule', 'Lilly']]
    result[rule_name] = result[rule_name].map(lambda x: len(eval(x)))
    machanism_name = [item + ' uncertainty' for item in
                      ['Colloidal aggregators', 'FLuc inhibitors', 'Blue fluorescence', 'Green fluorescence',
                       'Reactive compounds', 'Promiscuous compounds', 'Other assay interference', ]]
    for _, item in enumerate(machanism_name):
        result[item] = result[item].map(lambda x: 'High' if x < un_threhold[_] else 'Low')
    return result


def downloadCSV(request):
    if request.method == 'POST':
        taskid = request.POST.get('taskid')
        try:
            result = pd.DataFrame(
                list(chemfh.objects.filter(taskid=taskid).values()))
        except chemfh.DoesNotExist:
            result = None
        result = handleDownload(result)
        response = HttpResponse(content_type='text/csv')
        response['Content-Disposition'] = "attachment;filename='ChemFH_result.csv'"
        writer = csv.writer(response)
        writer.writerow(result.columns)
        writer.writerow(list(result.iloc[0]))
        return response
    else:
        return HttpResponseRedirect(reverse('home:index'))


def download(request):
    if request.method == 'POST':
        filename = request.POST.get('filename')
        result = pd.read_csv(os.path.join(SITE_ROOT, '../static/results') +
                             '/csv/' + filename + '.csv')
        result = handleDownload(result)
        response = HttpResponse(result.to_csv(), content_type='text/csv')
        response['Content-Disposition'] = "attachment;filename='ChemFH_result.csv'"
        return response
    else:
        return HttpResponseRedirect(reverse('home:index'))


def downloadPDF(request):
    if request.method == 'POST':
        taskid = request.POST.get('taskid')
        pdf_filepath = os.path.join(
            SITE_ROOT, '../static/results') + '/pdf/' + str(taskid) + '.pdf'
        try:
            result = pd.DataFrame(
                list(chemfh.objects.filter(taskid=taskid).values()))
        except chemfh.DoesNotExist:
            result = None
        result = result.rename(columns=name2dbName)
        gen_pdf(result, pdf_filepath, os.path.join(
            SITE_ROOT, '../static/results'), taskid)
        buffer = open(pdf_filepath, 'rb')
        return FileResponse(buffer, as_attachment=True, filename=pdf_filepath.split('/')[-1])
    else:
        return HttpResponseRedirect(reverse('home:index'))
