'''
Description: Generate and save the molecule object of SMARTS
Author: Kotori Y
Date: 2020-10-24 16:38:14
LastEditors: Kotori Y
LastEditTime: 2020-10-24 19:26:35
FilePath: \ChemFH\ChemFH\substructure_filter\load_pattern.py
AuthorMail: kotori@cbdd.me
'''

import csv
import os
import _pickle as cPickle
import gzip
from rdkit.Chem import AllChem as Chem

__all__ = ['loadpkl']
_dir = os.path.abspath(os.path.join(os.path.dirname("__file__")))
print(_dir)

def _Generatepkl(endpoint):
    """  
    *Internal Use Only*
    
    the pkl file in this package, storing the rdkit.Chem.rdchem.Mol object,
    was generated under the environment whose rdkit version is '2020.03.5'.
    Since, the file may can't be successfully loaded. This function is designedexit for
    re-generating a pkl file under this situation.
    
    :param endpoint: the name of file
    :type endpoint: str
        
    :return: None
    """
    file = os.path.join(_dir, 'static/media/rule/data/SMARTS', f'{endpoint}.txt')
    # try:
    #     os.mkdir(os.path.join(_dir, 'data/Pattern'))
    #     # print(os.path.join(_dir, 'data/Pattern'))
    # except FileExistsError:
    #     pass

    with open(file, 'r', encoding='utf-8') as f_obj:
        lines = csv.reader(f_obj, delimiter='\t')
        next(lines)
        lines = tuple(lines)
    f_obj.close()

    for line in lines:
        patt = Chem.MolFromSmarts(line[-1])
        line[0] = patt

    out = cPickle.dumps(lines, protocol=-1)
    outfile = os.path.join(_dir, 'static/media/rule/data/Pattern', f'{endpoint}.pkl.gz')
    with gzip.open(outfile, 'wb') as f_out:
        f_out.write(out)
    f_out.close()


def loadpkl(endpoint):
    """ 
    
    loading the specific pkl file which contain the 'Rejected' and 'Accepted' SMARTS
    in rdkit.Chem.rdchem.Mol object format to avoid repeated reading SMARTS by 'MolFromSmarts'
    
    :param endpoint: the endpoint of pkl file meant
    :type endpoint: str
    
    
    :return: whose element ia also a list with two elements: 
    0: the name of SMARTS, 1: Original SMARTS
    :rtype: list   
        
    """
    filename = os.path.join(_dir, 'static/media/rule/data/Pattern', f'{endpoint}.pkl.gz')
    try:
        pattl = cPickle.load(gzip.open(filename, 'rb'))
    except:
        _Generatepkl(endpoint)
        return loadpkl(endpoint)
    return pattl

# if '__main__' == __name__:
#     files = os.listdir('../data/SMARTS')
#     print(files)
#
#     for file in files:
#         endpoint = os.path.splitext(file)[0]
#         patt = loadpkl(endpoint)
#         print(patt)
#         # _Generatepkl(endpoint)
