'''
Description: 
Author: Kotori Y
Date: 2020-10-24 16:03:49
LastEditors: Kotori Y
LastEditTime: 2020-11-03 08:56:20
FilePath: \ChemFH\ChemFH\substructure_filter\__init__.py
AuthorMail: kotori@cbdd.me
'''

from collections.abc import Iterable
from functools import partial
import multiprocessing as mp

import os
import sys
from django.conf import settings

sys.path.append(os.path.join(settings.SITE_ROOT, 'static') + '/media/rule/substructure_filter')

try:
    from check_substructure import checkPattl
    from visualization import highlightAtoms
except Exception:
    from .check_substructure import checkPattl
    from .visualization import highlightAtoms


class FrequentHitterFilter:

    def __init__(self, nJobs=1):
        self.nJobs = nJobs if nJobs > 0 else None

    def screening(self, mols, endpoint):
        func = partial(checkPattl, endpoint=endpoint)
        mols = mols if isinstance(mols, Iterable) else (mols,)
        # pool = mp.Pool(self.nJobs)
        # res = pool.map_async(func, mols).get()
        # pool.close()
        # pool.join()
        res = list(func(mols))
        return res


# if '__main__' == __name__:
#     from rdkit import Chem
#     import pandas as pd
#
#     smis = [
#         'OC1=C[C-]2[OH+]C(c3ccc(O)c(O)c3)=C(O)C=C2C(O)=C1',
#         'O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12',
#         'OCC(O)C(O)C(O)C(O)CO',
#     ]
#     mols = [Chem.MolFromSmiles(smi) for smi in smis]
#
#     Filter = FrequentHitterFilter()
#     res = Filter.screening(mols, endpoint='PAINS')
#     print(res)
