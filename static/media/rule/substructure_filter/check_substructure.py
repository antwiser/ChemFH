'''
Description: the functions
Author: Kotori Y
Date: 2020-10-24 16:08:49
LastEditors: Kotori Y
LastEditTime: 2020-11-04 11:23:49
FilePath: \ChemFH\ChemFH\substructure_filter\check_substructure.py
AuthorMail: kotori@cbdd.me
'''

from collections import namedtuple
from functools import wraps

try:
    from load_pattern import loadpkl
except Exception:
    from .load_pattern import loadpkl


class InvalidInputError(Exception):
    pass


def withEndpoint(endpoint="PAINS"):
    def check(func):
        @wraps(func)
        def wrapper(mols, endpoint):
            pattl = loadpkl(endpoint)
            for mol in mols:
                res = func(mol=mol, pattl=pattl)
                yield res

        return wrapper

    return check


def checkValidMol(func):
    @wraps(func)
    def wrapper(mol, **kwgrs):
        if mol is not None:
            return func(mol, **kwgrs)
        else:
            return {}

    return wrapper


@withEndpoint()
@checkValidMol
def checkPattl(mol, pattl):
    checkRes = {}
    for patt, smart in pattl:
        atoms = mol.GetSubstructMatches(patt)
        if atoms:
            checkRes[smart] = atoms

    return checkRes

# if '__main__' == __name__:
#     from rdkit import Chem
#
#     smis = [
#         'OC1=C[C-]2[OH+]C(c3ccc(O)c(O)c3)=C(O)C=C2C(O)=C1',
#         'O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12',
#         'OCC(O)C(O)C(O)C(O)CO',
#         'O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12',
#     ]
#     mols = [Chem.MolFromSmiles(smi) for smi in smis]
#
#     print(list(checkPattl(mols=mols, endpoint="PAINS")))
