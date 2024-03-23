from rdkit import Chem
from .pretreat import StandardMol
import numpy as np


def check_SMILES_valid(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol if Chem.MolFromSmiles(smiles) and smiles != "" else False


def check_smiles_input(input):
    if np.sum([item != '' for item in input]) == 0:
        return False
    return input


def check_wash_mol(smiles):
    # smiles: list
    mol = check_SMILES_valid(smiles)
    if mol:
        mol = StandardMol(mol)
        newSmiles = Chem.MolToSmiles(mol)
        smiles = newSmiles if newSmiles else smiles
        return smiles
    else:
        return None


def check_input_smiles(data):
    smiles = [data.SMILES] if isinstance(data.SMILES, str) else data.SMILES
    # batch process, smiles is a smiles list
    smiles = check_smiles_input(smiles)
    if isinstance(smiles, bool) and not smiles:
        return {
            'success': 0,
            'info': 'Input Invalid!'
        }
    smiles = [check_wash_mol(item) for item in smiles]
    if (sum(item is not None for item in smiles)) == 0:
        return {
            'success': 0,
            'info': 'Molecule Invalid!'
        }
    smiles = list(filter(lambda s: s and s.strip(), smiles))
    return {
        'success': 1,
        'info': smiles
    }


def wash_input_mol(data, issmiles=True, invalidStr='', returnInvalidIdx=False):
    smiles = [data] if isinstance(data, str) else data
    result = []
    flagIdx = []
    result2 = []
    for item in smiles:
        try:
            mol = Chem.MolFromSmiles(item)
            mol = StandardMol(mol)
            if issmiles:
                returnItem = Chem.MolToSmiles(mol)
            else:
                returnItem = mol
            result.append(returnItem)
            flagIdx.append(True)
            result2.append(returnItem)
        except:
            flagIdx.append(False)
            result2.append(invalidStr)
    if returnInvalidIdx:
        return result, flagIdx
    if invalidStr:
        return result2
    return result
