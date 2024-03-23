# -*- coding: utf-8 -*-
"""
ChemFH was developed by Jiacai Yi et al. of CBDD GROUP of CSU China.
This project is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
Based on a work at http://home.scbdd.com.
Permissions beyond the scope of this license may be available at http://home.scbdd.com/. If you have any questions, please feel free to contact us.

# @Time    : 2020/10/23 下午1:51
# @Author  : Jiacai Yi
# @FileName: forms.py
# @E-mail  ：1076365758@qq.com
"""
from django import forms
from django.contrib.auth.models import User
import re

mechanism_choices = (
    ("0", "colloidal aggregators"),
    ("1", "Fluc inhibitors"),
    ("2", "blue fluorescence"),
    ("3", "green fluorescence"),
    ("4", "reactive compounds"),
    ("5", "other assay interferences"),
    ("6", "promiscuous compounds"),
)


class EvaluationSingleForm(forms.Form):
    smiles = forms.CharField(label="SMILES", required=True, widget=forms.TextInput(
        attrs={
            'class': 'form-control form-control-lg',
        }
    ), help_text='Please enter a valid SMILES!')
    method = forms.CharField(label='method', initial='1', widget=forms.HiddenInput(attrs={}))
    mechanism = forms.MultipleChoiceField(label="Mechanism", choices=mechanism_choices,
                                          widget=forms.CheckboxSelectMultiple(
                                              attrs={
                                                  'class': 'list-unstyled',
                                              }
                                          ))
