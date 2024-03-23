# -*- coding: utf-8 -*-
"""
ChemFH was developed by Jiacai Yi et al. of CBDD GROUP of CSU China.
This project is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
Based on a work at http://home.scbdd.com.
Permissions beyond the scope of this license may be available at http://home.scbdd.com/. If you have any questions, please feel free to contact us.

# @Time    : 2020/10/26 下午4:43
# @Author  : Jiacai Yi
# @FileName: index2value.py
# @E-mail  ：1076365758@qq.com
"""
from django.template import Library

# 将注册类实例化为register对象
register = Library()


@register.filter(name='index2value')
def index2value(array, arg):
    return array[arg]
