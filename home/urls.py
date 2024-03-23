# -*- coding: utf-8 -*-
"""
ChemFH was developed by Jiacai Yi et al. of CBDD GROUP of CSU China.
This project is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
Based on a work at http://home.scbdd.com.
Permissions beyond the scope of this license may be available at http://home.scbdd.com/. If you have any questions, please feel free to contact us.

# @Time    : 2020/10/22 下午6:46
# @Author  : Jiacai Yi
# @FileName: urls.py
# @E-mail  ：1076365758@qq.com
"""

from django.urls import path, re_path

from . import views
from django.views.generic import TemplateView

app_name = 'home'
urlpatterns = [
    path('', views.index, name='index'),
    path('contact/', views.contact, name="contact"),
    path('documentation/', views.documentation, name="documentation"),
    path('documentation/<str:path>/', views.documentationPath, name="documentationpath"),
    path('publication/', views.publication, name="publication"),
    path('apis/', views.apis, name="apis"),
    # path('apis/<str:path>/', views.apispath, name="apispath"),
]
