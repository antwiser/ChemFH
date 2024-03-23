# -*- coding: utf-8 -*-
"""
ChemFH was developed by Jiacai Yi et al. of CBDD GROUP of CSU China.
This project is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
Based on a work at http://home.scbdd.com.
Permissions beyond the scope of this license may be available at http://home.scbdd.com/. If you have any questions, please feel free to contact us.

# @Time    : 2020/10/22 下午10:14
# @Author  : Jiacai Yi
# @FileName: urls.py
# @E-mail  ：1076365758@qq.com
"""
from django.urls import path

from . import views

app_name = 'service'
urlpatterns = [
    path('evaluation/index/', views.evaluationIndex, name="evaluation"),
    path('evaluation/cal', views.evaluationCal, name="evaluationCal"),
    path('screening/index/', views.filterIndex, name="screening"),
    path('screening/cal', views.screeningCal, name="filterCal"),
    path('screening-result/', views.screening_result, name="screening_result"),
    path('screening/<str:filename>/<int:index>', views.screening_detail, name="screening_detail"),
    path('result-datasource/', views.datasource, name="result_datasource"),
    path('result/<str:filename>', views.result_file, name='result'),
    path('downpdf', views.downloadPDF, name="downpdf"),
    path('downcsv', views.downloadCSV, name="downcsv"),
    path('download', views.download, name="download"),
]


