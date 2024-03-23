from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import render, get_object_or_404
from django.urls import reverse
from django.views import generic
from django.utils import timezone
import pandas as pd


def index(request):
    return render(request, "home/index.html")


def contact(request):
    return render(request, 'contact/index.html', locals())


def documentation(request):
    return render(request, 'documentation/index.html', locals())


def documentationPath(request, path):
    return render(request, 'documentation/' + path + '/', locals())


def publication(request):
    from .static_data import publicationData, webData
    data = pd.DataFrame(publicationData)
    webData = pd.DataFrame(webData)
    return render(request, 'pub/index.html', {'datas': data, 'webData': webData})


def apis(request):
    return render(request, 'apis/ori_index.html', locals())


def apispath(request, path):
    return render(request, 'apis/' + path + '/', locals())


