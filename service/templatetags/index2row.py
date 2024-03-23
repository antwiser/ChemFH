from django.template import Library

# 将注册类实例化为register对象
register = Library()


@register.filter(name='index2row')
def index2row(datas, arg):
    return datas.iloc[arg].tolist()
