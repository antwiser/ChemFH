import hashlib
# import pyecharts.options as opts
# from pyecharts.charts import Radar


def getMD5(request: str) -> str:
    md5_obj = hashlib.md5()
    md5_obj.update(request.encode('utf-8'))
    return md5_obj.hexdigest()

#
# def getPlot():
#     v1 = [[4300, 10000, 28000, 35000, 50000, 19000]]
#     v2 = [[5000, 14000, 28000, 31000, 42000, 21000]]
#
#     radar_view = Radar(init_opts=opts.InitOpts(bg_color="#CCCCCC")).add_schema(
#         schema=[
#             opts.RadarIndicatorItem(name="销售（sales）", max_=6500),
#             opts.RadarIndicatorItem(name="管理（Administration）", max_=16000),
#             opts.RadarIndicatorItem(name="信息技术（Information Technology）", max_=30000),
#             opts.RadarIndicatorItem(name="客服（Customer Support）", max_=38000),
#             opts.RadarIndicatorItem(name="研发（Development）", max_=52000),
#             opts.RadarIndicatorItem(name="市场（Marketing）", max_=25000),
#         ],
#         splitarea_opt=opts.SplitAreaOpts(
#             is_show=True, areastyle_opts=opts.AreaStyleOpts(opacity=1)
#         ),
#         textstyle_opts=opts.TextStyleOpts(color="#fff"),
#     ).add(
#         series_name="预算分配（Allocated Budget）",
#         data=v1,
#         linestyle_opts=opts.LineStyleOpts(color="#CD0000"),
#     ).add(
#         series_name="实际开销（Actual Spending）",
#         data=v2,
#         linestyle_opts=opts.LineStyleOpts(color="#5CACEE"),
#     ).set_series_opts(label_opts=opts.LabelOpts(is_show=False)).set_global_opts(
#         title_opts=opts.TitleOpts(title="基础雷达图"), legend_opts=opts.LegendOpts()
#     ).render_embed()
#
#     return radar_view
