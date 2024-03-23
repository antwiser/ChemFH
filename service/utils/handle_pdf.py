from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, Image, TableStyle
from reportlab.platypus import Frame, ListFlowable, ListItem, BalancedColumns
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle, _baseFontNameB, _baseFontName
from reportlab.lib.pagesizes import A4
from reportlab.rl_config import defaultPageSize
from reportlab.lib.units import inch
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib import colors
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT, TA_CENTER, TA_RIGHT
from reportlab.pdfbase import pdfmetrics
from reportlab.lib.colors import Color, yellow, green, red, black, blue
from reportlab.graphics.shapes import Circle
from reportlab.lib.colors import tan, green
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.charts.barcharts import VerticalBarChart
import os
from django.conf import settings
import PIL.Image as Img
import numpy as np
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from rdkit import Chem
import cairosvg
import math


def HighlightAtoms(mol, highlightAtoms, figsize=[400, 200], kekulize=True):
    """This function is used for showing which part of fragment matched the SMARTS by the id of atoms.

    :param mol: The molecule to be visualized
    :type mol: rdkit.Chem.rdchem.Mol
    :param highlightAtoms: The atoms to be highlighted
    :type highlightAtoms: tuple
    :param figsize: The resolution ratio of figure
    :type figsize: list
    :return: a figure with highlighted molecule
    :rtype: IPython.core.display.SVG

    """

    def _revised(svg_words):
        """
        """
        svg_words = svg_words.replace(
            'stroke-width:2px', 'stroke-width:1.5px').replace(
            'fonts-size:17px', 'fonts-size:15px').replace(
            'stroke-linecap:butt', 'stroke-linecap:square').replace(
            'fill:#FFFFFF', 'fill:none').replace(
            'svg:', '')
        return svg_words

    mc = Chem.Mol(mol.ToBinary())

    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        rdDepictor.Compute2DCoords(mc)
    drawer = rdMolDraw2D.MolDraw2DSVG(*figsize)
    drawer.DrawMolecule(mc, highlightAtoms=highlightAtoms)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    # It seems that the svg renderer used doesn't quite hit the spec.
    # Here are some fixes to make it work in the notebook, although I think
    # the underlying issue needs to be resolved at the generation step
    return _revised(svg)


def get_Para(text, style):
    if type(text) == str:
        return Paragraph(text, style)
    elif type(text) == list:
        return [Paragraph(str(item), style) for item in text]
    # elif type(text) == np.ndarray:
    #     return np.fromfunction(lambda x: Paragraph(x, style), (text.shape[0], text.shape[1]))


pdfmetrics.registerFont(
    TTFont('Lobster', os.path.join(settings.SITE_ROOT, 'static') + '/home/fonts/Lobster-Regular.ttf'))

PAGE_HEIGHT = defaultPageSize[1]
PAGE_WIDTH = defaultPageSize[0]

A4_width, A4_height = A4

header_style = ParagraphStyle(
    name='header',
    fontSize=11,
    fontName=_baseFontNameB,
)

header_center_style = ParagraphStyle(
    name='header',
    fontSize=11,
    fontName=_baseFontNameB,
    alignment=TA_CENTER,
)

decision_style = ParagraphStyle(name='decision_style',
                                fontName=_baseFontName,
                                fontSize=10,
                                leading=12,
                                spaceBefore=6,
                                alignment=TA_CENTER)

justify_style = ParagraphStyle(name="justify_style",
                               justifyBreaks=1,
                               alignment=TA_JUSTIFY,
                               leading=16)

styles = getSampleStyleSheet()
pageinfo = ""


def myFirstPage(canvas, doc):
    canvas.saveState()
    canvas.drawImage(os.path.join(settings.SITE_ROOT, 'static') + '/icon.png', x=inch,
                     y=PAGE_HEIGHT - inch - 26,
                     width=160, height=48, mask='auto')
    # canvas.setFont('Lobster', 16)
    # canvas.drawCentredString(inch * 2 + 40, PAGE_HEIGHT - inch, doc.title)
    canvas.setFont('Times-Roman', 12)
    canvas.drawRightString(PAGE_WIDTH - inch, PAGE_HEIGHT - inch - 40, doc.smiles)
    canvas.setStrokeColor(Color(0, 0, 0, alpha=0.5))
    canvas.line(inch, PAGE_HEIGHT - 120, PAGE_WIDTH - inch, PAGE_HEIGHT - 120)
    canvas.setFont('Times-Roman', 9)
    canvas.drawString(inch, 0.75 * inch, "Page %d %s" % (doc.page, pageinfo))
    canvas.setFont('Times-Roman', 80)
    canvas.setFillAlpha(0.05)
    canvas.rotate(45)
    canvas.drawCentredString((PAGE_WIDTH / 2) * 1.5, 0, "ChemFH")
    # canvas.drawString(inch, 0.75 * inch, "First Page")
    canvas.restoreState()


def myLaterPages(canvas, doc):
    canvas.saveState()
    canvas.setFont('Times-Roman', 9)
    canvas.drawString(inch, 0.75 * inch, "Page %d %s" % (doc.page, pageinfo))
    canvas.setFont('Times-Roman', 80)
    canvas.setFillAlpha(0.05)
    canvas.rotate(45)
    canvas.drawCentredString((PAGE_WIDTH / 2) * 1.5, 0, "ChemFH")
    canvas.restoreState()


basic_property_comment = get_Para([
    'Contain hydrogen atoms. Optimal:100~600',
    'Van der Waals volume',
    'Density = MW / Volume',
    'Number of hydrogen bond acceptors. Optimal:0~12',
    'Number of hydrogen bond donors. Optimal:0~7',
    'Number of rotatable bonds. Optimal:0~11',
    'Number of rings. Optimal:0~6',
    'Number of atoms in the biggest ring. Optimal:0~18',
    'Number of heteroatoms. Optimal:1~15',
    'Formal charge. Optimal:-4 ~4',
    'Number of rigid bonds. Optimal:0~30',
    'Flexibility = nRot /nRig',
    'Optimal: ≤ 2',
    'Topological Polar Surface Area. Optimal:0~140',
], style=styles['BodyText'])
basic_property_property = get_Para([
    'Molecular Weight', 'Volume', 'Density', 'nHA', 'nHD', 'nRot', 'nRing', 'MaxRing', 'nHet', 'fChar', 'nRig',
    'Flexibility', 'Stereo Centers', 'TPSA'], style=styles['BodyText'])

value_style = ParagraphStyle(
    name='value_style',
    fontName=_baseFontName,
    spaceBefore=6,
    alignment=TA_CENTER
)

left_style = ParagraphStyle(
    name="left_style",
    alignment=TA_LEFT
)

un_threhold = [
    0.0010918042703688,
    0.004291940348845,
    0.0009872382265281,
    0.0023240092127062,
    0.0016870113071982,
    0.0061948686233965,
    9.781904685012926e-05,
]


def get_final_decision(list_datas):
    result = ['-'] * len(list_datas)
    for _, item in enumerate(list_datas):
        if item > 0.5:
            result[_] = Paragraph('<font color="red">●</font>', decision_style)
        else:
            result[_] = Paragraph('<font color="green">●</font>', decision_style)
    return result


def get_final_uncertainty(list_data):
    result = ['-'] * len(list_data)
    for _, item in enumerate(list_data):
        if item < un_threhold[_]:
            # 可信
            result[_] = Paragraph(
                '<font color="#79A1FB">High-confidence</font>',
                decision_style)
        else:
            result[_] = Paragraph('<font color="#F07875">Low-confidence</font>', decision_style)
    return result


IMAGES_FORMAT = ['.png', '.jpg']  # 图片格式


def image_compose(IMAGE_SIZE=200, IMAGE_ROW=1, IMAGE_COLUMN=3, IMAGE_SAVE_PATH='', IMAGES_PATH='', image_names=None):
    if image_names is None:
        image_names = []
    to_image = Img.new('RGBA', (IMAGE_COLUMN * IMAGE_SIZE, IMAGE_ROW * IMAGE_SIZE))  # 创建一个新图
    # 循环遍历，把每张图片按顺序粘贴到对应位置上
    for y in range(1, IMAGE_ROW + 1):
        for x in range(1, IMAGE_COLUMN + 1):
            from_image = Img.open(IMAGES_PATH + image_names[IMAGE_COLUMN * (y - 1) + x - 1]).resize(
                (IMAGE_SIZE, IMAGE_SIZE), Img.LANCZOS)
            to_image.paste(from_image, ((x - 1) * IMAGE_SIZE, (y - 1) * IMAGE_SIZE))
    return to_image.save(IMAGE_SAVE_PATH)  # 保存新图


def go(smiles, basic_property, model_prediction, filepath, row_data, filename, other_rules1, other_rules2):
    # doc = SimpleDocTemplate("result.pdf", pagesize=A4, topMargin=110)
    pass


def image_delete(image_names, images_path):
    for item in image_names:
        filename = images_path + item
        if os.path.exists(filename):
            os.remove(filename)


def gen_pdf(csv_file, filepath, dirpath, taskid):
    row_data = csv_file[:1]
    smiles = row_data['smiles'][0]
    model_prediction_header = np.array(
        [get_Para('Mechanism', header_style), get_Para('Pred. Score', header_center_style),
         get_Para('Decision', header_center_style), get_Para('Uncertainty', header_center_style)])
    mechanisms_name = ['Colloidal aggregators', 'FLuc inhibitors', 'Blue fluorescence', 'Green fluorescence',
                       'Reactive compounds', 'Promiscuous compounds',
                       'Other assay interference', ]
    uncertainty_name = [item + ' uncertainty' for item in mechanisms_name]
    rule_name = [item + '_index' for item in
                 ['Aggregators', 'Fluc', 'Blue_fluorescence', 'Green_fluorescence', 'Reactive', 'Promiscuous',
                  'Other_assay_interference']]
    mechanism = get_Para(mechanisms_name, style=styles['BodyText'])
    score_list = row_data[mechanisms_name].iloc[0].values.tolist()
    score = get_Para(["{:.3f}".format(item) for item in score_list], style=value_style)
    decision = get_final_decision(score_list)
    uncertainty = get_final_uncertainty(row_data[uncertainty_name].iloc[0].values.tolist())
    # uncertainty = get_Para(,style=value_style)
    atom_list = row_data[rule_name].iloc[0].values.tolist()
    cnt = get_Para([len(eval(item)) for item in atom_list], style=value_style)
    prediction = np.vstack((model_prediction_header, np.array([mechanism, score, decision, uncertainty]).T)).tolist()

    other_rule_name = ['ALARM_NMR', 'BMS', 'Chelator_Rule', 'GST_FHs_Rule', 'His_FHs_Rule',
                       'Luciferase_Inhibitor_Rule', 'NTD', 'PAINS',
                       'Potential_Electrophilic_Rule', 'Lilly']
    other_header = np.array(
        [get_Para('Mechanism', header_style), get_Para('Number of Rule', header_style),
         get_Para('Number of Match Structure', header_center_style), ])
    other_rule_result = row_data[[item + '_index' for item in other_rule_name]].iloc[0].values.tolist()
    rule_mechanism = get_Para(['ALARM NMR Rule', 'BMS Rule', 'Chelator Rule', 'GST FHs Rule', 'His FHs Rule',
                               'Luciferase Inhibitor Rule', 'NTD Rule', 'PAINS Rule',
                               'Potential Electrophilic Rule', 'Lilly Medchem Rules'], style=styles['Normal'])
    num_of_rule = get_Para([75, 176, 55, 34, 19, 3, 105, 480, 119, 273],
                           style=value_style)
    rule_cnt = get_Para([len(eval(item)) for item in other_rule_result], style=value_style)
    rule_smarts = get_Para([item for item in other_rule_result], style=value_style)
    other = np.vstack((other_header, np.array([rule_mechanism, num_of_rule, rule_cnt]).T)).tolist()

    mechanism_rule = np.vstack((
        np.array(
            [get_Para('Mechanism', header_style), get_Para('Number of alerts', header_center_style),
             get_Para('Number of matched alerts', header_center_style), ]),
        np.array([mechanism,
                  get_Para([26, 18, 26, 16, 3, 6, 7],
                           style=value_style),
                  cnt]).T
    )).tolist()

    doc = SimpleDocTemplate(filepath, pagesize=A4)
    doc.title = "ChemFH"
    doc.smiles = row_data['smiles'].values[0]
    story = []

    story.append(Spacer(width=0, height=50))

    story.append(Paragraph("1. Molecular structure", style=styles["Heading2"]))
    png_dir = dirpath + '/img/' + str(taskid)
    if not os.path.exists(png_dir):
        os.makedirs(png_dir)
    png_path = png_dir + '/' + 'smiles.png'
    svg = HighlightAtoms(Chem.MolFromSmiles(smiles), highlightAtoms=(), figsize=[400, 400])
    cairosvg.svg2png(bytestring=svg, write_to=png_path, dpi=300)
    story.append(Image(png_path, width=200, height=200, hAlign="CENTER"), )

    # story.append(Spacer(width=0, height=10))

    story.append(Paragraph("2. Results of ChemFH model", style=styles["Heading2"]))

    story.append(Spacer(width=0, height=10))

    story.append(Paragraph(
        "By harnessing the power of high-quality databases and advanced Multi-task DMPNN architectures, ChemFH consistently identifies seven distinct types of frequent hitters, spanning Colloidal aggregators, FLuc inhibitors, Blue/Green fluorescence, Reactive compounds, Promiscuous compounds, and Other assay interferences. This robust detection capability significantly elevates the efficiency of drug research and development processes. For a comprehensive understanding of these seven mechanisms, please explore the detailed information available at https://chemfh.scbdd.com/documentation/#/mechanism.",
        style=justify_style))
    # story.append(Paragraph(
    #     "<font color='green'>●</font>: Higher credibility.  <font color='yellow'>●</font>: Lower credibility",
    #     style=styles["Normal"]))
    story.append(Spacer(width=0, height=10))

    model_prediction_table = Table(prediction, spaceBefore=2,
                                   colWidths=[doc.width * 0.40, doc.width * 0.20, doc.width * 0.20, doc.width * 0.20])
    model_prediction_table.setStyle(TableStyle([
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('ALIGN', (1, 0), (1, -1), 'CENTER'),
        ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black),
        ('BOX', (0, 0), (-1, -1), 0.25, colors.black),
    ]))
    story.append(model_prediction_table)

    story.append(Spacer(width=0, height=10))

    story.append(Paragraph(
        "<strong>Tips</strong>: The corresponding relationships of the two labels "
        "are as "
        "follows: <font "
        "color='green'>● Pass</font> (Non-FH compound); <font color='red'>● Reject</font> (FH compound)",
        style=justify_style))

    story.append(Spacer(width=0, height=20))

    story.append(Paragraph("3. Results of rules", style=styles["Heading2"]))

    story.append(Spacer(width=0, height=10))

    story.append(Paragraph("3.1 Seven Mechanism Rules", style=styles["Heading3"]))

    mechanism_rule_table = Table(mechanism_rule, spaceBefore=2,
                                 colWidths=[doc.width * 0.35, doc.width * 0.25, doc.width * 0.40])
    mechanism_rule_table.setStyle(TableStyle([
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('ALIGN', (1, 0), (1, -1), 'CENTER'),
        ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black),
        ('BOX', (0, 0), (-1, -1), 0.25, colors.black),
    ]))
    story.append(mechanism_rule_table)

    # 下面是每个机制的规则匹配结果
    img_names = []
    all_rule_name = [item + '_index' for item in ['Aggregators', 'Fluc', 'Blue_fluorescence', 'Green_fluorescence',
                                                  'Reactive', 'Promiscuous', 'Other_assay_interference',
                                                  'ALARM_NMR', 'BMS', 'Chelator_Rule', 'GST_FHs_Rule',
                                                  'His_FHs_Rule', 'Luciferase_Inhibitor_Rule', 'NTD', 'PAINS',
                                                  'Potential_Electrophilic_Rule', 'Lilly']]
    display_name = ['Colloidal aggregators', 'FLuc inhibitors',
                    'Blue fluorescence', 'Green fluorescence', 'Reactive compounds', 'Promiscuous compounds',
                    'Other assay interference', 'ALARM NMR Rule', 'BMS Rule', 'Chelator Rule', 'GST FHs Rule',
                    'His FHs Rule', 'Luciferase Inhibitor Rule', 'NTD Rule', 'PAINS Rule',
                    'Potential Electrophilic Rule', 'Lilly Medchem Rules']

    if os.path.exists(png_dir):
        pass
    else:
        os.makedirs(png_dir)
    for _, item in (row_data[all_rule_name].iloc[0].items()):
        item = eval(item)
        if len(item) == 1:
            png_path = png_dir + '/' + str(_) + '.png'
            svg = HighlightAtoms(Chem.MolFromSmiles(smiles), highlightAtoms=item[0], figsize=[300, 300])
            cairosvg.svg2png(bytestring=svg, write_to=png_path, dpi=300)
        elif len(item) > 1:
            image_names = []
            for key, value in enumerate(item):
                sub_png_path = png_dir + '/' + str(_) + '#' + str(key) + '.png'
                image_names.append(str(_) + '#' + str(key) + '.png')
                svg = HighlightAtoms(Chem.MolFromSmiles(smiles), highlightAtoms=value, figsize=[300, 300])
                cairosvg.svg2png(bytestring=svg, write_to=sub_png_path, dpi=300)
            save_path = png_dir + '/' + str(_) + '.png'
            image_compose(IMAGE_ROW=math.ceil(len(item) / 3), IMAGE_COLUMN=3 if len(item) >= 3 else len(item),
                          IMAGE_SAVE_PATH=save_path,
                          IMAGES_PATH=png_dir + '/',
                          image_names=image_names)  # 调用函数
            # image_delete(image_names, image_path)
        else:
            continue
        img_names.append(str(_) + '.png')

    img_list = ['aggregator', 'luciferae inhibitor', 'Fluorescent', 'Fluorescent', 'reactive', 'Promiscuous', 'Other']
    introduction_list = [
        "Colloidal aggregation in High Throughput Screening (HTS) is driven by the critical aggregation concentration (CAC), "
        "forming 30–600 nm aggregators. They bind nonspecifically to proteins, destabilizing enzymes. Lower protein stability increases susceptibility to inhibition."
        " The strong binding affinity, often with picomolar KD values, poses challenges, especially for larger aggregators (>250 nm) with poorer absorption. "
        "Early drug discovery aims to filter out colloidal aggregators to prevent assay interference.",

        "Luciferase-based bioluminescence assays, renowned for high sensitivity, are widely used in High Throughput Screening (HTS). "
        "Firefly luciferase (FLuc) from Photinus pyralis is a common choice in HTS due to its intricate enzymatic mechanism. "
        "FLuc utilizes ATP and D-luciferin to produce a luciferyl-adenylate intermediate, which, upon nucleophilic attack by oxygen, "
        "generates green and red light. Despite FLuc's utility, interference can occur in HTS assays due to FLuc inhibitors, "
        "manifesting as specific inhibition or nonspecific interference through processes like enzyme denaturation or the inner-filter effect in "
        "fluorescent detection.",

        "Fluorescence involves a fluorophore absorbing light, elevating an electron's energy state. Widely used in assays, fluorophores pose interference in high-throughput screening (HTS) due to quenching and autofluorescence. Quenching diminishes emitted light, seen in molecular and luciferase assays. Autofluorescence, natural light emission, overlaps with detection fluorophores, causing interference. Managing these interferences is crucial for accurate HTS results in various applications.",

        "Fluorescence involves a fluorophore absorbing light, elevating an electron's energy state. Widely used in assays, fluorophores pose interference in high-throughput screening (HTS) due to quenching and autofluorescence. Quenching diminishes emitted light, seen in molecular and luciferase assays. Autofluorescence, natural light emission, overlaps with detection fluorophores, causing interference. Managing these interferences is crucial for accurate HTS results in various applications.",

        "Chemical reactive compounds (CRCs) modify proteins or assay reagents, complicating interference mechanisms. "
        "Some inert compounds transform into CRCs within cells, posing challenges in exploration. "
        "CRCs induce interference via redox cycling-induced hydrogen peroxide (H2O2) and electrophilic functionality. "
        "Compounds generating H2O2 in reducing agents may falsely affect targets by oxidizing cysteine residues, "
        "leading to shared concerns about false-positives. Illustratively, isoquinoline-1,3,4-trione derivatives inactivate caspase-3 "
        "by oxidizing its catalytic cysteine to sulfonic acid (–SO3H) through reactive oxygen species.",

        "The promiscuity of compounds is a broad concept. Here, promiscuous compounds refer to those binding specifically to different macromolecular targets. "
        "Unlikely for a molecule to occupy multiple, diverse binding sites through simple, noncovalent interactions typical in specific screening hits, "
        "a promiscuous activity profile usually implies a nonspecific mode of action. These interactions may involve unintended targets, "
        "triggering adverse reactions and safety issues.",

        "Homogeneous proximity assays, vital for high-throughput screening in drug discovery, include ALPHA, FRET, and TR-FRET. Despite their value, these technologies face compound-mediated interferences like signal attenuation (quenching, inner-filter effects, light scattering), signal emission (autofluorescence), and disruption of affinity capture components (affinity tags, antibodies). These interferences can yield false-positive results in assays."
    ]
    other_introduction_list = [
        "The ALARM NMR (Assay to Label Reactants in Mixtures Nuclear Magnetic Resonance) rule is introduced in the "
        "research [1]. This methodology serves as an efficient and reliable experimental approach for identifying "
        "reactive false positives in biochemical screening assays using Nuclear Magnetic Resonance (NMR). <strong>(75 "
        "substructures)</strong>",
        "The BMS rule, which stands for Biologically Meaningful Substructure, is an empirical process designed to "
        "enhance the effectiveness of high-throughput screening (HTS) deck filters. The approach is detailed in the "
        "research [2]. The BMS rule aims to identify and incorporate biologically relevant substructures into "
        "screening libraries, thereby improving the likelihood of discovering compounds with genuine biological "
        "activity. <strong>(176 substructures)</strong>",
        "The Chelator Rule is introduced in the research [3]. This rule serves as a guiding principle in the design "
        "and utilization of chelator-based fragment libraries for the purpose of targeting metalloproteinases. "
        "Metalloproteinases play crucial roles in various biological processes, and their dysregulation is associated with several diseases, "
        "making them attractive therapeutic targets. <strong>(55 substructures)</strong>",
        "The GST/GSH FH (Glutathione S-Transferase/Glutathione Frequent Hitters) filter rule is a methodology "
        "described in the research [4]. This rule is devised to identify small molecules that frequently interact "
        "with the glutathione S-transferase (GST) enzyme and its glutathione (GSH) substrate. GST is a critical "
        "enzyme involved in detoxification processes, and its interaction with GSH is essential for cellular "
        "function. <strong>(34 substructures)</strong>",
        "The His-tagged protein FH (Frequent Hitters) filter rule is a method outlined in the research [5]. This "
        "rule is designed to identify small molecules that frequently interact with proteins tagged with histidine ("
        "His-tagged proteins) in high-throughput screening assays conducted using AlphaScreen technology. <strong>(19 "
        "substructures)</strong>",
        "The Luciferase Inhibitor Rule is a concept presented in the research [6]. This rule aims to address "
        "challenges in high-throughput screening (HTS) assays that use luciferase as a reporter gene. "
        "Luciferase-based assays are widely employed in drug discovery to assess biological activity, "
        "and false-positive hits can be a significant concern. <strong>(3 substructures)</strong>",
        "The NTD (Neglected Tropical Diseases) Rule is introduced in the research [7]. This rule encapsulates "
        "insights and guidelines derived from the process of assembling screening libraries specifically tailored for "
        "drug discovery efforts targeting neglected tropical diseases. Neglected diseases, often afflicting populations "
        "in resource-limited regions, present unique challenges and necessitate targeted drug discovery approaches. "
        "<strong>(105 substructures)</strong>",
        "The PAINS rules (Pan Assay Interference Compounds) constitute a set of substructure filters designed to "
        "exclude interfering compounds from screening libraries. These rules aim to identify compounds that may lead "
        "to misleading results in bioassays, thereby enhancing the accuracy and reliability of screening experiments. "
        "The research [8] provides a comprehensive description of the development and optimization of PAINS rules, "
        "along with their application in the exclusion of potentially false-positive compounds from screening libraries. "
        "<strong>(480 substructures)</strong>",
        "The Potential Electrophilic Rule emerges as a groundbreaking paradigm, "
        "challenging conventional norms in the realm of chemical structural alerts. "
        "Unveiled in the research [9], this rule beckons a bold departure from the commonplace, "
        "ushering in a new era of avant-garde perspectives in toxicological assessments. "
        "<strong>(119 substructures)</strong>",
        "Lilly Medchem Rules, developed by Eli Lilly Company, "
        "represent a comprehensive set of guidelines designed to identify and eliminate compounds that may interfere with biological assays "
        "in the drug discovery process [10]. These rules were developed over an 18-year period and consist of 275 "
        "criteria aimed at refining "
        "screening sets and improving the efficiency of drug development. The rules cover various aspects, including chemical reactivity (e.g., acyl halides), "
        "potential interference with assay measurements (such as fluorescence, absorbance, or quenching), activities that may cause damage to proteins "
        "(oxidizers, detergents), instability issues (e.g., latent aldehydes), and a lack of druggability (compounds lacking both oxygen and nitrogen). "
        "<strong>(273 substructures)</strong>",
    ]
    for idx, (_, item) in enumerate((row_data[all_rule_name[:7]].iloc[0].items())):
        item = eval(item)
        story.append(Paragraph(str(idx + 1) + ') ' + display_name[idx], style=styles["Heading4"]))
        story.append(Spacer(width=0, height=10))
        F = [[Paragraph('Mechanism', style=header_center_style),
              Paragraph('Conceptual diagram', style=header_center_style)], [Paragraph(
            introduction_list[idx],
            style=justify_style),
                 Image(os.path.join(
                     settings.SITE_ROOT,
                     'static') + '/explanation/img/' +
                       img_list[idx] + '.png', width=200,
                       height=150)]]
        tmp = Table(F)
        tmp.setStyle(TableStyle([
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
            ('ALIGN', (1, 0), (1, -1), 'CENTER'),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black),
            ('BOX', (0, 0), (-1, -1), 0.25, colors.black),
        ]))
        story.append(tmp)
        story.append(Paragraph("<font color='red'>Match substructure: </font>", style=styles["Heading4"]))
        if len(item) != 0:
            filepath = png_dir + '/' + str(_) + '.png'
            story.append(Image(filepath, width=150 * len(item), height=150, hAlign='CENTER'))
        else:
            story.append(Paragraph("-", style=header_center_style))

    # 下面是10条其他规则的匹配结果

    story.append(Spacer(width=0, height=20))

    story.append(Paragraph("3.2 Other Rules", style=styles["Heading3"]))

    other_table = Table(other, spaceBefore=2,
                        colWidths=[doc.width * 0.35, doc.width * 0.25, doc.width * 0.40])
    other_table.setStyle(TableStyle([
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('ALIGN', (1, 0), (1, -1), 'CENTER'),
        ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black),
        ('BOX', (0, 0), (-1, -1), 0.25, colors.black),
    ]))
    story.append(other_table)

    # story.append(Spacer(width=0, height=30))
    #
    # story.append(Paragraph("4. Results of substructure matching", style=styles["Heading2"]))
    #
    # story.append(Spacer(width=0, height=10))
    #
    # story.append(Paragraph("4.1 Seven Mechanisms", style=styles["Heading3"]))

    # story.append(Spacer(width=0, height=10))
    # story.append(Paragraph("4.2 Other Rules", style=styles["Heading3"]))

    for idx, (_, item) in enumerate((row_data[all_rule_name[7:]].iloc[0].items())):
        item = eval(item)
        story.append(Paragraph(str(idx + 1) + ') ' + display_name[idx + 7], style=styles["Heading4"]))
        story.append(Spacer(width=0, height=10))
        story.append(Paragraph(other_introduction_list[idx], style=justify_style))
        story.append(Spacer(width=0, height=10))
        story.append(Paragraph("<font color='red'>Match substructure: </font>", style=styles["Heading4"]))
        if len(item) != 0:
            filepath = png_dir + '/' + str(_) + '.png'
            count = len(item)
            row = math.ceil(count / 3)
            col = 3
            width = 150 * count if count < 3 else 450
            height = 150 * row
            story.append(Image(filepath, width=width, height=height, hAlign='CENTER'))
        else:
            story.append(Paragraph("-", style=header_center_style))

    story.append(Spacer(width=0, height=30))

    story.append(Paragraph("References:", style=styles['Heading2']))

    references = [
        "Baell JB, Holloway GA. New Substructure Filters for Removal of Pan Assay Interference Compounds (PAINS) "
        "from Screening Libraries and for Their Exclusion in Bioassays, Journal of Medicinal Chemistry 2010;53:2719-2740.",
        "Pearce BC, Sofia MJ, Good AC et al. An Empirical Process for the Design of High-Throughput Screening Deck "
        "Filters, Journal of Chemical Information and Modeling 2006;46:1060-1068.",
        "Brenke JK, Salmina ES, Ringelstetter L et al. Identification of Small-Molecule Frequent Hitters of "
        "Glutathione S-Transferase–Glutathione Interaction, SLAS Discovery 2016;21:596-607.",
        "Schorpp K, Rothenaigner I, Salmina E et al. Identification of Small-Molecule Frequent Hitters from "
        "AlphaScreen High-Throughput Screens, Journal of Biomolecular Screening 2014;19:715-726.",
        "Huth JR, Mendoza R, Olejniczak ET et al. ALARM NMR: A Rapid and Robust Experimental Method To Detect "
        "Reactive False Positives in Biochemical Screens, Journal of the American Chemical Society 2005;127:217-224.",
        "Ghosh D, Koch U, Hadian K et al. Luciferase Advisor: High-Accuracy Model To Flag False Positive Hits in "
        "Luciferase HTS Assays, Journal of Chemical Information and Modeling 2018;58:933-942.",
        "Agrawal A, Johnson SL, Jacobsen JA et al. Chelator Fragment Libraries for Targeting Metalloproteinases, "
        "ChemMedChem 2010;5:195-199.",
        "Brenk R, Schipani A, James D et al. Lessons Learnt from Assembling Screening Libraries for Drug Discovery "
        "for Neglected Diseases, ChemMedChem 2008;3:435-444.",
        "Sushko I, Salmina E, Potemkin VA et al. ToxAlerts: A Web Server of Structural Alerts for Toxic Chemicals and Compounds with Potential Adverse Reactions, Journal of Chemical Information and Modeling 2012;52:2310-2316.",
        "Bruns RF, Watson IA. Rules for Identifying Potentially Reactive or Promiscuous Compounds, Journal of "
        "Medicinal Chemistry 2012;55:9763-9772.",
    ]
    story.append(ListFlowable([Paragraph(item, style=justify_style) for item in references], bulletType='1',
                              ListStyle=justify_style))

    doc.build(story, onFirstPage=myFirstPage, onLaterPages=myLaterPages)
