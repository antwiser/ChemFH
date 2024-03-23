from django.db import models


# Create your models here.

class chemfh(models.Model):
    taskid = models.CharField(db_column='taskId', max_length=255, primary_key=True)
    smiles = models.TextField(blank=True, null=True)
    aggregators = models.FloatField(db_column='Aggregators', blank=True, null=True)
    fluc = models.FloatField(db_column='FLuc inhibitors', blank=True, null=True)
    blue = models.FloatField(db_column='Blue fluorescence', blank=True, null=True)
    green = models.FloatField(db_column='Green fluorescence', blank=True, null=True)
    reactive = models.FloatField(db_column='Reactive compounds', blank=True, null=True)
    other = models.FloatField(db_column='Other assay interference', blank=True, null=True)
    promiscuous = models.FloatField(db_column='Promiscuous compounds', blank=True, null=True)
    aggregators_un = models.FloatField(db_column='Aggregators uncertainty', blank=True, null=True)
    fluc_un = models.FloatField(db_column='FLuc inhibitors uncertainty', blank=True, null=True)
    blue_un = models.FloatField(db_column='Blue fluorescence uncertainty', blank=True, null=True)
    green_un = models.FloatField(db_column='Green fluorescence uncertainty', blank=True, null=True)
    reactive_un = models.FloatField(db_column='Reactive compounds uncertainty', blank=True, null=True)
    other_un = models.FloatField(db_column='Other assay interference uncertainty', blank=True, null=True)
    promiscuous_un = models.FloatField(db_column='Promiscuous compounds uncertainty', blank=True, null=True)
    aggregators_idx = models.TextField(db_column='Aggregators_index', blank=True, null=True)
    fluc_idx = models.TextField(db_column='Fluc_index', blank=True, null=True)
    blue_idx = models.TextField(db_column='Blue_fluorescence_index', blank=True, null=True)
    green_idx = models.TextField(db_column='Green_fluorescence_index', blank=True, null=True)
    reactive_idx = models.TextField(db_column='Reactive_index', blank=True, null=True)
    other_idx = models.TextField(db_column='Other_assay_interference_index', blank=True, null=True)
    promiscuous_idx = models.TextField(db_column='Promiscuous_index', blank=True, null=True)
    alarm = models.TextField(db_column='ALARM_NMR_index', blank=True, null=True)
    bms = models.TextField(db_column='BMS_index', blank=True, null=True)
    chelator = models.TextField(db_column='Chelator_Rule_index', blank=True, null=True)
    gst = models.TextField(db_column='GST_FHs_Rule_index', blank=True, null=True)
    his = models.TextField(db_column='His_FHs_Rule_index', blank=True, null=True)
    luciferase = models.TextField(db_column='Luciferase_Inhibitor_Rule_index', blank=True, null=True)
    ntd = models.TextField(db_column='NTD_index', blank=True, null=True)
    pains = models.TextField(db_column='PAINS_index', blank=True, null=True)
    potential = models.TextField(db_column='Potential_Electrophilic_Rule_index', blank=True, null=True)
    lilly = models.TextField(db_column='Lilly_index', blank=True, null=True)

    class Meta:
        db_table = 'chemfh'
