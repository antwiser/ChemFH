# Generated by Django 4.1 on 2024-01-14 02:30

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='chemfh',
            fields=[
                ('taskid', models.CharField(db_column='taskId', max_length=255, primary_key=True, serialize=False)),
                ('smiles', models.TextField(blank=True, null=True)),
                ('aggregators', models.FloatField(blank=True, db_column='Aggregators', null=True)),
                ('fluc', models.FloatField(blank=True, db_column='FLuc inhibitors', null=True)),
                ('blue', models.FloatField(blank=True, db_column='Blue fluorescence', null=True)),
                ('green', models.FloatField(blank=True, db_column='Green fluorescence', null=True)),
                ('reactive', models.FloatField(blank=True, db_column='Reactive compounds', null=True)),
                ('other', models.FloatField(blank=True, db_column='Other assay interference', null=True)),
                ('promiscuous', models.FloatField(blank=True, db_column='Promiscuous compounds', null=True)),
                ('aggregators_un', models.FloatField(blank=True, db_column='Aggregators uncertainty', null=True)),
                ('fluc_un', models.FloatField(blank=True, db_column='FLuc inhibitors uncertainty', null=True)),
                ('blue_un', models.FloatField(blank=True, db_column='Blue fluorescence uncertainty', null=True)),
                ('green_un', models.FloatField(blank=True, db_column='Green fluorescence uncertainty', null=True)),
                ('reactive_un', models.FloatField(blank=True, db_column='Reactive compounds uncertainty', null=True)),
                ('other_un', models.FloatField(blank=True, db_column='Other assay interference uncertainty', null=True)),
                ('promiscuous_un', models.FloatField(blank=True, db_column='Promiscuous compounds uncertainty', null=True)),
                ('aggregators_idx', models.FloatField(blank=True, db_column='Aggregators_index', null=True)),
                ('fluc_idx', models.FloatField(blank=True, db_column='Fluc_index', null=True)),
                ('blue_idx', models.FloatField(blank=True, db_column='Blue_fluorescence_index', null=True)),
                ('green_idx', models.FloatField(blank=True, db_column='Green_fluorescence_index', null=True)),
                ('reactive_idx', models.FloatField(blank=True, db_column='Reactive_index', null=True)),
                ('other_idx', models.FloatField(blank=True, db_column='Other_assay_interference_index', null=True)),
                ('promiscuous_idx', models.FloatField(blank=True, db_column='Promiscuous_index', null=True)),
                ('alarm', models.FloatField(blank=True, db_column='ALARM_NMR_index', null=True)),
                ('bms', models.FloatField(blank=True, db_column='BMS_index', null=True)),
                ('chelator', models.FloatField(blank=True, db_column='Chelator_Rule_index', null=True)),
                ('gst', models.FloatField(blank=True, db_column='GST_FHs_Rule_index', null=True)),
                ('his', models.FloatField(blank=True, db_column='His_FHs_Rule_index', null=True)),
                ('luciferase', models.FloatField(blank=True, db_column='Luciferase_Inhibitor_Rule_index', null=True)),
                ('ntd', models.FloatField(blank=True, db_column='NTD_index', null=True)),
                ('pains', models.FloatField(blank=True, db_column='PAINS_index', null=True)),
                ('potential', models.FloatField(blank=True, db_column='Potential_Electrophilic_Rule_index', null=True)),
            ],
            options={
                'db_table': 'chemfh',
            },
        ),
    ]