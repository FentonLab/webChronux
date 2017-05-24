# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models, migrations
from django.conf import settings


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
    ]

    operations = [
        migrations.CreateModel(
            name='AnalysisDetail',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('name', models.CharField(max_length=255, blank=True, null=True)),
                ('description', models.CharField(max_length=255)),
                ('timeWindow', models.FloatField(blank=True, null=True)),
                ('frequencyBandWidth', models.FloatField(blank=True, null=True)),
                ('numTapers', models.IntegerField(blank=True, null=True)),
                ('stepSize', models.FloatField(blank=True, null=True)),
                ('padding', models.IntegerField(blank=True, null=True)),
                ('upperFrequency', models.FloatField(blank=True, null=True)),
                ('lowerFrequency', models.FloatField(blank=True, null=True)),
                ('analysisOutputPath', models.CharField(max_length=512, blank=True, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='AnalysisPlot',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('name', models.CharField(max_length=128)),
                ('plotPath', models.CharField(max_length=255)),
                ('plotFileName', models.CharField(max_length=128)),
                ('analysisDetail', models.ForeignKey(to='chronux.AnalysisDetail')),
            ],
        ),
        migrations.CreateModel(
            name='AnalysisResultFile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('name', models.CharField(max_length=128)),
                ('filePath', models.CharField(max_length=255)),
                ('resultFileName', models.CharField(max_length=128)),
                ('analysisDetail', models.ForeignKey(to='chronux.AnalysisDetail')),
            ],
        ),
        migrations.CreateModel(
            name='BipolarMontage',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
            ],
        ),
        migrations.CreateModel(
            name='Datafile',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('filePath', models.CharField(max_length=512)),
            ],
        ),
        migrations.CreateModel(
            name='Epoch',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('epochValue1', models.FloatField()),
                ('epochValue2', models.FloatField()),
                ('epochValue3', models.FloatField()),
                ('epochValue4', models.FloatField()),
                ('epochType', models.CharField(max_length=256)),
                ('datafile', models.ForeignKey(to='chronux.Datafile')),
            ],
        ),
        migrations.CreateModel(
            name='FileType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('name', models.CharField(max_length=256)),
                ('description', models.CharField(max_length=512)),
            ],
        ),
        migrations.CreateModel(
            name='JobStatusCode',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('code', models.CharField(max_length=10)),
                ('description', models.CharField(max_length=255, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='LeadList',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('name', models.CharField(max_length=50)),
                ('description', models.CharField(max_length=255, blank=True, null=True)),
            ],
        ),
        migrations.CreateModel(
            name='LeadListChannel',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('channelId', models.CharField(max_length=10)),
                ('baseChannelId', models.CharField(max_length=10)),
                ('channelOrder', models.IntegerField(blank=True, null=True)),
                ('xCoordinate', models.FloatField(blank=True, null=True)),
                ('yCoordinate', models.FloatField(blank=True, null=True)),
                ('zCoordinate', models.FloatField(blank=True, null=True)),
                ('isInterior', models.NullBooleanField()),
                ('leadList', models.ForeignKey(to='chronux.LeadList')),
            ],
        ),
        migrations.CreateModel(
            name='Project',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('name', models.CharField(max_length=256)),
                ('description', models.CharField(max_length=512, blank=True, null=True)),
                ('user', models.ForeignKey(blank=True, null=True, to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.CreateModel(
            name='SubmittedJob',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('name', models.CharField(max_length=512)),
                ('description', models.CharField(max_length=512, blank=True, null=True)),
                ('submittedOn', models.DateTimeField(blank=True, null=True)),
                ('downloadDataFileLink', models.CharField(max_length=512, blank=True, null=True)),
                ('completedTime', models.DateTimeField(blank=True, null=True)),
                ('analysisDetail', models.ForeignKey(blank=True, null=True, to='chronux.AnalysisDetail')),
                ('jobStatusCode', models.ForeignKey(to='chronux.JobStatusCode')),
                ('submittedBy', models.ForeignKey(to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.CreateModel(
            name='SubmittedJobType',
            fields=[
                ('id', models.AutoField(verbose_name='ID', primary_key=True, serialize=False, auto_created=True)),
                ('name', models.CharField(max_length=10)),
                ('description', models.CharField(max_length=255)),
            ],
        ),
        migrations.AddField(
            model_name='submittedjob',
            name='submittedJobType',
            field=models.ForeignKey(to='chronux.SubmittedJobType'),
        ),
        migrations.AddField(
            model_name='datafile',
            name='fileType',
            field=models.ForeignKey(blank=True, null=True, to='chronux.FileType'),
        ),
        migrations.AddField(
            model_name='datafile',
            name='project',
            field=models.ForeignKey(to='chronux.Project'),
        ),
        migrations.AddField(
            model_name='bipolarmontage',
            name='channelId1',
            field=models.ForeignKey(related_name='channelId1', to='chronux.LeadListChannel'),
        ),
        migrations.AddField(
            model_name='bipolarmontage',
            name='channelId2',
            field=models.ForeignKey(related_name='channelId2', to='chronux.LeadListChannel'),
        ),
        migrations.AddField(
            model_name='analysisdetail',
            name='project',
            field=models.ForeignKey(to='chronux.Project'),
        ),
    ]
