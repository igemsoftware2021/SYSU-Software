# Generated by Django 3.2.7 on 2021-10-20 07:39

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('api', '0005_alter_answer_ansto'),
    ]

    operations = [
        migrations.AlterField(
            model_name='answer',
            name='ansto',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='answers_question', to='api.answer'),
        ),
    ]
