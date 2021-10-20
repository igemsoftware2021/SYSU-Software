from django.db import models
from django.contrib.auth.models import User
import os
from django.conf import settings

class Topic(models.Model):
    # Create Time
    created = models.DateTimeField(auto_now_add=True)
    # Creater
    owner = models.ForeignKey(to=User, related_name='user_topic', on_delete=models.CASCADE)
    # Create Title
    title = models.CharField(max_length=50, unique=False)
    # Create Content
    content = models.TextField()
    # Type 0 Topic  1 Private Topic 2 Discuss
    type = models.IntegerField(default=0)

    @property
    def name(self):
        return "Topic"

    class Meta:
        unique_together=("title", "type")

class Comment(models.Model):
    # Create Time
    created = models.DateTimeField(auto_now_add=True)
    # Answer Question
    ansto = models.ForeignKey(to='self', null=True, blank=True, related_name='comment_comment', on_delete=models.CASCADE)
    # Commenter
    owner = models.ForeignKey(to=User, related_name='commenter', on_delete=models.CASCADE)
    # Comment
    content = models.TextField()
    # Comment To Topic
    topic = models.ForeignKey(to=Topic, related_name='comment_topic', on_delete=models.CASCADE)

    # Type of Comment Topic
    @property
    def type(self):
        return self.topic.type

    @property
    def name(self):
        return "Comment"


class Answer(models.Model):
    # Create Time
    created = models.DateTimeField(auto_now_add=True)
    # Answer Question
    ansto = models.ForeignKey(to='self', null=True, blank=True, related_name='answers_question', on_delete=models.CASCADE)
    # Answerer
    owner = models.ForeignKey(to=User, related_name='answerser', on_delete=models.CASCADE)
    # Content
    description = models.TextField()
    # Agree
    ansagree = models.IntegerField(default=0)
    # Keep unique
    keep = models.CharField(default='', max_length=100)

    # Type of Answer Topic
    @property
    def type(self):
        return self.ansto.topic.type

    @property
    def name(self):
        return "Answer"



class Attachment(models.Model):
    # Create Time
    created = models.DateTimeField(auto_now_add=True)
    # Owner
    owner = models.ForeignKey(to=User, related_name='attachmenter', on_delete=models.CASCADE)
    # Path
    data = models.FileField(default='', upload_to='attachments')
    # Attachment Topic
    topic = models.ForeignKey(to=Topic, related_name='attachment_topic', on_delete=models.CASCADE)

    @property
    def name(self):
        return "Attachment"

    # Type of Attachment Topic
    @property
    def type(self):
        return self.topic.type


class UserInfo(models.Model):
    # User Job
    job_title = models.CharField(max_length=20, choices=(('Researcher', '调查员'), ('Other Job', '其他')), default='male')
    # User image
    userimg = models.IntegerField(default=0)
    # User Info Owner
    owner = models.OneToOneField(to=User, related_name='info', on_delete=models.CASCADE)
    # User Company
    work_institution = models.CharField(max_length=50, blank = True)
    # User Researcher File
    research_field = models.CharField(max_length=50, blank = True)
    # User Nick Name
    nick_name = models.CharField(max_length=50, unique=True)
    
