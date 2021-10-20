# -*- coding=utf-8 -*-
from django.contrib.auth.models import User
from django.db.models.query import QuerySet
from rest_framework import serializers
from api.models import *


class UserSerializer(serializers.ModelSerializer):
    info = serializers.PrimaryKeyRelatedField(many=False, read_only=True)
    email = serializers.ReadOnlyField(source="username")
    
    class Meta:
        model = User
        fields = ('id', 'email', 'info')


class UserInfoSerializer(serializers.ModelSerializer):
    email = serializers.ReadOnlyField(source='owner.username')

    class Meta:
        model = UserInfo
        fields = ('id', 'email', 'nick_name', 'job_title', 'work_institution', 'research_field', 'userimg')


class SimUserInfoSerializer(serializers.ModelSerializer):

    class Meta:
        model = UserInfo
        fields = ('id', 'owner', 'nick_name', 'userimg')

class SimUserSerializer(serializers.ModelSerializer):
    info = SimUserInfoSerializer()
    email = serializers.ReadOnlyField(source="username")

    class Meta:
        model = User
        fields = ('id', 'email', 'info')


class ReturnTopicSerializer(serializers.ModelSerializer):
    owner = SimUserSerializer()

    class Meta:
        model = Topic
        fields = ('id', 'owner', 'title', 'created', 'type')


class TopicSerializer(serializers.ModelSerializer):
    owner = serializers.ReadOnlyField(source='owner.info.nick_name')

    class Meta:
        model = Topic
        fields = ('id', 'owner', 'title', 'content', 'created', 'type')


class SimTopicSerializer(serializers.ModelSerializer):
    class Meta:
        model = Topic
        fields = ('id', 'title')


class AnswerSerializer(serializers.ModelSerializer):
    owner = serializers.ReadOnlyField(source='owner.info.nick_name')

    class Meta:
        model = Answer
        fields = ('id', 'ansto', 'owner', 'description', 'created', 'ansagree', 'keep')

class ReturnAnswerSerializer(serializers.ModelSerializer):
    owner = SimUserSerializer()
    ansto = serializers.ReadOnlyField(source='ansto.id')

    class Meta:
        model = Answer
        fields = ('id', 'ansto', 'owner', 'description', 'created', 'ansagree', 'keep')

class CommentSerializer(serializers.ModelSerializer):
    owner = serializers.ReadOnlyField(source='owner.info.nick_name')
    topic = serializers.PrimaryKeyRelatedField(queryset=Topic.objects.all())
    ansto = serializers.ReadOnlyField(source='ansto.id')

    class Meta:
        model = Comment
        fields = ('id', 'ansto', 'owner', 'topic', 'content', 'created')

class ReturnCommentSerializer(serializers.ModelSerializer):
    owner = SimUserSerializer()
    topic = serializers.PrimaryKeyRelatedField(queryset=Topic.objects.all())
    ansto = serializers.ReadOnlyField(source='ansto.id')

    class Meta:
        model = Comment
        fields = ('id', 'ansto', 'owner', 'topic', 'content', 'created')

class AttachmentSerializer(serializers.ModelSerializer):
    owner = serializers.ReadOnlyField(source='owner.info.nick_name')
    topic = serializers.PrimaryKeyRelatedField(queryset=Topic.objects.all())

    class Meta:
        model = Attachment
        fields = ('id', 'owner', 'topic', 'data', 'created')

class ReturnAttachmentSerializer(serializers.ModelSerializer):
    owner = SimUserSerializer()
    topic = serializers.PrimaryKeyRelatedField(queryset=Topic.objects.all())

    class Meta:
        model = Attachment
        fields = ('id', 'owner', 'topic', 'data', 'created')
