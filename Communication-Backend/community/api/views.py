# -*- coding=utf-8 -*-
import re
from django.utils import tree
from rest_framework import permissions
from rest_framework import response
from rest_framework.response import Response
from django.contrib.auth.models import User
from rest_framework import viewsets, status
from rest_framework.serializers import Serializer
from api.models import *
from rest_framework.decorators import action
from django.contrib.auth import authenticate, login, logout
from rest_framework.views import APIView
from api.serializers import *
from api.permissions import IsResearcher
from api.utils import check_permission, random_str, check_owner
import logging
import json
from django.core.mail import send_mail

logger = logging.getLogger(__name__)

class ForgetPasswordViwe(APIView):

    def post(self, request):
        try:
            name=request.POST.get("name")
            password=request.POST.get("password")
            user=User.objects.get(username=name)
            code=request.POST.get("code") 
            if code==request.session["code"]:
                user.set_password(password)
                user.save()
                del request.session["code"] 
                return Response({"status": '0', 'msg': 'Password Reset！'}, status=status.HTTP_200_OK)
            else:
                return Response({"status": '1', 'msg': 'Wrong Code！'}, status=status.HTTP_400_BAD_REQUEST)
        except Exception as e:
            print(repr(e))
            return Response({"status": '1', 'msg': 'Wrong Param！'}, status=status.HTTP_400_BAD_REQUEST)

    def get(self, request):
        try:
            email = str(request.GET.get('email'))
            try:
                user = User.objects.get(username=email)
            except User.DoesNotExist:
                return Response({"status": '1', 'msg': 'Email NOT Exist！'}, status=status.HTTP_400_BAD_REQUEST)
            nick_name = request.GET.get('name')
            if str(user.info.nick_name) != str(nick_name):
                return Response({"status": '1', 'msg': "Email doesn't match username！"}, status=status.HTTP_400_BAD_REQUEST)
            code=random_str()
            request.session["code"]=code 
            email_title = "iGME Community Find Password"
            email_body = "Code：{0}".format(code)
            send_status = send_mail(email_title, email_body, None, [email,])
            return Response({"status": '0', 'msg': 'ok'}, status=status.HTTP_200_OK)
        except Exception as e:
            print(repr(e))
            return Response({"status": '1', 'msg': 'Wrong Param！'}, status=status.HTTP_400_BAD_REQUEST)

class Discuss(APIView):
    def get(self, request):
        if not request.user.is_authenticated:
            return Response({"status": '1', 'msg': 'No login！'}, status=status.HTTP_401_UNAUTHORIZED)
        if not os.path.exists("./discuss.txt"):
            return Response({"status": '0', 'msg': 'ok', 'data' : {}}, status=status.HTTP_200_OK)
        with open("discuss.txt", "r") as f:
            data = []
            for line in f.readlines():
                data.append(line)
            data = ''.join(data)
        return Response({"status": '0', 'msg': 'ok', 'data' : {data}}, status=status.HTTP_200_OK)
    
    def post(self, request):
        if not request.user.is_authenticated:
            return Response({"status": '1', 'msg': 'No login！'}, status=status.HTTP_401_UNAUTHORIZED)
        data = request.data['content']
        with open("discuss.txt", "w") as f:
            f.write(data)
        return Response({"status": '0', 'msg': 'ok', 'data' : {data}}, status=status.HTTP_200_OK)

class LogoutView(APIView):

    def get(self, request):
        user = request.user
        if user.is_authenticated:
            logout(request)  
            return Response({"status": '0', 'msg': 'Logout'}, status=status.HTTP_200_OK)
        return Response({"status": '1', 'msg': 'No login！'}, status=status.HTTP_400_BAD_REQUEST)

class LoginView(APIView):

    def post(self, request):
        if request.user.is_authenticated:
            return Response({"status": '1', 'msg': 'Now is Logined！'}, status=status.HTTP_406_NOT_ACCEPTABLE)
        try:
            email = request.data['email']
            pwd = request.data['password']
            user = authenticate(username=email, password=pwd)
            if user is not None:
                login(request, user)
                user_data = SimUserSerializer(user)
                data = (json.dumps(user_data.data))
                return Response({"status": '0', 'msg': 'Login Succeed！', 'data':{data}}, status=status.HTTP_200_OK)
            else:
                return Response({"status": '1', 'msg': 'Wrong email or password！'}, status=status.HTTP_401_UNAUTHORIZED)
        except Exception as e:
            print(repr(e))
            return Response({"status": '1', 'msg': 'Wrong Param！'}, status=status.HTTP_400_BAD_REQUEST)

class RegisterView(APIView):

    def post(self, request):
        if request.user.is_authenticated:
            return Response({"status": '1', 'msg': 'Logined，please logout out first.'}, status=status.HTTP_406_NOT_ACCEPTABLE)
        try:
            name = request.data['name']
            pwd = request.data['password']
            email = request.data['email']
            job = request.data['job']
            wi = request.data['work_institution']
            rf = request.data['research_field']
            usrimg = request.data['usrimg']
            try:
                user = User.objects.get(username=email)
                return Response({"status": '1', 'msg': 'Email Exist!'}, status=status.HTTP_406_NOT_ACCEPTABLE)
            except User.DoesNotExist:
                try:
                    userinfo = UserInfo.objects.get(nick_name = name)
                    return Response({"status": '1', 'msg': 'username Exist!'}, status=status.HTTP_406_NOT_ACCEPTABLE)
                except UserInfo.DoesNotExist as e:
                    user = User.objects.create_user(username=email, password=pwd, email=email)
                    userinfo = UserInfo.objects.create(owner=user, nick_name = name, job_title = job, userimg=int(usrimg), work_institution=wi, research_field = rf)
                user.save()
                userinfo.save()
                login(request, user)
                authenticate(username=email, password=pwd)
                user_data = UserSerializer(user)
                data = (json.dumps(user_data.data))
                return Response({"status": '0', 'msg': 'Registion Succeed！', 'data':{data}}, status=status.HTTP_200_OK)
        except Exception as e:
            print(repr(e))
            return Response({"status": '1', 'msg': 'Wrong Param！'}, status=status.HTTP_400_BAD_REQUEST)

class UserViewSet(viewsets.ReadOnlyModelViewSet):
    queryset = User.objects.all()
    serializer_class = UserSerializer
    permission_classes = (permissions.IsAuthenticated, )

    @action(detail=True)
    def my_info(self, request, pk=None):
        try:
            info = UserInfo.objects.get(pk=pk)
        except UserInfo.DoesNotExist:
            return Response({"status" : '0', 'msg': 'User Not Exist'}, status=status.HTTP_400_BAD_REQUEST)
        serializer = UserInfoSerializer(info)
        data = json.dumps(serializer.data)
        return Response({"status": '0', 'msg': 'ok', 'data':{data}}, status=status.HTTP_200_OK)

    @action(detail=False)
    def get_info(self, request):
        if not request.user.is_authenticated:
            return Response({"status": '1', 'msg': 'Please Login Firdt！'}, status=status.HTTP_400_BAD_REQUEST)
        try:
            info = UserInfo.objects.get(pk=request.user.info.id)
        except UserInfo.DoesNotExist:
            return Response({"status" : '0', 'msg': 'User Not Exist'}, status=status.HTTP_400_BAD_REQUEST)
        serializer = UserInfoSerializer(info)
        data = json.dumps(serializer.data)
        return Response({"status": '0', 'msg': 'ok', 'data':{data}}, status=status.HTTP_200_OK)

class TopicViewSet(viewsets.ModelViewSet):
    queryset = Topic.objects.all()
    permission_classes = (permissions.IsAuthenticated, IsResearcher,)

    def get_serializer_class(self):
        if self.action == 'list' or self.action == 'retrieve':
            return ReturnTopicSerializer
        else:
            return TopicSerializer

    def perform_create(self, serializer):
            serializer.save(owner=self.request.user)

    def create(self, request, *args, **kwargs):
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        type = int(self.request.data['type'])
        user_job = self.request.user.info.job_title
        if type == 1 and user_job != "Researcher":
            return Response({"status": "0", "msg" : "You have no permission"}, status=status.HTTP_401_UNAUTHORIZED)
        else:
            self.perform_create(serializer)
            headers = self.get_success_headers(serializer.data)
            return Response(serializer.data, status=status.HTTP_201_CREATED, headers=headers)

    @action(detail=True)
    def get_content(self, request, pk=None):
        try:
            topic = Topic.objects.get(pk=pk)
            if not check_permission(topic, request.user):
                return Response({"status": "0", "msg" : "You have no permission"}, status=status.HTTP_401_UNAUTHORIZED)
            # page = self.paginate_queryset(topic)
            # if page is not None:
            serializer = TopicSerializer(topic)
            data = json.dumps(serializer.data)
            return Response({"status": "0", "msg" : "ok", "data" : {data}}, status=status.HTTP_200_OK)
        except Topic.DoesNotExist:
            return Response({"status": "-1", "msg" : "Topic Not Exist！"}, status=status.HTTP_400_BAD_REQUEST)
        except Exception as e:
            print(repr(e))
            return Response({"status": "-1", "msg" : "Wrong Param！"}, status=status.HTTP_400_BAD_REQUEST)

    @action(detail=False)
    def search(self, request):
        title = request.query_params.get('title', None)
        if title is not None:
            lists1 = Topic.objects.filter(title__icontains=title)
            lists2 = Topic.objects.filter(content__icontains=title)
            if request.user.info.job_title != "Researcher":
                U = Topic.objects.filter(type__in=[0,2])
            else:
                U = Topic.objects.all()
            serializers = SimTopicSerializer(lists1.union(lists2).intersection(U), many=True)
            data = json.dumps(serializers.data)
            # page = self.paginate_queryset(lists1.union(lists2))
            # if page is not None:
            #     serializer = SimTopicSerializer(page, many=True)
            #     data = json.dumps(serializer.data)
            #     return self.get_paginated_response({"status": "0", "msg" : "ok", "data" : {data}}, status=status.HTTP_200_OK)
            # else:
            return Response({"status": "0", "msg" : "ok", "data" : {data}}, status=status.HTTP_200_OK)
        return Response({"status": '1', 'msg': 'Searching content cant be empty！'}, status=status.HTTP_400_BAD_REQUEST)

    @action(detail=False)
    def get_topic_list(self, request):
        if request.user.info.job_title != "Researcher":
            topics = Topic.objects.filter(type__in=[0,2])
        else:
            topics = Topic.objects.all()
        serializer = SimTopicSerializer(topics, many=True)
        return Response(serializer.data)
        
    # 获取话题下所有评论
    @action(detail=True)
    def get_comment(self, request, pk=None):
        topic = Topic.objects.get(pk=pk)
        if not check_permission(topic, request.user):
            return Response({"status": "0", "msg" : "You have no permission"}, status=status.HTTP_401_UNAUTHORIZED)
        comments = topic.comment_topic.all()
        # page = self.paginate_queryset(questions)
        # if page is not None:
        serializers = ReturnCommentSerializer(comments, many=True)
        data = json.dumps(serializers.data)
        return Response({'code' : 0, 'msg': 'ok', "data" : {data}}, status=status.HTTP_200_OK)

    # 获取话题下所有附件
    @action(detail=True)
    def get_attachment(self, request, pk=None):
        topic = Topic.objects.get(pk=pk)
        if not check_permission(topic, request.user):
            return Response({"status": "0", "msg" : "You have no permission"}, status=status.HTTP_401_UNAUTHORIZED)
        attachment = topic.attachment_topic.all()
        # page = self.paginate_queryset(attachment)
        # if page is not None:
        serializers = ReturnAttachmentSerializer(attachment, many=True)
        data = json.dumps(serializers.data)
        return Response({'code' : 0, 'msg': 'ok', "data" : {data}}, status=status.HTTP_200_OK)

    def destroy(self, request, *args, **kwargs):
        instance = self.get_object()
        if not check_owner(instance, self.request.user):
            return Response({"status": "0", "msg" : "You have no permission"}, status=status.HTTP_401_UNAUTHORIZED)
        else:
            self.perform_destroy(instance)
            return Response(status=status.HTTP_204_NO_CONTENT)


class CommentViewSet(viewsets.ModelViewSet):
    queryset = Comment.objects.all()
    serializer_class = CommentSerializer
    permission_classes = (permissions.IsAuthenticated,IsResearcher,)

    def perform_create(self, serializer):
        serializer.save(owner=self.request.user)

    def create(self, request, *args, **kwargs):
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        topic = Topic.objects.get(pk=int(self.request.data['topic']))
        if not check_permission(topic, self.request.user):
            return Response({"status": "0", "msg" : "You have no permission"}, status=status.HTTP_401_UNAUTHORIZED)
        else:
            self.perform_create(serializer)
            headers = self.get_success_headers(serializer.data)
            return Response(serializer.data, status=status.HTTP_201_CREATED, headers=headers)

    def destroy(self, request, *args, **kwargs):
        instance = self.get_object()
        if not check_owner(instance, self.request.user):
            return Response({"status": "0", "msg" : "You have no permission"}, status=status.HTTP_401_UNAUTHORIZED)
        else:
            self.perform_destroy(instance)
            return Response(status=status.HTTP_204_NO_CONTENT)

class AnswerViewSet(viewsets.ModelViewSet):
    queryset = Answer.objects.all()
    serializer_class = AnswerSerializer
    permission_classes = (permissions.IsAuthenticated,IsResearcher,)

    @action(detail=True)
    def agree_answer(self, request, pk=None):
        try:
            answer = Answer.objects.get(pk=pk)
            answer.ansagree = answer.ansagree + 1
            answer.save()
            return Response({"status": "0", "msg" : "ok"}, status=status.HTTP_200_OK)
        except Answer.DoesNotExist:
            return Response({"status": "-1", "msg" : "Suggwtion NOT Exist！"}, status=status.HTTP_400_BAD_REQUEST)
        except Exception as e:
            print(repr(e))
            return Response({"status": "-1", "msg" : "Wrong Param！"}, status=status.HTTP_400_BAD_REQUEST)

    def perform_create(self, serializer):
        serializer.save(owner=self.request.user)

    def destroy(self, request, *args, **kwargs):
        instance = self.get_object()
        if not check_owner(instance, self.request.user):
            return Response({"status": "0", "msg" : "You have no permission"}, status=status.HTTP_401_UNAUTHORIZED)
        else:
            self.perform_destroy(instance)
            return Response(status=status.HTTP_204_NO_CONTENT)

class AttachmentViewSet(viewsets.ModelViewSet):
    queryset = Attachment.objects.all()
    serializer_class = AttachmentSerializer
    permission_classes = (permissions.IsAuthenticated,IsResearcher,)

    def perform_create(self, serializer):
        serializer.save(owner=self.request.user)
        
    def create(self, request, *args, **kwargs):
        serializer = self.get_serializer(data=request.data)
        serializer.is_valid(raise_exception=True)
        topic = Topic.objects.get(pk=int(self.request.data['topic']))
        if not check_permission(topic, self.request.user):
            return Response({"status": "0", "msg" : "You have no permission"}, status=status.HTTP_401_UNAUTHORIZED)
        else:
            self.perform_create(serializer)
            headers = self.get_success_headers(serializer.data)
            return Response(serializer.data, status=status.HTTP_201_CREATED, headers=headers)

    def destroy(self, request, *args, **kwargs):
        instance = self.get_object()
        if not check_owner(instance, self.request.user):
            return Response({"status": "0", "msg" : "You have no permission"}, status=status.HTTP_401_UNAUTHORIZED)
        else:
            self.perform_destroy(instance)
            return Response(status=status.HTTP_204_NO_CONTENT)



