# -*- coding=utf-8 -*-
from rest_framework.routers import DefaultRouter
from api.views import *
from django.conf.urls import url, include

router = DefaultRouter()
router.register(r'user', UserViewSet)
router.register(r'answers', AnswerViewSet)
router.register(r'topics', TopicViewSet)
router.register(r'comments', CommentViewSet)
router.register(r'attachments', AttachmentViewSet)

urlpatterns = [
    url(r'login/', LoginView.as_view()),
    url(r'register/', RegisterView.as_view()),
    url(r'discuss/', Discuss.as_view()),
    url(r'logout/', LogoutView.as_view()),
    url(r'resetpwd/', ForgetPasswordViwe.as_view()),
    url(r'^', include(router.urls))
]