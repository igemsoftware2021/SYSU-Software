from rest_framework import permissions
from api.models import Topic


class IsResearcher(permissions.BasePermission):
    def has_object_permission(self, request, view, obj):
        if int(obj.type) == 1 and request.user.info.job_title != "Researcher":
            return False
        elif int(obj.type != 2) and obj.name == "Question":
            return False
        elif int(obj.type != 2) and obj.name == "Answer":
            return False
        else:
            return True