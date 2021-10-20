# -*- coding=utf-8 -*-
import random

def check_permission(topic, user):
    if topic.type != 1:
        return True
    if user.info.job_title == "Researcher":
        return True
    return False

def check_owner(obj, user):
    if obj.owner.id != user.id:
        return False
    return True

def random_str(randomlength=4):
    str = ''
    chars = 'abcdefghijklmnopqrstuvwsyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
    length = len(chars) - 1
    for i in range(randomlength):
        str += chars[random.randint(0, length)]
    return str