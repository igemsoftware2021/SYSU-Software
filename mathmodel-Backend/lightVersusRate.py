#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
from py_expression_eval import Parser
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import mean_squared_error
import datetime
import random

def lightVersusRateOfReaction(a1, d1, k1, a2, d2, a3, d3, k2, k3, s, etol, lkFuncStr, k1m, k2m, A1, A2, l12, l23, D1, D2):
    '''
    Parameters
    ----------
    a1 d1 k1 a2 d2 a3 d3 k2 k3 : number
        Corresonding to the parameters a_i，d_i，k_i, i=1,2,3
    s : number
        底物浓度，相对etol不能太大） 
    etol(单种酶（如酶1）的总浓度，假设了1:1:1混合，注意不是三种酶浓度的总和)
    lkFuncStr 光照-聚集函数（因变量是聚集程度比例k，自变量是光照强度，要求用户输入k关于光照强度的表达式，且该表达式中用‘*’表示乘号且不得省略，即使用计算机编程中数学表达式的通常记法。例如，“k=3l+1”应当表达为“3*l+1”）
    k1m（即k_1^M） k2m（即k_2^M）
    A1（酶1活性位点表面积） A2（酶2活性位点表面积）
    l12（即l_12） l23（即l_23） 
    D1 D2 (扩散系数都取作10**-9)
    '''
    
    parser = Parser()
    lkFunc = parser.parse(lkFuncStr)
    variables = lkFunc.variables()
    if len(variables) != 1:
        raise Exception("Incorrect number of independent variables in ", lkFuncStr)
    variable = variables[0]
    
    l = 0
    step = 0.1
    v = []
    x = np.arange(0, 100, step)
    for l in x:
        l12 *= 10
        l23 *= 10
        A1 *= 100
        A2 *= 100
        '''
        a1 *= 10 ** -3
        a2 *= 10 ** -3
        a3 *= 10 ** -3
        '''
        k = lkFunc.evaluate({variable : l})
        b=math.pow(12*math.pi*etol,-1/3)
        K1=k1*s*(1/l12-1/b)/(4*(k1m+s)*D1*math.pi)
        i1=(a1*k1*s*(k/2 - 1)*(d2 + k2)*(k - 2))/(a2*(2*d1*k2 + 2*k1*k2 - 2*d1*k*k2 - 2*a1*k1*s + 2*a1*k2*s - 2*k*k1*k2 + a1*k*k1*s - 2*a1*k*k2*s)) - (k*(K1*a2*d1*k2 - a1*k1*k2*s - a1*d2*k1*s + K1*a2*k1*k2 + K1*a1*a2*k2*s))/(2*a2*k2*(d1 + k1 + a1*s))
        K2=k2*i1/(k2m+i1)/D2/4/math.pi*(1/l23-1/b)
        i2=(2*a2*i1*k2*(k/2 - 1)*(d3 + k3)*(k - 1))/(a3*(2*d2*k3 + 2*k2*k3 - 2*a2*i1*k2 + 2*a2*i1*k3 - d2*k*k3 - k*k2*k3 + 2*a2*i1*k*k2 - a2*i1*k*k3)) - (k*(K2*a3*d2*k3 + K2*a3*k2*k3 - a2*d3*i1*k2 - a2*i1*k2*k3 + K2*a2*a3*i1*k3))/(2*a3*k3*(d2 + k2 + a2*i1))
        v.append((a3*etol*k3*(2*d3*i2 + 2*i2*k3 + 2*a3*i2**2 + 2*K2*a3*i2 + K2*d3*k + K2*k*k3))/(2*(d3 + k3 + a3*i2)*(d3 + k3 + K2*a3 + a3*i2)))
    fittedValues, fittedFunc = funcFit(x, v)
    graphFilename = drawGraph(x, fittedValues)
    return fittedFunc, graphFilename, fittedValues

def funcFit(x, v):
    rmses = []
    coefs = []
    values = []
    
    for i in range(2, 11):
        fit = np.polyfit(x, v, i)
        value = np.polyval(fit, x)
        values.append(value)
        coefs.append(fit)
        rmses.append(math.sqrt(mean_squared_error(v, value)))
    
    minValue = rmses[0]
    minIndex = 0
    i = 1
    while i < len(rmses):
        if rmses[i] < minValue:
            minValue = rmses[i]
            minIndex = i
        i += 1
        
    coef = coefs[minIndex]
    value = values[minIndex]
    
    return value, coef

def drawGraph(x, y):
    plt.figure(figsize = (8, 5))
    plt.plot(x, y)
    plt.xlabel('light energy density (W/m^2)')
    plt.ylabel('reaction rate (mM/s)')
    filename = datetime.datetime.now().strftime('%H.%M.%S.%f-') + str(random.randint(0, 1000000)) +'.png'
    plt.savefig(filename)
    return filename
