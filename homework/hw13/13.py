# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 14:28:13 2022

@author: YXY
"""

import numpy as np
import matplotlib.pyplot as plt
import math as ma


def sch_random(N, M = 1, a = 16807, b = 0, m = 2**31 - 1, seed = 1):
    #N为生成个数，M为生成间隔
    q, r = m // a, m % a    
    #得到p，r
    while N>0:
        N -= 1
        for j in range(M):
            seed = a * (seed % q) -r * (seed // q) 
            #进行schrage方法
            if seed < 0:
                seed += m
        yield seed/m
        #单位化schrage,现在支持一次性生成可用无穷次的迭代了
        

def MHM(fp, N, delta, seed = 114514, f_random = sch_random):#固定抽样范围（0，inf）
    i = 0#循环变量
    eff = 0#计算抽样效率
    res = np.zeros((N),dtype = 'double')
    x0 = 5 * delta
    random = f_random(ma.inf, seed = seed)
    temp1 = 0
    while i < N/100:#热化
        temp1 = (next(random) - 0.5) * delta
        a = fp(x0+temp1) / fp(x0)
        if x0 + temp1 < 0:
            continue
        if a > 1:
            x0 += temp1
        else:
            if a > next(random):
                x0 += temp1
            else:
                continue
        i += 1
    i = 0
    while i < N:#正式模拟
        temp1 = (next(random) - 0.5) * delta
        if x0 + temp1 < 0:
            continue
        a = fp(x0+temp1) / fp(x0)
        if a > 1:
            x0 += temp1
            res[i] = x0
            eff += 1
        else:
            eff += 2
            if a > next(random):
                x0 += temp1
                res[i] = x0
            else:
                continue
        i += 1
    return res, N / eff


def freq(c):#频数统计
    d = {}
    for i in c:
        if round(i,1) in d:
            d[round(i,1)] += 1
        else:
            d[round(i,1)] = 1
    x , y = np.zeros((len(d))),np.zeros((len(d)))
    j = 0
    for i in sorted(d.items(),key = lambda x:x[0]):
        x[j] = i[0]
        y[j] = i[1]
        j += 1
    return x,y,len(c)#返回长度计算抽样效率


def fp1(alpha,beta):
    def f(x):
        return x**(alpha-1.0)*ma.exp(-x/beta)
    return f


def fp2(alpha,beta):
    def f(x):
        return (x-alpha*beta)**2*x**(alpha-1.0)*ma.exp(-x/beta)
    return f
if __name__ == "__main__":
    intres1 = []
    effres1 = []
    intres2 = []
    effres2 = []
    M0 = 1000000
    alpha = 3.5
    beta = 1.5
    i = 0.1
    x = [y for y in range(10)]
    st1 = [alpha*beta*beta for y in range(10)]
    st2 = [0.5 for y in range(10)]
    for j in range(10):
        a, b = MHM(fp1(alpha,beta),M0,i)
        intres1 += [((a-alpha*beta)**2/M0).sum()]
        effres1 += [b]
        a, b = MHM(fp2(alpha,beta),M0,i)
        intres2 += [1/((1/(a-alpha*beta)**2/M0).sum())]
        effres2 += [b]
        i *= 2.0
    fig, ax = plt.subplots(2, 2, figsize = (12,12))
    print(intres1,intres2,effres1,effres2)
    ax[0][0].set(xlabel = r"$\log_2(10i)$",ylabel = 'integral result',title = '$p(x)=f(x)$,integral result')
    ax[0][0].plot( x, intres1,color = 'b', linestyle = '-', marker = '*')
    ax[0][0].plot( x, st1,color = 'r', linestyle = '-', marker = 'None')
    ax[0][1].set(xlabel = r"$\log_2(10i)$",ylabel = 'efficiency',title = '$p(x)=f(x)$,efficiency')
    ax[0][1].plot( x, effres1, color = 'b', linestyle = '-', marker = '*')
    ax[0][1].plot( x, st2, color = 'r', linestyle = '-', marker = 'None')
    ax[1][0].set(xlabel = r"$\log_2(10i)$",ylabel = 'integral result',title = '$p(x)=(x-\\alpha\\beta)^2f(x)$,integral result')
    ax[1][0].plot( x, intres2,color = 'b', linestyle = '-', marker = '*')
    ax[1][0].plot( x, st1,color = 'r', linestyle = '-', marker = 'None')
    ax[1][1].set(xlabel = r"$\log_2(10i)$",ylabel = 'efficiency',title = '$p(x)=(x-\\alpha\\beta)^2f(x)$,efficiency')
    ax[1][1].plot( x, effres2, color = 'b', linestyle = '-', marker = '*')
    ax[1][1].plot( x, st2,color = 'r', linestyle = '-', marker = 'None')
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    