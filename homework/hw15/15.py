# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 15:35:13 2022

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


def MHM_2d(fp, N, delta, beta, seed = 114514, f_random = sch_random):
    i = 0#循环变量
    eff = 0#计算抽样效率
    res = np.zeros((N,2),dtype = 'double')
    pos = np.zeros((2),dtype = 'double')
    temp1 = np.zeros((2),dtype = 'double')
    random = f_random(ma.inf, seed = seed)
    while i < N/100:#热化
        temp1[0] = delta*(next(random)-0.5)
        temp1[1] = delta*(next(random)-0.5)
        a = fp(pos) - fp(pos+temp1)
        if a > 0:
            pos += temp1
        else:
            if ma.exp(a) > next(random):
                pos += temp1
            else:
                continue
        i += 1
    i = 0
    while i < N:#正式模拟
        eff += 2
        temp1[0] = delta*(next(random)-0.5)
        temp1[1] = delta*(next(random)-0.5)
        a = fp(pos) - fp(pos+temp1)
        if a > 0:
            pos += temp1
            res[i][0] = pos[0]
            res[i][1] = pos[1]
        else:
            eff += 1
            if ma.exp(beta*a) > next(random):
                pos += temp1
                res[i][0] = pos[0]
                res[i][1] = pos[1]
            else:
                continue
        i += 1
    return res, N / eff


def fp(pos):
    return -2.0*(pos[0]**2+pos[1]**2)+0.5*(pos[0]**4+pos[1]**4)+0.5*(pos[0]-pos[1])**4

if __name__ == "__main__":
    N = 10000
    fig, ax = plt.subplots(4, 3, figsize = (24,32))
    res_r = np.zeros((3,3))
    res_x = np.zeros((3,3))
    reff = np.zeros((3,3))
    a = [0.2,1,5]
    for i in range(3):
        for j in range(3):
            res, eff = MHM_2d(fp, N, a[i], a[j], seed = 11451419)
            res = res.T
            ax[i][j].set(xlabel = 'x',ylabel = 'y',title = r'$\delta = {},\beta = {}$'.format(a[i], a[j]))
            ax[i][j].plot(res[0], res[1], color = 'b', linestyle = 'None', marker = '.')
            xsq = (res[0]**2).sum()/N
            ysq = (res[1]**2).sum()/N
            r = xsq+ysq
            res_r[j][i] = r 
            res_x[j][i] = xsq
            reff[j][i] = eff
    x0 = [-1,0,1]
    color = ['black','yellow','red']
    ax[3][0].set(xlabel = r'$\log_5(\delta)$',ylabel = r'$<r^2>$',title = r'$<r^2>-\delta$')
    ax[3][1].set(xlabel = r'$\log_5(\delta)$',ylabel = 'eff',title = r'$eff-\delta$')
    ax[3][2].set(xlabel = r'$\log_5(\delta)$',ylabel = r'$<x^2>$',title = r'$<x^2>-\delta$')
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    for i in range(3):
        temp1 += ax[3][0].plot(x0,res_r[i],color = color[i],linestyle = '-',marker = '*')
        temp2 += ax[3][1].plot(x0,reff[i],color = color[i],linestyle = '-',marker = '*')
        temp3 += ax[3][2].plot(x0,res_x[i],color = color[i],linestyle = '-',marker = '*')
        temp4 += [r"$\beta = {}$".format(a[i])]
    print(temp4)
    ax[3][0].legend(handles=temp1, labels=temp4)
    ax[3][1].legend(handles=temp2, labels=temp4)
    ax[3][2].legend(handles=temp3, labels=temp4)
        
        
    
    
    
    
    