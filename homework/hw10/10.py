# -*- coding: utf-8 -*-
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


def rw_sq(f_random = sch_random,n = 10000,seed = 114, r = [0,0],a = lambda x:[0.25,0.25,0.25,0.25]):
    res = []
    temp = [r[0],r[1]]
    res.append([r[0],r[1]])
    random = f_random(ma.inf,seed = seed)
    for i in range(n):
        array = a(i)
        ra = next(random)
        for i in range(4):
            ra -= array[i]
            if ra < 0:
                break
        if i<2:
            temp[i]+=1
        else:
            temp[i-2]-=1
        res.append([temp[0],temp[1]])
    return np.array(res).T
    
def f(omega = 0.001,e = 0.1):
    return lambda x:[0.25+e*ma.cos(omega*x),0.25,0.25-e*ma.cos(omega*x),0.25]
if __name__ == '__main__':
    fig, ax = plt.subplots(1, 1, figsize = (9,9))
    r = [0,0]
    path = rw_sq(sch_random,n = 10000,seed = 114, r = r,a = f(0.01,0.1))
    x1,y1 = path[0][::20],path[1][::20]
    for i in range(1,len(x1)):
        x1[-i] -= x1[-i-1]
        y1[-i] -= y1[-i-1]
    c = []
    d = []
    d += [0]
    c += [((x1*x1).sum()+(y1*y1).sum())/(len(x1))]
    for i in range(int(len(x1)/3)):
        d += [i+1]
        c += [((x1[i+1:]*x1[:-i-1]).sum()+(y1[i+1:]*y1[:-i-1]).sum())/(len(x1)-i-1)]
    ax.plot(d, c, color = 'b', linestyle = '-', marker = 'None')
    
    
    
    
    
    
    
    