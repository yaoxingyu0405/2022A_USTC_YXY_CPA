# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import math as m

def sch_random(N, M = 1, a = 16807, b = 0, m = 2**31 - 1, seed = 1):#N为生成个数，M为生成间隔
    q, r = m // a, m % a    #得到p，r
    for i in range(N):
        for j in range(M):
            seed = a * (seed % q) -r * (seed // q) #进行schrage方法
            if seed < 0:
                seed += m
        yield seed/m#单位化schrage
        

def fun_sel_sam(function,f_random,a,b,k=5,N=100000,seed1 = 114514,seed2 = 1919810):#抽样函数
    result = []
    temp = 0
    bmax = np.zeros((k+1),dtype = 'double')
    b0 = np.zeros((k),dtype = 'double')
    const = (b-a)/k/100
    for i in range(k):
        temp = 0
        for j in range(100):
            if temp < function((i*100+j+0.5)*const+a):
                temp = function((i*100+j+0.5)*const+a)
        b0[i] = temp
        bmax[i+1] = temp
        bmax[i+1] += bmax[i]
    bmax = bmax/bmax[-1]#产生阶梯分布函数
    
    
    brandom = f_random(N,seed = seed1)
    for i in f_random(N,seed = seed2):
        temp1,temp2 = 0,k
        while temp2 - temp1>1 :
            temp = (temp1+temp2)//2
            if bmax[temp]>i:
                temp2 = temp
            else:
                temp1 = temp
        #二分查找得到位置
        b_rand = next(brandom)*b0[temp1]
        i = ((i-bmax[temp1])/(bmax[temp2]-bmax[temp1])+temp1)*(b-a)/k+a
        if function(i)>b_rand:#判断是否保留
            result+=[i]
    return np.array(result)


def freq(c):#频数统计
    d = {}
    for i in c:
        if round(i,2) in d:
            d[round(i,2)] += 1
        else:
            d[round(i,2)] = 1
    x , y = np.zeros((len(d))),np.zeros((len(d)))
    j = 0
    for i in sorted(d.items(),key = lambda x:x[0]):
        x[j] = i[0]
        y[j] = i[1]
        j += 1
    return x,y,len(c)#返回长度计算抽样效率
    
    
if __name__ == '__main__':
    fig, ax = plt.subplots(1, 2, figsize = (8,2))
    
    x,y,eff= freq(fun_sel_sam(lambda x:m.exp(-2*x**2), sch_random, a=-5, b=5,k=10,N=100000))
    z = np.exp(-2*x**2)#产生标准曲线
    y = y/y.max()*1.1
    ax[0].plot(x, y, color = 'r', linestyle = '-', marker = 'None')
    ax[0].plot(x, z, color = 'b', linestyle = '-', marker = 'None')
    print('抽样效率:'+str(eff/100000))
    
    
    x,y,eff = freq(fun_sel_sam(lambda x:1/(2*x**4+1), sch_random, a=-5, b=5,k=10,N=100000))
    z = 1/(2*x**4+1)#产生标准曲线
    y = y/y.max()*1.1
    ax[1].plot(x, y, color = 'r', linestyle = '-', marker = 'None')
    ax[1].plot(x, z, color = 'b', linestyle = '-', marker = 'None')
    print('抽样效率:'+str(eff/100000))
    
    
    