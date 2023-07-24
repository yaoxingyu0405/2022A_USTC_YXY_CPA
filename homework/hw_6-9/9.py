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


def dis_array(f_dis,a,b):#由离散抽样函数生成抽样列
    i = 0
    res = [0]
    res = np.zeros((b-a+2),dtype = 'double')
    for i in range(b-a+1):
        res[i+1] = f_dis(i+a)
        res[i+1] += res[i]
    return res/res[i+1]

def con_sam_array(f_con,a,b,k = 10):#由连续函数生成舍选抽样区间列
    temp = 0
    bmax = np.zeros((k+1),dtype = 'double')
    b0 = np.zeros((k),dtype = 'double')
    const = (b-a)/k/100
    for i in range(k):
        temp = 0
        for j in range(100):
            if temp < f_con((i*100+j+0.5)*const+a):
                temp = f_con((i*100+j+0.5)*const+a)
        b0[i] = temp
        bmax[i+1] = temp
        bmax[i+1] += bmax[i]
    bmax = bmax/bmax[-1]
    return bmax,b0

def search(a,b):#二分查找
    temp1,temp2 = 0,len(b)
    while temp2 - temp1>1 :
        temp = (temp1+temp2)//2
        if b[temp] > a:
            temp2 = temp
        else:
            temp1 = temp
    return temp1
    
def sam(sam):
    a = sam['inf']
    b = sam['sup']
    fun = sam['function']
    N = sam['N']
    n = sam['n']
    if sam['ftype'] == 'con':
        k = sam['k']
        result = []
        bmax,b0 = con_sam_array(fun,a,b,k)
        random1 = sam['random'](N*n,seed = sam['seed1'])
        random2 = sam['random'](N*n,seed = sam['seed2'])
        for l in range(n):
            res = 0
            j = 0
            for mm in range(N):
                i = next(random1)
                num = search(i,bmax)
                b_rand = next(random2)*b0[num]
                i = ((i-bmax[num])/(bmax[num+1]-bmax[num])+num)*(b-a)/k+a
                if fun(i)>b_rand:
                    res += i
                    j += 1
            if j<100:
                continue
            result += [(res/j-sam['avg'])*j**0.5/sam['sigma']]
        return result
    if sam['ftype'] == 'dis':
        result = []
        bmax = dis_array(fun,a,b)
        random = sam['random'](N*n,seed = sam['seed1'])
        for l in range(n):
            res = 0
            for mm in range(N):
                i = next(random)
                num = search(i,bmax)
                res += num
            result += [(res/N-sam['avg'])*N**0.5/sam['sigma']]
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
    sam1 = {'ftype':'dis','sup':20,'inf':0,'function':lambda x:1/m.e**2/m.factorial(x)*2*x,
            'N':100,'n':10000,'random':sch_random,'avg':2,'sigma':2,'seed1':1145}#l=1的泊松分布
    # sam1 = {'ftype':'dis','sup':20,'inf':0,'function':lambda x:1/m.e**2/m.factorial(x)*2*x,
    #         'N':100,'n':10000,'random':sch_random,'avg':2,'sigma':2,'seed1':114514}#l=2的泊松分布
    fig, ax = plt.subplots(1, 1, figsize = (4,2))
    x,y,eff = freq(sam(sam1))
    y = y/y.max()
    ax.plot(x, y, color = 'r', linestyle = '-', marker = 'None')
            
        
        
        
        
        
    

















    
