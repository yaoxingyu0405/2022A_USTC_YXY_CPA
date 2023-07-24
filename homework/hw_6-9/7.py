# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
def data_read(f_name):
    try:
        f = open(f_name,'r')
    except FileNotFoundError:
        print('无法打开指定的文件!')
        return 0
    a = []
    b = []
    next(f)
    for line in f:
        tem = line.find('\t')
        a += [int(line[:tem])]
        b += [int(line[tem+1:-1])]
    return np.array(a), np.array(b)


def sch_random(N, M = 1, a = 16807, b = 0, m = 2**31 - 1, seed = 1):#N为生成个数，M为生成间隔
    q, r = m // a, m % a    #得到p，r
    for i in range(N):
        for j in range(M):
            seed = a * (seed % q) -r * (seed // q) #进行schrage方法
            if seed < 0:
                seed += m
        yield seed/m#单位化schrage


def direct_sam(f_name,f_random,N,seed = 114):#离散直接抽样
    res = []
    a,b = data_read(f_name)
    lenl = len(a)
    for i in range(1,lenl):
        b[i] += b[i-1]
    b =b / b[-1]
    for i in f_random(N,seed = seed):
        temp1 ,temp2= 0,lenl
        while temp2 - temp1!=1 :
            temp3 = int((temp1+temp2)/2)
            if b[temp3]>i:
                temp2 = temp3
            else:
                temp1 = temp3
        #二分查找
        res += [a[temp2]]
    return np.array(res)
    
    
def sel_sam(f_name,f_random,N,k=5,seed1 = 114,seed2 = 514):#k为分段数，离散舍选抽样
    res = []
    part = []
    bmax = []
    a,b = data_read(f_name)
    lenl = len(a)
    
    for i in range(k):
        part += [b[int(i*lenl/k):int((i+1)*lenl/k)]]
        bmax += [part[i].max()*(int((i+1)*lenl/k)-int(i*lenl/k))]
    bmax = np.array(bmax)
    for i in range(1,k):
        bmax[i] += bmax[i-1]
    bmax = bmax / bmax[-1]
    brandom = f_random(N,seed = seed1)
    for i in f_random(N,seed = seed2):
        temp1 ,temp2= 0,k
        while temp2 - temp1>1 :
            temp3 = (temp1+temp2)//2
            if bmax[temp3]>i:
                temp2 = temp3
            else:
                temp1 = temp3
        i = int((i-bmax[temp1])/(bmax[temp2]-bmax[temp1])*(int((temp2+1)*lenl/k)-int(temp2*lenl/k)))+int((temp2)*lenl/k)
        b_random = next(brandom)*part[temp2].max()
        if b_random<b[i]:
            res+=[a[i]]
    return np.array(res)
    

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

if __name__ =='__main__':
    fig, ax = plt.subplots(1, 2, figsize = (8,2)) 
    a,b = data_read("data.TXT")
    b = b/b.max()
    
    x,y,eff = freq(direct_sam('data.TXT',sch_random,N=100000))
    y = y/y.max()
    ax[0].plot(x, y, color = 'r', linestyle = '-', marker = 'None')
    ax[0].plot(a, b, color = 'b', linestyle = '--', marker = 'None')

    x,y,eff = freq(sel_sam('data.TXT',sch_random,N=100000))
    y = y/y.max()
    ax[1].plot(x, y, color = 'r', linestyle = '-', marker = 'None')
    ax[1].plot(a, b, color = 'b', linestyle = '--', marker = 'None')
    print("抽样效率："+str(eff/100000))