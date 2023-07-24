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
        
def DLA(f_random = sch_random,seed = 114,N = 10000):#生产有N个粒子的二维网格DLA算法，返回各点的坐标
    res = [[0,0]]
    r = 0
    pos = [0,0]
    bundary = [[0,1], [1,0], [0,-1], [-1,0]]
    random = f_random(ma.inf,seed = seed)
    for i in range(N):
        while True:
            flag = 0
            theta = 2 * ma.pi * next(random)
            pos[0] = int((r + 5.0) * ma.cos(theta))
            pos[1] = int((r + 5.0) * ma.sin(theta))
            while True:
                j = int(4 * next(random))
                if j < 2:
                    pos[j] += 1
                else:
                    pos[j-2] -= 1
                if pos in bundary:
                    res += [[pos[0],pos[1]]]
                    bundary.remove(pos)
                    for k in [[pos[0]+1,pos[1]],[pos[0],pos[1]+1]
                             ,[pos[0]-1,pos[1]],[pos[0],pos[1]-1]]:
                        if not (k in res):
                            if not (k in bundary):
                                bundary.append(k)
                    if(pos[0] ** 2 + pos[1] ** 2 > r ** 2):
                        r = ma.sqrt(pos[0] ** 2 + pos[1] ** 2)
                    flag = 1
                    break
                if(pos[0] ** 2 + pos[1] ** 2 > 1.1 * (r + 5) ** 2): 
                    break
            if flag == 1:
                break
    return np.array(res).T,np.array(bundary).T


def search(a,array):#一般路过二分查找
    temp1, temp2 = 0,len(array)
    while temp2-temp1>1:
        if array[int((temp1+temp2)/2)] > a:
            temp2 = int((temp1+temp2)/2)
        else:
            temp1 = int((temp1+temp2)/2)
    return temp1


def box_counting(res,bundary):#盒计数法计算分形维度
    xmax = 0
    ymax = 0
    xmin = 0
    ymin = 0
    rres = []
    width = 0
    length = 0
    tempw = 0
    templ = 0
    for pos in bundary[0]:
        if pos > xmax:
            xmax = pos
        if pos < xmin:
            xmin = pos
    for pos in bundary[1]:
        if pos > ymax:
            ymax = pos
        if pos < ymin:
            ymin = pos
    n = np.zeros((ymax-ymin-1,xmax-xmin-1),dtype = 'bool')
    for i in range(len(res[0])):
        n[res[1][i]-ymin-1][res[0][i]-xmin-1] = True
    length = ymax-ymin-1
    width = xmax-xmin-1
    rres += [len(res[0])]
    while True:
        temp = 0
        tempw = width%2
        templ = length%2
        width = width//2
        length = length//2
        if not width*length:
            break
        tempn = np.zeros((length+templ,width+tempw),dtype = 'bool')
        for i in range(length):
            for j in range(width):
                if (n[2*i][2*j] or n[2*i+1][2*j] or n[2*i+1][2*j+1] or n[2*i][2*j+1]):
                    tempn[i][j] = True
                    temp+=1
            if tempw:
                if(n[2*i][2*width] or n[2*i+1][2*width]):
                    tempn[i][width] = True
                    temp+=1
        if templ:
            for j in range(width):
                if (n[2*length][2*j] or n[2*length][2*j+1]):
                    tempn[i][j] = True
                    temp+=1
            if tempw:
                if n[2*length][2*width]:
                    tempn[length][width] = True
                    temp+=1
        n = tempn
        length+=templ 
        width+=tempw
        rres+=[temp]
    return np.array(rres)


def sandbox(res, bundary):
    xmax = 0
    ymax = 0
    xmin = 0
    ymin = 0
    count = 1
    result = []
    for pos in bundary[0]:
        if pos > xmax:
            xmax = pos
        if pos < xmin:
            xmin = pos
    for pos in bundary[1]:
        if pos > ymax:
            ymax = pos
        if pos < ymin:
            ymin = pos
    a = min([ymax,-ymin,xmax,-xmin])-1#最小边长，确定计数范围
    n = np.zeros((2*a+1,2*a+1),dtype = 'bool')
    for i in range(len(res[0])):
        try:
            n[res[1][i]+a][res[0][i]+a] = True
        except IndexError:
            continue
    for i in range(1,a+1):
        for j in range(a-i,a+i+1):
            if n[a+i][j]:
                count += 1
            if n[a-i][j]:
                count += 1
        for j in range(a-i+1,a+i):
            if n[j][a+i]:
                count += 1
            if n[j][a-i]:
                count += 1
        result += [count]
    return result
        
    
if __name__ == "__main__":
    seed = [114,514,1919,810]
    color = ['r','b','k','yellow']
    fig, ax = plt.subplots(1, 2, figsize = (18,9))
    ax[0].set(xlabel = r'$\log_2(b)$',ylabel = r'$\log_2(N)$',title = r'$\log_2(N)-\log_2(b)$,sandbox')
    ax[1].set(xlabel = r'$\log_2(b)$',ylabel = r'$\log_2(N)$',title = r'$\log_2(N)-\log_2(b)$,box-counting')
    temp1 = []
    temp2 = []
    temp3 = []
    temp4 = []
    for i in range(4):
        a, b = DLA(sch_random,seed = seed[i],N = 1000)
        
        res = sandbox(a, b)
        x1 = np.array([ma.log2(2*x+1) for x in  range(1,len(res)+1)])
        y = np.log2(res)
        temp1 += ax[0].plot(x1, y, color = color[i], linestyle = '-', marker = '*')
        avgx = (x1).sum()/len(x1)
        avgx2 = (x1**2).sum()/len(x1)
        avgy2 = (y**2).sum()/len(x1)
        avgy = (y).sum()/len(x1)
        avgxy = (x1*y).sum()/len(x1)
        r = (avgxy-avgx*avgy)/ma.sqrt(avgx2-avgx**2)/ma.sqrt(avgy2-avgy**2)
        aa = (avgxy-avgx*avgy)/(avgx2-avgx**2)
        ba = (avgx2*avgy-avgxy*avgx)/(avgx2-avgx**2)
        print('sandbox:a = '+str(aa)+'b = '+str(ba)+'r = '+str(r))
        temp3 += [r"$\nu = {}$".format(aa)]
        
        res = box_counting(a, b)
        x1 = np.array([x for x in  range(0,len(res))])
        y = np.log2(res)
        temp2 += ax[1].plot(x1, y, color = color[i], linestyle = '-', marker = '*')
        avgx = (x1).sum()/len(x1)
        avgx2 = (x1**2).sum()/len(x1)
        avgy2 = (y**2).sum()/len(x1)
        avgy = (y).sum()/len(x1)
        avgxy = (x1*y).sum()/len(x1)
        r = (avgxy-avgx*avgy)/ma.sqrt(avgx2-avgx**2)/ma.sqrt(avgy2-avgy**2)
        aa = (avgxy-avgx*avgy)/(avgx2-avgx**2)
        ba = (avgx2*avgy-avgxy*avgx)/(avgx2-avgx**2)
        print('box-counting:a = '+str(aa)+'b = '+str(ba)+'r = '+str(r))
        temp4 += [r"$\nu = {}$".format(-aa)]
    ax[0].legend(handles=temp1, labels=temp3)
    ax[1].legend(handles=temp2, labels=temp4)
        
        
    
    
    