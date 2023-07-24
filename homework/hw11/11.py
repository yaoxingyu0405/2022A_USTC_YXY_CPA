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
        
def DLA(f_random = sch_random,seed = 114,N = 1000):#生产有N个粒子的二维网格DLA算法，返回各点的坐标
    res = [[0,0]]
    r = 0
    pos = [0,0]
    bundary = [[0,1], [1,0], [0,-1], [-1,0]]
    random = f_random(ma.inf,seed = seed)
    for i in range(1,N):
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


def lighting(f_random,seed = 514, N = 1000, gamma = 1, M = 1000):#介电击穿模型
    res = [[0,0]]
    r = 0
    pos = [0,0]
    bundary_in = [[0,0]]
    bundary_out = [[0,-1],[0,1],[-1,0],[1,0]]
    random = f_random(ma.inf,seed = seed)
    for i in range(1,N):
        
        array = laplace(bundary_in, bundary_out, r, gamma, random, M)
        pos = bundary_out[search(next(random), array)]
        res += [[pos[0],pos[1]]]#pos在一块固定内存，故需要直接添加值
        bundary_out.remove(pos)#移除外边界中增加的点
        for k in [[pos[0]+1,pos[1]],[pos[0],pos[1]+1]
                 ,[pos[0]-1,pos[1]],[pos[0],pos[1]-1]]:
            if not (k in bundary_in):
                if not (k in bundary_out):
                    bundary_out.append(k)#更新外边界，之前不在外边界且未被添加的点将会被添加
                    
            if k in bundary_in:
                if not (([k[0]+1,k[1]] in bundary_out) or ([k[0]-1,k[1]] in bundary_out)
                    or ([k[0],k[1]+1] in bundary_out) or ([k[0],k[1]-1] in bundary_out)):
                    bundary_in.remove(k)#更新内边界，之前在内边界且周围一个外点也没有的点将会被移除
                    
        if (([pos[0]+1,pos[1]] in bundary_out) or ([pos[0]-1,pos[1]] in bundary_out)
            or ([pos[0],pos[1]+1] in bundary_out) or ([pos[0],pos[1]-1] in bundary_out)):
            bundary_in.append([pos[0],pos[1]])#更新内边界，新添加点若至少有一个外点即会被添加
        if(pos[0] ** 2 + pos[1] ** 2 > r ** 2):
            r = ma.sqrt(pos[0] ** 2 + pos[1] ** 2) #更新半径，确定行走外界
    return np.array(res).T,np.array(bundary_out).T
    

def laplace(bundary_in,bundary_out,r,gamma,random,M):
    k = 0
    res = np.zeros(len(bundary_out), dtype="double")
    for pos in bundary_out:
        temp = 0
        for i in range(M): 
            pos_tem = [pos[0], pos[1]]
            l = 0
            while True:
                l += 1
                j = int(4 * next(random))
                if j < 2:
                    pos_tem[j] += 1
                else:
                    pos_tem[j-2] -= 1
                if pos_tem in bundary_in:
                    temp += 1.0
                    break
                if(pos_tem[0] ** 2 + pos_tem[1] ** 2 > 2*(r + 10) ** 2): 
                    break
        res[k] = temp / M
        k += 1
    res = np.absolute(1-res) ** gamma
    result = np.zeros((len(bundary_out)+1),dtype="double")
    for i in range(len(bundary_out)):
        result[i + 1] = result[i] + res[i]
    return result/result[-1]
                
    
def search(a,array):
    temp1, temp2 = 0,len(array)
    while temp2-temp1>1:
        if array[int((temp1+temp2)/2)] > a:
            temp2 = int((temp1+temp2)/2)
        else:
            temp1 = int((temp1+temp2)/2)
    return temp1


if __name__ == "__main__":
    fig, ax = plt.subplots(1, 2, figsize = (18,9))
    res, bundary = DLA(sch_random,seed = 114,N = 1000)
    x1, y1 = res[0],res[1]
    x2, y2 = bundary[0],bundary[1]
    ax[0].plot(x2, y2, color = 'k', linestyle = 'None', marker = 'o')
    ax[0].plot(x1, y1, color = 'grey', linestyle = 'None', marker = 'o')
    res, bundary = lighting(sch_random,seed = 514,N = 300,M = 50)
    x1, y1 = res[0],res[1]
    x2, y2 = bundary[0],bundary[1]
    ax[1].plot(x2, y2, color = 'k', linestyle = 'None', marker = 'o')
    ax[1].plot(x1, y1, color = 'grey', linestyle = 'None', marker = 'o')
    
    
    
    