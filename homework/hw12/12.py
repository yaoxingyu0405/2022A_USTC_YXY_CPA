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
        
        
#————————————ver:bool矩阵，宽b高（b+1），垂直管道的打开与否，（i，j）是节点（i，j）下方的管道，
#                最外围的保护管道是为了方便编程而引入
#————————————par:bool矩阵，宽（b+1）高b，水平管道的打开与否，（i，j）是节点（i，j）左方的管道，
#                保护管道同上
#————————————site:b*b bool矩阵，节点与下端联通与否
#————————————上述各量均通过ndarray进行传址，方便在函数中进行修改
        
        
def map_roast(ver,par,b):#第一次传入网格进行烘焙，输入平行键与竖直键的状态，判断是否形成上下通路，并将通路进行标记
    site = np.zeros((b,b),dtype = "bool")
    flag = 0
    pio = []
    for i in range(b):
        if site[0][i]:
            continue
        site[0][i] = True
        pio.append([0,i])
        while len(pio):
            for j in pio:
                if j[0] == b-1:
                    flag = 1
                    site[j[0]][j[1]] = True
                    return site, flag
                if ver[j[0]][j[1]]:
                    if not site[j[0]-1][j[1]]:
                        pio.append([j[0]-1,j[1]])
                        site[j[0]-1][j[1]] = True
                if ver[j[0]+1][j[1]]:
                    if not site[j[0]+1][j[1]]:
                        pio.append([j[0]+1,j[1]])
                        site[j[0]+1][j[1]] = True
                if par[j[0]][j[1]]:
                    if not site[j[0]][j[1]-1]:
                        pio.append([j[0],j[1]-1])
                        site[j[0]][j[1]-1] = True
                if par[j[0]][j[1]+1]:
                    if not site[j[0]][j[1]+1]:
                        pio.append([j[0],j[1]+1])
                        site[j[0]][j[1]+1] = True
                pio.remove(j)
            if flag == 1:
                break
    return site, flag


def map_update(site,ver,par,new,b):#后续传入网格进行更新，输入平行键与竖直键的状态，判断是否形成上下通路，并将通路进行标记
    flag1 = False
    flag2 = False
    pio = []
    temp = new['coor'] 
    if new['type'] == 'ver':
         if site[temp[0]][temp[1]]^site[temp[0]-1][temp[1]]:
             if site[temp[0]][temp[1]]:
                 site[temp[0]-1][temp[1]] = True
                 pio.append([temp[0]-1,temp[1]])
             else:
                 site[temp[0]][temp[1]] = True
                 pio.append([temp[0],temp[1]])
         flag2 = True
    else:
        if site[temp[0]][temp[1]]^site[temp[0]][temp[1]-1]:
            if site[temp[0]][temp[1]]:
                site[temp[0]][temp[1]-1] = True
                pio.append([temp[0],temp[1]-1])
            else:
                site[temp[0]][temp[1]] = True
                pio.append([temp[0],temp[1]])
        flag2 = True
    if flag2:
        while len(pio):
            for j in pio:
                if j[0] == b-1:
                    flag1 = True
                    site[j[0]][j[1]] = True
                    return flag1
                if ver[j[0]][j[1]]:
                    if not site[j[0]-1][j[1]]:
                        pio.append([j[0]-1,j[1]])
                        site[j[0]-1][j[1]] = True
                if ver[j[0]+1][j[1]]:
                    if not site[j[0]+1][j[1]]:
                        pio.append([j[0]+1,j[1]])
                        site[j[0]+1][j[1]] = True
                if par[j[0]][j[1]]:
                    if not site[j[0]][j[1]-1]:
                        pio.append([j[0],j[1]-1])
                        site[j[0]][j[1]-1] = True
                if par[j[0]][j[1]+1]:
                    if not site[j[0]][j[1]+1]:
                        pio.append([j[0],j[1]+1])
                        site[j[0]][j[1]+1] = True
                pio.remove(j)
            if flag1:
                break
    return flag1
                 

def MC_normalize(random, b = 8, N = 1000,p = 0.5):#蒙特卡洛重整化系数计算
#输入重整化单元长度，计算次数，估计填充概率，返回从0.95p*(2*b*(b-1))-40开始计数的蒙特卡洛重整化
#返回概率的数组
    flag1 = False
    M = b*(b-1)
    const = int(1.9*p*M-40)
    new = {'coor':[0,0]}
    result = np.zeros((2*M-const+1),dtype='double')
    for i in range(N):
        par = np.zeros((b,b+1),dtype = "bool")
        ver = np.zeros((b+1,b),dtype = "bool")
        j = 0
        while j < (const):
            temp = int(2*M * next(random))
            if temp < b*(b-1):
                if not ver[temp//b+1][temp%b]:
                    ver[temp//b+1][temp%b] = True
                    j += 1
            else:
                if not par[temp%b][temp//b-b+2]:
                    par[temp%b][temp//b-b+2] = True
                    j+= 1
        site,flag = map_roast(ver,par,b)
        #首先，产生初始地图
        while not flag:
            temp = int(M * next(random))
            if temp < b*(b-1):
                if not ver[temp//b+1][temp%b]:
                    ver[temp//b+1][temp%b] = True
                    new['coor'] = [temp//b+1,temp%b]
                    new['type'] = 'ver'
                    flag1 = True
                    j += 1
            else:
                if not par[temp%b][temp//b-b+2]:
                    par[temp%b][temp//b-b+2] = True
                    new['coor'] = [temp%b,temp//b-b+2]
                    new['type'] = 'par'
                    flag1 = True
                    j += 1
            if flag1:
                    flag = map_update(site,ver,par,new,b)    
                    flag1 = False
        #然后对地图进行修改
        result[j-const] += 1
    for i in range(2*M-const):
        result[i+1] += result[i]
    result *= 1 / N
    return result


def p_star(array,n0,N,p):
    temp1 = 2*N*p*(1-p)
    sum = 0
    for i in range(N-n0+1):
        sum += array[i]*ma.exp(-(i+n0-p*N)**2/temp1)
    sum /= ma.sqrt(ma.pi*temp1)
    return sum - p
    
def div_p_star(array,n0,N,p):
    temp1 = 2*N*p*(1-p)
    sum = 0
    for i in range(N-n0+1):
        sum += array[i]*ma.exp(-(i+n0-p*N)**2/temp1)*((i+n0)/p-(N-i-n0)/(1-p))
    sum /= ma.sqrt(ma.pi*temp1)
    return sum


def solve_normalize(f_random = sch_random,b0 = 8, N = 100,p = 0.5,seed = 114):
    random = f_random(ma.inf,seed = seed)
    i = b0
    array = MC_normalize(random,b0,N,p)
    M = i*(i-1)
    n0 = int(1.9*p*M-40)
    f = p_star(array, n0, 2*M, p)
    while f<0:
        p += 0.001
        f = p_star(array, n0, 2*M, p)
        #定步长搜索
    return p,div_p_star(array,n0,2*M,p)
    #返回p与p的导数
        
      
if __name__ == "__main__":
    # 求解分型相关参数
    p = 0.3
    N = 5
    n = 3
    res1 = np.zeros((N))
    res2 = np.zeros((N))
    res3 = np.zeros((N))
    for i in range(n, n+N):
        res1[i-n] = i
        res2[i-n],res3[i-n] = solve_normalize(b0 = 2**i,p = p,seed = 114514)
    print(res2,np.log2(res3)/res1)
    st1 = [0.5 for y in range(N)]
    st2 = [1.3333 for y in range(N)]
    fig, ax = plt.subplots(1, 2, figsize = (12,6))
    ax[0].set(xlabel = r"$\log_2(b)$",ylabel = 'p',title = 'p-$\log_2(b)$')
    ax[0].plot( res1, res2,color = 'b', linestyle = '-', marker = '*')
    ax[0].plot( res1, st1,color = 'r', linestyle = '-', marker = 'None')
    ax[1].set(xlabel = r"$\log_2(b)$",ylabel = r'$\frac{\log_2(b)}{\log_2(dp^*/dp)}$',title = r'$\frac{\log_2(b)}{\log_2(dp^*/dp)}-\log_2(b)$')
    ax[1].plot( res1, res1/np.log2(res3),color = 'b', linestyle = '-', marker = '*')
    ax[1].plot( res1, st2,color = 'r', linestyle = '-', marker = 'None')