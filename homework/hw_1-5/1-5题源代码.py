# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


def sch_random(N, M = 1, a = 16807, b = 0, m = 2**31 - 1, seed = 1):#N为生成个数，M为生成间隔
    q, r = m // a, m % a    #得到p，r
    for i in range(N):
        for j in range(M):
            seed = a * (seed % q) -r * (seed // q) #进行schrage方法
            if seed < 0:
                seed += m
        yield seed  #产生可迭代对象，简化代码


def distribution(f_random, seed = 114, N = 100000, M = 1):#传入函数以减少耦合，下同
    l = list(f_random(2*N, M, seed = seed))#使用list将可迭代对象转换成list对象，下同
    x, y = l[::2], l[1::2] #进行切片，得到相邻的列，便于绘图
    fig, ax = plt.subplots(1, 1, figsize = (60,60))
    ax.plot(x, y, color = 'r', linestyle = 'None', marker = '.')
    fig.show()
    return True


def compare(f_random, seed = 1919, N = 100000):
    l = list(f_random(3*N, seed = seed))
    x, y, z = np.array(l[::3]), np.array(l[1::3]), np.array(l[2::3])#进行切片
    res = (x-z)/np.absolute(x-z)+(z-y)/np.absolute(z-y)#使用numpy自带的矩阵计算简化代码
    res = (res +np.absolute(res)).sum()/4#使用绝对值函数直接获得个数
    return res/N


def powered_k(f_random,seed=1919,N=100000,k=1):
    m = np.array(list(f_random(N,seed = seed)))
    m = m/m.max()
    return (m**k).sum()*(k+1)/N-1


def independence(f_random, seed =810, N=100000, l0 = 1):#C(l)检验
    l = list(f_random(N+l0, seed = seed))
    a, b = np.array(l[l0:])/N, np.array(l[:-l0])/N# 在此处先除以N是为了防止溢出
    avg1, avg2, avg3 = b.sum() , (b*b).sum()*N , (a*b).sum()*N 
    return (avg3 - avg1 ** 2) / (avg2 - avg1 ** 2)

#############抽样#############
def sph_sam(f_random,seed1 = 114514,seed2 = 1919810,N=1000):#上半球面抽样
    theta = np.arccos(np.array(list(f_random(N, seed = seed1)))/(2**31 - 1))#具体推导见文档捏
    phi = 2*3.1415926535*(np.array(list(f_random(N, seed = seed2)))/(2**31 - 1))
    x = np.cos(phi)*np.sin(theta)
    y = np.sin(phi)*np.sin(theta)
    z = np.cos(theta)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plt.rcParams['savefig.dpi']=1000
    plt.rcParams['figure.dpi']=1000
    ax.scatter(x, y, z, c='r', marker='*')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
    
    
def mar_sam(f_random,seed1=114514,seed2=1919810,N=1000):#Marsaglia抽样
    u = (np.array(list(f_random(N=int(1.27324*N), seed = seed1)))/(2**31 - 1))*2-1
    v = (np.array(list(f_random(N=int(1.27324*N), seed = seed2)))/(2**31 - 1))*2-1
    r2 = u**2+v**2
    t = (1-(r2-1)/np.absolute(r2-1))/2
    u,v,r2 = u*t,v*t,r2*t
    x,y = u*np.sqrt(1-r2),v*np.sqrt(1-r2)
    fig, ax = plt.subplots(1, 1, figsize = (10,10))
    ax.plot(x, y, color = 'r', linestyle = 'None', marker = '*')
    fig.show()
    
    
    
if __name__ == '__main__':
    distribution(sch_random)
    print(independence(sch_random, l=1))
    print(powered_k(sch_random, k=1))
    print(compare(sch_random))
    sph_sam(sch_random)
    mar_sam(sch_random)
    
    
    
    
    
    
    
    
    
    