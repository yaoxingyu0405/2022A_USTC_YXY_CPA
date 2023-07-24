# -*- coding: utf-8 -*-
import numpy as np
def sch_random(N, M = 1, a = 16807, b = 0, m = 2**31 - 1, seed = 1):#N为生成个数，M为生成间隔
    q, r = m // a, m % a    #得到p，r
    for i in range(N):
        for j in range(M):
            seed = a * (seed % q) -r * (seed // q) #进行schrage方法
            if seed < 0:
                seed += m
        yield seed/m#单位化schrage


def MC_int(f_array_int,f_random,N,range_array,seed=14514):#这里的mc方法只适用于立方体型积分区域
    seed0=list(f_random(len(range_array),M=191,seed=seed))#生成种子序列
    a = []
    for i in range(len(range_array)):
        a += [range_array[i][0]+(range_array[i][1]-range_array[i][0])*np.array(list(f_random(N,1,seed = seed0[i])))]#生成随机点坐标
    return f_array_int(a).sum()/N
    

if __name__ =='__main__':
    range_array1=[[0.0,0.7],[0.0,4/7],[0.0,0.9],[0.0,2.0],[0.0,13/11]]
    range_array2=[[0.0,5.0]]
    print(MC_int(lambda a:5+a[0]**2-a[1]**2+3*a[0]*a[1]-a[2]**2+a[3]**3-a[4]**3, sch_random, 100000, range_array1)-1)
    print(MC_int(lambda a:np.sqrt((a[0]**2+2*np.sqrt(a[0]))), sch_random, 100000, range_array2)*5)