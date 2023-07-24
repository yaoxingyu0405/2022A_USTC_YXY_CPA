# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 10:59:33 2022

@author: YXY
"""

import numpy as np
import matplotlib.pyplot as plt
import math as ma


def fun1():
    N1 = 1000
    N2 = 100
    L = 1.0
    fig, ax = plt.subplots(1, 1, figsize = (12,6))
    ax.set(xlabel = r"$\lambda$",ylabel = '$x_n$',title = '$x_n-\lambda$')
    for i in range(N1):
        y = [L/N1*(i+1+x/N2) for x in range(N2)]
        x = np.zeros((N2+100))
        x[0] = 0.3
        for j in range(N2+99):
            x[j+1] = ((L/N1*(i+1)))*ma.sin(ma.pi*x[j])
        ax.plot(y,x[100:],color = 'b',linestyle = 'None',marker = '.')
    #绘制图像，展示定值状态，倍周期分叉与混沌
    
def fun2():
    N1 = 10000
    N2 = 2000
    L = 1.0
    res0 = [[0,1]]
    for i in range(int(0.717*N1),N1):
        x = np.zeros((N2+100))
        x[0] = 0.3
        for j in range(N2+99):
            x[j+1] = ((L/N1*(i+1)))*ma.sin(ma.pi*x[j])
        res = np.abs(np.fft.fft(x[100:]))
        temp = 1
        for j in range(1,N2):
            if res[j] > 1:
                temp = N2/j
                break
        if res0[-1][1] == int(temp):
            res0[-1][0] = i
        else:
            res0 += [[i,int(temp)]]
        if len(res0)>5:
                break
    alpha = 0
    for i in range(5):
        print('n='+str(res0[i][1])+',lambda='+str(res0[i][0]/N1))
    for i in range(3):
        alpha += (res0[i][0]-res0[i+1][0]) / (res0[i+1][0]-res0[i+2][0])
        print('n='+str(res0[i][1])+",Feigenbaum = "+str((res0[i][0]-res0[i+1][0]) / (res0[i+1][0]-res0[i+2][0])))
    print("Feigenbaum = "+str(alpha/3))
    #找到分叉点，并计算Feigenbaum常数

if __name__ == "__main__":
    fun1()
    fun2()
    
