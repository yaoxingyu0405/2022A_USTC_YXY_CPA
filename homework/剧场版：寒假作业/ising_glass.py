# -*- coding: utf-8 -*-
"""
Created on Sun Jan 22 15:37:09 2023

@author: YXY
"""

import numpy as np
import math as ma
import matplotlib.pyplot as plt

inter_num = 20

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

class Ising_glass:
    def __init__(self,n,E1,E2,seed = 114514):
        self.scale = n
        self.N = n**2
        self.E1 = E1
        self.E2 = E2
        self.xy_data = np.zeros((n,n),dtype = 'float')
        self.ising_data = np.zeros((n,n),dtype = 'short')
        self.random = sch_random(ma.inf, seed = seed)
        for i in range(n):
            for j in range(n):
                self.xy_data[i][j] = next(self.random)
                self.ising_data[i][j] = 2 * int(next(self.random) + 0.5) - 1
        # 随机初始化
    def metropolis(self):
        energy = 0
        r2 = 0
        flag = True
        temp = np.zeros(4,dtype = 'float')
        while flag:
            r = int(2*next(self.random)*self.N)
            if r < self.N:
                i = r // self.scale
                j = r % self.scale
                if i == 0:
                    temp[0] = self.xy_data[1][j]
                    temp[2] = self.xy_data[self.scale-1][j]
                elif i == self.scale-1:
                    temp[0] = self.xy_data[0][j]
                    temp[2] = self.xy_data[self.scale-2][j]
                else:
                    temp[0] = self.xy_data[i+1][j]
                    temp[2] = self.xy_data[i-1][j]
                if j == 0:
                    temp[1] = self.xy_data[i][1]
                    temp[3] = self.xy_data[i][self.scale-1]
                elif j == self.scale-1:
                    temp[1] = self.xy_data[i][0]
                    temp[3] = self.xy_data[i][self.scale-2]
                else:
                    temp[1] = self.xy_data[i][j+1]
                    temp[3] = self.xy_data[i][j-1]
                r2 = next(self.random)
                energy = (np.cos(2*ma.pi*(self.xy_data[i][j] - temp)) - np.cos(2*ma.pi*(r2- temp))).sum()
                if energy < 0:
                    self.xy_data[i][j] = r2
                    break
                else:
                    r3 = next(self.random)
                    if r3 < ma.exp(-energy*self.E1):
                        self.xy_data[i][j] = r2
                        break    
            else:
                r -= self.N
                i = r // self.scale
                j = r % self.scale
                if i == 0:
                    temp[0] = self.ising_data[1][j]
                    temp[2] = self.ising_data[self.scale-1][j]
                elif i == self.scale-1:
                    temp[0] = self.ising_data[0][j]
                    temp[2] = self.ising_data[self.scale-2][j]
                else:
                    temp[0] = self.ising_data[i+1][j]
                    temp[2] = self.ising_data[i-1][j]
                if j == 0:
                    temp[1] = self.ising_data[i][1]
                    temp[3] = self.ising_data[i][self.scale-1]
                elif j == self.scale-1:
                    temp[1] = self.ising_data[i][0]
                    temp[3] = self.ising_data[i][self.scale-2]
                else:
                    temp[1] = self.ising_data[i][j+1]
                    temp[3] = self.ising_data[i][j-1]
                energy = 2*self.ising_data[i][j]*temp.sum()
                if energy < 0:
                    self.ising_data[i][j] *= -1
                    break
                else:
                    r1 = next(self.random)
                    if r1 < ma.exp(-energy*self.E2):
                        self.ising_data[i][j] *= -1
                        break 
    def Wolff(self,n):
        culster = []
        culster_mat = np.zeros((100,100),dtype = 'bool')
        if n == 0:    
            r = int(next(self.random)*self.N)
            i = r // self.scale
            j = r % self.scale
            data = self.ising_data[i][j]
            self.ising_data[i][j] = -data
            culster += [[i,j]]
            culster_mat[i][j] = True
            temp = [0,0,0,0]
            for k in culster:
                if k[0] == 0:
                    temp[0] = [self.scale-1,k[1]]
                    temp[1] = [1,k[1]]
                elif k[0] == self.scale-1:
                    temp[0] = [0,k[1]]
                    temp[1] = [self.scale-2,k[1]]
                else:
                    temp[0] = [k[0]+1,k[1]]
                    temp[1] = [k[0]-1,k[1]]
                if k[1] == 0:
                    temp[2] = [k[0],self.scale-1]
                    temp[3] = [k[0],1]
                elif k[1] == self.scale-1:
                    temp[2] = [k[0],0]
                    temp[3] = [k[0],self.scale-2]
                else:
                    temp[2] = [k[0],k[1]-1]
                    temp[3] = [k[0],k[1]+1]
                for l in temp:
                    if culster_mat[l[0]][l[1]]:
                        continue
                    if self.ising_data[l[0],l[1]] != data:
                        continue
                    r2 = next(self.random)
                    if r2 > ma.exp(-2*self.E2):
                        culster += [l.copy()]
                        self.ising_data[l[0]][l[1]] = -data  
                        culster_mat[l[0]][l[1]] = True
        else:
            r = int(next(self.random)*self.N)
            i = r // self.scale
            j = r % self.scale
            r1 = next(self.random)
            self.xy_data[i][j] =(-self.xy_data[i][j] + 2*r1) % 1
            culster += [[i,j]]
            temp = [0,0,0,0]
            culster_mat[i][j] = True
            for k in culster:
                if k[0] == 0:
                    temp[0] = [self.scale-1,k[1]]
                    temp[1] = [1,k[1]]
                elif k[0] == self.scale-1:
                    temp[0] = [0,k[1]]
                    temp[1] = [self.scale-2,k[1]]
                else:
                    temp[0] = [k[0]+1,k[1]]
                    temp[1] = [k[0]-1,k[1]]
                if k[1] == 0:
                    temp[2] = [k[0],self.scale-1]
                    temp[3] = [k[0],1]
                elif k[1] == self.scale-1:
                    temp[2] = [k[0],0]
                    temp[3] = [k[0],self.scale-2]
                else:
                    temp[2] = [k[0],k[1]-1]
                    temp[3] = [k[0],k[1]+1]
                for l in temp:
                    if culster_mat[l[0]][l[1]]:
                        continue
                    energy = (ma.cos(2*ma.pi*(self.xy_data[k[0]][k[1]]-
                                          self.xy_data[l[0]][l[1]]))-
                            ma.cos(2*ma.pi*(-self.xy_data[k[0]][k[1]]+2*r1-
                                         self.xy_data[l[0]][l[1]])))
                    if energy > 0:
                        continue
                    else:
                        r2 = next(self.random)
                        if r2 < ma.exp(self.E1*energy):
                            continue
                        else:
                            culster += [l.copy()]
                            self.xy_data[l[0]][l[1]] = (-self.xy_data[l[0]][l[1]] + 2*r1) % 1
                            culster_mat[l[0]][l[1]] = True
        return len(culster)
    def entropy(self):
        up = np.zeros(inter_num,dtype='float')
        down = np.zeros(inter_num,dtype = 'float')
        for i in range(self.scale):
            for j in range(self.scale):
                if self.ising_data[i][j] == 1:
                    up[int(inter_num*self.xy_data[i][j])] += 1 
                else:
                    down[int(inter_num*self.xy_data[i][j])] += 1
        up /= self.N
        down /= self.N
        sum = 0.0
        for i in range(inter_num):
            if up[i] != 0:
                sum -= up[i]*ma.log(up[i])
            if down[i] != 0:
                sum -= down[i]*ma.log(down[i])
        a = np.linspace(0, 2*ma.pi,inter_num+1)[:-1]
        res = ma.sqrt((np.sin(a)*(up+down)).sum()**2+(np.sin(a)*(up+down)).sum()**2)
        return sum, up.sum(), res
    def energy(self):
        res1 = 0.0
        res2 = 0.0
        for i in range(self.scale):
            for j in range(self.scale):
                if i != self.scale-1:
                    res1 += ma.cos(2*ma.pi*(self.xy_data[i][j]-self.xy_data[i+1][j]))
                    res2 += self.xy_data[i][j]*self.xy_data[i+1][j]
                else:
                    res1 += ma.cos(2*ma.pi*(self.xy_data[self.scale-1][j]-self.xy_data[0][j]))
                    res2 += self.xy_data[self.scale-1][j]*self.xy_data[0][j]
                if j != self.scale-1:
                    res1 += ma.cos(2*ma.pi*(self.xy_data[i][j]-self.xy_data[i][j+1]))
                    res2 += self.xy_data[i][j]*self.xy_data[i][j+1]
                else:
                    res1 += ma.cos(2*ma.pi*(self.xy_data[i][self.scale-1]-self.xy_data[i][0]))
                    res2 += self.xy_data[i][self.scale-1]*self.xy_data[i][0]
        return (res1+res2*(self.E2/self.E1))/self.N
    def dplot(self,figname):
        X, Y = np.meshgrid(np.linspace(1, self.scale, self.scale), 
                           np.linspace(1, self.scale, self.scale))
        U = np.cos(2 * ma.pi * self.xy_data)
        V = np.sin(2 * ma.pi * self.xy_data)
        fig, ax = plt.subplots(1,1,figsize = (12,12))
        ax.scatter(X,Y,s = 100,c = self.ising_data, marker = 's',cmap = 'Accent')
        ax.quiver(X, Y, U, V, color="k", angles='xy',
                  scale_units='xy', scale=1, width=0.0015)
        ax.set(title = 'E1 = {},E2 = {}'.format(self.E1, self.E2))
        fig.savefig(figname)

        
    
if __name__ == '__main__':
    a = Ising_glass(100,10, 0.1,seed = 514)
    flag = 0
    ana = 0
    
    for i in range(2000):
        ana += a.Wolff(1)
        if i % 10 == 9:
            print(ana/10,a.entropy()[2])
            
            ana = 0
    a.dplot('wolff0')
    
    
    
    
    
    
    