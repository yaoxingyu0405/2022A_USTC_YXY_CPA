# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 22:33:15 2023

@author: YXY
"""

from ising_glass import Ising_glass
from mpi4py import MPI
import numpy as np
import math as ma

k = 20
N = 58
e2 = 1.0

def mpi_ising(rank):
    if rank:
        while True:
            E1,E2 = comm.recv(source = 0)
            if E1<0 or E2<0:
                break
            simulation = Ising_glass(50, E1, E2,seed = 2004)
            entropy_res = np.zeros((100),dtype = 'float')
            energy_res = np.zeros((100),dtype = 'float')
            M1_res = np.zeros((100),dtype = 'float')
            M2_res = np.zeros((100),dtype = 'float')
            flag0 = True
            flag1 = True
            n1 = 0
            n2 = 0
            for n in range(2000):
                for i in range(k):
                    if flag0:
                        n1 += simulation.Wolff(0)
                    if flag1:
                        n2 += simulation.Wolff(1)
                ana1,ana2,ana3 = simulation.entropy()
                if ana2 > 0.99 or ana2 <0.01:
                    flag0 = False
                n1 /= simulation.N * k
                n2 /= simulation.N * k
                if n > 10:
                    if n1 > 0.99 or n1 < 0.01:
                        flag0 = False
                    if n2 > 0.99 or n2 < 0.01:
                        flag1 = False
                if not(flag0 or flag1):
                    break
            energy_res[0] = -simulation.energy()
            entropy_res[0] = ana1
            M1_res[0] = abs(2*ana2-1)
            M2_res[0] = ana3
            for i in range(99):
                for j in range(5):
                    simulation.Wolff(0)
                    simulation.Wolff(1)
                ana1,ana2,ana3 = simulation.entropy()
                energy_res[i+1] = -simulation.energy()
                entropy_res[i+1] = ana1
                M1_res[i+1] = abs(2*ana2-1)
                M2_res[i+1] = ana3
            msg = [1/E1, energy_res.sum()/100, entropy_res.sum()/100, M1_res.sum()/100,
                   M2_res.sum()/100, ma.sqrt(((energy_res-energy_res.sum()/100)**2).sum()/99),
                   ma.sqrt(((entropy_res-entropy_res.sum()/100)**2).sum()/99),
                   ma.sqrt(((M1_res-M1_res.sum()/100)**2).sum()/99),
                   ma.sqrt(((M2_res-M2_res.sum()/100)**2).sum()/99),rank]
            comm.send(msg,dest = 0)
    else:
        msg_list = 1/np.linspace(0.1, 3.0,N+1)
        res = []
        for i in range(1,16):
            comm.send((msg_list[i-1],e2*msg_list[i-1]),dest = i)
        for i in msg_list[15:]:
            temp = comm.recv()
            res += [temp[:-1]]
            comm.send((i,i*e2),dest = temp[-1])
        for i in range(15):
            temp = comm.recv()
            res += [temp[:-1]]
            comm.send((-1,-1),dest = temp[-1])
        res.sort(key = lambda x:x[0])
        f = open('C:\\Users\\YXY\\Desktop\\ising_glass\\resf2.txt','w')
        f.write('E_1 = {}E_2\nT\t\tE\t\tSigmaE\tS\t\tSigmaS\tM1\t\tSigmaM1\tM2\t\tSigmaM2\n'.format(e2))
        for i in res:
            f.write('{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}\n'
                    .format(i[0],i[1],i[5],i[2],i[6],i[3],i[7],i[4],i[8]))
        f.close()

            
if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    mpi_ising(rank)