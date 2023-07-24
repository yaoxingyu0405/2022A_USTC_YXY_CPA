# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 21:01:40 2022

@author: YXY
"""

import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 1, figsize = (6,4))
# a1 = [1,2]
# b1 = [-7.3,-22.3]
# a2 = [1,2,3]
# b2 = [-5.19,-15.98,-48.36]
# a3 = [1,2,3,4,5,6]
# b3 = [1.10,-0.35,-2.27,-5.52,-11.88,-24.61]
# a4 = [1,2,3,4,5,6,7,8,9]
# b4 = [4.23,2.99,1.72,0.35,-1.2,-3.76,-8.4,-17.64,-36.04]
# ax.plot(a1,b1,color = 'r',linestyle = '-',marker = '.',label = 'line1')
# ax.plot(a2,b2,color = 'b',linestyle = '-',marker = '.',label = 'line2')
b3 = [14.77,-0.42,-0.42,1.11,-1.29,-0.94,-2.89,-4.21,-7.91,-12.97,-21.72,-36.04]
a3 = [1,2,3,4,5,6,7,8,9,10,11,12]
b4 = [-0.41,-0.41,3.07,-0.42,-0.54,1.22,-0.90,-1.43,-2.47,-4.59,-7.96,-13.39,-22.19,-36.04]
a4 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
ax.plot(a3,b3,color = 'r',linestyle = '-',marker = '.',label = 'line3')
ax.plot(a4,b4,color = 'b',linestyle = '-',marker = '.',label = 'line4')
ax.legend()
ax.set(title = 'Secant Method',ylabel = 'ln f(x_n)',xlabel='n')

