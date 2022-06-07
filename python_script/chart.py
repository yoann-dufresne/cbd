#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 14:55:40 2020

@author: muratokutucu
"""

import matplotlib.pyplot as plt
import numpy as np


def getDataFromFile():
    ''' Return a tab composed of two lists :
        - first list = first column of the file (percentages)
        - second list = second column of the file (execution times) '''
    X, Y = [], []
    f = open("perf.txt", "r")
    for l in f:
        tab = l.split("\t")
        X.append(float(tab[0]))
        Y.append(float(tab[1])) 
    return [X, Y]
        

data = getDataFromFile()
X, Y = data[0], data[1]
#print(X)
#print(Y)
plt.axis([0, 100, 0, max(Y)*2])
plt.plot(X, Y)
plt.title("Evolution of taken time to accomplish 10000 calls to successors in function of percentage of k-mers\' presence")
plt.xlabel('Percentage of apparition (%)')
plt.ylabel('Time (microsecondes)')
plt.show()
