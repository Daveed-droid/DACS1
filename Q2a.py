"""
@Project ：DACS1
@File ：Q2a.py
@Author ：David Canosa Ybarra
@Date ：19/03/2024 11:45
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from LaminaClass import Lamina
from LaminateClass import Laminate
from AssignmentData import *
#import AssignmentData
#def Q2a(Lamina_mean = Lamina_mean, Laminate = Laminate):
def Q2a(Lamina_mean = Lamina_mean, Laminate = Laminate):
    LayUp = [0,90,45,-45]
    LayUp = np.append(LayUp, np.flip(LayUp))
    LayUp = np.append(LayUp, np.flip(LayUp))
    print(LayUp)
    Laminate = Laminate(LayUp, Lamina_mean)
    Laminate.calcABD()
    angle = np.arange(0,181,1)
    sig = np.array([Xt_mean,Yt_mean,Xc_mean,Yc_mean,S_mean])
    print(angle)
    """
    for i in angle:

        i = 0
        dn = 50 # Load increment
        while failure = False:

            dn = dn*(i+1)
            Ns = dn*np.sin(np.degtorad(angle(i)))
            Ny = dn*np.cos(np.degtorad(angle(i)))
            Load = np.array([0, Ny, Ns, 0, 0, 0]).T
            stresses = Laminate.calcPlyStresses(Load)
            
            Apply stress failure criteria
            """




    print(Laminate.ABD)
Q2a()

