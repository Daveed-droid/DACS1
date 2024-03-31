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

#def Q2a(Lamina_mean = Lamina_mean, Laminate = Laminate):
def Q2a(Lamina_mean = Lamina_mean, Laminate = Laminate):
    LayUp = [0,90,45,-45]
    LayUp = np.append(LayUp, np.flip(LayUp))
    LayUp = np.append(LayUp, np.flip(LayUp))

    #Laminate = Laminate(LayUp, Lamina_mean)
    #Laminate.calcABD()
    angle = np.arange(0,181,20)
    sig = np.array([Xt_mean,Yt_mean,Xc_mean,Yc_mean,S_mean])


    for i in angle:

        i = 0
        dn = 50 # Load increment
        lpf = False

        while lpf == False:

            dn = dn*(i+1)
            Ns = dn*np.sin(np.deg2rad(angle[i]))

            Ny = dn*np.cos(np.deg2rad(angle[i]))
            Load = np.array([0, Ny, Ns, 0, 0, 0]).T
            plys = Laminate(LayUp, Lamina_mean)
            stresses = plys.calcPlyStresses(Load)
            print(stresses)
            #Apply stress failure criteria
            for j in LayUp:

            # Puck Fibre Failure
                if stresses[i,0] > 0:
                    f_xff = stresses[i,0] / Xt_mean
                else:
                    f_xff = stresses[i, 0] / -Xc_mean
                if stresses[i,1] > 0:
                    f_yff = stresses[i,1] / Yt_mean
                else:
                    f_yff = stresses[i, 1] / -Yc_mean
                if stresses[i,2] > 0:
                    f_sff = stresses[i,2] / S_mean
                else:
                    f_sff = stresses[i, 2] / -S_mean
                if f_xff,f_yff







            i = i+1
            print(i)
            lpf = True



Q2a()

