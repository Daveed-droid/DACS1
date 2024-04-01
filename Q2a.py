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
    Nf = np.zeros([2,len(angle)])
    Nl = np.zeros([2, len(angle)])
    for i in range(0,len(angle)):

        failure = 0
        dn = 500000 # Load increment
        lpf = False
        fpf = False

        #Ns = dn * np.sin(np.deg2rad(angle[i]))

        #Ny = dn * np.cos(np.deg2rad(angle[i]))
        #Load = np.array([0, Ny, Ns, 0, 0, 0]).T

        plys = Laminate(LayUp, Lamina_mean)
        #stresses = plys.calcPlyStresses(Load)
        print(plys.ABD)
        while lpf == False:


            Ns = dn*np.sin(np.deg2rad(angle[i]))

            Ny = dn*np.cos(np.deg2rad(angle[i]))
            Load = np.array([0, Ny, Ns, 0, 0, 0]).T

            stresses = plys.calcPlyStresses(Load)
            #print(stresses)
            print(angle[i])
            print("test")
            #Apply stress failure criteria
            for j in range(0,len(LayUp)):

            # Puck Fibre Failure
                if stresses[0,j] > 0:
                    f_xff = stresses[0,j] / Xt_mean
                else:
                    f_xff = stresses[0,j] / -Xc_mean
                if stresses[1,j] > 0:
                    f_yff = stresses[1,j] / Yt_mean
                else:
                    f_yff = stresses[1,j] / -Yc_mean
                if stresses[2,j] > 0:
                    f_sff = stresses[2,j] / S_mean
                else:
                    f_sff = stresses[2,j] / -S_mean



                if angle[i] <90:   #Mode A
                    pass

                elif angle[i] >= 90:# and stresses[i,1] >= -?? Mode B
                    pass

                else:#Mode C
                    pass

                # Determine which failure mode, and apply deg rule
                if f_xff >= 1 or f_yff >= 1:
                    if fpf == False:
                        Nf[0,i] = Ny
                        Nf[1,i] = Ns
                        fpf = True
                    elif failure >= 3:
                        Nl[0, i] = Ny
                        Nl[1, i] = Ns
                        lpf = True
                    if LayUp[j] == 45 or LayUp[j] == -45:
                        plys.ABD[0:3,0:3] = plys.ABD[0:3,0:3] - 8*plys.QGlobalAr[j]*t
                    else:
                        plys.ABD[0:3, 0:3] = plys.ABD[0:3, 0:3] - 8 * plys.QGlobalAr[j]*t
                    failure = failure +1

                    break

            dn = dn * (i + 1)

    return [Nf, Nl]











Q2a()

print(Q2a()[0])

