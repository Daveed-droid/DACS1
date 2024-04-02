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
#from AssignmentData import *
from Q2TestData import *

def Q2a(Lamina_mean = Lamina_mean, Laminate = Laminate):
    #LayUp = [0,90,45,-45]
    #LayUp = np.append(LayUp, np.flip(LayUp))
    #LayUp = np.append(LayUp, np.flip(LayUp))
    LayUp = np.array([0,45,-45,90,90,-45,45,0]) # Test to compare to lecture
    angle = np.arange(0,10,20)
    sig = np.array([Xt_mean,Xc_mean,Yt_mean,Yc_mean,S_mean])
    strength = np.array([sig,sig,sig,sig])
    Nf = np.zeros([2,len(angle)])
    Nl = np.zeros([2, len(angle)])

    for i in range(0,len(angle)):
        failure = np.zeros(4)
        print("failure:", failure)
        dn = 50000 # Load increment
        lpf = False
        fpf = False
        Farg  = dn  # Initial length of load vector Ns-Ny
        #Ns = dn * np.sin(np.deg2rad(angle[i]))

        #Ny = dn * np.cos(np.deg2rad(angle[i]))
        #Load = np.array([0, Ny, Ns, 0, 0, 0]).T

        plys = Laminate(LayUp, Lamina_mean)
        #stresses = plys.calcPlyStresses(Load)
        #print("ABD", plys.ABD)
        Q = np.array([plys.QGlobalAr[0], plys.QGlobalAr[1], plys.QGlobalAr[2], plys.QGlobalAr[3]])
        while lpf == False:

            Ns = Farg*np.sin(np.deg2rad(angle[i]))
            print("Ns",Ns)
            Ny = Farg*np.cos(np.deg2rad(angle[i]))
            print("Ny", Ny)
            #Load = np.array([0, Ny, Ns, 0, 0, 0]).T
            Load = np.array([Farg, 0, 0, 0, 0, 0]).T    #Test
            stresses = plys.calcPlyStresses(Load)
            print("poission:", plys.calcEngConst()[2])
            print("inverse ABD:",np.linalg.inv(plys.ABD[0:3,0:3]))
            print("stress:", stresses[:,0:4])
            print("angle:", angle[i])
            print("test")
            #Apply max stress failure criteria
            for j in range(0,4):   # Max stress Failure
                print("ply:", LayUp[j])
                if stresses[0,j] > 0 and strength[j,0] != 0:
                    f_xm = stresses[0,j] / strength[j,0]#Xt_mean
                elif stresses[0,j] < 0 and strength[j,0] != 0:
                    f_xm = stresses[0,j] / -strength[j,1]#Xc_mean

                if stresses[1,j] > 0 and strength[j,0] != 0:
                    f_ym = stresses[1,j] / strength[j,2]#Yt_mean
                elif stresses[1,j] < 0 and strength[j,0] != 0:
                    f_ym = stresses[1,j] / -strength[j,3]#Yc_mean
                if stresses[2,j] > 0 and strength[j,0] != 0:
                    f_sm = stresses[2,j] / strength[j,4]#S_mean
                elif stresses[2,j] < 0 and strength[j,0] != 0:
                    f_sm = stresses[2,j] / -S_mean
                print("Failure value")
                print(f_xm,f_ym,f_sm)
                # Determine which failure mode, and apply deg rule
                if (f_xm >= 1 and failure[j] <= 1) or ((f_ym >= 1 or f_sm >= 1) and failure[j] == 1):
                    print("A")
                    failure[j] = failure[j] + 1
                    print("failure:", failure)
                    if fpf == False:
                        Nf[0,i] = Ny
                        Nf[1,i] = Ns
                        fpf = True
                        print("fpf")
                    elif np.count_nonzero(failure) == 4:
                        Nl[0, i] = Ny
                        Nl[1, i] = Ns
                        lpf = True
                        print("Last ply failure!!")
                        break
                    if LayUp[j] == 45 or LayUp[j] == -45:
                        failure[j] = failure[j]-1
                        #failure[2:4] = failure[2:4]+1
                        failure[1:3] = failure[1:3] + 1 #Test
                        #plys.ABD[0:3,0:3] = plys.ABD[0:3,0:3] - 8 * Q[j]*t
                        plys.ABD[0:3, 0:3] = plys.ABD[0:3, 0:3] - 4 * Q[j] * t  #Test
                        #Strength[2:4,:] = 0
                        strength[1:3, :] = 0    #Test
                    else:
                        #plys.ABD[0:3, 0:3] = plys.ABD[0:3, 0:3] - 4 * Q[j]*t
                        plys.ABD[0:3, 0:3] = plys.ABD[0:3, 0:3] - 2 * Q[j] * t  #Test
                        # Strength[j,:] = 0
                        strength[j, :] = 0  #Test
                    break

                elif (f_ym >= 1 or f_sm >= 1) and failure[j] == 0:
                    print("B")
                    failure[j] = failure[j] + 1
                    print("failure:", failure)
                    if fpf == False:
                        Nf[0,i] = Ny
                        Nf[1,i] = Ns
                        fpf = True
                        print("fpf")
                        """
                    elif np.count_nonzero(failure) == 4:
                        Nl[0, i] = Ny
                        Nl[1, i] = Ns
                        lpf = True
                        print("Last ply failure!!")
                        break
                        """
                    if LayUp[j] == 45 or LayUp[j] == -45:
                        failure[j] = failure[j] - 1
                        #failure[2:4] = failure[2:4] + 1
                        failure[1:3] = failure[1:3] + 1 #Test
                        Q[j] = Laminate([45], Lamina(t, E1_mean, E2_mean*0.1, v12_mean, G12_mean)).QGlobalAr[0]
                        #plys.ABD[0:3,0:3] = plys.ABD[0:3,0:3] - 8*plys.QGlobalAr[j] * t + 8 * Q[j] * t
                        plys.ABD[0:3, 0:3] = plys.ABD[0:3, 0:3] - 4 * plys.QGlobalAr[j] * t + 4 * Q[j] * t  #Test
                        # Strength[2:4,2:4] = 0
                        strength[1:3, 2:4] = 0.1*strength[1:3, 2:4]  # Test
                        print("failure update:",failure)
                    else:
                        Q[j] = Laminate([angle[i]], Lamina(t, E1_mean, E2_mean * 0.1, v12_mean, G12_mean)).QGlobalAr[0]
                        #plys.ABD[0:3,0:3] = plys.ABD[0:3,0:3] - 4*plys.QGlobalAr[j] * t + 4 * Q[j] * t
                        plys.ABD[0:3,0:3] = plys.ABD[0:3,0:3] - 2*plys.QGlobalAr[j] * t + 2 * Q[j] * t  #Test
                        # Strength[j,2:4] = 0
                        strength[j, 2:4] = 0.1*strength[j, 2:4]  # Test
                        print(strength[j,:])
                    break
                # When no failure detected, continue to next ply and increase
                elif j == 3:
                    print("C")
                    Farg = Farg + dn
                    print("Force:", Farg)
                    print("Failures:", failure)
                #lpf = True


    return [Nf, Nl]











Q2a()

print(Q2a()[0])

