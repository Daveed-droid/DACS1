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
#from Q2TestData import *

def Q2a(Lamina_mean = Lamina_mean, Laminate = Laminate):
    LayUp = [0,90,45,-45]
    #LayUp = [0,-45,30,80]
    LayUp = np.append(LayUp, np.flip(LayUp))
    LayUp = np.append(LayUp, np.flip(LayUp))
    #LayUp = np.array([0,45,-45,90,90,-45,45,0]) # Test to compare to lecture
    angle = np.arange(230,361,1000)
    print(angle)
    sig = np.array([Xt_mean,Xc_mean,Yt_mean,Yc_mean,S_mean])
    strength = np.array([sig,sig,sig,sig])
    Nf = np.zeros([4,len(angle)])
    Nl = np.zeros([4, len(angle)])
    D = 50000 # Main step size

    for i in range(0,len(angle)):

        failure = np.zeros(4)

        dn = D # Load increment
        lpf = False
        fpf = False
        Farg  = dn  # Initial length of load vector Ns-Ny
        plys = Laminate(LayUp, Lamina_mean)
        strength = np.array([sig, sig, sig, sig])
        #stresses = plys.calcPlyStresses(Load)
        print("ABD", plys.ABD)
        Q = np.array([plys.QGlobalAr[0], plys.QGlobalAr[1], plys.QGlobalAr[2], plys.QGlobalAr[3]])
        Qxyz = np.array([plys.QGlobalAr[0], plys.QGlobalAr[1], plys.QGlobalAr[2], plys.QGlobalAr[3]])
        Q123 = np.array([Lamina_mean.calcQMatrix(), Lamina_mean.calcQMatrix(), Lamina_mean.calcQMatrix(),
                         Lamina_mean.calcQMatrix()])

        while lpf == False:

            Ns = Farg*np.sin(np.deg2rad(angle[i]))

            Ny = Farg*np.cos(np.deg2rad(angle[i]))
            Load = np.array([0, Ny, Ns, 0, 0, 0]).T
            #Load = np.array([Ny, Ns, 0, 0, 0, 0]).T    #Test
            stresses = plys.calcPlyStresses2(Load)
            f_xm = 0
            f_ym = 0
            f_sm = 0
            #Apply max stress failure criteria
            for j in range(0,4):   # Max stress Failure
                if failure[j] == 2:
                    f_xm = 0
                    f_ym = 0
                    f_sm = 0
                    if j == 3:

                        Farg = Farg + dn

                    continue

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


                if (f_xm-1 >= 0.01 or f_ym-1 >= 0.01 or f_sm-1 >= 0.01) and fpf == False:

                    Farg = Farg - dn
                    dn = dn / 2
                    break

                # Determine which failure mode, and apply deg rule

                if (f_xm >= 1 and failure[j] <= 1) or ((f_ym >= 1 or f_sm >= 1) and failure[j] == 1):
                    failure[j] = 2
                    if fpf == False:
                        Nf[0,i] = Ny
                        Nf[1,i] = Ns
                        fpf = True



                    plys.ABD[0:3, 0:3] = plys.ABD[0:3, 0:3] - 4 * Q[j]*t
                    #plys.ABD[0:3, 0:3] = plys.ABD[0:3, 0:3] - 2 * Q[j] * t  #Test
                    strength[j,:] = 0
                    #strength[j, :] = 0  #Test

                    if np.count_nonzero(failure == 2) == 4:
                        Nl[0, i] = Ny
                        Nl[1, i] = Ns
                        lpf = True
                    dn = D
                    break

                elif (f_ym >= 1 or f_sm >= 1) and failure[j] == 0:
                    failure[j] = 1
                    if fpf == False:
                        Nf[0,i] = Ny
                        Nf[1,i] = Ns
                        fpf = True

                    Q[j] = Laminate([angle[i]], Lamina(t, E1_mean, E2_mean * 0.1, v12_mean, G12_mean)).QGlobalAr[0]
                    plys.ABD[0:3,0:3] = plys.ABD[0:3,0:3] - 4*plys.QGlobalAr[j] * t + 4 * Q[j] * t
                    #plys.ABD[0:3,0:3] = plys.ABD[0:3,0:3] - 2*plys.QGlobalAr[j] * t + 2 * Q[j] * t  #Test
                    strength[j,2:4] = 0.1*strength[j,2:4]
                    #strength[j, 2:4] = 0.1*strength[j, 2:4]  # Test

                    if np.count_nonzero(failure == 2) == 4:
                        Nl[0, i] = Ny
                        Nl[1, i] = Ns
                        lpf = True
                    dn = D
                    break
                # When no failure detected, continue back to first ply and increase load
                elif j == 3:
                    Farg = Farg + dn



        #Puck
        #print("---------------------------------------------------- PUCK -----------------------------------------------------")
        failure = np.zeros(4)


        lpf = False
        fpf = False
        dn = D
        Farg = dn
        Ns = Farg * np.sin(np.deg2rad(angle[i]))

        Ny = Farg * np.cos(np.deg2rad(angle[i]))

        plys = Laminate(LayUp, Lamina_mean)
        strength = np.array([sig, sig, sig, sig])
        # stresses = plys.calcPlyStresses(Load)
        StrainLst = []
        Qxyz = np.array([plys.QGlobalAr[0], plys.QGlobalAr[1], plys.QGlobalAr[2], plys.QGlobalAr[3]])
        Q123 = np.array([Lamina_mean.calcQMatrix(), Lamina_mean.calcQMatrix(), Lamina_mean.calcQMatrix(), Lamina_mean.calcQMatrix()])
        while lpf == False:
            Ns = Farg * np.sin(np.deg2rad(angle[i]))

            Ny = Farg * np.cos(np.deg2rad(angle[i]))

            Load = np.array([0, Ny, Ns, 0, 0, 0]).T


            # Apply max stress failure criteria
            Lamina_mean.setStrengths(Xt_mean, Yt_mean, Xc_mean, Yc_mean, S_mean)
            f_FF, f_IFF = plys.Puck(Load)
            for j in range(0, 4):  # Max stress Failure
                if failure[j] == 2:
                    f_FF[j] = 0
                    f_IFF[j] = 0
                    if j == 3:

                        Farg = Farg + dn

                    continue

                #Puck criteria
                # f_FF,f_IFF = plys.Puck(stresses[0:3,j], Ply_strength)
                #f_IFF = plys.Puck(stresses[0:3,j], Ply_strength)[1]


                if (f_FF[j]-1 >= 0.01 or f_IFF[j]-1 >= 0.01) and fpf == False:

                    Farg = Farg - dn
                    dn = dn / 2
                    break


                # Determine which failure mode, and apply deg rule
                if (f_FF[j] >= 1 and failure[j] <= 1) or (f_IFF[j] >= 1 and failure[j] == 1):

                    failure[j] = 2

                    if fpf == False:
                        Nf[2, i] = Ny
                        Nf[3, i] = Ns
                        fpf = True
                    plys.ABD[0:3, 0:3] = plys.ABD[0:3, 0:3] - 4 * Qxyz[j] * t

                    #strength[j, :] = 0

                    if np.count_nonzero(failure == 2) == 4:
                        Nl[2, i] = Ny
                        Nl[3, i] = Ns
                        lpf = True
                    dn = D
                    break
                #Test for inter fiber fracture
                elif f_IFF[j] >= 1 and failure[j] == 0:

                    failure[j] = 1

                    if fpf == False:
                        Nf[2, i] = Ny
                        Nf[3, i] = Ns
                        fpf = True

                    Q123[j] = Lamina(t, E1_mean, E2_mean * 0.1, v12_mean, G12_mean)
                    Qxyz[j] = Laminate([angle[i]], Lamina(t, E1_mean, E2_mean * 0.1, v12_mean, G12_mean)).QGlobalAr[0]
                    plys.ABD[0:3, 0:3] = plys.ABD[0:3, 0:3] - 4 * plys.QGlobalAr[j] * t + 4 * Qxyz[j] * t

                    #strength[j, 2:4] = 0.1 * strength[j, 2:4]


                    if np.count_nonzero(failure == 2) == 4:
                        Nl[2, i] = Ny
                        Nl[3, i] = Ns
                        lpf = True

                    dn = D
                    break
                # When no failure detected on the last ply, continue back to first ply and increase load
                elif j == 3:

                    Farg = Farg + dn



    plt.plot(Nf[0]*10**-3,Nf[1]*10**-3,marker="<", label="fpf Max")
    plt.plot(Nl[0]*10**-3, Nl[1]*10**-3,marker=">", label="lpf Max")
    plt.legend()
    plt.title("Stress Failure Envelope, Max")
    plt.xlabel("Ny[N/mm]")
    plt.ylabel("Ns[N/mm]")
    plt.savefig("Q2aMax.png")
    plt.show()

    plt.plot(Nf[2]*10**-3, Nf[3]*10**-3,marker="o", label="fpf Puck")
    plt.plot(Nl[2]*10**-3, Nl[3]*10**-3,marker="s", label="lpf Puck")
    plt.legend()
    plt.title("Stress Failure Envelope, Puck")
    plt.xlabel("Ny[N/mm]")
    plt.ylabel("Ns[N/mm]")
    plt.savefig("Q2aPuck.png")
    plt.show()
    #print(Nf)
    #print(Nl)

    return [Nf, Nl]
Q2a()

#print(Q2a())

