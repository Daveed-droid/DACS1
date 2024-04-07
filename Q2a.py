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
    #LayUp = [0,45,55,12]
    LayUp = np.append(LayUp, np.flip(LayUp))
    LayUp = np.append(LayUp, np.flip(LayUp))
    #LayUp = np.array([0,45,-45,90,90,-45,45,0]) # Test to compare to lecture
    angle = np.arange(0,361,10)
    print(angle)
    sig = np.array([Xt_mean,Xc_mean,Yt_mean,Yc_mean,S_mean])
    strength = np.array([sig,sig,sig,sig])
    Nf = np.zeros([4,len(angle)])
    Nl = np.zeros([4, len(angle)])
    D = 50000 # Main step size

    for i in range(0,len(angle)):

        failure = np.zeros(4)
        print("failure:", failure)
        dn = D # Load increment
        lpf = False
        fpf = False
        Farg  = dn  # Initial length of load vector Ns-Ny
        plys = Laminate(LayUp, Lamina_mean)
        strength = np.array([sig, sig, sig, sig])
        #stresses = plys.calcPlyStresses(Load)
        #print("ABD", plys.ABD)
        Q = np.array([plys.QGlobalAr[0], plys.QGlobalAr[1], plys.QGlobalAr[2], plys.QGlobalAr[3]])
        Qxyz = np.array([plys.QGlobalAr[0], plys.QGlobalAr[1], plys.QGlobalAr[2], plys.QGlobalAr[3]])
        Q123 = np.array([Lamina_mean.calcQMatrix(), Lamina_mean.calcQMatrix(), Lamina_mean.calcQMatrix(),
                         Lamina_mean.calcQMatrix()])

        while lpf == False:

            Ns = Farg*np.sin(np.deg2rad(angle[i]))
            print("Ns",Ns)
            Ny = Farg*np.cos(np.deg2rad(angle[i]))
            print("Ny", Ny)
            Load = np.array([0, Ny, Ns, 0, 0, 0]).T
            #Load = np.array([Farg, 0, 0, 0, 0, 0]).T    #Test
            stresses = plys.calcPlyStresses2(Load)
            print("poission:", plys.calcEngConst()[2])
            print("inverse ABD:",np.linalg.inv(plys.ABD[0:3,0:3]))
            print("stress:", stresses[:,0:4])
            print("angle:", angle[i])
            print("test")
            f_xm = 0
            f_ym = 0
            f_sm = 0
            #Apply max stress failure criteria
            for j in range(0,4):   # Max stress Failure
                if failure[j] == 2:
                    f_FF = 0
                    f_IFF = 0
                    if j == 3:
                        print("C")
                        Farg = Farg + dn
                        print("Force:", Farg)
                        print("Failures:", failure)
                    continue

                if (f_xm-1 >= 0.01 or f_ym-1 >= 0.01 or f_sm-1 >= 0.01) and fpf == False:

                    Farg = Farg - dn
                    dn = dn / 2
                    break

                print("ply:", LayUp[j])
                print(f_xm, f_ym, f_sm)
                print("Strength:", strength[j,:])
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
                    failure[j] = 2
                    print("failure:", failure)
                    if fpf == False:
                        Nf[0,i] = Ny
                        Nf[1,i] = Ns
                        fpf = True
                        print("fpf")



                    plys.ABD[0:3, 0:3] = plys.ABD[0:3, 0:3] - 4 * Q[j]*t
                    #plys.ABD[0:3, 0:3] = plys.ABD[0:3, 0:3] - 2 * Q[j] * t  #Test
                    strength[j,:] = 0
                    print(strength[j])
                    #strength[j, :] = 0  #Test

                    if np.count_nonzero(failure == 2) == 4:
                        Nl[0, i] = Ny
                        Nl[1, i] = Ns
                        lpf = True
                        print("Last ply failure!!")
                    dn = D
                    break

                elif (f_ym >= 1 or f_sm >= 1) and failure[j] == 0:
                    print("B")
                    failure[j] = 1
                    print("failure:", failure)
                    if fpf == False:
                        Nf[0,i] = Ny
                        Nf[1,i] = Ns
                        fpf = True
                        print("fpf")





                    Q[j] = Laminate([angle[i]], Lamina(t, E1_mean, E2_mean * 0.1, v12_mean, G12_mean)).QGlobalAr[0]
                    plys.ABD[0:3,0:3] = plys.ABD[0:3,0:3] - 4*plys.QGlobalAr[j] * t + 4 * Q[j] * t
                    #plys.ABD[0:3,0:3] = plys.ABD[0:3,0:3] - 2*plys.QGlobalAr[j] * t + 2 * Q[j] * t  #Test
                    strength[j,2:4] = 0.1*strength[j,2:4]
                    #strength[j, 2:4] = 0.1*strength[j, 2:4]  # Test
                    print(strength[j,:])

                    if np.count_nonzero(failure == 2) == 4:
                        Nl[0, i] = Ny
                        Nl[1, i] = Ns
                        lpf = True
                        print("Last ply failure!!")
                    dn = D
                    break
                # When no failure detected, continue back to first ply and increase load
                elif j == 3:
                    print("C")
                    Farg = Farg + dn
                    print("Force:", Farg)
                    print("Failures:", failure)
                #lpf = True

        



        #Puck
        print("---------------------------------------------------- PUCK -----------------------------------------------------")
        failure = np.zeros(4)
        #print("failure:", failure)

        lpf = False
        fpf = False
        dn = D
        Farg = dn# * 10
        Ns = Farg * np.sin(np.deg2rad(angle[i]))
        print("Ns", Ns)
        Ny = Farg * np.cos(np.deg2rad(angle[i]))
        print("Ny", Ny)
        plys = Laminate(LayUp, Lamina_mean)
        strength = np.array([sig, sig, sig, sig])
        # stresses = plys.calcPlyStresses(Load)
        # print("ABD", plys.ABD)
        Qxyz = np.array([plys.QGlobalAr[0], plys.QGlobalAr[1], plys.QGlobalAr[2], plys.QGlobalAr[3]])
        Q123 = np.array([Lamina_mean.calcQMatrix(), Lamina_mean.calcQMatrix(), Lamina_mean.calcQMatrix(), Lamina_mean.calcQMatrix()])
        while lpf == False:
            Ns = Farg * np.sin(np.deg2rad(angle[i]))
            #print("Ns", Ns)
            Ny = Farg * np.cos(np.deg2rad(angle[i]))
            #print("Ny", Ny)
            print(Farg)
            Load = np.array([0, Ny, Ns, 0, 0, 0]).T
            # Load = np.array([Farg, 0, 0, 0, 0, 0]).T    #Test
            stresses = plys.calcPlyStresses2(Load)
            
            #print("poission:", plys.calcEngConst()[2])
            #print("inverse ABD:", np.linalg.inv(plys.ABD[0:3, 0:3]))
            #print("stress:", stresses[:, 0:4])
            print("angle:", angle[i])
            #print("test")
            f_FF = 0
            f_IFF = 0
            # Apply max stress failure criteria

            for j in range(0, 4):  # Max stress Failure

                if failure[j] == 2:
                    f_FF = 0
                    f_IFF = 0
                    if j == 3:
                        print("C")
                        Farg = Farg + dn
                        print("Force:", Farg)
                        print("Failures:", failure)
                    continue
                print("ply:", LayUp[j])

                #print("Strength:", strength[j, :])
                #Puck criteria
                f_FF,f_IFF = plys.Puck(stresses[0:3,j], Ply_strength)
                #f_IFF = plys.Puck(stresses[0:3,j], Ply_strength)[1]


                print("Failure value")
                print(f_FF, f_IFF)

                if (f_FF-1 >= 0.01 or f_IFF-1 >= 0.01) and fpf == False:

                    Farg = Farg - dn
                    dn = dn / 2
                    break


                # Determine which failure mode, and apply deg rule
                if (f_FF >= 1 and failure[j] <= 1) or (f_IFF >= 1 and failure[j] == 1):
                    print("A")
                    failure[j] = 2
                    print("failure:", failure)
                    if fpf == False:
                        Nf[2, i] = Ny
                        Nf[3, i] = Ns
                        fpf = True
                        print("fpf")

                    plys.ABD[0:3, 0:3] = plys.ABD[0:3, 0:3] - 4 * Qxyz[j] * t

                    #strength[j, :] = 0

                    print("strength", strength[j])
                    print("stress", stresses[0:3,j])


                    if np.count_nonzero(failure == 2) == 4:
                        Nl[2, i] = Ny
                        Nl[3, i] = Ns
                        lpf = True
                        print("Last ply failure!!")
                    dn = D
                    break
                #Test for inter fiber fracture
                elif f_IFF >= 1 and failure[j] == 0:
                    print("B")
                    failure[j] = 1
                    print("failure:", failure)
                    if fpf == False:
                        Nf[2, i] = Ny
                        Nf[3, i] = Ns
                        fpf = True
                        print("fpf")
                    Q123[j] = Lamina(t, E1_mean, E2_mean * 0.1, v12_mean, G12_mean)
                    Qxyz[j] = Laminate([angle[i]], Lamina(t, E1_mean, E2_mean * 0.1, v12_mean, G12_mean)).QGlobalAr[0]
                    plys.ABD[0:3, 0:3] = plys.ABD[0:3, 0:3] - 4 * plys.QGlobalAr[j] * t + 4 * Qxyz[j] * t

                    #strength[j, 2:4] = 0.1 * strength[j, 2:4]

                    print(strength[j, :])

                    if np.count_nonzero(failure == 2) == 4:
                        Nl[2, i] = Ny
                        Nl[3, i] = Ns
                        lpf = True
                        print("Last ply failure!!")
                    dn = D
                    break
                # When no failure detected on the last ply, continue back to first ply and increase load
                elif j == 3:
                    print("C")
                    Farg = Farg + dn
                    print("Force:", Farg)
                    print("Failures:", failure)
                # lpf = True

    plt.plot(Nf[0],Nf[1],marker="<", label="fpf Max")
    plt.plot(Nl[0], Nl[1],marker=">", label="lpf Max")
    plt.plot(Nf[2], Nf[3],marker="o", label="fpf Puck")
    plt.plot(Nl[2], Nl[3],marker="s", label="lpf Puck")
    plt.legend()
    plt.show()
    print(Nf)
    print(Nl)
    return [Nf, Nl]
Q2a()

#print(Q2a()[0])

