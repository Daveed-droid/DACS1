"""
@Project ：DACS1
@File ：Q1a.py
@Author ：David Canosa Ybarra
@Date ：19/03/2024 11:45
"""
#from LaminaClass import Lamina
#from LaminateClass import Laminate
import numpy as np
import matplotlib.pyplot as plt


def Q1a(n):
    angle = np.arange(0,8.1,1)
    #print(np.flip(angle))
    for i in angle:
        laminates = []
        LayUp = [15, i, -i, 75, 75]

        for j in range(0, n):

            #LayUp.append(np.flip(LayUp))
            LayUp = np.append(LayUp, np.flip(LayUp))
        plys = Laminate(LayUp, Lamina)
        laminates.append[plys]
        laminates[i].calcEngConst()


    # Plot the engineering constants as a function of theta
    Ex = laminates[:].Ex
    Ey = laminates[:].Ey
    Gxy = laminates[:].Gxy
    vxy = laminates[:].vxy

    fig, axs = plt.subplots(2,2)
    axs[0,0].plot(angle, Ex)
    axs[0,0].set_title('Ex')
    axs[1, 0].plot(angle, Ey)
    axs[1, 0].set_title('Ey')
    axs[0, 1].plot(angle, Gxy)
    axs[0, 1].set_title('Gxy')
    axs[1, 1].plot(angle, vxy)
    axs[1, 1].set_title('vxy')
    for ax in axs.flat:
        ax.set(xlabel="angle")


    print(LayUp)


    #for i in angle:

        #ply = Laminate(LayUp, Lamina)


