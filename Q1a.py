"""
@Project ：DACS1
@File ：Q1a.py
@Author ：David Canosa Ybarra
@Date ：19/03/2024 11:45
"""
from LaminaClass import Lamina
from LaminateClass import Laminate
import numpy as np
import matplotlib.pyplot as plt


def Q1a(n):
    E1 = 140 * 10 ** 9
    E2 = 10 * 10 ** 9
    G12 = 5 * 10 ** 9
    v12 = 0.3
    t = 0.125

    v21 = v12 * E2 / E1

    Lamina_ = Lamina(t, E1, E2, v12, G12)
    angle = np.arange(0,90.1,5)
    #print(np.flip(angle))
    laminates = []
    k = 0
    Ex = []
    Ey = []
    vxy = []
    vyx = []
    Gxy = []
    for i in angle:

        LayUp = [15, i, -i, 75, 75]
        #LayUp = [20] #test, not correct compared to solutions using lec 2 eqs for poisson ratio
        for j in range(0, n):

            #LayUp.append(np.flip(LayUp))
            LayUp = np.append(LayUp, np.flip(LayUp))
        plys = Laminate(LayUp, Lamina_)
        laminates.append(plys)

        laminates[k].calcEngConst()
        #print(LayUp)
        #print(laminates[k].calcEngConst()[0]*10**-9, LayUp)
        Ex.append(laminates[k].calcEngConst()[0]*10**-9)
        Ey.append(laminates[k].calcEngConst()[1]*10**-9)
        vxy.append(laminates[k].calcEngConst()[2])
        vyx.append(laminates[k].calcEngConst()[3])
        Gxy.append(laminates[k].calcEngConst()[4]*10**-9)
        k = k + 1



    # Plot the engineering constants as a function of theta
    #print(vxy)
    #print(vyx)
    fig, axs = plt.subplots(2,2)
    axs[0, 0].plot(angle, Ex)
    axs[0, 0].set_title('Ex')
    axs[0, 0].set_ylabel('Ex [GPa]')
    axs[1, 0].plot(angle, Ey)
    axs[1, 0].set_title('Ey')
    axs[1, 0].set_ylabel('Ey [GPa]')
    axs[0, 1].plot(angle, Gxy)
    axs[0, 1].set_title('Gxy')
<<<<<<< Updated upstream
    axs[1, 1].plot(angle, vxy)
    axs[1, 1].plot(angle, vyx)

=======
    axs[0, 1].set_ylabel('Gxy [GPa]')
    axs[1, 1].plot(angle, vxy, label="vxy")
    axs[1, 1].plot(angle, vyx, label="vyx")
>>>>>>> Stashed changes
    axs[1, 1].set_title('vxy, vyx')
    axs[1, 1].set_ylabel('v [-]')
    axs[1, 1].legend()


    for ax in axs.flat:
        ax.set(xlabel="angle")
    plt.tight_layout()
<<<<<<< Updated upstream
    plt.savefig("Q1A")
    plt.show()
=======
    plt.savefig("Plots/engineeringConstants")
>>>>>>> Stashed changes

    print(LayUp)

Q1a(1)


