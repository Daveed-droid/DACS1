"""
@Project ：DACS1
@File ：Q3.py
@Author ：David Canosa Ybarra
@Date ：19/03/2024 11:46
"""
import pickle
import numpy as np
from AssignmentData import Lamina_props
from LaminaClass import Lamina
from LaminateClass import Laminate
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import norm
import multiprocessing


def Q3(LayUp, Load, n = 100):
    # Generate fpf values for prescribed lamina
    FPF = np.zeros((2,n), dtype = float)
    for n_i in tqdm(range(n)):
        # Initialize randomly picked lamina properties
        LaminaLst=[]
        for i in range(len(LayUp)):
            # Random pick
            props = []
            for [mean, std] in Lamina_props:
                props.append(np.random.normal(mean, std))
            E1, E2, v12, G12 = props[0:4]
            Lamina_ = Lamina(t, E1, E2, v12, G12)
            Xt, Xc, Yt, Yc, S = props[4:9]
            Lamina_.setStrengths(Xt, Yt, Xc, Yc, S)
            LaminaLst.append(Lamina_)
        Laminate_ = Laminate(LayUp, LaminaLst)
        LoadFPF= Laminate_.calcFailurePuck(Load[:,0])
        LoadNorm1 = np.linalg.norm(Load[:,0])
        LoadNorm2 = np.linalg.norm(Load[:,1])
        LoadFPFNorm = np.linalg.norm(LoadFPF)
        FPF[0, n_i] = LoadFPFNorm/LoadNorm1
        FPF[1, n_i] = LoadFPFNorm/LoadNorm2
    # Store Data
    # for i in range(len(Load[0])):
    #     with open(f'Q3_DataAnalysis/Data/reliability_analysis_Load{int(round(np.linalg.norm(Load[:,i])/1e3, 0))}_n{n}_FPF.pkl', 'wb') as f:  # open a text file
    #         pickle.dump(FPF[i, :], f) # serialize the list
    #     f.close()
    #     (mean, std) = norm.fit(FPF[i, :])
    #     with open(f'Q3_DataAnalysis/Data/reliability_analysis_Load{int(round(np.linalg.norm(Load[:,i])/1e3, 0))}_n{n}_meanStd.pkl', 'wb') as f:  # open a text file
    #         pickle.dump([mean, std], f) # serialize the list
    #     f.close()
    # x = np.linspace(mean - 8*std, mean + 8*std, 100)
    plt.hist(FPF[0, :])
    plt.show()
    plt.clf()
    plt.hist(FPF[1, :])
    plt.show()
    plt.clf()
    # plt.plot(x, norm.pdf(x, mean, std))
    # plt.show()
    # print("Failure Prob: ", norm.cdf(1, mean, std))
    # return mean, std, norm.cdf(1, mean, std)


if __name__ == '__main__':
    # Set seed
    np.random.seed(0)
    LayUp = np.array([0, 90, 45, -45])
    LayUp = np.append(LayUp, np.flip(LayUp))
    LayUp = np.append(LayUp, np.flip(LayUp))
    Load1 = (np.array([[np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30)), 0, 0, 0, 0]])*1200).T
    Load1 = Load1*1000 # Convert to N/m
    Load2 = (np.array([[np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30)), 0, 0, 0, 0]])*850).T
    Load2 = Load2*1000 # Convert to N/m
    print(Load1)
    Load = np.hstack((Load1, Load2))
    print()
    t = 0.125e-3

    # Q3(LayUp, Load, n=12800000) # ~ 4 hr
    # Q3(LayUp, Load, n=6400000) # ~ 2 hr
    # Q3(LayUp, Load, n=3200000) # ~ 1 hr ----- 7hr
    # Q3(LayUp, Load, n=1600000) # ~ 32 min
    # Q3(LayUp, Load, n=800000) # ~ 16 min
    # Q3(LayUp, Load, n=400000) # ~ 8 min
    # Q3(LayUp, Load, n=200000) # ~ 4 min
    # Q3(LayUp, Load, n=100000) # ~ 2 min
    # Q3(LayUp, Load, n=50000) # ~ 1 min
    Q3(LayUp, Load, n=25000) # ~ 30 sec
    # Q3(LayUp, Load, n=12500) # ~ 15 sec ----- 8hr
    # for i in range(len(Lamina_props)):
    #     print(fr"{Lamina_props[i,0]} & {Lamina_props[i,1]} \\ \hline")












    # n_start = 1000
    # meanB = 0
    # stdB = 0
    # PfB = 1
    # n = n_start
    # mean, std, Pf = Q3(LayUp, Load, n=n)
    # while not (np.isclose(meanB, mean, rtol = 0.01) and np.isclose(stdB, std, rtol = 0.01) and np.isclose(PfB, Pf, rtol = 0.01)):
    #     meanB = mean
    #     stdB = std
    #     PfB = Pf
    #     n = n*2
    #     mean, std, Pf = Q3(LayUp, Load, n=n)
    #     print(f"Pf B {PfB} Pf {Pf}")
    #     print(np.isclose(PfB, Pf, rtol = 0.01))

