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


def Q3(LayUp, Load, n = 100):
    FPF = []
    for n_i in tqdm(range(n)):
        LaminaLst=[]
        for i in range(len(LayUp)):
            # Random pick
            props = []
            for [mean, std] in Lamina_props:
                props.append(np.random.normal(mean, std))
            E1, E2, v12, G12 = props[0:4]
            Lamina_ = Lamina(t, E1, E2, v12, G12)
            Xt, Yt, Xc, Yc, S = props[4:9]
            Lamina_.setStrengths(Xt, Yt, Xc, Yc, S)
            LaminaLst.append(Lamina_)
        Laminate_ = Laminate(LayUp, LaminaLst)
        LoadFPF= Laminate_.calcFailurePuck(Load, dL_step = 5000)
        LoadNorm = np.linalg.norm(Load)
        LoadFPFNorm = np.linalg.norm(LoadFPF)
        FPF.append(LoadFPFNorm/LoadNorm)
    with open(f'Q3_DataAnalysis/Data/reliability_analysis_Load{round(np.linalg.norm(Load), 0)}_n{n}_FPF.pkl', 'wb') as f:  # open a text file
        pickle.dump(FPF, f) # serialize the list
    f.close()
    (mean, std) = norm.fit(FPF)
    with open(f'Q3_DataAnalysis/Data/reliability_analysis_Load{round(np.linalg.norm(Load), 0)}_n{n}_meanStd.pkl', 'wb') as f:  # open a text file
        pickle.dump([mean, std], f) # serialize the list
    f.close()
    x = np.linspace(mean - 8*std, mean + 8*std, 100)
    plt.hist(FPF)
    plt.show()
    plt.clf()
    plt.plot(x, norm.pdf(x, mean, std))
    plt.show()
    print("Failure Prob: ", norm.cdf(1, mean, std))


if __name__ == '__main__':
    np.random.seed(0)
    LayUp = np.array([0, 90, 45, -45])
    LayUp = np.append(LayUp, np.flip(LayUp))
    LayUp = np.append(LayUp, np.flip(LayUp))
    Load = (np.array([np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30)), 0, 0, 0, 0])*1200).T
    Load = Load*1000 # Convert to N/m
    t = 0.125e-3
    print(np.linalg.norm(Load))
    # Q3(LayUp, Load, n=102400)
    # Q3(LayUp, Load, n=51200)
    # Q3(LayUp, Load, n=25600)
    # Q3(LayUp, Load, n=12800)
    # Q3(LayUp, Load, n=6400)
    # Q3(LayUp, Load, n=3200)
    # Q3(LayUp, Load, n=1600)
    Q3(LayUp, Load, n=800)
    # Q3(LayUp, Load, n=400)
    # Q3(LayUp, Load, n=200)
    # Q3(LayUp, Load, n=100)