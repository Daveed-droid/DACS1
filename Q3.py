"""
@Project ：DACS1
@File ：Q3.py
@Author ：David Canosa Ybarra
@Date ：19/03/2024 11:46
"""
import numpy as np
from AssignmentData import Lamina_props
from LaminaClass import Lamina
from LaminateClass import Laminate
import matplotlib.pyplot as plt
from tqdm import tqdm


def Q3(LayUp, Load, n = 10000):
	LPF = []
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
		LoadFPF, LoadLPF = Laminate_.calcFailure(Load, dL_step = 20000)
		LoadNorm = np.linalg.norm(Load)
		LoadFPFNorm = np.linalg.norm(LoadFPF)
		LPF.append(LoadFPFNorm/LoadNorm)
	plt.hist(LPF)
	plt.show()


if __name__ == '__main__':
	np.random.seed(0)
	LayUp = np.array([0, 90, 45, -45])
	LayUp = np.append(LayUp, np.flip(LayUp))
	LayUp = np.append(LayUp, np.flip(LayUp))
	Load = (np.array([np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30)), 0, 0, 0, 0])*850).T
	Load = Load*1000 # Convert to N/m
	t = 0.125e-3
	Q3(LayUp, Load)