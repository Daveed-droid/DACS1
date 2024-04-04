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


def Q3(LayUp, Load, n = 1000):
	LPF = []
	for n in range(n):
		# Random pick
		props = []
		for [mean, std] in Lamina_props:
			props.append(np.random.normal(mean, std))
		E1, E2, v12, G12 = props[0:4]
		Lamina_ = Lamina(t, E1, E2, v12, G12)
		Xt, Yt, Xc, Yc, S = props[4:9]
		Lamina_.setStrengths(Xt, Yt, Xc, Yc, S)
		Laminate_ = Laminate(LayUp, Lamina_)
		LoadFPF, LoadLPF = Laminate_.calcFailure(Load)
		LoadNorm = np.linalg.norm(Load)
		LoadFPFNorm = np.linalg.norm(LoadFPF)
		LPF.append(LoadFPFNorm/LoadNorm)
	plt.hist(LPF)


if __name__ == '__main__':
	np.random.seed(0)
	LayUp = np.array([0, 90, 45, -45])
	LayUp = np.append(LayUp, np.flip(LayUp))
	LayUp = np.append(LayUp, np.flip(LayUp))
	Load = np.array([np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30)), 0, 0, 0, 0]).T
	t = 0.125
	Q3(LayUp, Load)