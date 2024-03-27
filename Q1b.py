"""
@Project ：DACS1
@File ：Q1b.py
@Author ：David Canosa Ybarra
@Date ：19/03/2024 11:45
"""
from AssignmentData import Lamina_mean
from LaminateClass import Laminate
import numpy as np
from matplotlib import pyplot as plt

def Q1b(Lamina_mean = Lamina_mean, Laminate = Laminate):
	LayUp = [0, 0, 90, 30, 90]
	Laminate = Laminate(LayUp, Lamina_mean)
	Load = np.array([0.2e2, 1.8e4, 0, 18e3, 0, 0]).T
	strains = Laminate.calcGloStrains(Load)
	stresses = Laminate.calcGloStresses(Load)
	print(strains)
	print(stresses)
	z = Laminate.zlst
	data = [[Laminate.zlst[0],stresses[0,1]]]+[[z, st] for z, st in zip(Laminate.zlst, stresses[0,:])]+[[Laminate.zlst[-1],stresses[0,-1]]]
	plt.plot(data)
	plt.show()


Q1b()
