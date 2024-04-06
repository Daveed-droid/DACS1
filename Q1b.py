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
	strains = Laminate.calcPlyStrains(Load)
	stresses = Laminate.calcPlyStresses(Load)
	print(strains)
	print(stresses)
	z = Laminate.zlst
	data = [[0,z[0]]]
	for k in range(len(LayUp)):
		point1 = [stresses[0, k], z[k]]
		point2 = [stresses[0, k], z[k+1]]
		data.append(point1)
		data.append(point2)
	data.append([0,z[-1]])
	data = np.asarray(data)
	print(data)
	plt.plot(data[:,0]/1e6,data[:,1]*1e3, color="red")
	plt.plot([0,0],[z[0]*1e3, z[-1]*1e3], linestyle="dashed", color="black", linewidth=1)
	plt.plot([np.min(stresses)/1e6,np.max(stresses)/1e6],[0, 0], linestyle="--", color="black", linewidth=1)
	plt.xlabel("Stress Fibre Direction [MPa]")
	plt.ylabel("Position From Midline [mm]")
	plt.show()


Q1b()
