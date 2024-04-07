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

def plotStresses(Laminate, stresses, stressi):
	title = {
		0: r"Stress In Direction 1 $\sigma_{1}$",
		1: r"Stress In Direction 2 $\sigma_{2}$",
		2: r"Shear Stress In Direction 12 $\tau_{12}$"
	}
	z = Laminate.zlst
	data = [[0,z[0]]]
	for k in range(len(Laminate.LayUp)):
		point1 = [stresses[stressi, k], z[k]]
		point2 = [stresses[stressi, k], z[k+1]]
		data.append(point1)
		data.append(point2)
	data.append([0,z[-1]])
	data = np.asarray(data)
	plt.title(title[stressi])
	plt.plot([0,0],[z[0]*1e3, z[-1]*1e3], linestyle="dashed", color="black", linewidth=1)
	min, max = np.min(stresses[stressi,:])/1e6, np.max(stresses[stressi,:])/1e6
	rangei = max-min
	for i in range(len(z)):
		plt.plot([min-0.2*rangei, max+0.05*rangei],[z[i]*1e3, z[i]*1e3], linestyle="-", color="black", linewidth=0.5)
	for i in range(len(z)-1):
		plt.text(min-0.19*rangei, (z[i]*1e3+z[i+1]*1e3)/2, f"Ply: {Laminate.LayUp[i]}$^\circ$")
	plt.plot([min, max],[0, 0], linestyle="--", color="black", linewidth=1)
	plt.plot(data[:,0]/1e6,data[:,1]*1e3, color="red")
	plt.xlim([min-0.2*rangei, max+0.05*rangei])
	plt.xlabel("Stress [MPa]")
	plt.ylabel("Position From Midline [mm]")
	plt.show()
	plt.savefig(f"Plots/Q1b/stress{stressi}")

def plotStrains(Laminate, strains, straini):
	title = {
		0: r"Strain In Direction 1 $\epsilon_{1}$",
		1: r"Strain In Direction 2 $\epsilon_{2}$",
		2: r"Shear Strain In Direction 12 $\epsilon_{12}$"
	}
	z = Laminate.zlst
	data = [[0,z[0]]]
	for k in range(len(Laminate.LayUp)):
		point1 = [strains[straini, k], z[k]]
		point2 = [strains[straini, k], z[k+1]]
		data.append(point1)
		data.append(point2)
	data.append([0,z[-1]])
	data = np.asarray(data)
	plt.title(title[straini])
	plt.plot([0,0],[z[0]*1e3, z[-1]*1e3], linestyle="dashed", color="black", linewidth=1)
	min, max = np.min(strains[straini,:]), np.max(strains[straini,:])
	rangei = max-min
	for i in range(len(z)):
		plt.plot([min-0.2*rangei, max+0.05*rangei],[z[i]*1e3, z[i]*1e3], linestyle="-", color="black", linewidth=0.5)
	for i in range(len(z)-1):
		plt.text(min-0.19*rangei, (z[i]*1e3+z[i+1]*1e3)/2, f"Ply: {Laminate.LayUp[i]}$^\circ$")
	plt.plot([min, max],[0, 0], linestyle="--", color="black", linewidth=1)
	plt.plot(data[:,0],data[:,1]*1e3, color="blue")
	plt.xlim([min-0.2*rangei, max+0.05*rangei])
	plt.xlabel("Strains [-]")
	plt.ylabel("Position From Midline [mm]")
	plt.show()
	plt.savefig(f"Plots/Q1b/strain{straini}")


def Q1b(Lamina_mean = Lamina_mean, Laminate = Laminate):
	LayUp = [0, 0, 90, 30, 90]
	Laminate = Laminate(LayUp, Lamina_mean)
	Load = np.array([0.2e2, 1.8e4, 0, 18e3, 0, 0]).T
	strains = Laminate.calcPlyStrains(Load)
	stresses = Laminate.calcPlyStresses(Load)
	print(strains)
	print(stresses)
	stresses2 = Laminate.calcGloStresses(Load)
	print(stresses2)
	stresses2 = Laminate.calcGloStrains(Load)
	print(stresses2)
	plotStresses(Laminate, stresses, 0)
	plotStresses(Laminate, stresses, 1)
	plotStresses(Laminate, stresses, 2)
	plotStrains(Laminate, strains, 0)
	plotStrains(Laminate, strains, 1)
	plotStrains(Laminate, strains, 2)

Q1b()
