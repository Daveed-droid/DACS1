import matplotlib.pyplot as plt
import numpy as np
from LaminateClass import Laminate
from LaminaClass import Lamina
from Assignment2.Assignment2Data import Metal, AssignmentLamina, AssignmentMetal



class Fuselage:
	def __init__(self, Material, Stiffeners: bool, dTheta = 20):
		dia = 6
		# Create fuselage nodes
		Theta = np.arange(-90, 271, dTheta)
		self.x = -dia/2 * np.cos(np.deg2rad(Theta))
		self.y = dia/2 * np.sin(np.deg2rad(Theta))
		self.element_pos = []
		self.element_ang = []
		for i in range(len(Theta)-1):
			self.element_pos.append([self.x[i], self.y[i], self.x[i+1], self.y[i+1]])
			self.element_ang.append(-(90+((Theta[i]+Theta[i+1])/2)))
		self.element_pos = np.asarray(self.element_pos, dtype = float)
		if type(Material) == Laminate:
			pass
		elif type(Material) == Metal:
			pass
		else:
			pass



	def PlotNodes(self) -> None:
		fig, ax = plt.subplots(figsize = (4, 4), layout = 'constrained')
		ax.plot([min(self.x), max(self.x)], [0, 0], color="black", linewidth=1, linestyle="-.", alpha=0.5)
		ax.plot([0, 0], [min(self.y), max(self.y)], color="black", linewidth=1, linestyle="-.", alpha=0.5)
		ax.plot(self.x, self.y, marker = "X", color = "red")
		for i in range(len(self.x)-1):
			plt.text((self.x[i]+self.x[i+1])*1.05/2, (self.y[i]+self.y[i+1])*1.05/2, f"{i+1}", fontsize=8, horizontalalignment='center', verticalalignment='center')
		plt.show()

	def Load(self, moment: float, shear: float) -> bool:
		# Global buckling
		# Yield
		
		return Failed


if __name__ == "__main__":
	Lam = Laminate([0,0,0,90], AssignmentLamina)
	a = Fuselage(Lam, False)
	a.PlotNodes()
