import matplotlib.pyplot as plt
import numpy as np
from LaminateClass import Laminate
from LaminaClass import Lamina
from Assignment2.Assignment2Data import Metal, AssignmentLamina, AssignmentMetal

def Rx(theta):
	R = np.zeros((3,3))
	R[0, 0] = 1
	R[1, 1] = np.cos(theta)
	R[2, 2] = np.cos(theta)
	R[2, 1] = np.sin(theta)
	R[1, 2] = -np.sin(theta)
	return R


class Fuselage:
	def __init__(self, Material, ratio = [1, 1, 1], Stiffeners: bool = False, Metal = False, dTheta = 20):
		dia = 6
		self.Metal = Metal
		# Create fuselage nodes
		Theta = np.arange(-90, 271, dTheta)
		self.nNodes = len(Theta)
		self.nElem = self.nNodes-1
		self.x = -dia/2 * np.cos(np.deg2rad(Theta))
		self.y = dia/2 * np.sin(np.deg2rad(Theta))
		self.element_pos = []
		self.element_ang = []
		for i in range(self.nNodes-1):
			self.element_pos.append([self.x[i], self.y[i], self.x[i+1], self.y[i+1]])
			self.element_ang.append(-(90+((Theta[i]+Theta[i+1])/2)))
		self.element_pos = np.asarray(self.element_pos, dtype = float)
		self.Material =  Material
		self.ratio = ratio
		if type(Material) == list:
			self.laminates = []
			assert self.nElem%2 == 0 # Must be even
			assert (self.nElem//2)%sum(ratio) == 0 # Must be divisible by the sum of the ratio
			halfnElem = self.nElem//2
			perRatio = halfnElem//sum(ratio)
			for i in range(len(self.Material)):
				self.laminates+=[self.Material[i]]*(perRatio*ratio[i])
			self.laminates = self.laminates + list(reversed(self.laminates))
			# Calc ABD
			self.ABD = np.zeros((6,6), dtype = float)
			self.bend_stiff = np.zeros(self.nElem)
			self.shear_stiff = np.zeros(self.nElem)
			axial_stiff = 0
			temp = 0
			for i in range(self.nElem):
				x1, y1, x2, y2 = self.element_pos[i, :]
				ang = self.element_ang[i]
				AMat = self.laminates[i].AMatrix
				h = self.laminates[i].h
				l = np.linalg.norm(np.array([x2, y2]).T - np.array([x1, y1]).T)
				A11 = AMat[0, 0]
				A33 = AMat[2, 2]
				# Adding steiner term
				self.bend_stiff[i] = A11*l*((y1+y2)/2)**2
				self.shear_stiff[i] = A33*l*abs(np.sin(np.deg2rad(ang)))
				# Find neutral axis
				axial_stiff += A11
				temp += A11*((y1+y2)/2)
			self.y_neutral_axis = temp/axial_stiff

		else:
			print("Invalid")



	def PlotNodes(self, Failed=None) -> None:
		fig, ax = plt.subplots(figsize = (6, 4), layout = 'constrained')
		ax.plot([min(self.x), max(self.x)], [0, 0], color="black", linewidth=1, linestyle="-.", alpha=0.5)
		ax.plot([0, 0], [min(self.y), max(self.y)], color="black", linewidth=1, linestyle="-.", alpha=0.5)
		cmap = plt.colormaps['cool']
		# Take colors at regular intervals spanning the colormap.
		colors = cmap(np.linspace(0, 1, len(self.Material)))
		temp = [int((self.nElem/2)*self.ratio[i]/sum(self.ratio)) for i in range(len(self.ratio))]
		temp = temp + list(reversed(temp))
		start = 0
		end = temp[0]
		imid = len(temp)//2
		start2 = sum(temp[0:imid])
		end2 = start2 + temp[imid]
		for i in range(len(self.Material)):
			ax.plot(self.x[0], self.y[0], color = colors[i], label=f"{self.Material[i].LayUp}") # Set legend
			ax.plot(self.x[start:end+1], self.y[start:end+1], marker = "X", color = colors[i])
			ax.plot(self.x[start2:end2+1], self.y[start2:end2+1], marker = "X", color = colors[(len(self.Material)-1)-i])
			if i == len(self.Material)-1:
				break
			i2 = i + imid
			start = end
			end = start+temp[i+1]
			start2 = end2
			end2 = start2+temp[i2+1]

		ax.plot([min(self.x), max(self.x)], [self.y_neutral_axis, self.y_neutral_axis], color="blue", linewidth=1, linestyle="-.", alpha=0.5)
		plt.text((min(self.x)+max(self.x))/2, self.y_neutral_axis + 0.2, f"N.A.", fontsize=8, color="blue", horizontalalignment='center', verticalalignment='center')
		for i in range(len(self.x)-1):
			plt.text((self.x[i]+self.x[i+1])*1.05/2, (self.y[i]+self.y[i+1])*1.05/2, f"{i+1}", fontsize=8, horizontalalignment='center', verticalalignment='center')
		if type(Failed)!=type(None):
			for i in range(len(self.x)-1):
				plt.text((self.x[i]+self.x[i+1])*0.90/2, (self.y[i]+self.y[i+1])*0.90/2, f"{round(Failed[i],2)}", fontsize=8, horizontalalignment='center', verticalalignment='center')

		xmin, xmax = ax.get_xlim()
		ax.set_xlim([xmin, xmax*2])
		plt.legend()

		plt.tight_layout()
		plt.show()

	def Load(self, moment: float, shear: float, plot_failure=False):
		# Global buckling
		# Material Failure
		EI = np.sum(self.bend_stiff)
		GA = np.sum(self.shear_stiff)
		y1 = self.element_pos[:, 1]
		y2 = self.element_pos[:, 3]
		y_avg = (y1+y2)/2
		ex = moment*(y_avg-self.y_neutral_axis)/EI
		exy = shear*np.sin(np.deg2rad(np.asarray(self.element_ang)))/GA
		LamStrains = np.zeros((3, self.nElem))
		LamStrains[0, :] = ex
		LamStrains[2, :] = exy
		Failed = np.zeros(self.nElem)
		for i in range(self.nElem):
			lam = self.laminates[i]
			if not self.Metal:
				f_FFp, f_IFFp = lam.PuckStrain(LamStrains[:, i])
				a = np.max(f_FFp)
				b = np.max(f_IFFp)
				d = max(a, b)
			elif self.Metal:
				f_xm, f_ym, f_sm = lam.MaxStrain(LamStrains[:, i])
				a = np.max(f_xm)
				b = np.max(f_ym)
				c = np.max(f_sm)
				temp = max(a, b)
				d = max(c, temp)
			Failed[i] = d
		if plot_failure:
			self.PlotNodes(Failed)
		return Failed


if __name__ == "__main__":
	CompLam = [0, 90, 90, 0]
	ShearLam = [45, -45, -45, 45]
	TensionLam = [0, 0, 0, 0]
	Lam1 = Laminate(CompLam, AssignmentLamina)
	Lam2 = Laminate(ShearLam, AssignmentLamina)
	Lam3 = Laminate(TensionLam, AssignmentLamina)
	# nelem should be divisible by 2 and divisible by sum of ratio
	nelem = 30
	a = Fuselage([Lam1, Lam2, Lam3], ratio = [1,1,1], Stiffeners = False, dTheta = 360//nelem)
	a.PlotNodes()
	Failed = a.Load(15e6, 1.5e6, plot_failure = True)
	print(Failed)

	t = 1.5e-3
	AssignmentMetalLamina = Lamina(t, 69e9, 69e9, 0.29, 26e9)
	AssignmentMetalLamina.setStrengths(410e6, 400e6, 430e6, 430e6, 230e6)
	Lam = Laminate([0], AssignmentMetalLamina)
	a = Fuselage([Lam], [1], Stiffeners = False, dTheta = 360 // nelem, Metal = True)
	a.PlotNodes()
	Failed = a.Load(15e6, 1.5e6, plot_failure = True)
	print(Failed)
