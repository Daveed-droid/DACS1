import matplotlib.pyplot as plt
import numpy as np

from Assignment2.StiffenerClass import Stiffener
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
	def __init__(self, Material, ratio = [1, 1, 1], Stiffeners: bool or list = False, Metal = False, dTheta = 20, rho = 0):
		dia = 6
		self.dia = dia
		self.Metal = Metal
		self.Stiffeners = Stiffeners
		self.dTheta = dTheta
		self.ratio = ratio
		# Create fuselage nodes
		Theta = np.arange(-90, 271, dTheta)
		self.nNodes = len(Theta)
		self.nElem = self.nNodes-1
		self.x = -dia/2 * np.cos(np.deg2rad(Theta))
		self.y = dia/2 * np.sin(np.deg2rad(Theta))
		self.element_pos = []
		self.element_ang = []
		self.thickness = []
		self.stiffenerElem = []
		self.Material = Material
		if self.Stiffeners:

			self.stiffener_pos = []
			for i in range(len(Stiffeners)):
				x_stiff = -dia/2 * np.cos(np.deg2rad(self.Stiffeners[i][0]))
				y_stiff = dia/2 * np.sin(np.deg2rad(self.Stiffeners[i][0]))
				stiff_L = self.Stiffeners[i][1].set_angle(-self.Stiffeners[i][0])
				stiff_R = self.Stiffeners[i][1].set_angle(self.Stiffeners[i][0])
				self.stiffenerElem.append(stiff_L)
				self.stiffenerElem.append(stiff_R)
				stiffner_pos_left = [x_stiff+stiff_L.x_cg, y_stiff+stiff_L.y_cg]
				stiffner_pos_right = [-x_stiff+stiff_R.x_cg, y_stiff+stiff_R.y_cg]
				self.stiffener_pos.append(stiffner_pos_right)
				self.stiffener_pos.append(stiffner_pos_left)
			self.stiffener_pos = np.asarray(self.stiffener_pos)


		for i in range((self.nNodes-1)):
			self.element_pos.append([self.x[i], self.y[i], self.x[i+1], self.y[i+1]])
			self.element_ang.append(-(90+((Theta[i]+Theta[i+1])/2)))
		self.element_pos = np.asarray(self.element_pos, dtype = float)
		for i in range(int((self.nNodes-1)/2)):
			if i/self.nNodes < self.ratio[0] /np.sum(self.ratio)/2:
				#print(i, self.ratio,self.nNodes,self.thickness,self.Material[0].h)
				self.thickness.append(self.Material[0].h)
			elif len(ratio) > 1 and self.ratio[0] /np.sum(self.ratio)/2  <= i/self.nNodes <  (self.ratio[1]+self.ratio[0]) /np.sum(self.ratio)/2:
				self.thickness.append(self.Material[1].h)
			elif len(ratio) > 2:
				self.thickness.append(self.Material[2].h)
		self.thickness = np.append(self.thickness, np.flip(self.thickness))



		A = 0
		for i in range(len(self.ratio)):
			t = self.Material[i].h
			A += np.pi*((dia/2 + t/2)**2-(dia/2 - t/2)**2)*self.ratio[i]/sum(self.ratio)
			if Stiffeners:
				for i in range(len(self.stiffenerElem)):
					A += self.stiffenerElem[i].A
		self.mass = A*rho
		if type(Material) == list:
			self.laminates = []
			assert self.nElem%2 == 0  # Must be even
			assert (self.nElem//2)%sum(ratio) == 0  # Must be divisible by the sum of the ratio
			halfnElem = self.nElem//2
			perRatio = halfnElem//sum(ratio)
			for i in range(len(self.Material)):
				self.laminates += [self.Material[i]]*(perRatio*ratio[i])
			self.laminates = self.laminates + list(reversed(self.laminates))
			if Stiffeners:
				self.bend_stiff = np.zeros(self.nElem+len(self.stiffenerElem))
				self.shear_stiff = np.zeros(self.nElem+len(self.stiffenerElem))
			else:
				self.bend_stiff = np.zeros(self.nElem)
				self.shear_stiff = np.zeros(self.nElem)

			axial_stiff = 0
			temp = 0
			# Finding neutral axis and setting shear stiffness
			for i in range(self.nElem):
				x1, y1, x2, y2 = self.element_pos[i, :]
				Ex, Ey, vxy, vyx, Gxy = self.laminates[i].calcEngConst()
				ang = self.element_ang[i]
				h = self.laminates[i].h
				l = np.linalg.norm(np.array([x2, y2]).T - np.array([x1, y1]).T)
				self.shear_stiff[i] = Gxy*h*l*abs(np.sin(np.deg2rad(ang)))
				# Find neutral axis
				axial_stiff += Ex*h*l
				temp += Ex*h*l*((y1+y2)/2)
			if Stiffeners:
				for j in range(len(self.stiffenerElem)):
					index = i+j
					stiff = self.stiffenerElem[j]
					x, y = self.stiffener_pos[j]
					axial_stiff += stiff.EA
					temp += stiff.EA*(y)
					self.shear_stiff[index] = stiff.GA


			self.y_neutral_axis = temp/axial_stiff
			self.Inertia = np.zeros(self.nElem)
			# Setting bending stiffness
			for i in range(self.nElem):
				x1, y1, x2, y2 = self.element_pos[i, :]
				Ex, Ey, vxy, vyx, Gxy = self.laminates[i].calcEngConst()
				h = self.laminates[i].h
				l = np.linalg.norm(np.array([x2, y2]).T - np.array([x1, y1]).T)
				# Adding steiner term
				self.bend_stiff[i] = Ex*h*l*(((y1+y2)/2)-self.y_neutral_axis)**2
				self.Inertia[i] = self.bend_stiff[i] / Ex

			if Stiffeners:
				for j in range(len(self.stiffenerElem)):
					index = i+j
					stiff = self.stiffenerElem[j]
					x, y = self.stiffener_pos[j]
					self.bend_stiff[index] = stiff.EA * (y - self.y_neutral_axis) ** 2

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
		t = self.thickness
		# Material Failure
		EI = np.sum(self.bend_stiff)
		I = np.sum(self.Inertia)
		GA = np.sum(self.shear_stiff)
		y1 = self.element_pos[:, 1]
		y2 = self.element_pos[:, 3]
		y_avg = (y1+y2)/2
		ex = moment*(y_avg-self.y_neutral_axis)/EI
		exy = shear*np.sin(np.deg2rad(np.asarray(self.element_ang)))/GA
		LamStrains = np.zeros((3, self.nElem))
		LamStrains[0, :] = ex
		LamStrains[2, :] = exy
		Stress = np.zeros([3,self.nElem])
		Stress[0,:] = ex *EI/I
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
		if self.Stiffeners:
			stif_pos = np.asarray(self.stiffener_pos)
			y = stif_pos[:, 1]
			ex = moment*(y-self.y_neutral_axis)/EI
			exy = shear/GA
			self.FailedStiffners = []
			for i in range(len(self.stiffenerElem)):
				LamStrains = np.zeros((3, 1))
				LamStrains[0, :] = ex[i]
				LamStrains[2, :] = exy
				FailLam = self.stiffenerElem[i].get_fpf(LamStrains)
				self.FailedStiffners.append(max(FailLam))
		for i in range(self.nElem):

			#	Buckling with stiffners
			if self.Stiffeners == True:
				#Check for stiffener crippling
				b = 2 * np.pi / self.nNodes * self.dia / 2  # Widt of element
				As = [] # Insert Aerea of the stiffeners!!!
				if Stress[0,i] >= StiffenerCrippling()[1]/As[i]:
					print("Failure Stiffener Crippling", i, Stress[0,i]/ (StiffenerCrippling()[1]/(t[i]*b)))
				#Check for skin buckling due to compression
				if Stress[0,i] >= SkinBuckling()[0]/(t[i]*b):
					print("Failure Skin Buckling Composite stiffener", i, Stress[0,i]/ (SkinBuckling()[0]/(t[i]*b)))
				# Check for skin buckling due to shear
				if ShearFlow() >= SkinBuckling()[1]:
					print("Failure Skin due to shear Composite stiffener", i, ShearFlow()[i] / SkinBuckling()[1])
			else:
				b = self.dia*np.pi*3/16 #Length of element without stiffeners
				#Check for skin buckling due to compression
				if Stress[0,i] >= self.SkinBuckling()[0]/(t[i]*b):
					print("Failure Skin Buckling", i, Stress[0,i], (self.SkinBuckling()[0]/(t[i]*b)) )

				# Check for skin buckling due to shear
				if self.ShearFlow()[i] >= self.SkinBuckling()[1]/t[i]:
					print("Failure Skin due to shear", i, self.ShearFlow()[i] / (self.SkinBuckling()[1]/t[i]))



		if plot_failure:
			self.PlotNodes(Failed)
		return Failed

	def ShearFlow(self):
		sy = 1.5*10**6
		if self.Stiffeners == True:

			b = 2*np.pi/self.nNodes*D/2
		else:
			b = np.pi*self.dia*3/16
		# Structural Idealization
		angle = np.arange(0,360,self.dTheta)
		t = self.thickness
		B = []
		qb = np.zeros(self.nElem)
		qs = np.zeros(self.nElem)
		qm = np.zeros(self.nElem)
		As = [0,0,0] # Stiffener area
		if self.Metal == True:

			t = self.Material[0].h
			for i in range(self.nElem):
				B.append(t*b*(2+self.y[i-1]/self.y[i])/6 + t*b*(2+self.y[i+1]/self.y[i])/6)
		elif self.Stiffeners == False:

			for i in range(self.nElem):
				B.append(t[i]*b*(2+self.y[i-1]/self.y[i])/6 + t[i]*b*(2+self.y[i+1]/self.y[i])/6)

		else:

			for i in range(self.nElem):
				B.append(As[i] + t[i]*b*(2+self.y[i-1]/self.y[i])/6 + t[i]*b*(2+self.y[i+1]/self.y[i])/6)


		# print(B)
		D = self.dia
		if self.Metal == True:
			Ixx = (D**4-(D-2*t)**4)*3.1415/64
		elif self.Stiffeners == False:
			Ixx = np.sum(self.y[0:-1]**2*(self.thickness*b))

		else:
			Ixx = np.sum(self.y[0:-1] ** 2*(self.thickness*(b)+As))
			pass

		for i in range(self.nElem):
			qb[i] = B[i]+self.y[i]*-sy*Ixx+qb[i-1]

			qm[i] = -qb[i]*(self.x[i+1]-self.x[i]*self.y[i]+qb[i]*(self.y[i+1]-self.y[i])*self.x[i])
			#print(qb[i],qm[i],self.x[i],self.y[i])

		qs0 = 0	#np.sum(qm) / (2*np.pi*(self.dia/2)**2)
		qs = qb + qs0

		#print(np.sum(qm),qs0, qb,self.y[0])




		return qs
	def SkinBuckling(self):
		#Compression
		a = 1 #Length
		if self.Stiffeners == True:

			b = 2*np.pi/self.nNodes*D/2
		else:
			b = np.pi*self.dia*3/16
		AR = a/b
		m = 2 #???
		if self.Metal == True:

			D11 = self.Material[0].ABD[3,3]
			D66 = self.Material[0].ABD[5,5]
			D12 = self.Material[0].ABD[3,4]
			D22 = self.Material[0].ABD[4,4]
		else:
			D11 = self.Material[0].ABD[3, 3]
			D66 = self.Material[0].ABD[5, 5]
			D12 = self.Material[0].ABD[3, 4]
			D22 = self.Material[0].ABD[4, 4]


		NCom = (np.pi**2*(D11*m**4+2*(D12+2*D66)*m**2*AR**2+D22*AR**4))/(a**2*m**2)

		#Shear
		if self.Metal == True:

			D11 = self.Material[0].ABD[3,3]
			D66 = self.Material[0].ABD[5,5]
			D12 = self.Material[0].ABD[3,4]
			D22 = self.Material[0].ABD[4,4]
		else:
			D11 = self.Material[1].ABD[3, 3]
			D66 = self.Material[1].ABD[5, 5]
			D12 = self.Material[1].ABD[3, 4]
			D22 = self.Material[1].ABD[4, 4]
		A = -0.27 + 0.185*(D12+2*D66)/(D11*D22)**0.5
		B = 0.82+0.46*(D12+2*D66)/(D11*D22)-0.2*((D12+2*D66)/(D11*D66)**0.5)**2
		beta = (D11/D22)**0.25
		K = 8.2 + 5*(D12+2*D66)/((D11*D22)**0.5*(A/beta+B*beta))
		NShear = 4*(D11*D22**3)**0.25*K/b**2
		return NCom, NShear

	def StiffenerCrippling(self, Stiff):
		b = Stiff.buckProp[1]
		StiffCond = Stiff.buckProp[0]
		LamStiffElem = Stiff.buckProp[2]
		D66 = LamStiffElem.ABD[5,5]

		if StiffCond == "OEF":	#OEF
			Sig_OEFr = 1.63/((b/t)**0.717)	#Ratio of crippling strength to compressive strength of stiffner
			Sig_stif = Sig_OEFr * XcStif
			Nxstif = 12*D66/b**2
		elif StiffCond == "NEF":	#NEF
			Sig_stifr = 11/((b/t)**1.124)	#Ratio of crippling strength to compressive strength of stiffner
			Sig_stif = Sig_NEFr * XCStif
			Nxstif = 2*3.1415**2/b**2*((D11*D22)**0.5+D12+2*D66)
		else:
			raise NotImplementedError

		return Sig_stif, Nxstif

if __name__ == "__main__":
	# CompLam = [0, 0, 0]
	# ShearLam = [45, -45, 45, -45]
	# TensionLam = [0, 0, 0]
	# Lam1 = Laminate(CompLam, AssignmentLamina)
	# Lam2 = Laminate(ShearLam, AssignmentLamina)
	# Lam3 = Laminate(TensionLam, AssignmentLamina)
	# # nelem should be divisible by 2 and divisible by sum of ratio
	# nelem = 30
	# a = Fuselage([Lam1, Lam2, Lam3], ratio = [1,1,1], Stiffeners = False, dTheta = 360//nelem, rho = 1610)
	# a.PlotNodes()
	# a.Load(15e6, 1.5e6, plot_failure = True)
	# print(a.mass)

	nelem = 30
	t = 1.43e-3
	AssignmentMetalLamina = Lamina(t, 69e9, 69e9, 0.29, 26e9)
	AssignmentMetalLamina.setStrengths(410e6, 400e6, 430e6, 430e6, 230e6)
	Lam = Laminate([0], AssignmentMetalLamina)
	a = Fuselage([Lam], [1], Stiffeners = False, dTheta = 360 // nelem, Metal = True, rho = 2770)
	a.PlotNodes()
	a.Load(15e6, 1.5e6, plot_failure = True)
	print(a.mass)

	CompLam = [0, 0, 0]
	ShearLam = [45, -45, 45, -45]
	TensionLam = [0, 0, 0]
	Lam1 = Laminate(CompLam, AssignmentLamina)
	Lam2 = Laminate(ShearLam, AssignmentLamina)
	Lam3 = Laminate(TensionLam, AssignmentLamina)
	LamF = Laminate(CompLam, AssignmentLamina)
	LamW = Laminate(ShearLam, AssignmentLamina)
	StringH = Stiffener(LamW, LamF, "hat")
	StringL = Stiffener(LamW, LamF, "L")
	StringI = Stiffener(LamW, LamF, "I")

	# nelem should be divisible by 2 and divisible by sum of ratio
	nelem = 30
	A = nelem//6
	B = [[(360//nelem)*i, StringH] for i in range(A)]
	i = len(B)
	C = [[(360//nelem)*(j+i), StringI] for j in range(1,A)]
	j = len(C)
	D = [[(360//nelem)*(k+j), StringL] for k in range(1,A)]
	E = B + C + D
	print(E)
	a = Fuselage([Lam1, Lam2, Lam3], ratio = [1,1,1], Stiffeners = E, dTheta = 360//nelem, rho = 1610)
	a.PlotNodes()
	a.Load(15e6, 1.5e6, plot_failure = True)
	print(a.mass)