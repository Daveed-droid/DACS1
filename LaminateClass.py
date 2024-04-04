"""
@Project ：DACS1
@File ：LaminateClass.py
@Author ：David Canosa Ybarra
@Date ：19/03/2024 11:44
"""
from LaminaClass import Lamina
import numpy as np

def StressGloTOPly(theta):
	"""
	Calculate the rotation matrix to go from stresses in the laminate ref frame to stresses in the lamina ref frame
	:param theta: Angle of the ply in Radians
	:return: Rotation matrix
	"""
	n = np.sin(np.deg2rad(theta))
	m = np.cos(np.deg2rad(theta))
	RotMat = np.array([[m**2, n**2, 2*m*n],
					   [n**2, m**2, -2*m*n],
					   [-m*n, m*n, m**2-n**2]])
	return RotMat

def StrainGloTOPly(theta):
	"""
	Calculate the rotation matrix to go from strains in the laminate ref frame to strains in the lamina ref frame
	:param theta: Angle of the ply in Radians
	:return: Rotation matrix
	"""
	n = np.sin(np.deg2rad(theta))
	m = np.cos(np.deg2rad(theta))
	RotMat = np.array([[m**2, n**2, m*n],
					   [n**2, m**2, -m*n],
					   [-2*m*n, 2*m*n, m**2-n**2]])
	return RotMat


class Laminate():
	"""
	The laminate class will be used to do operations on the layup
	"""
	def __init__(self, LayUp: list, Lamina: Lamina | list):
		"""
		Initializes the laminate class, this class takes in a layup and lamina
		:param LayUp: A list of angles [deg] starting from the bottom ply
		:param Lamina: A lamina object
		:param t: The lamina thickness
		"""
		self.LayUp = LayUp
		if type(Lamina) == list:
			self.Lamina = Lamina
			constantUD = False
		elif type(Lamina) == Lamina:
			self.Lamina = [Lamina]*len(LayUp)
			constantUD = True
		self.t = self.Lamina.t
		self.h = self.t * len(LayUp)
		self.zlst = np.linspace(-self.h / 2, self.h / 2, len(self.LayUp) + 1, endpoint = True)
		self.calcQGlobalLaminas()
		self.calcABD()

	def calcQGlobalLaminas(self):
		"""
		Calculates the Q matrix of all the laminas when placed at an angle
		:return: None
		"""
		self.QGlobalAr = list(range(len(self.LayUp)))

		for i, theta in enumerate(self.LayUp):
			Q11 = self.Lamina[i].QMatrix[0, 0]
			Q12 = self.Lamina[i].QMatrix[0, 1]
			Q22 = self.Lamina[i].QMatrix[1, 1]
			Q66 = self.Lamina[i].QMatrix[2, 2]
			n = np.sin(np.deg2rad(theta))
			m = np.cos(np.deg2rad(theta))
			Qxx = Q11 * m ** 4 + 2 * (Q12 + 2 * Q66) * m ** 2 * n ** 2 + Q22 * n ** 4
			Qxy = (Q11 + Q22 - 4 * Q66) * m ** 2 * n ** 2 + Q12 * (m ** 4 + n ** 4)
			Qyy = Q11 * n ** 4 + 2 * (Q12 + 2 * Q66) * m ** 2 * n ** 2 + Q22 * m ** 4
			Qxs = (Q11 - Q12 - 2 * Q66) * n * m ** 3 + (Q12 - Q22 + 2 * Q66) * n ** 3 * m
			Qys = (Q11 - Q12 - 2 * Q66) * m * n ** 3 + (Q12 - Q22 + 2 * Q66) * m ** 3 * n
			Qss = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * n ** 2 * m ** 2 + Q66 * (m ** 4 + n ** 4)
			QMatrix_glo = np.array([[Qxx, Qxy, Qxs],
									[Qxy, Qyy, Qys],
									[Qxs, Qys, Qss]])
			self.QGlobalAr[i] = QMatrix_glo
	def calcABD(self):
		"""
		Calculates the ABD matrix of the laminate
		:return:
		"""
		self.AMatrix = np.zeros((3, 3))
		self.BMatrix = np.zeros((3, 3))
		self.DMatrix = np.zeros((3, 3))
		# i and j are matrix components, k is the ply in the layup
		for k in range(len(self.LayUp)):
			Q = self.QGlobalAr[k]
			# calculating the A Matrix
			self.AMatrix += Q * (self.zlst[k + 1] - self.zlst[k])
			# calculating the B Matrix
			self.BMatrix += 0.5 * Q * (self.zlst[k + 1] ** 2 - self.zlst[k] ** 2)
			# calculating the D Matrix
			self.DMatrix += 3 ** -1 * Q * (self.zlst[k + 1] ** 3 - self.zlst[k] ** 3)
		# place the matrices together into the ABD matrix
		ABD_top = np.hstack((self.AMatrix, self.BMatrix))
		ABD_bottom = np.hstack((self.BMatrix, self.DMatrix))
		self.ABD = np.vstack((ABD_top, ABD_bottom))

	def calcGloStrainsNoCurve(self, Load):
		"""
		Calculates the strain and curvature of the laminate at a prescribed load
		:param Load: A 6x1 numpy array column vector with the loads applied [Nx, Ny, Nz, Mx, My, Mz]
		:return: The deflections of the laminate
		"""
		return np.linalg.inv(self.ABD) @ Load

	def calcPlyStrains(self, Load):
		GloStrains = self.calcGloStrains(Load)
		PlyStrains = np.zeros((3, len(self.LayUp)))
		for k, theta in enumerate(self.LayUp):
			theta = np.deg2rad(theta)
			PlyStrains[:, k] = StrainGloTOPly(theta)@GloStrains[:, k]
		return PlyStrains

	def calcGloStrains(self, Load):
		FlatStrains = np.linalg.inv(self.ABD)@Load
		GloStrains = np.zeros((3, len(self.LayUp)))
		zAvg = [self.zlst[k + 1] + self.zlst[k] for k in range(len(self.zlst) - 1)]
		for k in range(len(self.LayUp)):
			GloStrains[:, k] = (FlatStrains[0:3] + zAvg[k] * FlatStrains[3:6]).T
		return GloStrains



	def calcPlyStresses(self, Load):
		GloStresses = self.calcGloStresses(Load)
		PlyStresses = np.zeros((3, len(self.LayUp)))
		for k, theta in enumerate(self.LayUp):
			PlyStresses[:, k] = StressGloTOPly(theta)@GloStresses[:, k]
		return PlyStresses

	def calcGloStresses(self, Load):
		plyStrains = self.calcPlyStrains(Load)
		plyStresses = np.zeros((3, len(self.LayUp)))
		for k in range(len(self.LayUp)):
			plyStresses[:, k] = self.QGlobalAr[k] @ plyStrains[:, k]
		return plyStresses
	def calcEngConst(self):
		"""
		Calculates the engineering constants of the laminate
		:return: [Ex, Ey, vxy, vyx, Gxy]
		"""
		Axx = self.ABD[0, 0]
		Ayy = self.ABD[1, 1]
		Axy = self.ABD[0, 1]
		Ass = self.ABD[2, 2]

		Ex = (Axx * Ayy - Axy ** 2) / (self.h * Ayy)
		Ey = (Axx * Ayy - Axy ** 2) / (self.h * Axx)
		vxy = Axy / Ayy
		vyx = Axy / Axx
		Gxy = Ass / self.h
		return [Ex, Ey, vxy, vyx, Gxy]

	def calcStressEnvelope(self):
		pass

	def Puck(self, Load, Strength):		# Strength is list of ply properties: [Xt_mean,Xc_mean,Yt_mean,Yc_mean,S_mean]

		pt12 = 0.3
		pc12 = 0.25
		pc11 = 0.225
		Xt = Strength[0]
		Xc = Strength[1]
		Yt = Strength[2]
		Yc = Strength[3]
		S = Strength[4]
		Force = np.array([Load[0], Load[1], Load[2], 0, 0, 0]).T
		stresses = self.calcPlyStresses(Force)
		N1 = stresses[0]
		N2 = stresses[1]
		N12 = stresses[2]
		N12c = S*(1+2*pc11)**0.5

		Ra = pc11*S/pc12
		if N2 >= 0:	#Mode A

			f_IFFp = (((1/Yt-pt12/S)*N2)**2+(N12/S)**2)**0.5+pt12*S/N2
		elif abs(N2/N12) <= abs(Ra/N12c):#Mode B
			f_IFFp = ((N12/S)**2+(pc12*N2/S)**2)**0.5+pc12*N2/S
		else:#Mode C
			f_IFFp = ((N12/(2*(1+pc11)*S))**2+(N2/Yc)**2)*Yc/-N2
		if N1 >= 0: #FF tension
			f_FFp = N1/Xt
		else: #FF compression
			f_FFp = -N1/Xc

		return f_FFp, f_IFFp	#  result of analysis, if f_p is below 1 lamina did not fail, if it is 1 or higher lamina has failed
	def calcFailure(self, Load):
		failure = np.zeros(4)
		print("failure:", failure)
		dn = 50000 # Load increment
		lpf = False
		LoadFPF = np.zeros_like(Load)
		fpf = False
		LoadLPF = np.zeros_like(Load)
		Farg  = dn*10  # Initial length of load vector Ns-Ny
		strength = []
		for j in range(len(self.LayUp)):
			strength.append([self.Lamina[j].Xt, self.Lamina[j].Xc, self.Lamina[j].Yt, self.Lamina[j].Yc, self.Lamina[j].S])
		strength = np.asarray(strength, dtype=float)
		strength = strength.reshape((-1, 5))
		while lpf == False:
			stresses = self.calcPlyStresses(Load)
			f_xm = 0
			f_ym = 0
			f_sm = 0

			#Apply max stress failure criteria
			for j in range(0,len(self.LayUp)):   # Max stress Failure
				print("ply:", self.LayUp[j])
				print(f_xm, f_ym, f_sm)

				print("Strength:", strength[j,:])
				if stresses[0,j] > 0 and strength[j,0] != 0:
					f_xm = stresses[0,j] / strength[j,0]#Xt_mean
				elif stresses[0,j] < 0 and strength[j,0] != 0:
					f_xm = stresses[0,j] / -strength[j,1]#Xc_mean

				if stresses[1,j] > 0 and strength[j,0] != 0:
					f_ym = stresses[1,j] / strength[j,2]#Yt_mean
				elif stresses[1,j] < 0 and strength[j,0] != 0:
					f_ym = stresses[1,j] / -strength[j,3]#Yc_mean
				if stresses[2,j] > 0 and strength[j,0] != 0:
					f_sm = stresses[2,j] / strength[j,4]#S_mean
				elif stresses[2,j] < 0 and strength[j,0] != 0:
					f_sm = stresses[2,j] / -strength[j,4]
				print("Failure value")
				print(f_xm,f_ym,f_sm)
				# Determine which failure mode, and apply deg rule
				if (f_xm >= 1 and failure[j] <= 1) or ((f_ym >= 1 or f_sm >= 1) and failure[j] == 1):
					print("A")
					failure[j] = failure[j] + 1
					print("failure:", failure)
					if fpf == False:
						LoadFPF = Load
						fpf = True
						print("fpf")

					if self.LayUp[j] == 45 or self.LayUp[j] == -45:
						failure[j] = failure[j]-1
						failure[2:4] = failure[2:4]+1
						self.ABD[0:3,0:3] = self.ABD[0:3,0:3] - 8 * Q[j]*t
						strength[2:4,:] = 0
					else:
						self.ABD[0:3, 0:3] = self.ABD[0:3, 0:3] - 4 * Q[j]*t
						strength[j,:] = 0
						print(strength[j])

					if np.count_nonzero(failure) == 4:
						LoadLPF = Load
						lpf = True
						print("Last ply failure!!")

					break

				elif (f_ym >= 1 or f_sm >= 1) and failure[j] == 0:
					print("B")
					failure[j] = failure[j] + 1
					print("failure:", failure)
					if fpf == False:
						LoadFPF = Load
						fpf = True
						print("fpf")



					if self.LayUp[j] == 45 or self.LayUp[j] == -45:
						failure[j] = failure[j] - 1
						failure[2:4] = failure[2:4] + 1
						Q[j] = Laminate([45], Lamina(t, E1_mean, E2_mean*0.1, v12_mean, G12_mean)).QGlobalAr[0]
						self.ABD[0:3,0:3] = self.ABD[0:3,0:3] - 8*self.QGlobalAr[j] * t + 8 * Q[j] * t
						strength[2:4,2:4] = 0.1*strength[2:4,2:4]
						print("failure update:",failure)
					else:
						# Degrade transverse elastic properties
						t, E1, E2, v12, G12 = self.Lamina[j].t, self.Lamina[j].E1, self.Lamina[j].E2, self.Lamina[j].v12, self.Lamina[j].G12
						self.Lamina[j] = Lamina(t, E1, E2*0.1, v12, G12)
						self.calcQGlobalLaminas()
						self.calcABD()

						strength[j,2:4] = 0.1*strength[j,2:4]
						print(strength[j,:])

					if np.count_nonzero(failure) == 4:
						LoadLPF = Load
						lpf = True
						print("Last ply failure!!")

					break
				# When no failure detected, continue to next ply and increase
				elif j == 3:
					print("C")
					Farg = Farg + dn
					print("Force:", Farg)
					print("Failures:", failure)
				#lpf = True

	def __repr__(self):
		return f"Laminate of layup {self.LayUp}"

if __name__ == "__main__":
	# Simple test to see if a lamina turned 90 deg gives the same tensile stiffness
	E1 = 140 * 10 ** 9
	E2 = 10 * 10 ** 9
	G12 = 5 * 10 ** 9
	v12 = 0.3
	t = 0.125e-3

	v21 = v12 * E2 / E1
	Q = 1 - v12 * v21
	Q11 = E1 / Q
	Q22 = E2 / Q
	Q12 = v12 * E2 / Q
	Q66 = G12
	Lamina_ = Lamina(t, E1, E2, v12, G12)
	Laminate_1 = Laminate([15, 15, 15, 15], Lamina_)
	Laminate_2 = Laminate([105, 105, 105, 105], Lamina_)
	print(Laminate_1.ABD)
	print(Laminate_2.ABD)
	print(Laminate_1.calcEngConst())
	print(Laminate_2.calcEngConst())
	_, Ey, vxy, _, _ = Laminate_1.calcEngConst()
	Ex, _, _, _, _ = Laminate_2.calcEngConst()
	assert np.isclose(Ex, Ey)
	Laminate_3 = Laminate([15, 0, 0, 75], Lamina_)
	print(vxy)