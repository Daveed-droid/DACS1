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
	def __init__(self, LayUp: list, Lamina: Lamina):
		"""
		Initializes the laminate class, this class takes in a layup and lamina
		:param LayUp: A list of angles [deg] starting from the bottom ply
		:param Lamina: A lamina object
		:param t: The lamina thickness
		"""
		self.LayUp = LayUp
		self.Lamina = Lamina
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
		Q11 = self.Lamina.QMatrix[0, 0]
		Q12 = self.Lamina.QMatrix[0, 1]
		Q22 = self.Lamina.QMatrix[1, 1]
		Q66 = self.Lamina.QMatrix[2, 2]
		for i, theta in enumerate(self.LayUp):
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


		f_p = 0.5	# Placeholder

		return f_p	# f_p is result of analysis, if f_p is below 1 lamina did not fail, if it is 1 or higher lamina has failed

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