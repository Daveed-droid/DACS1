"""
@Project ：DACS1
@File ：LaminateClass.py
@Author ：David Canosa Ybarra
@Date ：19/03/2024 11:44
"""
from LaminaClass import Lamina
import numpy as np

class Laminate():
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
		self.h = t * len(LayUp)
		self.zlst = np.linspace(-self.h / 2, self.h / 2, len(self.LayUp), endpoint = True)
		self.calcQGlobalLaminas()
		self.calcABD()

	def calcQGlobalLaminas(self):
		"""
		Calculates the Q matrix of all the laminas when placed at an angle
		:return: None
		"""
		self.QGlobalar = list(range(len(self.LayUp)))
		Q11 = self.Lamina.QMatrix[0, 0]
		Q12 = self.Lamina.QMatrix[0, 1]
		Q22 = self.Lamina.QMatrix[1, 1]
		Q66 = self.Lamina.QMatrix[2, 2]
		for i, theta in enumerate(self.LayUp):
			m = np.sin(np.deg2rad(theta))
			n = np.cos(np.deg2rad(theta))
			Qxx = Q11 ** 4 + 2 * (Q12 + 2 * Q66) * m ** 2 * n ** 2 + Q22 * n ** 4
			Qxy = (Q11 + Q22 - 4 * Q66) * m ** 2 * n ** 2 + Q12 * (m ** 4 + n ** 4)
			Qyy = Q11 * n ** 4 + 2 * (Q12 + 2 * Q66) * m ** 2 * n ** 2 + Q22 * m ** 4
			Qxs = (Q11 - Q12 - 2 * Q66) * n * m ** 3 + (Q12 - Q22 + 2 * Q66) * n ** 3 * m
			Qys = (Q11 - Q12 - 2 * Q66) * m * n ** 3 + (Q12 - Q22 + 2 * Q66) * m ** 3 * n
			Qss = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * n ** 2 * m ** 2 + Q66 * (m ** 4 + n ** 4)
			QMatrix_glo = np.array([[Qxx, Qxy, Qxs],
									[Qxy, Qyy, Qys],
									[Qxs, Qys, Qss]])
			self.QGlobalar[i] = QMatrix_glo
	def calcABD(self):
		"""
		Calculates the ABD matrix of the laminate
		:return:
		"""
		self.AMatrix = np.zeros((3, 3))
		self.BMatrix = np.zeros((3, 3))
		self.DMatrix = np.zeros((3, 3))
		# calculating the A Matrix
		for i in range(3):
			for j in range(3):
				for k in range(len(self.LayUp) - 1):
					Q = self.QGlobalar[k]
					self.AMatrix[i, j] += Q[i, j] * (self.zlst[k + 1] - self.zlst[k])
		# calculating the B Matrix
		for i in range(3):
			for j in range(3):
				for k in range(len(self.LayUp) - 1):
					Q = self.QGlobalar[k]
					self.BMatrix[i, j] += 0.5 * Q[i, j] * (self.zlst[k + 1] ** 2 - self.zlst[k] ** 2)
		# calculating the D Matrix
		for i in range(3):
			for j in range(3):
				for k in range(len(self.LayUp) - 1):
					Q = self.QGlobalar[k]
					self.DMatrix[i, j] += 3 ** -1 * Q[i, j] * (self.zlst[k + 1] ** 3 - self.zlst[k] ** 3)
		ABD_top = np.hstack((self.AMatrix, self.BMatrix))
		ABD_bottom = np.hstack((self.BMatrix, self.DMatrix))
		self.ABD = np.vstack((ABD_top, ABD_bottom))
	def calcStrains(self, Load):
		"""
		Calculates the deflection of the laminate at a prescribed load
		:param Load: A 6x1 numpy array column vector with the loads applied
		:return: The deflections of the laminate
		"""
		return np.linalg.inv(self.ABD) @ Load

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


if __name__ == "__main__":
	E1 = 140 * 10 ** 9
	E2 = 10 * 10 ** 9
	G12 = 5 * 10 ** 9
	v12 = 0.3
	t = 0.125

	v21 = v12 * E2 / E1
	Q = 1 - v12 * v21
	Q11 = E1 / Q
	Q22 = E2 / Q
	Q12 = v12 * E2 / Q
	Q66 = G12
	Lamina_ = Lamina(t, E1, E2, G12, v12)
	Laminate_1 = Laminate([0, 0, 0, 0], Lamina_)
	Laminate_2 = Laminate([90, 90, 90, 90], Lamina_)
	_, Ey, _, _, _ = Laminate_1.calcEngConst()
	Ex, _, _, _, _ = Laminate_2.calcEngConst()
	assert np.isclose(Ex, Ey)
