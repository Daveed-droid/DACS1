"""
@Project ：DACS1
@File ：LaminateClass.py
@Author ：David Canosa Ybarra
@Date ：19/03/2024 11:44
"""
from LaminaClass import Lamina
import numpy as np

class Laminate():
	def __init__(self, LayUp: list, Lamina: Lamina, t: float):
		self.LayUp = LayUp
		self.Lamina = Lamina
		self.t = t
		self.h = t * len(LayUp)
		self.zlst = np.linspace(-self.h / 2, self.h / 2, len(self.LayUp), endpoint = True)
		self.calcQGlobalLaminas()
		self.calcABD()

	def calcQGlobalLaminas(self):
		self.QGlobalar = np.zeros((1, len(self.LayUp)))
		Q11 = self.Lamina.QMatrix[0, 0]
		Q12 = self.Lamina.QMatrix[0, 1]
		Q22 = self.Lamina.QMatrix[1, 1]
		Q66 = self.Lamina.QMatrix[2, 2]
		for i, theta in enumerate(self.LayUp):
			m = np.sin(np.rad2deg(theta))
			n = np.cos(np.rad2deg(theta))
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
		self.AMatrix = np.zeros((3, 3))
		self.BMatrix = np.zeros((3, 3))
		self.DMatrix = np.zeros((3, 3))
		# calculating the A Matrix
		for i in range(3):
			for j in range(3):
				for k in range(len(self.LayUp)):
					Q = self.QGlobalar[k]
					self.AMatrix[i, j] += Q[i, j] * (self.zlst[k + 1] - self.zlst[k])
		# calculating the B Matrix
		for i in range(3):
			for j in range(3):
				for k in range(len(self.LayUp)):
					Q = self.QGlobalar[k]
					self.BMatrix[i, j] += 0.5 * Q[i, j] * (self.zlst[k + 1] ** 2 - self.zlst[k] ** 2)
		# calculating the D Matrix
		for i in range(3):
			for j in range(3):
				for k in range(len(self.LayUp)):
					Q = self.QGlobalar[k]
					self.DMatrix[i, j] += 3 ** -1 * Q[i, j] * (self.zlst[k + 1] ** 3 - self.zlst[k] ** 3)
		ABD_top = np.hstack((self.AMatrix, self.BMatrix))
		ABD_bottom = np.hstack((self.BMatrix, self.DMatrix))
		self.ABD = np.vstack((ABD_top, ABD_bottom))
	def calcStrains(self, Load):
		return np.linalg.inv(self.ABD) @ Load

	def calcEngConst(self):
		pass


def Ex(theta, E1, E2, G12, v12):
	m = np.sin(theta)
	n = np.cos(theta)
	_a = m ** 4 / E1
	_b = (G12 ** -1 - 2 * v12 / E1) * m ** 2 * n ** 2
	_c = n ** 4 / E2
	_Ex = (_a + _b + _c) ** -1
	return _Ex


def vxy(theta, Ex, E1, E2, G12, v12):
	m = np.sin(theta)
	n = np.cos(theta)
	_a = Ex(theta, E1, E2, G12, v12)
	_b = v12 / E1 * (m ** 4 + n ** 4)
	_c = (1 / E1 + 1 / E2 - 1 / G12) * m ** 2 * n ** 2
	_vxy = _a * (_b - _c)
	return _vxy


def Ey(theta, E1, E2, G12, v12):
	m = np.sin(theta)
	n = np.cos(theta)
	_a = n ** 4 / E1
	_b = (G12 ** -1 - 2 * v12 / E1) * m ** 2 * n ** 2
	_c = m ** 4 / E2
	_Ey = (_a + _b + _c) ** -1
	return _Ey


def Gxy(theta, E1, E2, G12, v12):
	m = np.sin(theta)
	n = np.cos(theta)
	_a = 2 * (2 / E1 + 2 / E2 + 4 * v12 / E1 - 1 / G12) * m ** 2 * n ** 2
	_b = 1 / G12 * (m ** 4 + n ** 4)
	_Gxy = (_a + _b) ** -1
	return _Gxy


def etaxs(theta, Ex, E1, E2, G12, v12):
	m = np.sin(theta)
	n = np.cos(theta)
	_a = Ex(theta, E1, E2, G12, v12)
	_b = (2 / E1 + 2 * v12 / E1 - 1 / G12) * m ** 3 * n
	_c = (2 / E2 + 2 * v12 / E1 - 1 / G12) * n ** 3 * m
	_etaxs = _a * (_b - _c)
	return _etaxs


def etays(theta, Ex, E1, E2, G12, v12):
	m = np.sin(theta)
	n = np.cos(theta)
	_a = Ey(theta, E1, E2, G12, v12)
	_b = (2 / E1 + 2 * v12 / E1 - 1 / G12) * n ** 3 * m
	_c = (2 / E2 + 2 * v12 / E1 - 1 / G12) * m ** 3 * n
	_etays = _a * (_b - _c)
	return _etays
