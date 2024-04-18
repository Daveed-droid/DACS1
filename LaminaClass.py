"""
@Project ：DACS1
@File ：LaminaClass.py
@Author ：David Canosa Ybarra
@Date ：19/03/2024 11:44
"""
import numpy as np

class Lamina():
	"""
	This class serves as a material data input to the Laminate class
	"""
	def __init__(self, t, E1, E2, v12, G12):
		self.t = t
		self.E1 = E1
		self.E2 = E2
		self.v12 = v12
		self.G12 = G12
		self.calcQMatrix()

	def setStrengths(self, Xt, Yt, Xc, Yc, S):
		self.Xt = Xt
		self.Yt = Yt
		self.Xc = Xc
		self.Yc = Yc
		self.S = S

	def calcQMatrix(self):
		self.QMatrix = np.zeros((3, 3))
		self.v21 = self.v12 * self.E2 / self.E1
		self.Q = 1 - self.v12 * self.v21
		self.QMatrix[0, 0] = self.E1 / self.Q
		self.QMatrix[1, 1] = self.E2 / self.Q
		self.QMatrix[1, 0] = (self.v12 * self.E2) / self.Q
		self.QMatrix[0, 1] = self.QMatrix[1, 0]
		self.QMatrix[2, 2] = self.G12



if __name__ == '__main__':
	# Testing function to match lecture values
	E1 = 140 * 10 ** 9
	E2 = 10 * 10 ** 9
	G12 = 5 * 10 ** 9
	v12 = 0.3
	t = 0.125e-3
	Lamina = Lamina(t, E1, E2, v12, G12)
	print(Lamina.QMatrix)
