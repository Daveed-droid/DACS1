"""
@Project ：DACS1
@File ：LaminaClass.py
@Author ：David Canosa Ybarra
@Date ：19/03/2024 11:44
"""


class Lamina():
	def __init__(self, t, E1, E2, v12, G12):
		self.t = t
		self.E1 = E1
		self.E2 = E2
		self.v12 = v12
		self.G12 = G12

	def setStrengths(self, Xt, Yt, Xc, Yc, S):
		self.Xt = Xt
		self.Yt = Yt
		self.Xc = Xc
		self.Yc = Yc
		self.S = S
