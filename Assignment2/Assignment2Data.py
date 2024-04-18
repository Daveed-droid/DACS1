from LaminaClass import Lamina

Ex = 142e9
Ey = 11.2e9
Gxy = 5e9
vxy = 0.3
Xt = 2200e6
Xc = 1800e6
Yt = 70e6
Yc = 300e6
S = 100e6
t = 0.135e-3
P = 1610

AssignmentLamina = Lamina(t, Ex, Ey, vxy, Gxy)
AssignmentLamina.setStrengths(Xt, Yt, Xc, Yc, S)
class Metal:
	def __init__(self, E, G, v):
		self.E, self.G, self.v = E, G, v
	def setStrengths(self, Xt, Yt, Xc, Yc, S):
		self.Xt, self.Xc = Xt, Xc
		self.Yt, self.Yc = Yt, Yc
		self.S = S

AssignmentMetal = Metal(69e9, 26e9, 0.29)
AssignmentMetal.setStrengths(410e6, 400e6, 430e6, 430e6, 230e6)

AssignmentMetalLamina = Lamina(t, 69e9, 69e9, 0.29, 26e9)
AssignmentMetalLamina.setStrengths(410e6, 400e6, 430e6, 430e6, 230e6)
