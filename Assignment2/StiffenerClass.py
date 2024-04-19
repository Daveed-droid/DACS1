import matplotlib.pyplot as plt
import numpy as np
from Assignment2Data import AssignmentLamina
from LaminateClass import Laminate

def Rx(theta):
	R = np.zeros((2,2))
	R[0, 0] = np.cos(theta)
	R[1, 1] = np.cos(theta)
	R[1, 0] = np.sin(theta)
	R[0, 1] = -np.sin(theta)
	return R
def get_length(v):
	v1, v2 = v[0, :], v[1, :]
	l = np.linalg.norm(v2-v1)
	return l

def get_center(v):
	v1, v2 = v[0, :], v[1, :]
	d = v2-v1
	center = d/2 + v1
	return center

def get_angle(v):
	v1, v2 = v[0, :], v[1, :]
	d = v2-v1
	angle = np.arctan2(d[1], d[0])
	return angle

class Stiffener():
	def __init__(self, LaminateWeb, LaminateFlange, CrossSection, angle = 0):
		tw, tf = LaminateWeb.h/2, LaminateFlange.h/2
		self.LaminateWeb, self.LaminateFlange = LaminateWeb, LaminateFlange
		self.buckProp = 0
		if CrossSection == "L":
			self.Flange = [True, False]
			self.midlines = np.zeros((4, 2))
			self.midlines[:, 0] = np.array([0, 0, 0.04, 0])
			self.midlines[:, 1] = np.array([0, 0+tf, 0, 0.04])
			self.midlines += tf
			self.buckProp = [["OEF", 0.04-tf, self.LaminateWeb]]
		elif CrossSection == "I":
			self.Flange = [True, False, True]
			self.midlines = np.zeros((4, 3))
			self.midlines[:, 0] = np.array([-0.02, 0, 0.02, 0])
			self.midlines[:, 1] = np.array([0, 0+tf, 0, 0.02-tf])
			self.midlines[:, 2] = np.array([-0.02, 0.02, 0.02, 0.02])
			self.midlines += tf
			self.buckProp = [["NEF", 0.02-2*tf, self.LaminateWeb],
							 ["OEF", 0.02-2*tw, self.LaminateFlange]]
		elif CrossSection == "hat":
			self.Flange = [True, False, False, True]
			self.midlines = np.zeros((4, 4))
			self.midlines[:, 0] = np.array([-0.03, 0, 0.03, 0])
			self.midlines[:, 1] = np.array([-0.02+tf, 0+tf, -0.01-tf, 0.02-tf])
			self.midlines[:, 2] = np.array([0.02-tf, 0+tf, 0.01+tf, 0.02-tf])
			self.midlines[:, 3] = np.array([-0.01, 0.02, 0.01, 0.02])
			self.midlines += tf
			self.buckProp = [["NEF", 0.02-2*tw, self.LaminateFlange],
							 ["NEF", (0.01**2+0.02**2)**0.5-2*tf, self.LaminateWeb]]
		else:
			raise TypeError


		# Rotate
		for i, line in enumerate(self.midlines.T):
			self.midlines[0:2, i] = (Rx(np.deg2rad(angle))@line[0:2].reshape((-1, 1))).reshape(-1)
			self.midlines[2:4, i] = (Rx(np.deg2rad(angle))@line[2:4].reshape((-1, 1))).reshape(-1)

		self.EAx = 0
		self.EAy = 0
		self.EA = 0
		self.GA = 0
		self.A = 0
		# Setting axial stiffness and shear stiffness
		for i, line in enumerate(self.midlines.T):
			line = line.reshape((2,2))
			l = get_length(line)
			cen = get_center(line)
			ang = get_angle(line)
			if self.Flange[i]:
				Ex, Ey, vxy, vyx, Gxy = LaminateFlange.calcEngConst()
				h = LaminateFlange.h
			elif self.Flange[i]:
				Ex, Ey, vxy, vyx, Gxy = LaminateWeb.calcEngConst()
				h = LaminateWeb.h
			self.EAx += Ex*l*h*cen[0]
			self.EAy += Ex*l*h*cen[1]
			self.EA += Ex*l*h
			self.GA += Gxy*l*h*np.sin(ang)
			self.A += l*h
		self.x_cg = self.EAx/self.EA
		self.y_cg = self.EAy/self.EA

	def set_angle(self, angle):
		# Rotate
		for i, line in enumerate(self.midlines.T):
			self.midlines[0:2, i] = (Rx(np.deg2rad(angle))@line[0:2].reshape((-1, 1))).reshape(-1)
			self.midlines[2:4, i] = (Rx(np.deg2rad(angle))@line[2:4].reshape((-1, 1))).reshape(-1)

		self.EAx = 0
		self.EAy = 0
		self.EA = 0
		self.GA = 0
		# Setting axial stiffness and shear stiffness
		for i, line in enumerate(self.midlines.T):
			line = line.reshape((2,2))
			l = get_length(line)
			cen = get_center(line)
			ang = get_angle(line)
			if self.Flange[i]:
				Ex, Ey, vxy, vyx, Gxy = self.LaminateFlange.calcEngConst()
				h = self.LaminateFlange.h
			elif self.Flange[i]:
				Ex, Ey, vxy, vyx, Gxy = self.LaminateWeb.calcEngConst()
				h = self.LaminateWeb.h
			self.EAx += Ex*l*h*cen[0]
			self.EAy += Ex*l*h*cen[1]
			self.EA += Ex*l*h
			self.GA += Gxy*l*h*np.sin(ang)
		self.x_cg = self.EAx/self.EA
		self.y_cg = self.EAy/self.EA
		return self

	def get_fpf(self, strains_fuselage):
		Fail = []
		for i, line in enumerate(self.midlines.T):
			if self.Flange[i]:
				Lam = self.LaminateFlange
			elif self.Flange[i]:
				Lam = self.LaminateWeb
			line = line.reshape((2,2))
			l = get_length(line)
			cen = get_center(line)
			ang = get_angle(line)
			ex_f = strains_fuselage[0]
			exy_f = strains_fuselage[2]
			ex = ex_f
			exy = np.sin(np.deg2rad(ang))*exy_f
			LamStrains = np.zeros((3, 1))
			LamStrains[0, :] = ex
			LamStrains[2, :] = exy
			f_FFp, f_IFFp = Lam.PuckStrain(LamStrains.reshape(-1))
			a = np.max(f_FFp)
			b = np.max(f_IFFp)
			d = max(a, b)
			Fail.append(d)
		return Fail


if __name__=="__main__":
	Lamw = Laminate([0, 0, 0], AssignmentLamina)
	Lamf = Laminate([0, 0, 0], AssignmentLamina)
	Stif = Stiffener(Lamw, Lamf, "I", angle = 0)
	mid = Stif.midlines
	fig, ax = plt.subplots()
	for i in range(len(mid[0])):
		v = mid[:, i].reshape((2, 2))
		plt.plot(v[:, 0], v[:, 1], color="blue", linewidth=10)
		center = get_center(v)
		plt.plot(center[0], center[1], color="red", marker="o", linewidth=10)

	ax.set_aspect('equal', 'box')
	plt.show()