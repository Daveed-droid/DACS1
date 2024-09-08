import numpy as np

from Assignment2.Assignment2Data import AssignmentLamina
from LaminateClass import Laminate
def sym(layup, n):
	for i in range(n):
		layup += list(np.flip(layup))
	return layup
CompLam = sym([0, 90, 45, -45, -45, 45, 45, -45, -45, 45], 1)*3
Lam1 = Laminate(CompLam, AssignmentLamina)

LamTest = Lam1

#Compression
a = 1 #Length

b = 1
AR = a/b
if AR<1.4:
	m = 1
elif AR<2.45:
	m=2
else:
	m=3

D11 = LamTest.ABD[3, 3]
D66 = LamTest.ABD[5, 5]
D12 = LamTest.ABD[3, 4]
D22 = LamTest.ABD[4, 4]


NCom = (np.pi**2*(D11*m**2+2*(D12+2*D66)*AR**2+D22*AR**4/m**2))/(a**2)

#Shear

D11 = LamTest.ABD[3,3]
D66 = LamTest.ABD[5,5]
D12 = LamTest.ABD[3,4]
D22 = LamTest.ABD[4,4]
A = -0.27 + 0.185*(D12+2*D66)/(D11*D22)**0.5
B = 0.82+0.46*(D12+2*D66)/(D11*D22)-0.2*((D12+2*D66)/(D11*D66)**0.5)**2
beta = (D11/D22)**0.25
K = 8.2 + 5*(D12+2*D66)/((D11*D22)**0.5*(A/beta+B*beta))
NShear = 4*(D11*D22**3)**0.25*K/b**2
Nyx = (9*np.pi**4*b)/(32*a**3) * (D11 + 2*(D12 + 2*D66)*a**2/b**2 + D22*a**4/b**4)

print(NCom)
print(NShear)
print(Nyx)