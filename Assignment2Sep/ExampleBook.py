from PlateClass import *
from LaminaClass import Lamina

Ex = 69e9
Ey = 69e9
Gxy = 5.1e9
vxy = 0.05
Xt = 2200e6
Xc = 1800e6
Yt = 70e6
Yc = 300e6
S = 100e6
t = 0.19e-3

BookLamina = Lamina(t, Ex, Ey, vxy, Gxy)
BookLamina.setStrengths(Xt, Yt, Xc, Yc, S)

A = Plate(0.2, 2152)
A.ABD(Laminate([45,0,0,0,45], BookLamina))
A.calcLoads()
A.calcFailureAll(verbose = True)
print(Laminate([45,0,0,0,45], BookLamina).ABD)
print(A.Pcr)
print(A.A11)
print(A.w11)
print(A.K02)
print(A.K20)