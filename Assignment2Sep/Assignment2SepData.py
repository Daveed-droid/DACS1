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

AssignmentLamina = Lamina(t, Ex, Ey, vxy, Gxy)
AssignmentLamina.setStrengths(Xt, Yt, Xc, Yc, S)

