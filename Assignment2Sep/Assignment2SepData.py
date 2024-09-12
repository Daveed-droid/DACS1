from LaminaClass import Lamina

KnockDownFactor = 0.8 * 0.65 * 0.8

Ex = 142e9
Ey = 11.2e9
Gxy = 5e9
vxy = 0.3
Xt = 2200e6 * KnockDownFactor
Xc = 1800e6 * KnockDownFactor
Yt = 70e6 * KnockDownFactor
Yc = 300e6 * KnockDownFactor
S = 100e6 * KnockDownFactor
t = 0.135e-3

AssignmentLamina = Lamina(t, Ex, Ey, vxy, Gxy)
AssignmentLamina.setStrengths(Xt, Yt, Xc, Yc, S)

