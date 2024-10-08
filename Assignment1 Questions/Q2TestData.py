"""
@Project ：DACS1
@File ：AssignmentData.py
@Author ：David Canosa Ybarra
@Date ：25/03/2024 15:43
"""
from LaminaClass import Lamina

t = 0.125e-3  # mm
E1_mean = 140e9  # Pa
E1_std = 3.28e9  # Pa
E2_mean = 10e9  # Pa
E2_std = 1.28e9  # Pa
v12_mean = 0.30  # -
v12_std = 0.018  # -
G12_mean = 5e9  # Pa
G12_std = 0.83e9  # Pa
Xt_mean = 1500e6  # Pa
Xt_std = 128.3e6  # Pa
Yt_mean = 50e6  # Pa
Yt_std = 8.2e6  # Pa
Xc_mean = 1200e6  # Pa
Yc_mean = 250e6  # Pa
S_mean = 70e6  # Pa
S_std = 6.21e6  # Pa

Lamina_mean = Lamina(t, E1_mean, E2_mean, v12_mean, G12_mean)
