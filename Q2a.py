"""
@Project ：DACS1
@File ：Q2a.py
@Author ：David Canosa Ybarra
@Date ：19/03/2024 11:45
"""
import numpy as np
import math
import matplotlib.pyplot as plt
from LaminaClass import Lamina
from LaminateClass import Laminate
#def Q2a(Lamina_mean = Lamina_mean, Laminate = Laminate):
def Q2a():
    LayUp = [0,90,45,-45]
    LayUp = np.append(LayUp, np.flip(LayUp))
    LayUp = np.append(LayUp, np.flip(LayUp))
    print(LayUp)


Q2a()