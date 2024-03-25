"""
@Project ：DACS1
@File ：Q1b.py
@Author ：David Canosa Ybarra
@Date ：19/03/2024 11:45
"""
from AssignmentData import Lamina_mean
from LaminateClass import Laminate
import numpy as np


def Q1b(Lamina_mean = Lamina_mean, Laminate = Laminate):
	LayUp = [0, 0, 90, 30, 90]
	Laminate = Laminate(LayUp, Lamina_mean)
	Load = np.array([0.2e2, 1.8e4, 0, 18e3, 0, 0]).T
	print(Laminate.calcPlyStrains(Load))
	print(Laminate.calcStresses(Load))


Q1b()
