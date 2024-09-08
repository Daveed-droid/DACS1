import numpy as np
import matplotlib.pyplot as plt
from Assignment2SepData import AssignmentLamina
from LaminateClass import Laminate


class Plate:
    def __init__(self, Laminate, a, Px):
        self.Pcr = None
        self.Px = Px
        self.w11 = None
        self.K02 = None
        self.K20 = None
        self.Laminate = Laminate
        self.a = a
        self.A11 = Laminate.ABD[1, 1]
        self.A12 = Laminate.ABD[1, 2]
        self.A22 = Laminate.ABD[2, 2]
        self.D11 = Laminate.ABD[3, 3]
        self.D66 = Laminate.ABD[5, 5]
        self.D12 = Laminate.ABD[3, 4]
        self.D22 = Laminate.ABD[4, 4]

    def calcLoads(self):
        # Step 1: Calculate buckling load Pcr
        self.Pcr = np.pi ** 2 / self.a * (self.D11 + 2 * (self.D12 + 2 * self.D66) + self.D22) / (
                    1 + self.A12 / self.A11)
        # Step 2: Find deflection w_11
        self.w11 = np.sqrt((16*self.A11*self.A22*(self.D11 + 2*(self.D12 + 2*self.D66) + self.D22))/((self.A11*self.A22 - self.A12**2)*(self.A11 + 3*self.A22))*(self.Px/self.Pcr - 1))
        # Step 3: Find K values
        self.K02 = ((self.A11*self.A22 - self.A12**2)/self.A22) * self.w11**2 / 32
        self.K20 = ((self.A11*self.A22 - self.A12**2)/self.A11) * self.w11**2 / 32
        # Step 4: Find w and its derivatives
        self.Py = self.Px*self.A12/self.A11 - self.w11**2 * (np.pi**2/(8*self.a))*(self.A11*self.A22 - self.A12**2)/self.A11
        x = np.linspace(0, self.a, 50, endpoint = True)
        y = np.linspace(0, self.a, 50, endpoint = True)
        self.xv, self.yv = np.meshgrid(x, y)
        self.w = self.w11*np.sin(np.pi*self.xv/self.a)*np.sin(np.pi*self.yv/self.a)
        self.dwdxx = -(np.pi/self.a)**2*self.w
        self.dwdyy = -(np.pi/self.a)**2*self.w
        self.dwdxy = self.w11*(np.pi/self.a)**2*np.cos(np.pi*self.xv/self.a)*np.cos(np.pi*self.yv/self.a)
        # Step 5: Find loads
        self.Nx = -(self.Px/self.a + 4*np.pi**2/self.a**2 * self.K02 * np.cos(2*np.pi*self.yv/self.a))
        self.Ny = -(self.Py / self.a + 4 * np.pi ** 2 / self.a ** 2 * self.K20 * np.cos(2 * np.pi * self.xv / self.a))
        self.Nxy = self.xv*0
        self.Mx = -self.D11*self.dwdxx-self.D12*self.dwdyy
        self.My = -self.D12 * self.dwdxx - self.D22 * self.dwdyy
        self.Mxy = -2*self.D66*self.dwdxy

        fig, ax = plt.subplots()
        c = ax.pcolormesh(self.xv, self.yv, self.w, cmap='RdBu')
        ax.set_title('Plot')
        fig.colorbar(c, ax=ax)
        plt.show()

    def calcFailureNx(self):
        Nx, Ny, Nxy, Mx, My, Mxy = np.reshape(self.Nx, -1), np.reshape(self.Ny, -1), np.reshape(self.Nxy, -1), np.reshape(self.Mx, -1), np.reshape(self.My, -1), np.reshape(self.Mxy, -1)
        Failed = np.reshape(np.zeros_like(self.xv), -1)
        for i in range(len(np.reshape(self.xv, -1))):
            Load = np.array([Nx[i], 0, 0, 0, 0, 0])
            f_FFp, f_IFFp = self.Laminate.Puck(Load)
            Failed[i] = 1.0 if 1.0 < f_FFp.max() or 1.0 < f_IFFp.max() else 0.0
        Failed = np.reshape(Failed, (self.xv.shape[0], self.xv.shape[1]))
        fig, ax = plt.subplots()
        c = ax.pcolormesh(self.xv, self.yv, Failed, cmap='RdBu')
        ax.set_title('Plot')
        fig.colorbar(c, ax=ax)
        plt.show()
    def calcFailureAll(self):
        Nx, Ny, Nxy, Mx, My, Mxy = np.reshape(self.Nx, -1), np.reshape(self.Ny, -1), np.reshape(self.Nxy, -1), np.reshape(self.Mx, -1), np.reshape(self.My, -1), np.reshape(self.Mxy, -1)
        Failed = np.reshape(np.zeros_like(self.xv), -1)
        for i in range(len(np.reshape(self.xv, -1))):
            Load = np.array([Nx[i], Ny[i], Nxy[i], Mx[i], My[i], Mxy[i]])
            f_FFp, f_IFFp = self.Laminate.Puck(Load)
            Failed[i] = 1.0 if 1.0 < f_FFp.max() or 1.0 < f_IFFp.max() else 0.0
        Failed = np.reshape(Failed, (self.xv.shape[0], self.xv.shape[1]))
        fig, ax = plt.subplots()
        c = ax.pcolormesh(self.xv, self.yv, Failed, cmap='RdBu')
        ax.set_title('Plot')
        fig.colorbar(c, ax=ax)
        plt.show()

if __name__=="__main__":
    Lam = Laminate([0,20,30,50], AssignmentLamina)
    A = Plate(Lam, 1, 150000)
    A.calcLoads()
    A.calcFailureNx()
