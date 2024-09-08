import numpy as np


class Plate:
    def __init__(self, Laminate, a):
        self.Pcr = None
        self.Px = None
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
        self.w = self.w11*np.sin(np.pi*x/self.a)*np.sin(np.pi*y/self.a)
        self.dwdxx = -(np.pi/self.a)**2*self.w
        self.dwdyy = -(np.pi/self.a)**2*self.w
        self.dwdxy = self.w11*(np.pi/self.a)**2*np.cos(np.pi*x/self.a)*np.cos(np.pi*y/self.a)
        # Step 5: Find loads
        y = np.linspace(0, self.a, 50, endpoint = True)
        self.Nx = -(self.Px/self.a + 4*np.pi**2/self.a**2 * self.K02 * np.cos(2*np.pi*y/self.a))
        self.Ny = -(self.Py / self.a + 4 * np.pi ** 2 / self.a ** 2 * self.K20 * np.cos(2 * np.pi * x / self.a))
        self.Nxy = 0
        self.Mx = -self.D11*self.dwdxx-self.D12*self.dwdyy
        self.My = -self.D12 * self.dwdxx - self.D22 * self.dwdyy
        self.Mxy = -2*self.D66*self.dwdxy
    def calcFailureNx(self):
        pass

    def calcFailureAll(self):
        pass


def SkinBuckling(Laminate):
    # Compression
    a = 1  # Length
    AR = a / a  # is square
    if AR < 1.4:
        m = 1
    elif AR < 2.45:
        m = 2
    else:
        m = 3

    D11 = Laminate.ABD[3, 3]
    D66 = Laminate.ABD[5, 5]
    D12 = Laminate.ABD[3, 4]
    D22 = Laminate.ABD[4, 4]

    # Simply supported plate
    Pcr = (np.pi ** 2 * (D11 * m ** 2 + 2 * (D12 + 2 * D66) * AR ** 2 + D22 * AR ** 4 / m ** 2)) / (a ** 2)
    beff = a * (2 * (1 + 2 * (1 + A12 / A11) * (1 - Pcr / Px) * A11 / (A11 + 3 * A22))) ** -1
    return Pc
