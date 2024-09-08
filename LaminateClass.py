"""
@Project ：DACS1
@File ：LaminateClass.py
@Author ：David Canosa Ybarra
@Date ：19/03/2024 11:44
"""
from LaminaClass import Lamina
import numpy as np


def StressGloTOPly(theta):
    """
	Calculate the rotation matrix to go from stresses in the laminate ref frame to stresses in the lamina ref frame
	:param theta: Angle of the ply in Radians
	:return: Rotation matrix
	"""
    n = np.sin(np.deg2rad(theta))
    m = np.cos(np.deg2rad(theta))
    RotMat = np.array([[m ** 2, n ** 2, 2 * m * n],
                       [n ** 2, m ** 2, -2 * m * n],
                       [-m * n, m * n, m ** 2 - n ** 2]])
    return RotMat


def StrainGloTOPly(theta):
    """
	Calculate the rotation matrix to go from strains in the laminate ref frame to strains in the lamina ref frame
	:param theta: Angle of the ply in Radians
	:return: Rotation matrix
	"""
    n = np.sin(np.deg2rad(theta))
    m = np.cos(np.deg2rad(theta))
    RotMat = np.array([[m ** 2, n ** 2, m * n],
                       [n ** 2, m ** 2, -m * n],
                       [-2 * m * n, 2 * m * n, m ** 2 - n ** 2]])
    return RotMat


class Laminate:
    """
	The laminate class will be used to do operations on the layup
	"""

    def __init__(self, LayUp: list, Lamina: Lamina | list):
        """
		Initializes the laminate class, this class takes in a layup and lamina
		:param LayUp: A list of angles [deg] starting from the bottom ply
		:param Lamina: A lamina object
		:param t: The lamina thickness
		"""
        self.LayUp = LayUp
        if type(Lamina) == list:
            self.Lamina = Lamina
            constantUD = False
        else:
            self.Lamina = [Lamina] * len(LayUp)
            constantUD = True
        self.t = [self.Lamina[i].t for i in range(len(LayUp))]
        self.h = np.sum(self.t)
        # TODO: Make zlist allow varying t
        self.zlst = np.linspace(-self.h / 2, self.h / 2, len(self.LayUp) + 1, endpoint = True)
        self.calcQGlobalLaminas()
        self.calcABD()

    def calcQGlobalLaminas(self):
        """
		Calculates the Q matrix of all the laminas when placed at an angle
		:return: None
		"""
        self.QGlobalAr = list(range(len(self.LayUp)))

        for i, theta in enumerate(self.LayUp):
            Q11 = self.Lamina[i].QMatrix[0, 0]
            Q12 = self.Lamina[i].QMatrix[0, 1]
            Q22 = self.Lamina[i].QMatrix[1, 1]
            Q66 = self.Lamina[i].QMatrix[2, 2]
            n = np.sin(np.deg2rad(theta))
            m = np.cos(np.deg2rad(theta))
            Qxx = Q11 * m ** 4 + 2 * (Q12 + 2 * Q66) * m ** 2 * n ** 2 + Q22 * n ** 4
            Qxy = (Q11 + Q22 - 4 * Q66) * m ** 2 * n ** 2 + Q12 * (m ** 4 + n ** 4)
            Qyy = Q11 * n ** 4 + 2 * (Q12 + 2 * Q66) * m ** 2 * n ** 2 + Q22 * m ** 4
            Qxs = (Q11 - Q12 - 2 * Q66) * n * m ** 3 + (Q12 - Q22 + 2 * Q66) * n ** 3 * m
            Qys = (Q11 - Q12 - 2 * Q66) * m * n ** 3 + (Q12 - Q22 + 2 * Q66) * m ** 3 * n
            Qss = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * n ** 2 * m ** 2 + Q66 * (m ** 4 + n ** 4)
            QMatrix_glo = np.array([[Qxx, Qxy, Qxs],
                                    [Qxy, Qyy, Qys],
                                    [Qxs, Qys, Qss]])
            self.QGlobalAr[i] = QMatrix_glo

    def calcABD(self):
        """
		Calculates the ABD matrix of the laminate
		:return:
		"""
        self.AMatrix = np.zeros((3, 3))
        self.BMatrix = np.zeros((3, 3))
        self.DMatrix = np.zeros((3, 3))
        # i and j are matrix components, k is the ply in the layup
        for k in range(len(self.LayUp)):
            Q = self.QGlobalAr[k]
            # calculating the A Matrix
            self.AMatrix += Q * (self.zlst[k + 1] - self.zlst[k])
            # calculating the B Matrix
            self.BMatrix += 0.5 * Q * (self.zlst[k + 1] ** 2 - self.zlst[k] ** 2)
            # calculating the D Matrix
            self.DMatrix += 3 ** -1 * Q * (self.zlst[k + 1] ** 3 - self.zlst[k] ** 3)
        # place the matrices together into the ABD matrix
        ABD_top = np.hstack((self.AMatrix, self.BMatrix))
        ABD_bottom = np.hstack((self.BMatrix, self.DMatrix))
        self.ABD = np.vstack((ABD_top, ABD_bottom))

    def calcGloStrainsNoCurve(self, Load):
        """
		Calculates the strain and curvature of the laminate at a prescribed load
		:param Load: A 6x1 numpy array column vector with the loads applied [Nx, Ny, Nz, Mx, My, Mz]
		:return: The deflections of the laminate
		"""
        return np.linalg.inv(self.ABD) @ Load

    def calcPlyStrains(self, Load):
        GloStrains = self.calcGloStrains(Load)
        PlyStrains = np.zeros((3, len(self.LayUp)))
        for k, theta in enumerate(self.LayUp):
            theta = np.deg2rad(theta)
            PlyStrains[:, k] = StrainGloTOPly(theta) @ GloStrains[:, k]
        return PlyStrains

    def calcGloStrains(self, Load):
        """
		Calculates global strains from a laminate and a load
		"""
        FlatStrains = np.linalg.inv(self.ABD) @ Load
        GloStrains = np.zeros((3, len(self.LayUp)))
        zAvg = [(self.zlst[k + 1] + self.zlst[k]) / 2 for k in range(len(self.zlst) - 1)]
        for k in range(len(self.LayUp)):
            GloStrains[:, k] = (FlatStrains[0:3] + zAvg[k] * FlatStrains[3:6]).T
        return GloStrains

    def calcPlyStresses(self, Load):
        """
		Calculates ply stresses from load
		"""
        gloStresses = self.calcGloStresses(Load)
        plyStresses = np.zeros((3, len(self.LayUp)))
        for k, theta in enumerate(self.LayUp):
            plyStresses[:, k] = StressGloTOPly(theta) @ gloStresses[:, k]
        return plyStresses

    def calcPlyStressesFromStrain(self, Strain):
        """
		Calculates ply stresses from load
		"""
        gloStresses = self.calcGloStressesFromStrain(Strain)
        plyStresses = np.zeros((3, len(self.LayUp)))
        for k, theta in enumerate(self.LayUp):
            plyStresses[:, k] = StressGloTOPly(theta) @ gloStresses[:, k]
        return plyStresses

    def calcGloStresses(self, Load):
        """
		Calculates global stresses from load
		"""
        gloStrains = self.calcGloStrains(Load)
        gloStresses = np.zeros((3, len(self.LayUp)))
        for k in range(len(self.LayUp)):
            gloStresses[:, k] = self.QGlobalAr[k] @ gloStrains[:, k]
        return gloStresses

    def calcGloStressesFromStrain(self, Strain):
        """
		Calculates global stresses from load
		"""
        gloStrains = Strain
        gloStresses = np.zeros((3, len(self.LayUp)))
        for k in range(len(self.LayUp)):
            gloStresses[:, k] = self.QGlobalAr[k] @ gloStrains[:]
        return gloStresses

    def calcPlyStresses2(self, Load):
        strain = self.calcPlyStrains(Load)
        plyStrain = np.zeros([3, len(self.LayUp)])
        plyStresses = np.zeros([3, len(self.LayUp)])
        for k in range(0, 4):
            plyStrain[:, k] = np.matmul(StrainGloTOPly(self.LayUp[k]), strain[:, k])
            plyStresses[:, k] = np.matmul(self.Lamina[k].QMatrix, plyStrain[:, k])

        return plyStresses

    def calcEngConst(self):
        """
		Calculates the engineering constants of the laminate
		:return: [Ex, Ey, vxy, vyx, Gxy]
		"""
        Axx = self.ABD[0, 0]
        Ayy = self.ABD[1, 1]
        Axy = self.ABD[0, 1]
        Ass = self.ABD[2, 2]

        Ex = (Axx * Ayy - Axy ** 2) / (self.h * Ayy)
        Ey = (Axx * Ayy - Axy ** 2) / (self.h * Axx)
        vxy = Axy / Ayy
        vyx = Axy / Axx
        Gxy = Ass / self.h
        return [Ex, Ey, vxy, vyx, Gxy]

    def calcStressEnvelope(self):
        pass

    def Puck(self, Load):  # Strength is list of ply properties: [Xt_mean,Xc_mean,Yt_mean,Yc_mean,S_mean]
        """
		Vectorized puck failure criteria calculator for given laminate and load
		"""
        Stress = self.calcPlyStresses(Load)
        pt12 = 0.3
        pc12 = 0.25
        pc11 = 0.225
        Xt = np.asarray([self.Lamina[k].Xt for k in range(len(self.LayUp))])
        Xc = np.asarray([self.Lamina[k].Xc for k in range(len(self.LayUp))])
        Yt = np.asarray([self.Lamina[k].Yt for k in range(len(self.LayUp))])
        Yc = np.asarray([self.Lamina[k].Yc for k in range(len(self.LayUp))])
        S = np.asarray([self.Lamina[k].S for k in range(len(self.LayUp))])
        N1 = Stress[0, :].T
        N2 = Stress[1, :].T
        N12 = Stress[2, :].T
        N12c = S * (1 + 2 * pc11) ** 0.5

        Ra = pc11 * S / pc12
        f_FFp, f_IFFp = np.zeros_like(Xt), np.zeros_like(Xt)
        # Mode C
        f_IFFp = ((N12 / (2 * (1 + pc11) * S)) ** 2 + (N2 / Yc) ** 2) * Yc / -N2
        mask = np.positive(N2 / N12) <= np.positive(Ra / N12c)  # Mode B
        f_IFFp[mask] = ((N12[mask] / S[mask]) ** 2 + (pc12 * N2[mask] / S[mask]) ** 2) ** 0.5 + pc12 * N2[mask] / S[
            mask]
        mask = N2 >= 0  # Mode A
        f_IFFp[mask] = (((1 / Yt[mask] - pt12 / S[mask]) * N2[mask]) ** 2 + (N12[mask] / S[mask]) ** 2) ** 0.5 + pt12 * \
                       N2[mask] / S[mask]
        f_FFp = -N1 / Xc
        mask = N1 >= 0  # FF Tension
        f_FFp[mask] = N1[mask] / Xt[mask]
        return f_FFp, f_IFFp  # result of analysis, if f_p is below 1 lamina did not fail, if it is 1 or higher lamina has failed

    def PuckStrain(self, Strain):  # Strength is list of ply properties: [Xt_mean,Xc_mean,Yt_mean,Yc_mean,S_mean]
        """
		Vectorized puck failure criteria calculator for given laminate and load
		"""
        Stress = self.calcPlyStressesFromStrain(Strain)
        pt12 = 0.3
        pc12 = 0.25
        pc11 = 0.225
        Xt = np.asarray([self.Lamina[k].Xt for k in range(len(self.LayUp))])
        Xc = np.asarray([self.Lamina[k].Xc for k in range(len(self.LayUp))])
        Yt = np.asarray([self.Lamina[k].Yt for k in range(len(self.LayUp))])
        Yc = np.asarray([self.Lamina[k].Yc for k in range(len(self.LayUp))])
        S = np.asarray([self.Lamina[k].S for k in range(len(self.LayUp))])
        N1 = Stress[0, :].T
        N2 = Stress[1, :].T
        N12 = Stress[2, :].T
        N12c = S * (1 + 2 * pc11) ** 0.5

        Ra = pc11 * S / pc12
        f_FFp, f_IFFp = np.zeros_like(Xt), np.zeros_like(Xt)
        # Mode C
        f_IFFp = ((N12 / (2 * (1 + pc11) * S)) ** 2 + (N2 / Yc) ** 2) * Yc / -N2
        mask = np.positive(N2 / N12) <= np.positive(Ra / N12c)  # Mode B
        f_IFFp[mask] = ((N12[mask] / S[mask]) ** 2 + (pc12 * N2[mask] / S[mask]) ** 2) ** 0.5 + pc12 * N2[mask] / S[
            mask]
        mask = N2 >= 0  # Mode A
        f_IFFp[mask] = (((1 / Yt[mask] - pt12 / S[mask]) * N2[mask]) ** 2 + (N12[mask] / S[mask]) ** 2) ** 0.5 + pt12 * \
                       N2[mask] / S[mask]
        f_FFp = -N1 / Xc
        mask = N1 >= 0  # FF Tension
        f_FFp[mask] = N1[mask] / Xt[mask]
        return f_FFp, f_IFFp  # result of analysis, if f_p is below 1 lamina did not fail, if it is 1 or higher lamina has failed

    def MaxStrain(self, Strain):  # Strength is list of ply properties: [Xt_mean,Xc_mean,Yt_mean,Yc_mean,S_mean]
        """
		Vectorized max stress failure criteria calculator for given laminate and load
		"""
        Stress = self.calcPlyStressesFromStrain(Strain)
        Xt = np.asarray([self.Lamina[k].Xt for k in range(len(self.LayUp))])
        Xc = np.asarray([self.Lamina[k].Xc for k in range(len(self.LayUp))])
        Yt = np.asarray([self.Lamina[k].Yt for k in range(len(self.LayUp))])
        Yc = np.asarray([self.Lamina[k].Yc for k in range(len(self.LayUp))])
        S = np.asarray([self.Lamina[k].S for k in range(len(self.LayUp))])
        N1 = Stress[0, :].T
        N2 = Stress[1, :].T
        N12 = Stress[2, :].T

        f_xm, f_ym, f_sm = np.zeros_like(Xt), np.zeros_like(Xt), np.zeros_like(Xt)
        # x max stress
        T, C = N1 >= 0, N1 < 0
        f_xm[T] = N1[T] / Xt[T]
        f_xm[C] = -N1[C] / Xc[C]
        # y max stress
        T, C = N2 >= 0, N2 < 0
        f_ym[T] = N2[T] / Yt[T]
        f_ym[C] = -N2[C] / Yc[C]
        # x max stress
        T, C = N12 >= 0, N12 < 0
        f_sm[T] = N12[T] / S[T]
        f_sm[C] = -N12[C] / S[C]

        return f_xm, f_ym, f_sm

    def calcFailure(self, Load, dL_step = 5000):
        """
		Failure load calculator, gives fpf and lpf load.
		"""
        failure = np.zeros(len(self.LayUp), dtype = int)
        lpf = False
        LoadFPF = np.zeros_like(Load)
        fpf = False
        LoadLPF = np.zeros_like(Load)
        initialLoad = np.zeros_like(Load)
        dL = dL_step * Load / np.linalg.norm(Load)
        Load = initialLoad
        strength = []
        for j in range(len(self.LayUp)):
            strength.append(
                [self.Lamina[j].Xt, self.Lamina[j].Xc, self.Lamina[j].Yt, self.Lamina[j].Yc, self.Lamina[j].S])
        strength = np.asarray(strength, dtype = float)
        strength = strength.reshape((-1, 5))
        while lpf == False:
            stresses = self.calcPlyStresses(Load)
            f_xm = 0
            f_ym = 0
            f_sm = 0

            # Apply max stress failure criteria
            for j in range(0, len(self.LayUp)):  # Max stress Failure
                if stresses[0, j] > 0 and strength[j, 0] != 0:
                    f_xm = stresses[0, j] / strength[j, 0]  # Xt_mean
                elif stresses[0, j] < 0 and strength[j, 0] != 0:
                    f_xm = stresses[0, j] / -strength[j, 1]  # Xc_mean

                if stresses[1, j] > 0 and strength[j, 0] != 0:
                    f_ym = stresses[1, j] / strength[j, 2]  # Yt_mean
                elif stresses[1, j] < 0 and strength[j, 0] != 0:
                    f_ym = stresses[1, j] / -strength[j, 3]  # Yc_mean
                if stresses[2, j] > 0 and strength[j, 0] != 0:
                    f_sm = stresses[2, j] / strength[j, 4]  # S_mean
                elif stresses[2, j] < 0 and strength[j, 0] != 0:
                    f_sm = stresses[2, j] / -strength[j, 4]

                # Determine which failure mode, and apply deg rule
                if (f_xm >= 1 and failure[j] <= 1) or ((f_ym >= 1 or f_sm >= 1) and failure[j] == 1):
                    failure[j] = failure[j] + 1
                    if fpf == False and np.sum(failure) >= 1:
                        LoadFPF = Load
                        fpf = True

                    self.Lamina[j] = Lamina(self.Lamina[j].t, 1e-8, 1e-8, 1e-8, 1e-8)
                    self.calcQGlobalLaminas()
                    self.calcABD()

                    if np.count_nonzero(failure) == len(self.LayUp):
                        LoadLPF = Load
                        lpf = True

                    break

                elif (f_ym >= 1 or f_sm >= 1) and failure[j] == 0:
                    failure[j] = failure[j] + 1
                    if fpf == False and np.sum(failure) >= 1:
                        LoadFPF = Load
                        fpf = True
                    # Degrade transverse elastic properties
                    t, E1, E2, v12, G12 = self.Lamina[j].t, self.Lamina[j].E1, self.Lamina[j].E2, self.Lamina[j].v12, \
                    self.Lamina[j].G12
                    self.Lamina[j] = Lamina(t, E1, E2 * 0.1, v12, G12)
                    self.calcQGlobalLaminas()
                    self.calcABD()
                    # Degrade strengths
                    strength[j, 2:4] = 0.1 * strength[j, 2:4]

                    if np.count_nonzero(failure) == len(self.LayUp):
                        LoadLPF = Load
                        lpf = True

                    break
                # When no failure detected, continue to next ply and increase
                elif j == len(self.LayUp) - 1:
                    Load = Load + dL
        return LoadFPF, LoadLPF

    def calcFailurePuck(self, Load, dL_step = 100000000):
        """
		Calculates the failure load just using puck's failure criteria
		"""
        Failed = False
        LoadI = np.zeros_like(Load)
        dL = Load * dL_step / np.linalg.norm(Load)

        while not Failed:
            f_FFp, f_IFFp = self.Puck(LoadI)
            if np.any(f_FFp > 1) or np.any(f_IFFp > 1):
                a = np.max([np.max(f_FFp), np.max(f_IFFp)])
                while not np.isclose(a, 1, atol = 0.01):
                    LoadI = LoadI / a
                    f_FFp, f_IFFp = self.Puck(LoadI)
                    a = np.max([np.max(f_FFp), np.max(f_IFFp)])
                Failed = True
                FailLoad = LoadI
                break
            LoadI += dL
        return FailLoad

    def __repr__(self):
        """
		String representation of the laminate class
		"""
        return f"Laminate of layup {self.LayUp}"


if __name__ == "__main__":
    # Simple test to see if a lamina turned 90 deg gives the same tensile stiffness
    E1 = 140 * 10 ** 9
    E2 = 10 * 10 ** 9
    G12 = 5 * 10 ** 9
    v12 = 0.3
    t = 0.125e-3

    v21 = v12 * E2 / E1
    Q = 1 - v12 * v21
    Q11 = E1 / Q
    Q22 = E2 / Q
    Q12 = v12 * E2 / Q
    Q66 = G12
    Lamina_ = Lamina(t, E1, E2, v12, G12)
    Laminate_1 = Laminate([15, 15, 15, 15], Lamina_)
    Laminate_2 = Laminate([105, 105, 105, 105], Lamina_)
    print(Laminate_1.ABD)
    print(Laminate_2.ABD)
    print(Laminate_1.calcEngConst())
    print(Laminate_2.calcEngConst())
    _, Ey, vxy, _, _ = Laminate_1.calcEngConst()
    Ex, _, _, _, _ = Laminate_2.calcEngConst()
    assert np.isclose(Ex, Ey)
    from AssignmentData import Xt_mean, Yt_mean, Xc_mean, Yc_mean, S_mean

    Lamina_.setStrengths(Xt_mean, Yt_mean, Xc_mean, Yc_mean, S_mean)

    LayUp = np.array([0, 90, 45, -45])
    LayUp = np.append(LayUp, np.flip(LayUp))
    LayUp = np.append(LayUp, np.flip(LayUp))
    Laminate_3 = Laminate(LayUp, Lamina_)

    LoadFPF, LoadLPF = Laminate_3.calcFailure(
        (np.array([np.cos(np.deg2rad(30)), np.sin(np.deg2rad(30)), 0, 0, 0, 0]) * 850).T)
    print("FPF", LoadFPF)
    print("LPF", LoadLPF)
