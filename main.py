import numpy as np

def crossMat(r):
    rx = np.zeros(3, 3)
    rx[0, 1] = -r[2]
    rx[0, 2] = r[1]
    rx[1, 2] = -r[0]
    
    rx[1, 0] = r[2]
    rx[2, 0] = -r[1]
    rx[2, 1] = r[0]

    return rx

def transf(E, r):
    X = np.zeros(6, 6)
    X[:3, :3] = E
    X[3:, 3:] = E
    X[3:, :3] = -np.dot(E, crossMat(r))

    return X

def dtransf(E, r):
    X = np.zeros(6, 6)
    X[:3, :3] = E
    X[3:, 3:] = E
    X[:3, 3:] = -np.dot(E, crossMat(r))

    return X

def invtransf(E, r):
    X = np.zeros(6, 6)
    X[:3, :3] = E.T
    X[3:, 3:] = E.T
    X[3:, :3] = np.dot(crossMat(r), E.T)

    return X

def dinvtransf(E, r):
    X = np.zeros(6, 6)
    X[:3, :3] = E
    X[3:, 3:] = E
    X[:3, 3:] = np.dot(crossMat(r), E.T)

    return X

def rot(E):
    X = np.zeros(6, 6)
    X[:3, :3] = E
    X[3:. 3:] = E

    return X

def xlt(r):
    X = np.eye(6)
    X[3:, :3] = -crossMat(r)

def spatialCross(v):
    X = np.zeros(6, 6)
    X[:3, :3] = crossMat(v[:3])
    X[3:, :3] = crossMat(v[3:])
    X[3:, 3:] = crossMat(v[:3])

    return X

def spatialdCross(v):
    X = np.zeros(6, 6)
    X[:3, :3] = crossMat(v[:3])
    X[:3, 3:] = crossMat(v[3:])
    X[3:, 3:] = crossMat(v[:3])

    return X


class Body():
    def __init__(self):
        self.parentJoint = Joint()
        self.parentBody = Body()

        self.I = np.zeros(6, 6)

class Joint():
    def __init__(self):
        pass

class skeleton():
    def __init__(self):
        self.bodies = []
        self.joints = []
        self.lamb = []

    def dof(self):
        return 1

    def calcBiasForces(self, q, dq):
        bodies = self.bodies
        joints = self.joints
        a_vp_0 = -0_a_g

        for i in range(1, len(bodies)):
            X_J, S_i, v_J, c_J = jcalc(jtype(i), q_i, dq_i)
            i_X_p = np.dot(X_J, X_T_i)
            if self.lamb[i] != 0:
                i_X_0 = np.dot(i_X_p, p_X_0)
            v_i = np.dot(i_X_p, v_p) + v_J
            a_vp_i = np.dot(i_X_p, a_vp_p) + c_J + np.dot(spatialCross(v_i), v_J)
            f_i = np.dot(I_i, a_vp_i) + np.dot(spatialdCross(v_i), np.dot(I_i, v_i)) - np.dot(i_Xd_0, 0_f_x_i)
        f_0 = np.dot(I_0, a_vp_0) + np.dot(spatialdCross(v_0), np.dot(I_0, v_0)) - 0_f_x_0
        for i = range(len(bodies), 0, -1):
            C_i = np.dot(S_i.T, f_i)
            f_p = f_p + np.dot(p_Xd_i, f_i)
        p_c_0 = f_0

    def calcMassMatrix(self):
        H = np.zeros(dof, dof)
        for i in range(len(bodies)):
            I_c_i = I_i
        for i in range(len(bodies), 0, -1):
            I_c_p = I_c_p + np.dot(p_Xd_i, np.dot(I_c_i, i_X_p))
            F_i = np.dot(I_c_i, S_i)
            H_ii = np.dot(S_i.T, F_i)
            j=i
            while self.lamb[j] != 0:
                F_i = np.dot(p_Xd_j, F_i)
                j = self.lamb[j]
                H_ij = np.dot(F_i.T, S_j)
                H_ji = H_ij.T
            F_i = np.dot(0_Xd_j, F_i)
