import numpy as np
import numpy.linalg as npl
from enum import Enum
import math
from queue import Queue

def crossMat(r):
    rx = np.zeros((3, 3))
    rx[0, 1] = -r[2]
    rx[0, 2] = r[1]
    rx[1, 2] = -r[0]

    rx[1, 0] = r[2]
    rx[2, 0] = -r[1]
    rx[2, 1] = r[0]

    return rx

def angleAxisToRotation(angle, axis):
    cos = math.cos(angle)
    sin = math.sin(angle)
    axis_mat = np.reshape(axis, (3, 1))
    return cos * np.eye(3) + sin*crossMat(axis) + (1.-cos)*np.dot(axis_mat, axis_mat.T)

def dual(_X):
    X = np.zeros((6, 6))
    X[:3, :3] = _X[:3, :3]
    X[:3, 3:] = _X[3:, :3]
    X[3:, :3] = _X[:3, 3:]
    X[3:, 3:] = _X[3:, 3:]

    return X

def inv(_X):
    X = np.zeros((6, 6))
    X[:3, :3] = _X[:3, :3].T
    X[:3, 3:] = _X[:3, 3:].T
    X[3:, :3] = _X[3:, :3].T
    X[3:, 3:] = _X[3:, 3:].T

    return X

def transf(E, r):
    X = np.zeros((6, 6))
    X[:3, :3] = E
    X[3:, 3:] = E
    X[3:, :3] = -np.dot(E, crossMat(r))

    return X

def dtransf(E, r):
    X = np.zeros((6, 6))
    X[:3, :3] = E
    X[3:, 3:] = E
    X[:3, 3:] = -np.dot(E, crossMat(r))

    return X

def invtransf(E, r):
    X = np.zeros((6, 6))
    X[:3, :3] = E.T
    X[3:, 3:] = E.T
    X[3:, :3] = np.dot(crossMat(r), E.T)

    return X

def dinvtransf(E, r):
    X = np.zeros((6, 6))
    X[:3, :3] = E
    X[3:, 3:] = E
    X[:3, 3:] = np.dot(crossMat(r), E.T)

    return X

def rot(E):
    X = np.zeros((6, 6))
    X[:3, :3] = E
    X[3:, 3:] = E

    return X

def xlt(r):
    X = np.eye(6)
    X[3:, :3] = -crossMat(r)

def spatialCross(v):
    X = np.zeros((6, 6))
    X[:3, :3] = crossMat(v[:3])
    X[3:, :3] = crossMat(v[3:])
    X[3:, 3:] = crossMat(v[:3])

    return X

def spatialdCross(v):
    X = np.zeros((6, 6))
    X[:3, :3] = crossMat(v[:3])
    X[:3, 3:] = crossMat(v[3:])
    X[3:, 3:] = crossMat(v[:3])

    return X

def BoxInertia(density, size):
    mass = 8.0 * density * size[0] * size[1] * size[2]
    ix = mass * (size[1] * size[1] + size[2] * size[2]) / 3.
    iy = mass * (size[0] * size[0] + size[2] * size[2]) / 3.
    iz = mass * (size[0] * size[0] + size[1] * size[1]) / 3.

    I = np.zeros((6, 6))
    I[:3, :3] = np.diag(np.array((ix, iy, iz)))
    I[3:, 3:] = mass * np.eye(3)
    return I

def SphereInertia(density, rad):
    rad_2 = rad * rad
    mass = density * math.pi * rad_2
    i = 0.4 * mass * rad_2
    I = np.zeros((6, 6))
    I[:3, :3] = np.diag(np.array((i, i, i)))
    I[3:, 3:] = mass * np.eye(3)
    return I

def CylinderInertia(density, rad, height):
    rad_2 = rad * rad
    mass = density * math.pi * rad_2 * height
    ix = mass * height * height  / 12.0 + .25 * mass * rad_2
    iy = ix
    iz = .5* mass * rad_2

    I = np.zeros((6, 6))
    I[:3, :3] = np.diag(np.array((ix, iy, iz)))
    I[3:, 3:] = mass * np.eye(3)
    return I

def TorusInertia(density, ring_rad, tube_rad):
    mass = density * 2. * (math.pi*math.pi) * ring_rad * tube_rad * tube_rad
    ix = mass * (0.625 * tube_rad * tube_rad + .5 * ring_rad + ring_rad)
    iy = ix
    iz = mass * (0.75 * tube_rad * tube_rad + ring_rad + ring_rad)

    I = np.zeros((6, 6))
    I[:3, :3] = np.diag(np.array((ix, iy, iz)))
    I[3:, 3:] = mass * np.eye(3)
    return I

class JointType(Enum):
    NullJoint = 0
    FixedJoint = 1
    RevoluteJoint = 2
    UniversalJoint = 3
    BallJoint = 4
    PrismaticJoint = 5
    FreeJoint = 6
    PlanarJoint = 7

    Length = 8


class BoxGeom():
    def __init__(self, x, y, z):
        self.size = np.array([x, y, z])


class Body():
    def __init__(self):
        self.parentJoint = None # type : Joint
        self.parentBody = None  # type : Body
        self.children = []      # type : list[Body]
        self.geom = None        # type : BoxGeom

        self.I = np.zeros((6, 6))
        self.X_T = np.zeros((6, 6))  # parent body to parent joint on parent body coordinates

        self.initE = np.eye(3)
        self.initr = np.zeros(3)
        self.com = np.zeros(3)

        # internal variables
        self._X_p = np.zeros((6, 6))  # parent body to this body on parent body coordinates
        self._X_0 = np.zeros((6, 6))   # global transform
        self._a_vp = np.zeros(6)
        self._v = np.zeros(6)
        self._f = np.zeros(6)
        self._C = np.zeros(6)

        self._F = np.zeros(6)
        self._I_c = np.zeros((6, 6))

    def setInitPosition(self, r):
        self.initr = r

    def setInitRotation(self, E):
        self.initE = E

    def setInertia(self, I):
        self.I = I


class Joint():
    def __init__(self, jointType, axis=list()):
        # assert jointType != JointType.NullJoint and jointType < JointType.Length
        assert JointType.NullJoint <= jointType < JointType.Length
        self.jointType = jointType
        self.S = None  # type: np.array
        self.dof = 0

        self.initE = np.eye(3)

        if jointType == JointType.FixedJoint:
            raise NotImplementedError
        elif jointType == JointType.PlanarJoint:
            raise NotImplementedError
        elif jointType == JointType.UniversalJoint:
            raise NotImplementedError
        elif jointType == JointType.RevoluteJoint:
            self.S = np.zeros((6, 1))
            self.S[:3, 0] = axis
            self.dof = 1
        elif jointType == JointType.BallJoint:
            self.S = np.zeros((6, 3))
            self.S[:3, :3] = np.eye(3)
            self.dof = 3
        elif jointType == JointType.PrismaticJoint:
            self.S = np.zeros((6, 1))
            self.S[3:, 0] = axis
            self.dof = 1
        elif jointType == JointType.FreeJoint:
            self.S = np.eye(6)
            self.dof = 6

    def setInitRotation(self, _E):
        self.initE = _E

    def getTransform(self, q):
        jointType = self.jointType
        E = np.eye(3)
        r = np.zeros(3)
        if jointType == JointType.FixedJoint:
            raise NotImplementedError
        elif jointType == JointType.PlanarJoint:
            raise NotImplementedError
        elif jointType == JointType.UniversalJoint:
            raise NotImplementedError
        elif jointType == JointType.RevoluteJoint:
            E = angleAxisToRotation(q, self.S[:3, 0].reshape(3))
        elif jointType == JointType.BallJoint:
            norm = npl.norm(q)
            E = angleAxisToRotation(norm, q/norm)
        elif jointType == JointType.PrismaticJoint:
            r = q * self.S[3:, 0].reshape(3)
        elif jointType == JointType.FreeJoint:
            norm = npl.norm(q[:3])
            E = angleAxisToRotation(norm, q[:3]/norm)
            r = q[3:]

        return E, r

    def getdSdq(self):
        pass

class Skeleton():
    def __init__(self):
        self.bodies = []    # type: list[Body]
        self.joints = []    # type: list[Joint]
        self.jointqidx = [] # type: list[int]
        self.lamb = []      # type: list[int]
        self.q = np.zeros(0)
        self.dq = np.zeros(0)

        self.gravity = [0., -9.81, 0., 0., 0., 0.]

        self.dof = 1

    def dof(self):
        return self.dof

    def contstruct(self, rootBody):
        '''
        :type rootBody: Body
        :return:
        '''
        Q = Queue()
        Q.put(rootBody)
        qIdx = 0
        while not Q.empty():
            body = Q.get()
            self.bodies.append(body)
            self.joints.append(body.parentJoint)
            self.jointqidx.append(qIdx)
            qIdx += body.parentJoint.dof
            for i in range(len(body.children)):
                Q.put(body.children[i])

        self.dof = qIdx

        self.lamb.append(0)
        for i in range(1, len(self.bodies)):
            parentIdx = self.bodies.index(self.bodies[i].parentBody)
            self.lamb.append(parentIdx)

        self.q = np.zeros(qIdx)
        self.dq = np.zeros(qIdx)
        self.ddq = np.zeros(qIdx)

        self.bodies[0].X_T = np.eye(6)

    def jcalc(self, jointIdx, q, dq):
        # return X_J, S_i, v_J, c_J
        joint = self.joints[jointIdx]
        X_J = transf(*joint.getTransform(q))
        return X_J, joint.S, np.dot(joint.S, dq), np.zeros(6)

    def jcalcPos(self, jointIdx, q):
        # return X_J, S_i
        joint = self.joints[jointIdx]
        X_J = transf(*joint.getTransform(q))
        return X_J, joint.S

    def calcBiasForces(self, q, dq, f_ext):
        bodies = self.bodies
        joints = self.joints

        '''
        a_vp_0 = -0_a_g

        for i in range(1, len(bodies)):
            X_J, S_i, v_J, c_J = jcalc(i, q_i, dq_i)
            i_X_p = np.dot(X_J, X_T_i)
            if self.lamb[i] != 0:
                i_X_0 = np.dot(i_X_p, p_X_0)
            v_i = np.dot(i_X_p, v_p) + v_J
            a_vp_i = np.dot(i_X_p, a_vp_p) + c_J + np.dot(spatialCross(v_i), v_J)
            f_i = np.dot(I_i, a_vp_i) + np.dot(spatialdCross(v_i), np.dot(I_i, v_i)) - np.dot(i_Xd_0, 0_f_x_i)
        f_0 = np.dot(I_0, a_vp_0) + np.dot(spatialdCross(v_0), np.dot(I_0, v_0)) - 0_f_x_0
        for i in range(len(bodies), 0, -1):
            C_i = np.dot(S_i.T, f_i)
            f_p = f_p + np.dot(p_Xd_i, f_i)
        p_c_0 = f_0
        '''
        self.bodies[0]._a_vp = self.gravity.copy()

        for i in range(1, len(bodies)):
            body = self.bodies[i]
            joint = self.joints[i]
            jointDof = joint.dof
            jointqIdx = self.jointqidx[i]

            X_J, S_i, v_J, c_J = self.jcalc(i, q[jointqIdx:jointqIdx+jointDof], dq[jointqIdx:jointqIdx+jointDof])
            body._X_p = np.dot(X_J, body.X_T)
            if self.lamb[i] != 0:
                body.X_0 = np.dot(body._X_p, body.parentBody._X_0)
            body._v = np.dot(body._X_p, body._v) + v_J
            body._a_vp = np.dot(body._X_p, body.parentBody._a_vp) + c_J + np.dot(spatialCross(body._v), v_J)
            body._f = np.dot(body.I, body._a_vp) + np.dot(spatialdCross(body._v), np.dot(body.I, body._v)) - np.dot(dual(body.X_0), f_ext[6*i:6*i+6])
        body = self.bodies[0]
        body._f = np.dot(body.I, body._a_vp) + np.dot(spatialdCross(body._v), np.dot(body.I, body._v)) - f_ext[:6]
        for i in range(len(bodies), 0, -1):
            body = self.bodies[i]
            body._C = np.dot(body.parentJoint.S, body._f)
            body.parentBody._f = body.parentBody._f + np.dot(dual(inv(body._X_p)), body._f)
        p_c_0 = self.bodies[0]._f

    def calcMassMatrix(self, q):
        H = np.zeros(self.dof, self.dof)
        '''
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
        '''

        for i in range(len(self.bodies)):
            self.bodies[i]._I_c = self.bodies[i].I.copy()

            body = self.bodies[i]
            joint = self.joints[i]
            dof = joint.dof
            qIdx = self.jointqidx[i]
            X_J, S_i = self.jcalcPos(i, q[qIdx:qIdx+dof])
            body._X_p = np.dot(X_J, body.X_T)
            if self.lamb[i] != 0:
                body.X_0 = np.dot(body._X_p, body.parentBody._X_0)
            else:
                body.X_0 = body._X_p.copy()

        for i in range(len(self.bodies), 0, -1):
            body = self.bodies[i]
            joint = self.joints[i]
            dof = joint.dof
            qIdx = self.jointqidx[i]

            body.parentBody._I_c = body.parentBody._I_c + np.dot(dual(inv(body._X_p)), np.dot(body._I_c, body._X_p))
            body._F = np.dot(body._I_c, joint.S)
            H[qIdx:qIdx+dof, qIdx:qIdx+dof] = np.dot(joint.S.T, body._F)
            j = i
            while self.lamb[j] != 0:
                bodyj = self.bodies[j]
                jointj = self.joints[j]
                body._F = np.dot(dual(inv(bodyj._X_p)), body._F)
                j = self.lamb[j]
                jdof = joint.dof
                jqIdx = self.jointqidx[i]
                H[qIdx:qIdx+dof, jqIdx:jqIdx+jdof] = np.dot(body._F.T, jointj.S)
                H[jqIdx:jqIdx+jdof, qIdx:qIdx+dof] = H[qIdx:qIdx+dof, jqIdx:jqIdx+jdof].T
            bodyj = self.bodies[j]
            body._F = np.dot(dual(inv(bodyj._X_0)), body._F)
