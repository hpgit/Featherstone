import numpy as np
import math


def rotz(rad):
    E = np.eye(3)
    E[:2, :2] = np.array([[math.cos(rad), -math.sin(rad)],[math.sin(rad), math.cos(rad)]])
    return E


def rotx(rad):
    E = np.eye(3)
    E[1:, 1:] = np.array([[math.cos(rad), -math.sin(rad)],[math.sin(rad), math.cos(rad)]])
    return E


def roty(rad):
    E = np.eye(3)
    E[0, 0] = math.cos(rad)
    E[0, 2] = math.sin(rad)
    E[2, 0] = -math.sin(rad)
    E[2, 2] = math.cos(rad)

    return E

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
    X[:3, :3] = E.T
    X[3:, 3:] = E.T
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
