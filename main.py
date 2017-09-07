from hpPhysics import JointType, Body, Joint, Skeleton, BoxGeom, BoxInertia
from Bvh import Bvh, BvhJoint
from SpatialMath import *
import numpy.linalg as npl

np.set_printoptions(precision=6, suppress=True, threshold=np.inf)

DEBUG = False

def importVectorFromFile(filename):
    file = open(filename)
    s = file.read().split()
    file.close()
    return np.array([float(d) for d in s])


bodyIdxToBvhJoint = []  # type: list[BvhJoint]


def buildHumanoidFromBvh(bvh):
    '''

    :type bvh: Bvh
    :return:
    '''

    bone_width = .15
    boxes_index = []
    num_boxes = bvh.getJointNum()

    bodies = []  # type: list[Body]
    geoms = []   # type: list[BoxGeom]

    offset_z = np.array((0., 0., 0.1))


    for i in range(num_boxes):
        joint = bvh.joints[i]
        children_size = len(joint.children)
        bodies.append(Body(joint.name))

        if children_size == 1:
            pos = bvh.getJointPosition(-1, joint)
            pos_child = bvh.getJointPosition(-1, joint.children[0])
            bone_vec = pos_child - pos
            bone_pos = .5*(pos_child + pos)
            bone_pos += offset_z

            bone_length = npl.norm(bone_vec) * .9

            if DEBUG:
                print('name: {}, length: {}'.format(joint.name, bone_length))
                print('jointpos: '+str(pos))
                print('jointchildpos: '+str(pos_child))

            bone_size = np.array((bone_width, bone_width, bone_length))
            bone_mass = bone_width * bone_width * bone_length * 1000.

            bodies[i].setInitPosition(bone_pos)

            unitz = np.array((0., 0., 1.))
            axis = np.cross(unitz, bone_vec)
            unit_axis = axis/npl.norm(axis)

            angle = math.acos(np.dot(unitz, bone_vec)/npl.norm(bone_vec))
            bodies[i].setInitRotation(angleAxisToRotation(angle, unit_axis))
            bodies[i].setInertia(BoxInertia(1000., bone_size/2.))
            geoms.append(BoxGeom(bone_size[0], bone_size[1], bone_size[2]))
            bodies[i].setGeom(geoms[i])

        else:
            joint_pos_ave = np.zeros(3)
            joint_pos_max = np.zeros(3)
            joint_pos_min = np.zeros(3)

            for j in range(children_size):
                joint_pos = bvh.getJointPosition(-1, joint.children[j])
                joint_pos_ave += joint_pos
                if j == 0:
                    joint_pos_max = joint_pos.copy()
                    joint_pos_min = joint_pos.copy()
                else:
                    for k in range(3):
                        joint_pos_max[k] = max([joint_pos_max[k], joint_pos[k]])
                        joint_pos_min[k] = min([joint_pos_min[k], joint_pos[k]])

            bone_size = joint_pos_max - joint_pos_min
            for k in range(3):
                bone_size[k] = max([bone_size[k], bone_width])

            joint_pos_ave /= children_size
            joint_pos_ave += offset_z

            bodies[i].setInitPosition(joint_pos_ave)
            bodies[i].setInitRotation(np.eye(3))
            bodies[i].setInertia(BoxInertia(1000., bone_size/2.))
            geoms.append(BoxGeom(bone_size[0], bone_size[1], bone_size[2]))
            bodies[i].setGeom(geoms[i])

        if DEBUG:
            print('name: ' + bodies[i].name +', size: '+str(bodies[i].geom.size) + ',  inertia: \n'+str(bodies[i].I))
            print('pos: ' + str(bodies[i].initr))
            print('rot: ' + str(bodies[i].initE))

        if i == 0:
            bodies[i].parentJoint = Joint(JointType.FreeJoint)
            # bodies[i].parentJoint.setInitPosition(bvh.getJointPosition(-1, bvh.joints[0]))
            bodies[i].parentJoint.setInitPosition(bodies[i].initr)
            bodies[i].parentJoint.parentBody = None
            bodies[i].parentJoint.childBody = bodies[i]
            bodies[i].parentJoint.name = bvh.joints[0].name

        else:
            joint_pos = bvh.getJointPosition(-1, bvh.joints[i]) + offset_z
            bodies[i].parentJoint = Joint(JointType.BallJoint)
            bodies[i].parentJoint.childBody = bodies[i]
            bodies[i].parentJoint.setInitPosition(joint_pos)
            bodies[i].parentJoint.setInitRotation(np.eye(3))
            bodies[i].parentJoint.name = joint.name
            bodies[i].parentBody = bodies[bvh.getJointIndex(joint.parent.name)]
            bodies[i].parentBody.children.append(bodies[i])
            bodies[i].parentJoint.parentBody = bodies[i].parentBody

    humanoid = Skeleton()
    humanoid.contstruct(bodies[0])

    for i in range(num_boxes):
        for j in range(num_boxes):
            if bvh.joints[j].name == humanoid.bodies[i].name:
                bodyIdxToBvhJoint.append(bvh.joints[j])
                break
    return humanoid

bvh = Bvh('test.bvh', 0.02)
bvh.rot = rotx(math.pi/2.)

humanoid = buildHumanoidFromBvh(bvh)
dof = humanoid.dof
numBody = len(humanoid.bodies)


qidx = humanoid.getJointQidxByName('LeftToes')

dq = np.zeros(dof)
dq[qidx+1] = 100.
q = np.zeros(dof)
q[0] = 0.

C = humanoid.calcBiasForces(q=q, dq=dq, f_ext=np.zeros(numBody*6)).round(decimals=6)
Ctrue = importVectorFromFile('vector.txt')

nameSet = [
    'Hips',
    'LeftUpLeg',
    'Spine',
    'RightUpLeg',
    'RightLeg',
    'RightFoot',
    'RightToes',
    'Spine1',
    'RightShoulder',
    'LeftShoulder1',
    'Head',
    'LeftArm',
    'LeftForeArm',
    'LeftHand',
    'RightArm',
    'RightForeArm',
    'RightHand',
    'LeftLeg',
    'LeftFoot',
    'LeftToes']

for i in range(len(C)):
    if (i >= 6 and i % 3 == 0) or i==0:
        nameIdx = int(i/3)
        if i>=6:
            nameIdx -= 1
        print('\n', nameSet[nameIdx])

    print(C[i], Ctrue[i])

print('\n\n', npl.norm(C-Ctrue))

# humanoid.calcMassMatrix(q=np.zeros(dof))


