from pyVirtualPhysics import *
from Bvh import *
from SpatialMath import angleAxisToRotation

bodyIdxToBvhJoint = []  # type: list[BvhJoint]

def buildVpFromBvh(bvh):
    '''

    :type bvh: Bvh
    :return:
    '''

    bone_width = .15
    boxes_index = []
    num_boxes = bvh.getJointNum()
    bodies = []  #type: list[vpBody]

    offset_z = np.array((0., 0., 0.1))


    for i in range(num_boxes):
        joint = bvh.joints[i]
        children_size = len(joint.children)
        bodies.append(vpBody())

        if children_size == 1:
            pos = bvh.getJointPosition(-1, joint)
            pos_child = bvh.getJointPosition(-1, joint.children[0])
            bone_vec = pos_child - pos
            bone_pos = .5*(pos_child + pos)
            bone_pos += offset_z
            bone_length = npl.norm(bone_vec) * .9

            bone_size = np.array((bone_width, bone_width, bone_length))
            bone_mass = bone_size.prod() * 1000.

            unitz = np.array((0., 0., 1.))
            axis = np.cross(unitz, bone_vec)
            unit_axis = axis/npl.norm(axis)

            angle = math.acos(np.dot(unitz, bone_vec)/npl.norm(bone_vec))
            # bodies[i].setInitPosition(bone_pos)
            # bodies[i].setInitRotation(angleAxisToRotation(angle, unit_axis))
            bone_size_Vec3 = Vec3(bone_size[0]/2., bone_size[1]/2., bone_size[2]/2.)
            bodies[i].SetInertia(BoxInertia(1000., bone_size_Vec3))

            bodies[i].m_szName = joint.name
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

            # bodies[i].setInitPosition(joint_pos_ave)
            # bodies[i].setInitRotation(np.eye(3))
            bone_size_Vec3 = Vec3(bone_size[0]/2., bone_size[1]/2., bone_size[2]/2.)
            bodies[i].SetInertia(BoxInertia(1000., bone_size_Vec3))

            bodies[i].m_szName = joint.name

        if i != 0:
            joint_pos = bvh.getJointPosition(-1, bvh.joints[i]) + offset_z
            vpjoint = vpBJoint()
            bodies[i].SetJoint(vpjoint)
            bodies[bvh.getJointIndex(joint.parent.name)].SetJoint(vpjoint)

            # bodies[i].parentJoint.setInitPosition(joint_pos)
            # bodies[i].parentJoint.setInitRotation(np.eye(3))

    humanoid = vpWorld()
    humanoid.AddBody(bodies[0])
    humanoid.SetGravity(Vec3(0., -9.81, 0.))
    humanoid.Initialize()
    print(humanoid.GetNumBody())

    for i in range(num_boxes):
        for j in range(num_boxes):
            if bvh.joints[j].name == humanoid.GetBody(i).m_szName:
                bodyIdxToBvhJoint.append(bvh.joints[j])
                break
    return humanoid



bvh = Bvh('test.bvh', 0.02)
bvh.rot = rotx(math.pi/2.)

world = buildVpFromBvh(bvh)
world.InverseDynamics()

