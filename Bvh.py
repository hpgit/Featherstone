import numpy as np

class BvhJoint():
    def __init__(self, name):
        self.name = name
        self.hasEndChild = False
        self.offsetFromParent = np.zeros(3)
        self.dof = 0
        self.parent = None  # type: BvhJoint
        self.children = []  # type: list[BvhJoint]

class Bvh():
    def __init__(self, filename=None):
        self.rot = np.eye(3)
        if filename is not None:
            self.constructFromBvhFile(filename)

        self.root = None  # type: BvhJoint
        self.joints = []  # type: list[BvhJoint]
        self.jointNames = []  # type: list[str]
        self.jointIdxToChannelStartIndex  = []  # type: list[int]
        self.channelIndexToJoint = []   # type: list[BvhJoint]
        self.channels = [] #type: list[list[float]]

        self.numFrames = 0
        self.frameTime = 0.

    def getJointNum(self):
        return len(self.joints)

    def getJoint(self, idx):
        return self.joints[idx]

    def constructFromBvhFile(self, filename, ratio=1.):
        file = open(filename)
        text = file.read()
        file.close()

        s = text.split()

        channelBeginIdx = 0
        stack = []  # type: list[BvhJoint]

        # HIERARCHY
        s.pop(0)

        while len(s) > 0 :
            temp = s.pop(0)
            if temp == 'ROOT' or temp == 'JOINT':
                isRoot = (temp == 'ROOT')

                # get joint name
                joint = BvhJoint(s.pop(0))
                if not isRoot:
                    joint.parent = stack[-1]
                    stack[-1].children.append(joint)
                else:
                    joint.parent = None
                    self.root = joint

                stack.append(joint)
                self.joints.append(joint)
                self.jointNames.append(joint.name)
                self.jointIdxToChannelStartIndex.append(channelBeginIdx)

                # offset
                s.pop(0) # {
                s.pop(0) # OFFSET
                joint.offsetFromParent = ratio * np.array([float(s.pop(0)), float(s.pop(0)), float(s.pop(0))])

                # channels
                s.pop(0) # CHANNELS
                joint.dof = int(s.pop(0)) # dof
                for i in range(joint.dof):
                    s.pop(0)

                channelBeginIdx += joint.dof
                for i in range(joint.dof):
                    self.channelIndexToJoint.append(joint)

            elif temp == '}':
                stack.pop()

            elif temp == 'End':
                s.pop(0) # Site
                s.pop(0) # {
                joint = BvhJoint('End_'+stack[-1].name)

                joint.parent = stack[-1]
                joint.parent.children.append(joint)
                joint.parent.hasEndChild = True

                s.pop(0) # OFFSET
                joint.offsetFromParent = ratio * np.array([float(s.pop(0)), float(s.pop(0)), float(s.pop(0))])

                s.pop(0) # }

            elif temp == 'MOTION':
                break

            else:
                print('BVH File Error!')

        s.pop(0)  # Frames:
        self.numFrames = int(s.pop(0))  # # of frame

        s.pop(0)  # Frame
        s.pop(0)  # Time:
        self.frameTime = float(s.pop(0))  # per-frame time

        for i in range(self.numFrames):
            self.channels.append([])
            for j in range(3):
                self.channels[i].append(ratio * float(s.pop(0)))
            for j in range(3, channelBeginIdx):
                self.channels[i].append(ratio * float(s.pop(0)))

    def getJointIndex(self, name):
        for i in range(len(self.joints)):
            if name == self.joints[i].name:
                return i

        return -1

    def getJointChannelIndex(self, joint):
        return self.jointIdxToChannelStartIndex[self.joints.index(joint)]

    def getJointPosition(self, frame, joint):
        '''

        :type frame: int
        :type joint: BvhJoint
        :return:
        '''
        pos = np.zeros(3)
        if joint.parent is None:
            pos = joint.offsetFromParent.copy()
            if frame >= 0:
                pos += np.array(self.channels[frame][:3])
            return np.dot(self.rot, pos)

        else:
            parent_rotation = np.eye(3)
            if frame >= 0:
                parent_rotation = self.getJointGlobalRotationAfter(frame, joint.parent)

            parent_pos = self.getJointPosition(frame, joint.parent)
            pos_temp = np.dot(self.rot, joint.offsetFromParent)
            pos_temp1 = np.dot(parent_rotation, pos_temp)
            return pos_temp1 + parent_pos


    def getJointGlobalRotationAfter(self, frame, joint):
        return np.dot(self.rot, np.dot(self.getOriginalJointGlobalRotationAfter(frame, joint), self.rot.T))

    def getBonePosition(self, frame, jointIdx):
        joint = self.joints[jointIdx]

        pass

    def getBoneLength(self, jointIdx):
        pass

    def getJointLocalRotationAtPredecessor(self, frame, joint):
        pass

    def roate(self, rotation):
        self.rot = rotation

    def getOriginalJointLocalRotationAtPredecessor(self, frame, joint):
        pass

    def getOriginalJointGlobalRotationAfter(self, frame, joint):
        pass