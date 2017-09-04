#include "myBvh.h"

Bvh::~Bvh()
{
	channelIndexToJoint.clear();
	while (joints.size() > 0)
	{
		BvhJoint *joint = joints.back();
		joints.pop_back();
		delete joint;
	}
}

void Bvh::constructFromBvhFile(const char *filename, sm::real ratio)
{
	FILE *fp = fopen(filename, "r");
	char temp[50];

	// Bvh bvh;
	std::vector<BvhJoint*> stack;
	
	// bvh parsing begin
	fscanf_s(fp, "%s", temp, sizeof(temp)); // HIERARCHY

	int channelBeginIdx = 0;

	while (fscanf_s(fp, "%s", temp, sizeof(temp)))
	{
		if (!strcmp(temp, "ROOT") || !strcmp(temp, "JOINT"))
		{
			bool isRoot = !strcmp(temp, "ROOT");

			// get joint info
			fscanf_s(fp, "%s", temp, sizeof(temp));  // joint name
			BvhJoint *joint = new BvhJoint(temp);
			if (!isRoot)
			{
				joint->parent = stack.back();
				stack.back()->children.push_back(joint);
			}
			else
			{
				joint->parent = nullptr;
				this->root = joint;
			}
			
			stack.push_back(joint);
			this->joints.push_back(joint);
			this->jointNames.push_back(joint->name);
			this->JointIdxToChannelStartIndex.push_back(channelBeginIdx);


			fscanf_s(fp, "%s", temp, sizeof(temp));  // {

			// offset
			fscanf_s(fp, "%s", temp, sizeof(temp));  // OFFSET
			fscanf_s(fp, "%s", temp, sizeof(temp));  // value of offset_x
			joint->offsetFromParent[0] = atof(temp)*ratio;
			fscanf_s(fp, "%s", temp, sizeof(temp));  // value of offset_y
			joint->offsetFromParent[1] = atof(temp)*ratio;
			fscanf_s(fp, "%s", temp, sizeof(temp));  // value of offset_z
			joint->offsetFromParent[2] = atof(temp)*ratio;

			// channels
			fscanf_s(fp, "%s", temp, sizeof(temp));  // CHANNELS
			fscanf_s(fp, "%s", temp, sizeof(temp));  // value of channels ( dof )
			joint->dof = atoi(temp);
			for (auto i = 0; i < joint->dof; i++)
				fscanf_s(fp, "%s", temp, sizeof(temp));  // position or rotation
			
			channelBeginIdx += joint->dof;

			for (auto i = 0; i < joint->dof; i++)
				this->channelIndexToJoint.push_back(joint);
		}
		else if (!strcmp(temp, "}"))
		{
			// parent node end
			stack.pop_back();
		}
		else if (!strcmp(temp, "End"))
		{
			// leaf node reached
			fscanf_s(fp, "%s", temp, sizeof(temp));  // Site
			fscanf_s(fp, "%s", temp, sizeof(temp));  // {
			char endname[50] = "End_";
			for (auto i = 0; i < stack.back()->name.length(); i++)
				endname[4 + i] = stack.back()->name.c_str()[i];
			endname[4 + stack.back()->name.length()] = NULL;
			BvhJoint *joint = new BvhJoint(endname);

			joint->parent = stack.back();
			stack.back()->children.push_back(joint);
			stack.back()->hasEndChild = true;

			// offset
			fscanf_s(fp, "%s", temp, sizeof(temp));  // OFFSET
			fscanf_s(fp, "%s", temp, sizeof(temp));  // value of offset_x
			joint->offsetFromParent[0] = atof(temp)*ratio;
			fscanf_s(fp, "%s", temp, sizeof(temp));  // value of offset_y
			joint->offsetFromParent[1] = atof(temp)*ratio;
			fscanf_s(fp, "%s", temp, sizeof(temp));  // value of offset_z
			joint->offsetFromParent[2] = atof(temp)*ratio;
			fscanf_s(fp, "%s", temp, sizeof(temp));  // }
		}
		else if (!strcmp(temp, "MOTION"))
		{
			break;
		}
		else
		{
			printf("BVH File Error!");
		}
	}

	fscanf_s(fp, "%s", temp, sizeof(temp));  // Frames:
	fscanf_s(fp, "%s", temp, sizeof(temp));  // # of frame
	this->numFrames = atoi(temp);

	fscanf_s(fp, "%s", temp, sizeof(temp));  // Frame
	fscanf_s(fp, "%s", temp, sizeof(temp));  // Time:
	fscanf_s(fp, "%s", temp, sizeof(temp));  // per-frame time
	this->frameTime = atof(temp);

	for (auto frame = 0; frame < this->numFrames; frame++)
	{
		std::vector<double> channel;
		for (auto chIdx = 0; chIdx < 3; chIdx++)
		{
			fscanf_s(fp, "%s", temp, sizeof(temp));
			channel.push_back(atof(temp)/100.);
		}
		for (auto chIdx = 3; chIdx < channelBeginIdx; chIdx++)
		{
			fscanf_s(fp, "%s", temp, sizeof(temp));
			channel.push_back(atof(temp));
		}
		this->channels.push_back(channel);
	}

	fclose(fp);
}

inline int Bvh::getJointIndex(const char* name)
{
	for (auto i = 0; i < jointNames.size(); i++)
	{
		if (!strcmp(name, jointNames.at(i).c_str()))
			return i;
	}
	return -1;
}

inline int Bvh::getJointIndex(const BvhJoint * joint)
{
	for (auto i = 0; i < joints.size(); i++)
	{
		if (joints[i] == joint)
			return i;
	}
	return -1;
}

inline int Bvh::getJointChannelIndex(const char* name)
{
	int idx = getJointIndex(name);
	assert(idx >= 0);
	return JointIdxToChannelStartIndex[idx];
}

inline int Bvh::getJointChannelIndex(const int jointIdx)
{
	assert(jointIdx < JointIdxToChannelStartIndex.size());
	return JointIdxToChannelStartIndex[jointIdx];
}

inline int Bvh::getJointChannelIndex(const BvhJoint * joint)
{
	int idx = getJointIndex(joint);
	assert(idx >= 0);
	return JointIdxToChannelStartIndex[idx];
}

void Bvh::getJointPosition(sm::real * pos, int frame, int jointIdx)
{
	getJointPosition(pos, frame, joints[jointIdx]);
}

void Bvh::getJointPosition(sm::real *pos, int frame, const BvhJoint* joint)
{
	assert(frame >= -1 && frame < numFrames);
	if (joint->parent == nullptr)
	{
		sm::RVector3 pos_temp;
		VCOPY3(pos_temp, joint->offsetFromParent);
		if (frame >= 0)
		{
			sm::RVector3 pos_temp2;
			VADD3(pos_temp2, pos_temp, channels[frame].data());
			VCOPY3(pos_temp, pos_temp2);
		}
		MMULT331_X(pos, rot, pos_temp);
	}
	else
	{
		sm::RMatrix33 parent_rotation;
		sm::RVector3 parent_pos, pos_temp, pos_temp1;
		MIDEN33(parent_rotation);
		if (frame >= 0)
			getJointGlobalRotationAfter(parent_rotation, frame, joint->parent);

		getJointPosition(parent_pos, frame, joint->parent);
		
		MMULT331_X(pos_temp, rot, joint->offsetFromParent);
		MMULT331_X(pos_temp1, parent_rotation, pos_temp);
		VADD3(pos, pos_temp1, parent_pos);
	}
}

void Bvh::getJointGlobalRotationAfter(sm::real *rotation, int frame, const BvhJoint * joint)
{
	assert(frame >= 0 && frame < numFrames);
	sm::RMatrix33 ori_rotation, rotation_temp;
	getOriginalJointGlobalRotationAfter(rotation, frame, joint);
	MMULT333_XX(rotation_temp, rot, ori_rotation);
	MMULT333_XT(rotation, rotation_temp, rot);
}

void Bvh::getBonePosition(sm::real * pos, int frame, int jointIdx)
{
	BvhJoint *joint = joints[jointIdx];
	if (joint->children.size() > 1)
	{
		VZERO3(pos);
		sm::real tempPos[3];
		for (auto i = 0; i < joint->children.size(); i++)
		{
			getJointPosition(tempPos, frame, joint->children[i]);
			VADD3(pos, pos, tempPos);
		}
		VDIV3S(pos, pos, joint->children.size());
	}
	else
	{
		sm::real pos1[3], pos2[3];
		getJointPosition(pos1, frame, jointIdx);
		getJointPosition(pos2, frame, joints[jointIdx]->children[0]);
		VADD3(pos, pos1, pos2);
		VMULT3S(pos, .5, pos);
	}
}

sm::real Bvh::getBoneLength(int jointIdx)
{
	const BvhJoint *joint = joints[jointIdx];
	if (joint->children.size() > 0)
	{
		return VLEN3(joint->offsetFromParent);
	}
	
	else
	{
		return -1.;
	}
}

void Bvh::getJointLocalRotationAtPredecessor(
	sm::real * rotation,
	int frame, 
	const BvhJoint * joint)
{
	assert(frame >= 0 && frame < numFrames);
	sm::RMatrix33 ori_rotation, rotation_temp;
	getOriginalJointLocalRotationAtPredecessor(ori_rotation, frame, joint);
	MMULT333_XX(rotation_temp, rot, ori_rotation);
	MMULT333_XT(rotation, rotation_temp, rot);
}

void Bvh::rotate(sm::real * rotation)
{
	MCOPY33(rot, rotation);
}

void Bvh::getOriginalJointLocalRotationAtPredecessor(sm::real * rotation, int frame, const BvhJoint * joint)
{
	assert(frame >= 0 && frame < numFrames);
	const int channelBeginIdx = getJointChannelIndex(joint);
	sm::RVector3 angles;
	sm::RMatrix33 rotz, rotx, roty;

	// get angles of joint at each channel
	if (joint->parent == nullptr)
	{
		angles[0] = channels.at(frame).at(channelBeginIdx + 3);
		angles[1] = channels.at(frame).at(channelBeginIdx + 4);
		angles[2] = channels.at(frame).at(channelBeginIdx + 5);
	}
	else
	{
		angles[0] = channels.at(frame).at(channelBeginIdx + 0);
		angles[1] = channels.at(frame).at(channelBeginIdx + 1);
		angles[2] = channels.at(frame).at(channelBeginIdx + 2);
	}

	// degree to radian
	for (auto i = 0; i < 3; i++)
		angles[i] = TORAD(angles[i]);

	sm::RMatrix33 rotation_temp;
	// get rotation matrix about single channel respectively
	MROTZ33(rotz, cos(angles[0]), sin(angles[0]));
	MROTX33(rotx, cos(angles[1]), sin(angles[1]));
	MROTY33(roty, cos(angles[2]), sin(angles[2]));

	// calculate rotation matrix
	MMULT333_XX(rotation_temp, rotz, rotx);
	MMULT333_XX(rotation, rotation_temp, roty);
}

void Bvh::getOriginalJointGlobalRotationAfter(sm::real * rotation, int frame, const BvhJoint * joint)
{
	if (joint->parent == nullptr)
	{
		getOriginalJointLocalRotationAtPredecessor(rotation, frame, joint);
	}
	else
	{
		sm::RMatrix33 prerotation, localRotation;
		getOriginalJointGlobalRotationAfter(prerotation, frame, joint->parent);
		getOriginalJointLocalRotationAtPredecessor(localRotation, frame, joint);
		MMULT333_XX(rotation, prerotation, localRotation);
	}
}

