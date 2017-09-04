#ifndef _MYBVH_H_
#define _MYBVH_H_

#pragma once

#include <sml.h>
#include <vector>
#include <string>

class BvhJoint
{
public:
	BvhJoint(char *_name) : name(_name)
	{
		hasEndChild = false;
	}
	~BvhJoint()
	{
		children.clear();
	}

	std::string name;
	sm::RVector3 offsetFromParent;
	int dof;

	bool hasEndChild;

	BvhJoint *parent;
	std::vector<BvhJoint*> children;
};

class Bvh
{
public:
	Bvh() { MIDEN33(rot); }
	Bvh(const char * filename) 
	{ 
		MIDEN33(rot); 
		constructFromBvhFile(filename); 
	}

	~Bvh();
	BvhJoint *root;
	int numFrames;
	double frameTime;
	
	std::vector<BvhJoint*> joints;
	std::vector<std::string> jointNames;
	std::vector<int> JointIdxToChannelStartIndex;
	std::vector<BvhJoint*> channelIndexToJoint;
	std::vector<std::vector<double>> channels;

	sm::RMatrix33 rot;

	void constructFromBvhFile(const char *filename, sm::real ratio=1.);

	int getJointNum() { return joints.size(); }
	BvhJoint* getJoint(int idx) { return joints[idx]; }
	
	int getJointIndex(const char* name);
	int getJointIndex(const BvhJoint* joint);
	int getJointChannelIndex(const char* name);
	int getJointChannelIndex(const int jointIdx);
	int getJointChannelIndex(const BvhJoint* joint);

	void getJointPosition(sm::real *pos, int frame, const int jointIdx);
	void getJointPosition(sm::real *pos, int frame, const char* name);
	void getJointPosition(sm::real *pos, int frame, const BvhJoint* joint);
	void getJointGlobalRotationAfter(sm::real *rotation, int frame, const int jointIdx);
	void getJointGlobalRotationAfter(sm::real *rotation, int frame, const char* name);
	void getJointGlobalRotationAfter(sm::real *rotation, int frame, const BvhJoint* joint);

	void getBonePosition(sm::real *pos, int frame, const int jointIdx);
	void getBonePosition(sm::real *pos, int frame, const char* name);
	void getBonePosition(sm::real *pos, int frame, const BvhJoint* joint);
	sm::real getBoneLength(int jointIdx);
	sm::real getBoneLength(const char *name);
	sm::real getBoneLength(const BvhJoint* joint);

	void getJointLocalRotationAtPredecessor(sm::real *rotation, int frame, const int jointIdx);
	void getJointLocalRotationAtPredecessor(sm::real *rotation, int frame, const char *name);
	void getJointLocalRotationAtPredecessor(sm::real *rotation, int frame, const BvhJoint* joint);

	void rotate(sm::real *rotation);

private:
	void getOriginalJointLocalRotationAtPredecessor(sm::real *rotation, int frame, const BvhJoint* joint);
	void getOriginalJointGlobalRotationAfter(sm::real *rotation, int frame, const BvhJoint* joint);
};

#endif // _MYBVH_H_