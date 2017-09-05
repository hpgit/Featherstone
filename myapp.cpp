#include <sg/sg.h>
#include <sml.h>
#include "myapp.h"
#include "util.h"
#include "mycontact.h"
#include "pgs.h"

using namespace sm;
using namespace sml;
using namespace sml::rose;
using namespace sml::util;

using sml::mocap::BVHCharacter;

// Common instances
SGGraphicsInterface			gi;
RigidBodySystem				robot;			// articulated robot
RigidBodySystem				humanoid;		// humanoid robot
GeomBox						ground;			// ground geometry
GeomVectorSpace				robot_cspace;	// collision space for robot
GeomVectorSpace				world_cspace;	// collision space for world
MyContact					contact;
ElapsedTimeChecker<real>	sim_timer;

bool g_bCaptured	= false;
int  g_body_idx		= 0;
int  g_x_prev		= 0;
int  g_y_prev		= 0;

// Input parameters
real			h;		// time step size
RMatrix			B;		// next-step velocity damping matrix
RMatrix			Ma;		// armature inertia of motors
RMatrix			R;		// regularizer of contact impulse

// bvh
Bvh myBvh;
int frame = 0;
static std::vector<BvhJoint *> bodyIdxToBvhJoint;

sml::mocap::BVHParser bvhParser;
sml::mocap::BVHCharacter *bvhCharacter;

void MyApp::start()
{
	sm::RMatrix33 rot;
	MROTX33(rot, 0., 1.);

	//mybvh
	myBvh.constructFromBvhFile("../motions/wd2_WalkForwardNormal00.bvh", 0.02);
	myBvh.rotate(rot);

	// hp
	bvhParser.parse("../motions/wd2_WalkForwardNormal00.bvh");
	bvhParser.scaleOffsets(0, 0.02);
	bvhParser.scaleRootPosition(0, 0.02);
	bvhCharacter = new sml::mocap::BVHCharacter(&bvhParser, 0);
	bvhCharacter->rotateRoot(rot);
	bvhCharacter->translateRoot(0., 0., -2.);

#if 0
	size_t jointNum = bvhCharacter->getJointSize();
	for (auto jj = 0; jj < jointNum; jj++)
	{
		sml::mocap::JOINT_PTR joint = bvhCharacter->getJoint(jj);
		printf("[%d] %s\n", jj, bvhCharacter->getJointName(joint));
	}
	printf("\nBones: \n");
	size_t bonenum = bvhCharacter->getBoneSize();
	for (auto jj = 0; jj < jointNum; jj++)
	{
		BVHCharacter::BVHBone *bone = (BVHCharacter::BVHBone *)bvhCharacter->getBone(jj);
		printf("[%d] %s\n", jj, bvhCharacter->getBoneName((sml::mocap::BONE_PTR *)bone));
	}
#endif
	const sm::real *rootPos = bvhCharacter->getJointPosition(bvhCharacter->getJoint(0));
	
	//buildHumanoidFromBvh(bvhCharacter);
	buildHumanoidFromMyBvh(myBvh);


	// Build a serial-link robot.
	// buildSerialLinks(robot, 4);
	humanoid.getBodyGeomObjects(&robot_cspace);
	world_cspace.add(&robot_cspace);

	// Make a ground geometry.
	ground.setSize(100, 100, 10);
	ground.setPosition(0, 0, -5);
	ground.updateBV();
	world_cspace.add(&ground);

	// Set up simulation parameters.
	int dof = humanoid.dof();
	h	= REAL(0.005);
	B	= REAL(0.001) * eye(dof);
	Ma	= zeros(dof);

	contact.update_contacts(world_cspace);
	sim_timer.setInterval(h);
	sim_timer.start(get_high_res_time());

	sgSetViewpointD(sm::RVector("-5 0 2"), sm::RVector("0 0 0"));
}

void MyApp::stop()
{
}

void MyApp::step(int pause)
{
	if( !pause && sim_timer.elapsed(get_high_res_time()) )
	{
		// Given parameters
		int dof		= humanoid.dof();
		int rqdim	= humanoid.rqdim();
		real* q		= humanoid.getData(ROSE_RBS_DATA_q);
		real* qd	= humanoid.getData(ROSE_RBS_DATA_qd);
		RVector v(qd, dof);
		RVector u(dof);

		//***** Phase I: Compute M, c, p
		RMatrix M0		= sysmass(&humanoid);
		RMatrix M		= M0 + Ma + h*B;
		RMatrix c		= biasforce(&humanoid, humanoid.getGravAcc());
		RVector	p		= -B*v;

		//***** Phase II: Impulse dynamics
		RMatrix M_hat		= M;
		RMatrix M_hat_inv	= pinv(M_hat);	// it can be done efficiently by exploiting branch-induced sparcity.
		//RMatrix M_hat_inv	= invsysmass(&robot);
		RVector v_hat		= v + h*M_hat_inv*(p+u-c);

		// Compute v_prime = v_hat + M_hat_inv * J' * f.
		RVector v_prime		= v_hat;

		if( contact.get_noc() > 0 )
		{
			int noc = contact.get_noc();

			// contact dynamics (v_plus = v_minus + A*f)
			RMatrix J			= contact.get_jacob(humanoid);
			RMatrix A			= J*M_hat_inv*J.tp();
			RVector v_minus		= J*v_hat;
			RVector v_star		= contact.get_v_star(humanoid);

			// solve min 0.5 * f'(A + R)f + f'(v_minus - v*) for f
			int n = J.rows();
			std::vector<real> cone_coeff(n);
			std::vector<int> cone_begin(n);
			std::vector<int> cone_end(n);
			std::vector<int> cone_normal(n);
			std::vector<int> idx_ubc(2*n/3);
			for(int i = 0; i < n; i+=3)
			{
				cone_coeff[i+0] = 0.8;
				cone_coeff[i+1] = 0.8;
				cone_coeff[i+2] = 0;

				cone_begin[i+0] = i+0;
				cone_begin[i+1] = i+0;
				cone_begin[i+2] = 0;

				cone_end  [i+0] = i+2;
				cone_end  [i+1] = i+2;
				cone_end  [i+2] = 0;

				cone_normal[i+0] = i+2;
				cone_normal[i+1] = i+2;
				cone_normal[i+2] = 0;

				idx_ubc[2*i/3+0] = i+0;
				idx_ubc[2*i/3+1] = i+1;
			}

			RVector f(n);
			R = REAL(1.0e-5) * eye(n);
			pgssolve(f, A+R, v_minus - v_star, n, f, 50, &idx_ubc[0], (int)idx_ubc.size(),
				&cone_coeff[0], &cone_begin[0], &cone_end[0], &cone_normal[0]);

			v_prime += M_hat_inv * J.tp() * f;

			for(int i = 0; i < noc; ++i) {
				contact.set_contact_force(i, f.getData()+i*3);
			}
		}

		// ***** Phase III: Numerical integration
		implicit_euler_integration(humanoid, v_prime, h);
		contact.update_contacts(world_cspace);
	}

	contact.draw(&gi);
	humanoid.draw(&gi);
	humanoid.resetExternalForce();
	humanoid.resetJointTorque();

	if( false && g_bCaptured )
	{
		// Apply a drag force to a selected body towards a user-specified point on the y-z plane.
		RVector e(3);

		// find the user-specified point.
		gi.castRayOnPlane(e, g_x_prev, g_y_prev, RVector("1 0 0 0"));

		RVector pb(robot.getBody(g_body_idx)->getPosition(), 3);
		RVector vb(robot.getBody(g_body_idx)->getLinearVel(), 3);
		real kp = 100;
		real kd = 2*sqrtr(kp);
		robot.addBodyForce(g_body_idx, kp * (e - pb) - kd * vb);

		gi.setColor(1,0,0);
		gi.drawSphere(e, NULL, 0.05);
		gi.drawLine(pb, e);
	}
}

void MyApp::command(int cmd)
{
}

int MyApp::mouseButton(int cur_button, int button_all, int state,
					   int x, int y, bool ctrl, bool shift)
{
	int num_bodies = robot.nob();	// Get the number of bodies.
	RigidBody* body;
	GeomBox* box;
	real p0[3], ray[3];

	if( !g_bCaptured && state == 1 )	// mouse button pressed
	{
		sgWndToWrdD(x, y, 0, p0);	// Coordinate conversion (Window to World)
			
		for(int i = 0; i < num_bodies; ++i)
		{
			body		= robot.getBody(i);
			box			= (GeomBox*) body->getGeomObject();
			sgWndToRayD(x, y, ray);
			if( vray_hit_box(ray, p0, box->getPosition(), box->getRotation(), box->getSize()) )
			{
				g_x_prev = x;
				g_y_prev = y;
				g_body_idx = i;
				box->setColor(1, 0, 0);
				g_bCaptured = true;
				break;
			}
		}
	}
	else
	{
		if( g_bCaptured )
		{
			box = (GeomBox*) robot.getBody(g_body_idx)->getGeomObject();
			box->setColor(0,1,0, 0.5);
			g_bCaptured = false;
			return 1;
		}
	}
	return 0;
}

int MyApp::mouseMotion(int button_all, int x, int y, 
	int lastx, int lasty, bool ctrl, bool shift)
{
	if( g_bCaptured )
	{
		g_x_prev = x;
		g_y_prev = y;
		return 1;
	}

	return 0;
}

void buildHumanoidFromMyBvh(Bvh &bvh)
{
	const sm::real bone_width = 0.15;
	std::vector<int> boxes_index;
	const int num_boxes = bvh.getJointNum();
	
	MemoryManager<RigidBody>		mmrb;
	MemoryManager<GSphericalJoint>	mmrj;
	MemoryManager<GeomBox>			mmgb;

	RigidBody*			box = mmrb.malloc_arr(num_boxes);
	GFreeMotionJoint	fjoint;
	GSphericalJoint*	rjoint = mmrj.malloc_arr(num_boxes);
	GeomBox*			geom = mmgb.malloc_arr(num_boxes);
	RMatrix33 I;
	RMatrix33 R;

	std::vector<RigidBody*>bodies;

	MIDEN33(I);
	MIDEN33(R);

	sm::RVector offset_z("0 0 0.1");

	for (int i = 0; i < num_boxes; ++i)
	{
		const BvhJoint *joint = bvh.joints[i];
		const int children_size = joint->children.size();
		if (children_size == 1)
		{
			// printf("%d\n", i);
			sm::real pos[3];
			bvh.getJointPosition(pos, -1, i);
			sm::real pos_child[3];
			bvh.getJointPosition(pos_child, -1, joint->children[0]);

			sm::real bone_vec[3];
			VSUB3(bone_vec, pos_child, pos);

			sm::real bone_pos[3];
			VADD3(bone_pos, pos_child, pos);
			VMULT3S(bone_pos, .5, bone_pos);

			VADD3(bone_pos, bone_pos, offset_z);

			const sm::real bone_length = VLEN3(bone_vec) * 0.9;
#if 0
			printf("name: %s, length: %f\n", joint->name.c_str(), bone_length);
			vprint3("jointpos: \n", "%f\t", pos);
			vprint3("jointchildpos: \n", "%f\t", pos_child);
#endif
			const sm::real bone_size[3] = { bone_width, bone_width, bone_length };
			const sm::real bone_mass = bone_width*bone_width*bone_length * 10000.;
			box[i].setInitPosition(bone_pos);		// box position

			sm::real unitz[3] = { 0., 0., 1. };
			sm::real axis[3];
			sm::real unit_axis[3];
			VCROSS3(axis, unitz, bone_vec);
			vunit3(unit_axis, axis);
			sm::real angle = acos(VDOT3(unitz, bone_vec) / VLEN3(bone_vec));
			sm::RQuaternion q;
			qset_au(q, angle, unit_axis);
			sm::RMatrix33 rm;
			q2rm(rm, q);

			box[i].setInitRotation(rm);										// box rotation
			box[i].getInertia()->setInertiaBoxTotal(bone_mass, bone_size);	// box inertia
			geom[i].setSize(bone_size);										// box geometry
			geom[i].setColor(0.2, 0.7, 0.7, 1.0);
			box[i].setGeomObject(geom + i);
			box[i].setName(joint->name.c_str());
		}
		else
		{
			sm::real joint_pos_ave[3] = { 0., 0., 0. };
			sm::real joint_pos_max[3] = { 0., 0., 0. };
			sm::real joint_pos_min[3] = { 0., 0., 0. };

			for (auto j = 0; j < children_size; j++)
			{
				sm::real joint_pos[3];
				bvh.getJointPosition(joint_pos, -1, joint->children[j]);
				VADD3(joint_pos_ave, joint_pos_ave, joint_pos);

				if (j == 0)
				{
					VCOPY3(joint_pos_max, joint_pos);
					VCOPY3(joint_pos_min, joint_pos);
				}
				else
				{
					for (auto k = 0; k < 3; k++)
					{
						if (joint_pos_max[k] < joint_pos[k]) joint_pos_max[k] = joint_pos[k];
						if (joint_pos_min[k] > joint_pos[k]) joint_pos_min[k] = joint_pos[k];
					}
				}
			}
			sm::real bone_size[3];
			VSUB3(bone_size, joint_pos_max, joint_pos_min);
			for (auto k = 0; k < 2; k++)
				if (bone_size[k] < bone_width)
					bone_size[k] = bone_width;
			const sm::real bone_mass = bone_size[0] * bone_size[1] * bone_size[2] * 10000;

			VDIV3S(joint_pos_ave, joint_pos_ave, children_size);

			VADD3(joint_pos_ave, joint_pos_ave, offset_z);

			box[i].setInitPosition(joint_pos_ave);
			box[i].setInitRotation(I);										// box rotation
			box[i].getInertia()->setInertiaBoxTotal(bone_mass, bone_size);	// box inertia
			geom[i].setSize(bone_size);										// box geometry
			geom[i].setColor(0.2, 0.7, 0.7, 1.0);
			box[i].setGeomObject(geom + i);
			box[i].setName(joint->name.c_str());
		}

		if (i == 0) {										// body joint
			fjoint.setInitPosition(box[i].getInitPosition());
			box[i].setJoint(&fjoint);
		}
		else {
			sm::real joint_pos[3];
			bvh.getJointPosition(joint_pos, -1, i);
			VADD3(joint_pos, joint_pos, offset_z);
			rjoint[i].setInitPosition(joint_pos);
			rjoint[i].setInitRotation(R);
			rjoint[i].setName(joint->name.c_str());
			box[i].setJoint(rjoint + i);
			box[i].setParent(&box[bvh.getJointIndex(joint->parent)]);					// parent body
		}
	}
	humanoid.constructSystem(box);
	humanoid.setGravAcc(0, 0, REAL(-9.81));

	for (auto i = 0; i < num_boxes; ++i)
	{
		for (auto j = 0; j < bvh.getJointNum(); j++)
		{
			if (!strcmp(bvh.joints[j]->name.c_str(), humanoid.getBody(i)->getName()))
			{
				bodyIdxToBvhJoint.push_back(bvh.joints[j]);
				break;
			}
		}
	}

	///*
	// init pose setting
	sm::real *_rq = humanoid.getData(sml::rose::ROSE_RBS_DATA_rq);
	sm::real *_q = humanoid.getData(sml::rose::ROSE_RBS_DATA_q);
	// VZERON(_q, humanoid.getDataSize(sml::rose::ROSE_RBS_DATA_q));
	// VSET3(_q, bvh.channels[0][0], bvh.channels[0][1], bvh.channels[0][2]);
	for (int i = 0; i < num_boxes; ++i)
	{
		sm::RQuaternion quat;
		sm::RMatrix33 rm;
		bvh.getJointLocalRotationAtPredecessor(rm, 0, bodyIdxToBvhJoint[i]);
		rm2q(quat, rm);
		QCOPY(&_rq[humanoid.rqidx(i)], quat);
	}
	humanoid.refreshState();
    printf("dof: %d\n", humanoid.dof());
	//*/
    humanoid.export_to_xml("humanoid.xml");
    
}

void buildHumanoidFromBvh(sml::mocap::BVHCharacter *character)
{
	const sm::real bone_width = 0.15;
	std::vector<int> boxes_index;
	const int num_bones = character->getBoneSize();
	for (auto i = 0; i < num_bones; i++)
		if (character->getBoneChildrenSize(character->getBone(i)) > 0)
			boxes_index.push_back(i);

	const int num_boxes = boxes_index.size();

	MemoryManager<RigidBody>		mmrb;
	MemoryManager<GSphericalJoint>	mmrj;
	MemoryManager<GeomBox>			mmgb;

	RigidBody*			box = mmrb.malloc_arr(num_boxes);
	GFreeMotionJoint	fjoint;
	GSphericalJoint*	rjoint = mmrj.malloc_arr(num_boxes);
	GeomBox*			geom = mmgb.malloc_arr(num_boxes);
	RMatrix33 I;
	RMatrix33 R;

	MIDEN33(I);
	MIDEN33(R);
	character->seekFrame(0);
	for (int i = 0; i < num_boxes; ++i)
	{
		const int bone_idx = boxes_index[i];
		BVHCharacter::BVHBone *bone_parent = (BVHCharacter::BVHBone *)character->getBone(bone_idx);
		const int bone_children_size = character->getBoneChildrenSize(bone_parent);
		if (bone_children_size == 1)
		{
			BVHCharacter::BVHBone *bone = (BVHCharacter::BVHBone *)character->getBoneChild(bone_parent, 0);
			const sm::real * pos = character->getBonePosition(bone);
			//sm::real pos[3];
			//getTPoseBonePosition(character, bone, pos);
			const sm::real * rot = character->getBoneRotation(bone);
			const sm::real bone_length = character->getBoneLength(bone) * 0.9;
			const sm::real bone_size[3] = { bone_width, bone_width, bone_length };
			const sm::real bone_mass = bone_width*bone_width*bone_length * 1000.;
			box[i].setInitPosition(pos);		// box position
			box[i].setInitRotation(rot);										// box rotation
			box[i].getInertia()->setInertiaBoxTotal(bone_mass, bone_size);	// box inertia
			geom[i].setSize(bone_size);										// box geometry
			geom[i].setColor(0.2, 0.7, 0.7, 1.0);
			box[i].setGeomObject(geom + i);
		}
		else
		{
			sm::real joint_pos_ave[3] = { 0., 0., 0. };
			sm::real joint_pos_max[3] = { 0., 0., 0. };
			sm::real joint_pos_min[3] = { 0., 0., 0. };
			for (auto j = 0; j < bone_children_size; j++)
			{
				const sm::real *joint_pos = character->getJointPosition(
					character->getJoint(
						character->getBoneIndex(
							character->getBoneChild(bone_parent, j))));
				VADD3(joint_pos_ave, joint_pos_ave, joint_pos);
				
				if (j == 0)
				{
					VCOPY3(joint_pos_max, joint_pos);
					VCOPY3(joint_pos_min, joint_pos);
				}
				else
				{
					for (auto k = 0; k < 3; k++)
					{
						if (joint_pos_max[k] < joint_pos[k]) joint_pos_max[k] = joint_pos[k];
						if (joint_pos_min[k] > joint_pos[k]) joint_pos_min[k] = joint_pos[k];
					}
				}
			}
			sm::real bone_size[3];
			VSUB3(bone_size, joint_pos_max, joint_pos_min);
			for (auto k = 0; k < 3; k++)
				if (bone_size[k] < bone_width)
					bone_size[k] = bone_width;

			VDIV3S(joint_pos_ave, joint_pos_ave, bone_children_size);

			box[i].setInitPosition(joint_pos_ave);
			box[i].setInitRotation(I);										// box rotation
			box[i].getInertia()->setInertiaBoxTotal(0.1, bone_size);	// box inertia
			geom[i].setSize(bone_size);										// box geometry
			geom[i].setColor(0.2, 0.7, 0.7, 1.0);
			box[i].setGeomObject(geom + i);
		}

		if (i == 0) {										// body joint
			fjoint.setInitPosition(box[i].getInitPosition());
			box[i].setJoint(&fjoint);
		}
		else {
			const sm::real *parent_joint_pos = character->getJointPosition(character->getJoint(character->getBoneIndex(bone_parent)));
			const sm::real *parent_joint_rot = character->getJointRotation(character->getJoint(character->getBoneIndex(bone_parent)));
			rjoint[i].setInitPosition(parent_joint_pos);
			rjoint[i].setInitRotation(parent_joint_rot);
			box[i].setJoint(rjoint + i);
			for (auto j = 0; j < boxes_index.size(); j++)
			{
				if (boxes_index[j] == character->getBoneIndex(character->getBoneParent(bone_parent)))
				{
					box[i].setParent(&box[j]);					// parent body
					break;
				}
			}
		}
	}
	humanoid.constructSystem(box);
	humanoid.setGravAcc(0, 0, REAL(-9.81));
}