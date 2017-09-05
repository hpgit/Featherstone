#include "util.h"

using namespace sm;
using namespace sml;
using namespace sml::rose;

RMatrix toDigMatrix(int n, ...)
{
	RMatrix D(n,n);

	int i;
	va_list vl;
	va_start(vl, n);
	for(i = 0; i < n; ++i) {
		D(i,i) = va_arg(vl,real);
	}
	va_end(vl);

	return D;
}

sm::RMatrix toDigMatrix(sm::RVector v)
{
	int n = v.dim();
	RMatrix D(n,n);
	for(int i = 0; i < n; ++i) {
		D(i,i) = v[i];
	}
	return D;
}

sm::RMatrix toDigMatrix(const char* str)
{
	RVector v(str);
	return toDigMatrix(v);
}

void buildSerialLinks(sml::rose::RigidBodySystem& robot, int num_of_links)
{
	const int num_boxes = num_of_links;
	const real box_size[] = {REAL(0.375), REAL(0.375), REAL(0.5)};
	const real box_mass = 5;
	const int z_offset = 1;

	MemoryManager<RigidBody>		mmrb;
	// MemoryManager<RevoluteJoint>	mmrj;
	MemoryManager<GSphericalJoint>	mmrj;
	MemoryManager<GeomBox>			mmgb;

	RigidBody*			box		= mmrb.malloc_arr(num_boxes);
	GFreeMotionJoint	fjoint;
	// RevoluteJoint*		rjoint	= mmrj.malloc_arr(num_boxes);
	GSphericalJoint*	rjoint = mmrj.malloc_arr(num_boxes);
	GeomBox*			geom	= mmgb.malloc_arr(num_boxes);
	RMatrix33 I;
	RMatrix33 R;

	MIDEN33(I);

	// Initialize initial rotation of joints.
#if 1
	mroty33(R, TORAD(-90));
#else
	MIDEN33(R);
#endif

	for(int i = 0; i < num_boxes; ++i)
	{
		box[i].setInitPosition(
			0, 0, box_size[2]*((num_boxes-i-1+z_offset)+REAL(0.5)));	// box position
		box[i].setInitRotation(I);										// box rotation
		box[i].getInertia()->setInertiaBoxTotal(box_mass, box_size);	// box inertia
		geom[i].setSize(box_size);										// box geometry
		geom[i].setColor(0,1,0, 0.5);
		box[i].setGeomObject(geom+i);

		if( i == 0 ) {										// body joint
			fjoint.setInitPosition(box[i].getInitPosition());
			box[i].setJoint(&fjoint);
		}
		else {
			rjoint[i].setInitPosition(
				0, 0, box_size[2]*((num_boxes-i-1+z_offset)+REAL(1.0)));
			rjoint[i].setInitRotation(R);
			box[i].setJoint(rjoint+i);
			box[i].setParent(&box[i-1]);					// parent body
		}
	}

	robot.constructSystem(box);
	robot.setGravAcc(0, 0, REAL(-9.81));
}

void implicit_euler_integration(sml::rose::RigidBodySystem& robot, sml::real* qd_next, sml::real h)
{
	int dof		= robot.dof();
	int rqdim	= robot.rqdim();
	real* qd	= robot.getData(ROSE_RBS_DATA_qd);

	vcopyN(qd, qd_next, dof);	// Update joint velocities first.
	robot.update_rqd();			// Then, update redundant joint velocities using qd.
	
	RVector x(robot.xdim());
	RVector xd(robot.xdim());
	robot.getState(x);
	robot.getStateDerivative(xd);

	// vprintN("qd: ", "%f\t", qd_next, robot.dof());
	// vprintN("x : ", "%f\t", x, robot.xdim());
	// vprintN("xd: ", "%f\t", xd, robot.xdim());

	// Update joint positions with next-step joint velocities.
	// (we don't need to update joint velocities because it has already been done above)
	for(int i = 0; i < dof; ++i)	x[i]		+= h*xd[i];
	for(int i = 0; i < rqdim; ++i)	x[dof*2+i]	+= h*xd[dof*2+i];

	robot.setState(x);
}