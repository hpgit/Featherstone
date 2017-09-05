#include "mycontact.h"

using namespace sm;
using namespace sml;
using namespace sml::rose;
using namespace sml::util;


static int nearCallback(const GeomRealizer* pGeom1, const GeomRealizer* pGeom2, void* pData)
{
	extern GeomVectorSpace robot_cspace;

	MyContact* mycontact = (MyContact*)pData;

	if( pGeom1->isSpace() || pGeom2->isSpace() ) {
		testGeomCollision(pGeom1, pGeom2, nearCallback, pData);
		return 0;
	}

	// If both geometries are in the same articulated body system, just return.
	if( robot_cspace.contains(pGeom1) && robot_cspace.contains(pGeom2) )
		return 0;

	int noc = mycontact->get_noc();
	mycontact->set_noc(
		noc + detectGeomObjectContacts((const GeomObject*)pGeom1, (const GeomObject*)pGeom2,
			mycontact->get_contact_data()+noc, sizeof(MY_CONTACT_DATA), MyContact::MAX_NOC-noc)
	);

	return 0;
}

MyContact::MyContact()
{
	m_noc = 0;
}

int MyContact::update_contacts(sml::rose::GeomSpace& space)
{
	set_noc(0);
	testGeomSpaceCollision(&space, nearCallback, this);
	return m_noc;
}

void MyContact::set_contact_force(int i, const sml::real* f)
{
	assert( 0<=i && i<get_noc()*3 );
	vcopy3(m_contact[i].f, f);
}

void MyContact::draw(sml::GraphicsInterface* gi) const
{
	extern real h;

	gi->setColor(1, 1, 0, REAL(0.5));
	for(int i = 0; i < m_noc; ++i) {
		gi->drawSphere(m_contact[i].p, NULL, REAL(0.05));

		real fp[3];
		real scale = 100;
		vdiv3s(fp, m_contact[i].f, h*scale);
		vadd3(fp, m_contact[i].p, fp);
		gi->drawLine(m_contact[i].p, fp);
	}
}

sm::RMatrix MyContact::get_jacob(sml::rose::RigidBodySystem& robot)
{
	extern GeomBox ground;

	int dof = robot.dof();
	int noc = get_noc();

	RMatrix M(3*noc, dof);
	RMatrix J(6, dof);
	for(int i = 0; i < noc; ++i) {
		const GeomRealizer* geom = m_contact[i].g1;
		if( geom == &ground )
			geom = m_contact[i].g2;
 
		robot.jacob(J, gr2bidx(geom), m_contact[i].p, 0);
		M(3*i, 0, J+3*dof, 3, dof);
	}

	return M;
}

sm::RVector MyContact::get_v_star(
	sml::rose::RigidBodySystem& robot,
	sml::real kappa/*=REAL(1.0e-2)*/,
	sml::real epsilon/*=REAL(1.0e-2)*/
	)
{
	extern real h;
	extern GeomBox ground;

	int noc = get_noc();

	RVector v_star(noc*3);

	real kp =   (1+epsilon) / pow2(kappa);
	real kd = 2*(1+epsilon) / kappa;
	for(int i = 0; i < noc; ++i)
	{
		const GeomRealizer* geom = m_contact[i].g1;
		if( geom == &ground )
			geom = m_contact[i].g2;

		int body_idx = gr2bidx(geom);
		RigidBody* body = robot.getBody(body_idx);

		real vc[3];
		body->getLinearVelAt(vc, m_contact[i].p);

		v_star[i*3+0] = 0;
		v_star[i*3+1] = 0;
		v_star[i*3+2] = vc[2] + h*(kp * (-0 + m_contact[i].depth) + kd * (0 - vc[2]));
	}

	return v_star;
}