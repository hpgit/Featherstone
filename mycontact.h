#pragma once

#include <sml.h>

//
struct MY_CONTACT_DATA : public sml::rose::CONTACT_DATA
{
	sml::real f[3];
	sml::real mu;	// friction coefficient
	sml::real e;	// restitution coefficient
};


class MyContact
{
public:
	enum { MAX_NOC = 100 };

public:
	MyContact();

	int update_contacts(sml::rose::GeomSpace& space);

	MY_CONTACT_DATA* get_contact_data() { return m_contact; }
	const MY_CONTACT_DATA* get_contact_data() const { return m_contact; }

	void set_noc(int noc) { assert( 0<=noc && noc<MAX_NOC ); m_noc = noc; }
	int get_noc() const { return m_noc; }

	void set_contact_force(int i, const sml::real* f);

	void draw(sml::GraphicsInterface* gi) const;

	sm::RMatrix get_jacob(sml::rose::RigidBodySystem& robot);
	sm::RVector get_v_star(sml::rose::RigidBodySystem& robot,
		sml::real kappa=REAL(1.0e-2), sml::real epsilon=REAL(1.0e-2));

protected:
	int m_noc;	// number of contacts
	MY_CONTACT_DATA	m_contact[MAX_NOC];
};