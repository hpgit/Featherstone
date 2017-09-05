#pragma once

#include <sml.h>
#include <stdarg.h>

void buildSerialLinks(sml::rose::RigidBodySystem& robot, int num_of_links);
void implicit_euler_integration(sml::rose::RigidBodySystem& robot, sml::real* qd_next, sml::real h);

sm::RMatrix toDigMatrix(int n, ...);
sm::RMatrix toDigMatrix(sm::RVector v);
sm::RMatrix toDigMatrix(const char* str);