#pragma once

#include "myBvh.h"

class MyApp
{
public:
	static void start();
	static void stop();
	static void step(int pause);
	static void command(int cmd);

	static int mouseButton(int cur_button, int button_all, int state,
		int x, int y, bool ctrl, bool shift);
	static int mouseMotion(int button_all, int x, int y, 
		int lastx, int lasty, bool ctrl, bool shift);
};

void buildHumanoidFromBvh(sml::mocap::BVHCharacter *character);
void buildHumanoidFromMyBvh(Bvh &bvh);