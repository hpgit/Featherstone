#include <sg/sg.h>
#include <sml.h>
#include "myapp.h"
#include "myBvh.h"

int main()
{
	const int argc = 3;
	char* argv[] = {"", "-pause", "-notex"};

	sgFunctions sgFunc;

	sgFunc.start		= MyApp::start;
	sgFunc.stop			= MyApp::stop;
	sgFunc.step			= MyApp::step;
	sgFunc.command		= MyApp::command;
	sgFunc.mouse_button = MyApp::mouseButton;
	sgFunc.mouse_motion = MyApp::mouseMotion;
	sgFunc.path_to_textures	= "../../sg/textures";

	sgSimulationLoop(argc, argv, 1024, 768, &sgFunc);
	
	return 0;
}