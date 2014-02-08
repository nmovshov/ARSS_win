///////////////////////////////////////////////////////////////////////////////
// Source file for project BoxOfRox.
//
// Author: Viranga
///////////////////////////////////////////////////////////////////////////////

#include "ARSS.h"
#include "BoxOfRox.h"

// Request lint level warnings with the LINT macro on Microsoft compilers
#ifdef LINT
#pragma warning(push,4)
#pragma warning(disable:4100) // unreferenced formal parameter
#endif

int main(int argc, char** argv)
{

#ifdef _DEBUG
	cout << "***DEBUG BUILD***" << endl;
#endif

	ncc::whoami();

	InitRun(argc,argv);
	InitGlut(argc,argv);
	InitHUD();
	InitDevIL();
	InitPhysX();
	InitCUDA();
	InitExperiment();
	glutMainLoop(); // enter event processing

	return 0; // never actually reached.
}

void CustomizeRun(int,char**)
{
	// Read experiment-specific program options from file.
	if (!ConfigExperimentOptions())
	{ncc__warning("Could not find ARSS config file. Attempting to continue with all default options.\a");}
}
bool ConfigExperimentOptions()
{
	// First check that an options file exists
	ifstream fp(gRun.iniFile.c_str());
	bool success = fp.good();
	fp.close();

	// Read in parameters by group, file.ini style
	char buf[MAX_CHARS_PER_NAME];

	// Experiment Type
	ncc::GetStrPropertyFromINIFile("experiment","experiment_type","",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if		(strcmp(buf,"fill_box")==0)
		rox::eExperimentType = rox::eFILL_BOX;
	else
		rox::eExperimentType = rox::eBAD_EXPERIMENT_TYPE;
	
	// Parameters of the grain size distribution
	ncc::GetStrPropertyFromINIFile("experiment","gsd_type","uniform",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if		(strcmp(buf,"uniform")==0)
		rox::gsd.type = rox::gsd.eGSD_UNIFORM;
	else
		rox::gsd.type = rox::gsd.eBAD_GSD_TYPE;

	return success;
}
void CustomizeScene(PxSceneDesc &sceneDesc)
{
	
}
void CustomizeGLUT()
{
	
}
void CustomizeHUD()
{

}
void RefreshCustomHUDElements()
{

}
void FireAction()
{
	
}
void LogExperiment()
{

}
void PrintDebug()
{
}
void ApplyCustomInteractions()
{
	
}
void RenderOtherStuff()
{

}
void CreateExperiment()
{
	gDebug.bXYGridOn=true;
	
	
	// Move the camera to where you can see TODO: replace with a function call
	FindExtremers();
	if (gExp.VIPs.extremers.outmost)
		gCamera.pos.z = gExp.VIPs.extremers.outmost->getGlobalPose().p.z + 10*gExp.defGrainSize;

	// Start the action
	gSim.isRunning=true;
	gSim.bPause=false;
	gCUDA.cudaCapable=false; // TODO: remove when cuda gravity is implemented
	RefreshHUD();
}
void RebootExperiment()
{
	gSim.isRunning=false;
	DestroyPhysX();
	InitPhysX();
	CreateExperiment();
}
void LoadExperiment()
{

}
void UpArrowAction()
{

}
void DownArrowAction()
{

}
void LeftArrowAction()
{

}
void RightArrowAction()
{

}

// BoxOfRox namespace functions

// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif
