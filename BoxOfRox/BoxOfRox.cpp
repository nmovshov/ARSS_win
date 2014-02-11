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
	
	// Parameters for the box
	ncc::GetStrPropertyFromINIFile("box","box_height","100",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	rox::params.boxHeight = atof(buf);
	ncc::GetStrPropertyFromINIFile("box","box_length","100",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	rox::params.boxLength = atof(buf);
	ncc::GetStrPropertyFromINIFile("box","box_width","100",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	rox::params.boxWidth = atof(buf);

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
	
	switch (rox::eExperimentType)
	{
	case rox::eFILL_BOX:
		rox::CreateFillBoxExperiment();
		break;
	case rox::eBAD_EXPERIMENT_TYPE: // intentional fall through
	default:
		ncc__error("Unknown experiment type. Experiment aborted.\a");
	}

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
void rox::CreateContainment()
{
	// Define height, length, and width
	PxReal H = rox::params.boxHeight;
	PxReal L = rox::params.boxLength;
	PxReal W = rox::params.boxWidth;

	// Now place the walls
	PxActor* anActor;

	// Left wall
	anActor = PxCreatePlane(*gPhysX.mPhysics,PxPlane(PxVec3(1,0,0),L/2),*gPhysX.mDefaultMaterial);
	if (!anActor)
		ncc__error("Containment wall creation failed. Experiment aborted.\a");
	anActor->setName("~leftwall"); // ~name is a quick way to request invisibility
	gPhysX.mScene->addActor(*anActor);

	// Right wall
	anActor = PxCreatePlane(*gPhysX.mPhysics,PxPlane(PxVec3(-1,0,0),L/2),*gPhysX.mDefaultMaterial);
	if (!anActor)
		ncc__error("Containment wall creation failed. Experiment aborted.\a");
	anActor->setName("~rightwall"); // ~name is a quick way to request invisibility
	gPhysX.mScene->addActor(*anActor);

	// Back wall
	anActor = PxCreatePlane(*gPhysX.mPhysics,PxPlane(PxVec3(0,0,-1),W/2),*gPhysX.mDefaultMaterial);
	if (!anActor)
		ncc__error("Containment wall creation failed. Experiment aborted.\a");
	anActor->setName("~backwall"); // ~name is a quick way to request invisibility
	gPhysX.mScene->addActor(*anActor);

	// Front wall
	anActor = PxCreatePlane(*gPhysX.mPhysics,PxPlane(PxVec3(0,0,1),W/2),*gPhysX.mDefaultMaterial);
	if (!anActor)
		ncc__error("Containment wall creation failed. Experiment aborted.\a");
	anActor->setName("~frontwall"); // ~name is a quick way to request invisibility
	gPhysX.mScene->addActor(*anActor);

	// Top wall
	anActor = PxCreatePlane(*gPhysX.mPhysics,PxPlane(PxVec3(0,1,0),H/2),*gPhysX.mDefaultMaterial);
	if (!anActor)
		ncc__error("Containment wall creation failed. Experiment aborted.\a");
	anActor->setName("~topwall"); // ~name is a quick way to request invisibility
	gPhysX.mScene->addActor(*anActor);

	// Bottom wall
	anActor = PxCreatePlane(*gPhysX.mPhysics,PxPlane(PxVec3(0,-1,0),H/2),*gPhysX.mDefaultMaterial);
	if (!anActor)
		ncc__error("Containment wall creation failed. Experiment aborted.\a");
	anActor->setName("~bottomwall"); // ~name is a quick way to request invisibility
	gPhysX.mScene->addActor(*anActor);

	// For convenience, place the camera for a nice viewpoint
	gCamera.pos.y = H;
	gCamera.pos.z = 10*H;
}

void rox::CreateFillBoxExperiment()
{
	rox::CreateContainment();
}

// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif
