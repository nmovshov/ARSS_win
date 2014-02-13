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
	ncc::GetStrPropertyFromINIFile("box","box_unit_size","100",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	rox::params.boxSize = atof(buf);

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
	// Define wall dimension and thickness
	PxReal U = rox::params.boxSize;
	PxReal t = 0.02*U;

	// We'll make the containment with a kinematic actor
	PxRigidDynamic* theBox = gPhysX.mPhysics->createRigidDynamic(PxTransform(PxVec3(0)));
	if (!theBox)
		ncc__error("actor creation failed!");
	theBox->setRigidDynamicFlag(PxRigidDynamicFlag::eKINEMATIC, true);

	// Define sides
	PxBoxGeometry box_side(U/2,t/2,U/2);
	PxMaterial* defmat=gPhysX.mDefaultMaterial;

	// Attach the sides
	theBox->createShape(box_side,*defmat); // the bottom
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(0,U,0))); // the top
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(-U/2,U/2,0),PxQuat(PxPi/2,PxVec3(0,0,1)))); // left wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(U/2,U/2,0),PxQuat(PxPi/2,PxVec3(0,0,1)))); // right wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(0,U/2,-U/2),PxQuat(PxPi/2,PxVec3(1,0,0)))); // back wall
	PxShape* fwall = theBox->createShape(box_side,*defmat,PxTransform(PxVec3(0,U/2,U/2),PxQuat(PxPi/2,PxVec3(1,0,0)))); // front wall

	// Make the front wall transparent
	fwall->setName("~fwall");

	// Register the box
	gPhysX.mScene->addActor(*theBox);
	rox::VIPs.theBox = theBox;
}

void rox::CreateFillBoxExperiment()
{
	rox::CreateContainment();
}

// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif
