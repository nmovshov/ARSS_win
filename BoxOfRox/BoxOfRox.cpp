///////////////////////////////////////////////////////////////////////////////
// Source file for project BoxOfRox.
//
// Authors: Viranga and Naor
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

	// Experiment parameters
	ncc::GetStrPropertyFromINIFile("experiment","experiment_type","",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if		(strcmp(buf,"fill_box")==0)
		rox::eExperimentType = rox::eFILL_BOX;
	else if (strcmp(buf,"shake_box")==0)
		rox::eExperimentType = rox::eSHAKE_BOX;
	else
		rox::eExperimentType = rox::eBAD_EXPERIMENT_TYPE;
	
	ncc::GetStrPropertyFromINIFile("experiment","shake_magnitude","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	rox::params.shakeMagnitude = atof(buf);

	// Parameters for the box
	ncc::GetStrPropertyFromINIFile("box","box_unit_size","100",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	rox::params.boxSize = atof(buf);
	ncc::GetStrPropertyFromINIFile("box","box_design","",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if (strcmp(buf,"AOSAT1")==0)
		rox::eBoxDesign = rox::eAOSAT1;
	else
		rox::eBoxDesign = rox::eBAD_BOX_DESIGN;

	// Parameters for the regolith
	ncc::GetStrPropertyFromINIFile("regolith","regolith_type","uniform",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if		(strcmp(buf,"uniform")==0)
		rox::regolith.type = rox::regolith.eREGOLITH_UNIFORM;
	else if	(strcmp(buf,"bimodal")==0)
		rox::regolith.type = rox::regolith.eREGOLITH_BIMODAL;
	else
		rox::regolith.type = rox::regolith.eBAD_REGOLITH_TYPE;

	ncc::GetStrPropertyFromINIFile("regolith","grain_size1","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	rox::regolith.size1 = atof(buf);
	ncc::GetStrPropertyFromINIFile("regolith","grain_size2","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	rox::regolith.size2 = atof(buf);

	rox::regolith.totalNumber = ncc::GetIntPropertyFromINIFile("regolith","grain_total_number",0,gRun.iniFile.c_str());
	
	// Code units and scaling
	ncc::GetStrPropertyFromINIFile("units","big_g","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	rox::units.bigG = atof(buf);

	return success;
}
void CustomizeScene(PxSceneDesc &sceneDesc)
{
	sceneDesc.gravity = PxVec3(0);
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
	rox::GravitateSelf();
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
	rox::OpenLeftDoor();
}
void RightArrowAction()
{
	rox::OpenRightDoor();
}

// BoxOfRox namespace functions
void rox::CreateTheBox()
{
	switch (rox::eBoxDesign)
	{
	case rox::eAOSAT1:
		rox::CreateAOSAT1();
		break;
	case rox::eBAD_BOX_DESIGN: // intentional fall through
	default:
		ncc__error("Unknown box design\a");
	}
}
void rox::FillTheBox()
{
	CreateRubbleGrain(PxVec3(-0.2,rox::params.boxSize/2,0),gExp.defGrainType,gExp.defGrainSize,*gPhysX.mDefaultMaterial,gExp.defGrainDensity);
	CreateRubbleGrain(PxVec3(+0.2,rox::params.boxSize/2,0),gExp.defGrainType,gExp.defGrainSize,*gPhysX.mDefaultMaterial,gExp.defGrainDensity);
}
void rox::OpenLeftDoor()
{
	if (rox::VIPs.ldoor) {
		rox::VIPs.ldoor->release();
		rox::VIPs.ldoor = NULL;
	}
}
void rox::OpenRightDoor()
{
	if (rox::VIPs.rdoor) {
		rox::VIPs.rdoor->release();
		rox::VIPs.rdoor = NULL;
	}
}
void rox::CreateFillBoxExperiment()
{
	rox::CreateTheBox();
	gCamera.pos = PxVec3(0,rox::params.boxSize,3*rox::params.boxSize);
	rox::FillTheBox();
}
void rox::GravitateSelf()
{
	if (gCUDA.cudaCapable)
		rox::GravitateOnDevice();
	else
		rox::GravitateOnHost();
}
void rox::GravitateOnDevice()
{
	rox::GravitateOnHost(); //TODO: include device implementation when needed, for speed
}
void rox::GravitateOnHost()
/*
 * This is a semi-optimized all-pairs gravity calculation in a single host thread. Not
 * usually the best option, but a necessary fall back option.
*/
{
	// Prepare
	PxU32 nbActors = gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,MAX_ACTORS_PER_SCENE);
	if (nbActors<2) return;
	float *bodies = new float[4*nbActors];
	float *forces = new float[3*nbActors];
	if (bodies==NULL || forces==NULL)
		ncc__error("Error allocating memory for body/forces arrays.");

	// Linearize actor positions
	for (PxU32 k=0; k<nbActors; k++)
	{
		PxRigidDynamic* actor = gPhysX.cast[k]->isRigidDynamic();
		PxTransform pose = actor->getGlobalPose().transform(actor->getCMassLocalPose());
		PxVec3 pos = pose.p;
		bodies[4*k+0] = actor->getMass();
		bodies[4*k+1] = pos.x;
		bodies[4*k+2] = pos.y;
		bodies[4*k+3] = pos.z;
	}

	// And now, the N^2 double loop.
	for (PxU32 j=0; j<nbActors; j++)
	{
		forces[3*j+0] = forces[3*j+1] = forces[3*j+2] = 0.0f;
		for (PxU32 k=0; k<j; k++)
		{
			PxReal x = bodies[4*k+1] - bodies[4*j+1];
			PxReal y = bodies[4*k+2] - bodies[4*j+2];
			PxReal z = bodies[4*k+3] - bodies[4*j+3];

			float distSqr = x*x + y*y + z*z;
			float distSix = distSqr*distSqr*distSqr;
			float invDistCube = 1.0f/sqrtf(distSix);

			float s = rox::units.bigG * bodies[4*j+0] * bodies[4*k+0] * invDistCube;

			forces[3*j+0] += x*s;
			forces[3*j+1] += y*s;
			forces[3*j+2] += z*s;
			forces[3*k+0] -= x*s;
			forces[3*k+1] -= y*s;
			forces[3*k+2] -= z*s;
		}
	}

	// Add accumulated forces to actors
	for (PxU32 k=0; k<nbActors; k++)
	{
		PxRigidDynamic* actor = gPhysX.cast[k]->isRigidDynamic();
		PxRigidDynamicFlags flags = actor->getRigidDynamicFlags();
		if (flags & PxRigidDynamicFlag::eKINEMATIC) continue;
		PxVec3 F(forces[3*k+0],forces[3*k+1],forces[3*k+2]);
		actor->addForce(F);
	}

	// Clean up
	delete [] bodies;
	delete [] forces;

	return;
}
void rox::CreateAOSAT1()
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
	PxBoxGeometry minibox_side(U/10,U/2,t/2);
	PxBoxGeometry minibox_door(t/2,U/2,U/10);
	PxMaterial* defmat=gPhysX.mDefaultMaterial;

	// Attach the sides
	// Middle Chamber
	theBox->createShape(box_side,*defmat); // middle chamber bottom wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(0,U,0))); // middle chamber top wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(-U/2,U/2,0),PxQuat(PxPi/2,PxVec3(0,0,1)))); // middle chamber left wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(U/2,U/2,0),PxQuat(PxPi/2,PxVec3(0,0,1)))); // middle chamber right wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(0,U/2,-U/2),PxQuat(PxPi/2,PxVec3(1,0,0)))); // middle chamber back wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(0,U/2,U/2),PxQuat(PxPi/2,PxVec3(1,0,0)))); // middle chamber front wall

	// Left Chamber
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(-U,0,0))); // left chamber bottom wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(-U,U,0))); // left chamber top wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(-U,U/2,-U/2),PxQuat(PxPi/2,PxVec3(1,0,0)))); // left chamber back wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(-U,U/2,U/2),PxQuat(PxPi/2,PxVec3(1,0,0)))); // left chamber front wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(-(3*U/2),U/2,0),PxQuat(PxPi/2,PxVec3(0,0,1)))); // left chamber left wall

	// Right Chamber
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(U,0,0))); // right chamber bottom wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(U,U,0))); // right chamber top wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(U,U/2,-U/2),PxQuat(PxPi/2,PxVec3(1,0,0)))); // right chamber back wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3(U,U/2,U/2),PxQuat(PxPi/2,PxVec3(1,0,0)))); // right chamber front wall
	theBox->createShape(box_side,*defmat,PxTransform(PxVec3((3*U/2),U/2,0),PxQuat(PxPi/2,PxVec3(0,0,1)))); // right chamber right wall

	// Mini Boxes (regolith containment)
	// PxShape* mini = theBox->createShape(minibox_side,*defmat,PxTransform(PxVec3(-(7*U/10),U/2,-U/10))); // left chamber minibox
	// gColors.colorBucket.push_back(vector<GLubyte>(3));
	// gColors.colorBucket.back()[0] = ncc::rgb::rRed[0];
	// gColors.colorBucket.back()[1] = ncc::rgb::rRed[1];
	// gColors.colorBucket.back()[2] = ncc::rgb::rRed[2];
	// mini->userData = &(gColors.colorBucket.back()[0]);
	// mini = theBox->createShape(minibox_side,*defmat,PxTransform(PxVec3((7*U/10),U/2,-U/10))); // right chamber minibox
	// mini->userData = &(gColors.colorBucket.back()[0]);

	// Mini Boxes (regolith containment)
	PxShape* mini = theBox->createShape(minibox_side,*defmat,PxTransform(PxVec3(-(3*U/5),U/2,3*U/10))); // left chamber minibox
	gColors.colorBucket.push_back(vector<GLubyte>(3));
	gColors.colorBucket.back()[0] = ncc::rgb::rRed[0];
	gColors.colorBucket.back()[1] = ncc::rgb::rRed[1];
	gColors.colorBucket.back()[2] = ncc::rgb::rRed[2];
	mini->userData = &(gColors.colorBucket.back()[0]);
	mini = theBox->createShape(minibox_side,*defmat,PxTransform(PxVec3((3*U/5),U/2,3*U/10))); // right chamber minibox
	mini->userData = &(gColors.colorBucket.back()[0]);

	// Doors
	rox::VIPs.ldoor = theBox->createShape(minibox_door,*defmat,PxTransform(PxVec3(-(7*U/10),U/2,2*U/5)));
	rox::VIPs.rdoor = theBox->createShape(minibox_door,*defmat,PxTransform(PxVec3(+(7*U/10),U/2,2*U/5)));
	gColors.colorBucket.push_back(vector<GLubyte>(3));
	gColors.colorBucket.back()[0] = ncc::rgb::bAqua[0];
	gColors.colorBucket.back()[1] = ncc::rgb::bAqua[1];
	gColors.colorBucket.back()[2] = ncc::rgb::bAqua[2];
	rox::VIPs.ldoor->userData = &(gColors.colorBucket.back()[0]);
	rox::VIPs.rdoor->userData = &(gColors.colorBucket.back()[0]);

	// Name, color, and register the box
	theBox->setName("#the_box");
	gColors.colorBucket.push_back(vector<GLubyte>(3));
	gColors.colorBucket.back()[0] = ncc::rgb::yLightYellow[0];
	gColors.colorBucket.back()[1] = ncc::rgb::yLightYellow[1];
	gColors.colorBucket.back()[2] = ncc::rgb::yLightYellow[2];
	theBox->userData = &(gColors.colorBucket.back()[0]);
	gPhysX.mScene->addActor(*theBox);
	rox::VIPs.theBox = theBox;
}

// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif
