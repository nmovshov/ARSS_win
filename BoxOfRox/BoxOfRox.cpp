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
	ncc::GetStrPropertyFromINIFile("regolith","grain_type","uniform",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if		(strcmp(buf,"uniform")==0)
		rox::regolith.type = rox::regolith.eGRAIN_UNIFORM;
	else if	(strcmp(buf,"bimodal")==0)
		rox::regolith.type = rox::regolith.eGRAIN_BIMODAL;
	else
		rox::regolith.type = rox::regolith.eBAD_GRAIN_TYPE;

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

	// Move the camera to where you can see TODO: replace with a function call
	FindExtremers();
	if (gExp.VIPs.extremers.outmost)
		gCamera.pos.z = gExp.VIPs.extremers.outmost->getGlobalPose().p.z + 10*gExp.defGrainSize;

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

}
void RightArrowAction()
{

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
void rox::CreateFillBoxExperiment()
{
	rox::CreateTheBox();


	//// Make a random convex mesh to be referenced by all grains
	//vector<PxVec3> verts = MakeRandomVertexList();
	//PxConvexMesh* theMesh = MakePxMeshFromVertexList(verts);

	//// Calculate placement positions
	//vector<PxVec3> positions(rox::grain.totalNumber);
	//PxReal safeDL = PxMax(rox::grain.size1,rox::grain.size2) * 2.06;

	//// Place actors
	//for (PxU32 k=0; k<positions.size(); k++)
	//{
	//	// select a size for the next grain
	//	PxReal grainScale = rox::grain.size1;
	//	//bool isSize2 = (k % (rox::grain.numberRatio)) == 0;
	//	bool isSize2 = (rox::grain.type==rox::grain.eGRAIN_BIMODAL && k < rox::grain.numberRatio);
	//	if (rox::grain.type==rox::grain.eGRAIN_BIMODAL && isSize2)
	//		grainScale = rox::grain.size2;
	//	
	//	// create a convex geometry for the next grain
	//	PxMeshScale meshScale;
	//	meshScale.scale = PxVec3(grainScale);
	//	PxMat33 M = meshScale.toMat33();
	//	PxConvexMeshGeometry meshGeometry(theMesh,meshScale);

	//	// create a convex actor from this geometry
	//	PxRigidDynamic* aGrain = PxCreateDynamic(*gPhysX.mPhysics,PxTransform(positions[k]),meshGeometry,*gPhysX.mDefaultMaterial,gExp.defGrainDensity);
	//	if (aGrain)
	//	{
	//		if (gPhysX.props.sleepThreshold > -1) aGrain->setSleepThreshold(0.5*gPhysX.props.sleepThreshold*gPhysX.props.sleepThreshold);
	//		if (gPhysX.props.angularDamping > -1) aGrain->setAngularDamping(gPhysX.props.angularDamping);
	//		if (gPhysX.props.linearDamping  > -1) aGrain->setLinearDamping(gPhysX.props.linearDamping);

	//		if (isSize2)
	//		{
	//			aGrain->setName("size2");
	//			ColorActor(aGrain,ncc::rgb::oDarkOrange);
	//		}
	//		else
	//		{
	//			aGrain->setName("size1");
	//		}

	//		RandOrientActor(aGrain);

	//		gPhysX.mScene->addActor(*aGrain);
	//	}
	//}
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
	PxBoxGeometry minibox_side(U/10,t/4,U/2);
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

	// Mini Boxes
	theBox->createShape(minibox_side,*defmat,PxTransform(PxVec3(-(3*U/5),U/2,-(3*U/10)),PxQuat(PxPi/2,PxVec3(1,0,0)))); // left chamber minibox
	theBox->createShape(minibox_side,*defmat,PxTransform(PxVec3((3*U/5),U/2,-(3*U/10)),PxQuat(PxPi/2,PxVec3(1,0,0)))); // right chamber minibox

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
