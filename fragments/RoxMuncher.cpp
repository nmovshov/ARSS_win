/////////////////////////////////////////////////////////////////////////////////
// Source file for project RoxMuncher, an asteroid mining simulator.
//
// Authors: Viranga and Naor
/////////////////////////////////////////////////////////////////////////////////

#include "ARSS.h"
#include "RoxMuncher.h"

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

// Experiment functions
void CustomizeRun(int argc, char** argv)
{
	// Read experiment-specific program options from file.
	if (!ConfigExperimentOptions())
	{ncc__warning("Could not find project config file. Attempting to continue with all default options.\a");}
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
	ncc::GetStrPropertyFromINIFile("experiment","Design","",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if (strcmp(buf,"RoxMuncher1")==0)
		munch::eClawDesign = munch::eRoxMuncher1;
	else
		munch::eClawDesign = munch::eBAD_CLAW_DESIGN;

	ncc::GetStrPropertyFromINIFile("experiment","experiment_type","",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if		(strcmp(buf,"fill_box")==0)
		munch::eExperimentType = munch::eFILL_BOX;
	else if (strcmp(buf,"shake_box")==0)
		munch::eExperimentType = munch::eSHAKE_BOX;
	else
		munch::eExperimentType = munch::eBAD_EXPERIMENT_TYPE;

	// Parameters for the bed
	ncc::GetStrPropertyFromINIFile("bed","bed_height","100",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	munch::params.bedHeight = atof(buf);
	ncc::GetStrPropertyFromINIFile("bed","bed_width","100",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	munch::params.bedWidth = atof(buf);

	// Parameters for the muncher
	ncc::GetStrPropertyFromINIFile("muncher","muncher_height","100",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	munch::params.muncherHeight = atof(buf);
	ncc::GetStrPropertyFromINIFile("muncher","muncher_width","100",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	munch::params.muncherWidth = atof(buf);
	ncc::GetStrPropertyFromINIFile("muncher","muncher_length","100",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	munch::params.muncherLength = atof(buf);

	// Parameters for the regolith
	ncc::GetStrPropertyFromINIFile("regolith","regolith_type","uniform",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if		(strcmp(buf,"uniform")==0)
		munch::regolith.type = munch::regolith.eREGOLITH_UNIFORM;
	else if	(strcmp(buf,"bimodal")==0)
		munch::regolith.type = munch::regolith.eREGOLITH_BIMODAL;
	else
		munch::regolith.type = munch::regolith.eBAD_REGOLITH_TYPE;

	ncc::GetStrPropertyFromINIFile("regolith","grain_size1","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	munch::regolith.size1 = atof(buf);
	ncc::GetStrPropertyFromINIFile("regolith","grain_size2","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	munch::regolith.size2 = atof(buf);

	munch::regolith.totalNumber = ncc::GetIntPropertyFromINIFile("regolith","grain_total_number",0,gRun.iniFile.c_str());
	
	// Code units and scaling
	ncc::GetStrPropertyFromINIFile("units","big_g","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	munch::units.bigG = atof(buf);

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
	char buf[MAX_CHARS_PER_NAME];
	int ch2px = 18; // hud uses 18 pt font
	float px2width = 1.0/glutGet(GLUT_WINDOW_WIDTH);
	float scrPos = glutGet(GLUT_WINDOW_WIDTH);
	unsigned int eid;
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
	munch::GravitateSelf();
}
void CreateExperiment()
{
	gDebug.bXYGridOn=true;
	switch (munch::eExperimentType)
	{
	case munch::eFILL_BOX:
		munch::CreateFillBoxExperiment();
		break;
	case munch::eBAD_EXPERIMENT_TYPE: // intentional fall through
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

// BoxOfRox specific functions
void munch::CreateTheBox()
{
	switch (munch::eClawDesign)
	{
	case munch::eRoxMuncher1:
		munch::CreateRoxMuncher();
		break;
	case munch::eBAD_CLAW_DESIGN: // intentional fall through
	default:
		ncc__error("Unknown claw design\a");
	}
}
void munch::FillTheBox()
{
	// Create a few randomized convex meshes
	int nbMeshes = 10;
	vector<PxConvexMesh*> mesh_list(nbMeshes);
	for (unsigned int k=0; k<mesh_list.size(); k++)
	{
		vector<PxVec3> verts = MakeRandomVertexList();
		mesh_list[k] = MakePxMeshFromVertexList(verts);
	}

	// Calculate a position grid that fits in the regolith bed
	PxBoxGeometry geometry;
	PxReal X_max = munch::params.bedWidth;
	PxReal Z_max = munch::params.bedWidth;	
	PxReal Y_max = munch::params.bedHeight;
	
	PxReal dl = munch::regolith.size1*2.06;
	vector<PxVec3> positions;
	unsigned int nbGrains = munch::regolith.totalNumber;
	PxReal x, y, z;
	PxReal clearance = dl + geometry.halfExtents.minElement();
	x = y = z = clearance;
	while (nbGrains--)
	{
		positions.push_back(PxVec3(x,y,z));
		x += dl;
		if (x > X_max - clearance) {x=clearance; z+=dl;}
		if (z > Z_max - clearance) {x=clearance; z=clearance; y+=dl;}
		if (y > Y_max - clearance) ncc__error("Can't fit grains in container, experiment aborted.");
	}

	// Place the grains in the bed
	for (unsigned int k=0; k<positions.size(); k++)
	{
		// Find the corner of the bed
		munch::VIPs.bedBottom->getBoxGeometry(geometry);
		PxVec3 center = munch::VIPs.bedBottom->getLocalPose().p;
		PxVec3 corner = center - geometry.halfExtents;
		// Create a mesh geometry from one of the meshes
		PxMeshScale meshScale;
		meshScale.scale = PxVec3(munch::regolith.size1);
		unsigned int m_i = k%mesh_list.size();
		PxConvexMeshGeometry meshGeometry(mesh_list[m_i],meshScale);

		// Create an actor with this geometry
		PxRigidDynamic* aGrain = PxCreateDynamic(*gPhysX.mPhysics,PxTransform(positions[k]+corner),meshGeometry,*gPhysX.mDefaultMaterial,gExp.defGrainDensity);
		if (aGrain)
		{
			if (gPhysX.props.sleepThreshold > -1) aGrain->setSleepThreshold(0.5*gPhysX.props.sleepThreshold*gPhysX.props.sleepThreshold);
			if (gPhysX.props.angularDamping > -1) aGrain->setAngularDamping(gPhysX.props.angularDamping);
			if (gPhysX.props.linearDamping  > -1) aGrain->setLinearDamping(gPhysX.props.linearDamping);
			aGrain->setName("regolith");
			RandOrientActor(aGrain);
			gPhysX.mScene->addActor(*aGrain);
		}
	}
}
void CreateFillBoxExperiment()
{
	munch::CreateTheBox();
	gCamera.pos = PxVec3(0,munch::params.bedWidth,3*munch::params.bedWidth);
	munch::FillTheBox();
}
void munch::GravitateSelf()
{
	if (gCUDA.cudaCapable)
		munch::GravitateOnDevice();
	else
		munch::GravitateOnHost();
}
void munch::GravitateOnDevice()
{
	munch::GravitateOnHost(); //TODO: include device implementation when needed, for speed
}
void munch::GravitateOnHost()
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

			float s = munch::units.bigG * bodies[4*j+0] * bodies[4*k+0] * invDistCube;

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
void munch::CreateRoxMuncher()
{
	// Define wall dimension and thickness
	PxReal BedWidth = munch::params.bedWidth;
	PxReal BedHeight = munch::params.bedHeight;
	PxReal muncherHeight = munch::params.muncherHeight;
	PxReal muncherWidth = munch::params.muncherWidth;
	PxReal muncherLength = munch::params.muncherLength;
	PxReal thickness = 0.005*BedWidth;
	PxReal muncherStartY = BedWidth/2;
	PxReal startMuncherAngle = PxPi/6;

	// We'll make the containment with a kinematic actor
	PxRigidDynamic* theBed = gPhysX.mPhysics->createRigidDynamic(PxTransform(PxVec3(0)));
	if (!theBed)
		ncc__error("actor creation failed!");
	theBed->setRigidDynamicFlag(PxRigidDynamicFlag::eKINEMATIC, true);

	PxRigidDynamic* MuncherLeft = gPhysX.mPhysics->createRigidDynamic(PxTransform(PxVec3(0)));
	if (!MuncherLeft)
		ncc__error("actor creation failed!");
	MuncherLeft->setRigidDynamicFlag(PxRigidDynamicFlag::eKINEMATIC, true);

	PxRigidDynamic* MuncherRight = gPhysX.mPhysics->createRigidDynamic(PxTransform(PxVec3(0)));
	if (!MuncherRight)
		ncc__error("actor creation failed!");
	MuncherRight->setRigidDynamicFlag(PxRigidDynamicFlag::eKINEMATIC, true);

	// Define sides
	PxBoxGeometry bed_bottom(BedWidth/2,thickness/2,BedWidth/2);
	PxBoxGeometry bed_side(BedWidth/2,thickness/2,BedHeight/2);
	PxBoxGeometry bed_side2(BedHeight/2,thickness/2,BedWidth/2);
	PxBoxGeometry muncher_side1(thickness/2,muncherWidth/2,muncherLength/2);
	PxBoxGeometry muncher_side2(muncherHeight/2,muncherWidth/2,thickness/2);
	PxBoxGeometry muncher_side3(muncherHeight/2,thickness/2,muncherLength/2);
	PxMaterial* defmat=gPhysX.mDefaultMaterial;

	// Attach the sides
	// Regolith Bed
	munch::VIPs.bedBottom = theBed->createShape(bed_bottom,*defmat); // bottom
	theBed->createShape(bed_side,*defmat,PxTransform(PxVec3(0,BedHeight/2,-BedWidth/2),PxQuat(PxPi/2,PxVec3(1,0,0)))); // back wall
	theBed->createShape(bed_side,*defmat,PxTransform(PxVec3(0,BedHeight/2,BedWidth/2),PxQuat(PxPi/2,PxVec3(1,0,0)))); // front wall
	theBed->createShape(bed_side2,*defmat,PxTransform(PxVec3(-BedWidth/2,BedHeight/2,0),PxQuat(PxPi/2,PxVec3(0,0,1)))); // left wall
	theBed->createShape(bed_side2,*defmat,PxTransform(PxVec3(BedWidth/2,BedHeight/2,0),PxQuat(PxPi/2,PxVec3(0,0,1)))); // right wall
	
	// Left Muncher
	MuncherLeft->createShape(muncher_side1,*defmat,PxTransform(PxVec3(-muncherHeight,muncherStartY,0),PxQuat(0,PxVec3(0,0,1)))); // left wall
	MuncherLeft->createShape(muncher_side2,*defmat,PxTransform(PxVec3(-muncherHeight/2,muncherStartY,muncherLength/2),PxQuat(0,PxVec3(0,0,1)))); // front wall
	MuncherLeft->createShape(muncher_side2,*defmat,PxTransform(PxVec3(-muncherHeight/2,muncherStartY,-muncherLength/2),PxQuat(0,PxVec3(0,0,1)))); // back wall
	MuncherLeft->createShape(muncher_side3,*defmat,PxTransform(PxVec3(-muncherHeight/2,muncherStartY-(muncherWidth/2),0),PxQuat(0,PxVec3(0,0,1)))); // bottom wall
	MuncherLeft->createShape(muncher_side3,*defmat,PxTransform(PxVec3(-muncherHeight/2,muncherStartY+(muncherWidth/2),0),PxQuat(0,PxVec3(0,0,1)))); // top wall
	
	// Right Muncher
	MuncherRight->createShape(muncher_side1,*defmat,PxTransform(PxVec3(muncherHeight,muncherStartY,0),PxQuat(0,PxVec3(0,0,1)))); // right wall
	MuncherRight->createShape(muncher_side2,*defmat,PxTransform(PxVec3(muncherHeight/2,muncherStartY,muncherLength/2),PxQuat(0,PxVec3(0,0,1)))); // front wall
	MuncherRight->createShape(muncher_side2,*defmat,PxTransform(PxVec3(muncherHeight/2,muncherStartY,-muncherLength/2),PxQuat(0,PxVec3(0,0,1)))); // back wall
	MuncherRight->createShape(muncher_side3,*defmat,PxTransform(PxVec3(muncherHeight/2,muncherStartY-(muncherWidth/2),0),PxQuat(0,PxVec3(0,0,1)))); // bottom wall
	MuncherRight->createShape(muncher_side3,*defmat,PxTransform(PxVec3(muncherHeight/2,muncherStartY+(muncherWidth/2),0),PxQuat(0,PxVec3(0,0,1)))); // top wall

	// Name, color, and register the bed
	theBed->setName("the_bed");
	gColors.colorBucket.push_back(vector<GLubyte>(3));
	gColors.colorBucket.back()[0] = ncc::rgb::rRed[0];
	gColors.colorBucket.back()[1] = ncc::rgb::rRed[1];
	gColors.colorBucket.back()[2] = ncc::rgb::rRed[2];
	theBed->userData = &(gColors.colorBucket.back()[0]);
	gPhysX.mScene->addActor(*theBed);
	munch::VIPs.theBed = theBed;

	// Name, color, and register MuncherLeft
	MuncherLeft->setName("#MuncherLeft");
	gColors.colorBucket.push_back(vector<GLubyte>(3));
	gColors.colorBucket.back()[0] = ncc::rgb::bDodgerBlue[0];
	gColors.colorBucket.back()[1] = ncc::rgb::bDodgerBlue[1];
	gColors.colorBucket.back()[2] = ncc::rgb::bDodgerBlue[2];
	MuncherLeft->userData = &(gColors.colorBucket.back()[0]);
	gPhysX.mScene->addActor(*MuncherLeft);
	munch::VIPs.MuncherLeft = MuncherLeft;

	// Name, color, and register MuncherRight
	MuncherRight->setName("#MuncherRight");
	gColors.colorBucket.push_back(vector<GLubyte>(3));
	gColors.colorBucket.back()[0] = ncc::rgb::bDodgerBlue[0];
	gColors.colorBucket.back()[1] = ncc::rgb::bDodgerBlue[1];
	gColors.colorBucket.back()[2] = ncc::rgb::bDodgerBlue[2];
	MuncherRight->userData = &(gColors.colorBucket.back()[0]);
	gPhysX.mScene->addActor(*MuncherRight);
	munch::VIPs.MuncherRight = MuncherRight;


	// TO DO: Make a proper joint hinge 
	// Going to just open the muncher with rotations and translations for now...
	
	PxTransform leftNewPose;
	leftNewPose.p = PxVec3(-0.55,0,0);
	leftNewPose.q = PxQuat(startMuncherAngle,PxVec3(0.0f,0.0f,-1.0f));
	munch::VIPs.MuncherLeft->setKinematicTarget(leftNewPose);

	PxTransform rightNewPose;
	rightNewPose.p = PxVec3(0.55,0,0);
	rightNewPose.q = PxQuat(startMuncherAngle,PxVec3(0.0f,0.0f,1.0f));
	munch::VIPs.MuncherRight->setKinematicTarget(rightNewPose);

}

// Unused functions

void UpArrowAction()
{
	Reveille();
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


void RenderOtherStuff()
{

}


// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif