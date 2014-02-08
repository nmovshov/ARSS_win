///////////////////////////////////////////////////////////////////////////////
// Source file for project Cliff. This project implements modeling of several variants
// of the cliff collapse problem. A review of the problem description and the expected
// relevant model parameters can be found in Holsapple (2012) and references therein.
//
// Author: Me (Naor)
///////////////////////////////////////////////////////////////////////////////

#include "ARSS.h" // all supporting global, non project-specific entities
#include "Cliff.h" // global project-specific entities

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
	if (!fp.good()) return false;
	fp.close();

	// Read in parameters by group, file.ini style
	char buf[MAX_CHARS_PER_NAME];

	// Experiment Type
	ncc::GetStrPropertyFromINIFile("experiment","experiment_type","",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if      (strcmp(buf,"rect_smooth")==0)
		cliff::eExperimentType=cliff::eRECT_SMOOTH_BASE_COLLAPSE;
	else if (strcmp(buf,"rect_rough")==0)
		cliff::eExperimentType=cliff::eRECT_ROUGH_BASE_COLLAPSE;
	else if (strcmp(buf,"rect_fill")==0)
		cliff::eExperimentType=cliff::eRECT_FILL;
	else
		cliff::eExperimentType=cliff::eBAD_EXPERIMENT_TYPE;

	// Experiment parameters
	cliff::params.H0	= ncc::GetIntPropertyFromINIFile("experiment","initial_height",0,gRun.iniFile.c_str());
	cliff::params.L0	= ncc::GetIntPropertyFromINIFile("experiment","initial_length",0,gRun.iniFile.c_str());
	cliff::params.NMax	= ncc::GetIntPropertyFromINIFile("experiment","max_nb_grains", 0,gRun.iniFile.c_str());

	return true;
}
void CustomizeScene(PxSceneDesc &sceneDesc)
{

}
void CustomizeGLUT()
{
	
}
void CustomizeHUD()
{
	cliff::hudMsgs.systemDiag = gHUD.hud.AddElement("",0.9,0.1);
}
void RefreshCustomHUDElements()
{
	char buf[MAX_CHARS_PER_NAME];
	int   ch2px		= 18; // hud uses 18pt font
	float px2width	= 1.0/glutGet(GLUT_WINDOW_WIDTH);
	float scrPos	= glutGet(GLUT_WINDOW_WIDTH);

	switch (cliff::eExperimentType)
	{
	case cliff::eRECT_SMOOTH_BASE_COLLAPSE: // intentional fall through
	case cliff::eRECT_FILL:
		int nbGrains = gPhysX.mScene->getNbActors(gPhysX.roles.dynamics);
		int nbSleepers	= CountSleepers();
		sprintf_s(buf,MAX_CHARS_PER_NAME,"Grains = %d (%d sleeping)",nbGrains,nbSleepers);
		scrPos = 1.0 - strlen(buf)*ch2px*px2width*0.56;
		gHUD.hud.SetElement(cliff::hudMsgs.systemDiag,buf,scrPos,0.1);
		break;
	}

}
void FireAction()
{
	if (cliff::VIPs.floodGate)
	{
		gPhysX.mScene->removeActor(*cliff::VIPs.floodGate);
	}
}
void PrintDebug()
{

}
void ApplyCustomInteractions()
{
	static PxReal elapsed = 0.0f;
	switch (cliff::eExperimentType)
	{
	case cliff::eRECT_FILL:
		if ((gSim.codeTime - elapsed) > 1.0f)
		{
			if (gPhysX.mScene->getNbActors(gPhysX.roles.dynamics) < cliff::params.NMax)
			{
				FindExtremers();
				PxReal y = 2*gExp.defGrainSize;
				if (gExp.VIPs.extremers.upmost)
					y += gExp.VIPs.extremers.upmost->getGlobalPose().p.y;
				PxVec3 spout(0,y,0);
				PxRigidDynamic* grain = CreateRubbleGrain(spout,gExp.defGrainType,gExp.defGrainSize,*gPhysX.mDefaultMaterial,gExp.defGrainDensity,gExp.bUniformRubble,gExp.defGrainSizeScatter);
				if (grain)
				{
					RandOrientActor(grain);
					RandSpinActor(grain,2);
					RandLaunchActor(grain,6);
				}
				else
				{
					ncc__error("Grain creation failed, experiment aborted.\a");
				}
			} 
			if (CountSleepers() == cliff::params.NMax)
			{
				SaveSceneToRepXDump();
				gSim.isRunning = false;
				DestroyPhysX();
				exit(0);
			}

			elapsed = gSim.codeTime;
		}
		break;
	}
}
void RenderOtherStuff()
{

}
void CreateExperiment()
{
	switch (cliff::eExperimentType)
	{
	case cliff::eRECT_FILL:
		CreateGroundPlane();
		cliff::CreateContainment();
		break; // the "fill" action will happen through ApplyCustomInteractions()
	case cliff::eRECT_ROUGH_BASE_COLLAPSE: // intentional fall through, to be eventually replaced
	case cliff::eRECT_SMOOTH_BASE_COLLAPSE:
		if (!LoadSceneFromFile(gRun.loadSceneFromFile))
		{
			ncc__error("Scene could not be loaded. Check that repx file exists. Experiment aborted.\a");
		}
		else
		{
			// find the floodgate
			PxU32 walls = gPhysX.mScene->getActors(gPhysX.roles.statics,gPhysX.cast,MAX_ACTORS_PER_SCENE);
			while (walls--)
			{
				const char* buf = gPhysX.cast[walls]->getName();
				if (strcmp(buf,"~floodgate")==0)
					cliff::VIPs.floodGate = gPhysX.cast[walls]->isRigidStatic();
			}

			// place the camera
			FindExtremers();
			if (gExp.VIPs.extremers.inmost) // or whoever
			{
				PxReal H = gExp.VIPs.extremers.upmost->getGlobalPose().p.y;
				gCamera.pos.y = H;
				gCamera.pos.z = 3*H;
			}
		}
		break;
	case cliff::eBAD_EXPERIMENT_TYPE: // intentional fall through
	default:
		ncc__error("Unknown experiment type. Experiment aborted.\a");
	}

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

// cliff namespace functions
void cliff::CreateContainment()
{
	// First, define height and length scales based on grain size
	PxReal H=cliff::params.H0 * gExp.defGrainSize * 2;
	PxReal L=cliff::params.L0 * gExp.defGrainSize * 2;
	PxReal W=PxMax(L,H);
	if (gExp.defGrainType == eBOX_GRAIN || gExp.defGrainType == ePYRAMID_GRAIN)
	{
		H *= 0.5;
		L *= 0.5;
		W *= 0.5;
	}

	// Estimate the required number of grains
	PxRigidDynamic* tmp = CreateRubbleGrain(PxVec3(0),gExp.defGrainType,gExp.defGrainSize,*gPhysX.mDefaultMaterial,gExp.defGrainDensity,gExp.bUniformRubble,gExp.defGrainSizeScatter);
	if (tmp)
	{
		PxReal V = tmp->getMass()/gExp.defGrainDensity;
		PxU32 nb = (H*L*W)/V;
		if (nb > cliff::params.NMax || !gExp.bUniformRubble)
		{
			ncc__warning("Requested number of grains may not be sufficient to fill requested dimensions.");
		} 
		else
		{
			cliff::params.NMax = nb;
		}
		gPhysX.mScene->removeActor(*tmp);
	}

	// Now place the walls
	PxActor* anActor;

	// Left wall
	anActor = PxCreatePlane(*gPhysX.mPhysics,PxPlane(PxVec3(1,0,0),L/2),*gPhysX.mDefaultMaterial);
	if (!anActor)
		ncc__error("Containment wall creation failed. Experiment aborted.\a");
	anActor->setName("~leftwall"); // ~name is a quick way to request invisibility
	gPhysX.mScene->addActor(*anActor);

	// Back wall
	anActor = PxCreatePlane(*gPhysX.mPhysics,PxPlane(PxVec3(0,0,1),W/2),*gPhysX.mDefaultMaterial);
	if (!anActor)
		ncc__error("Containment wall creation failed. Experiment aborted.\a");
	anActor->setName("~backwall"); // ~name is a quick way to request invisibility
	gPhysX.mScene->addActor(*anActor);

	// Front wall
	anActor = PxCreatePlane(*gPhysX.mPhysics,PxPlane(PxVec3(0,0,-1),W/2),*gPhysX.mDefaultMaterial);
	if (!anActor)
		ncc__error("Containment wall creation failed. Experiment aborted.\a");
	anActor->setName("~frontwall"); // ~name is a quick way to request invisibility
	gPhysX.mScene->addActor(*anActor);

	// Finally, the flood gate
	anActor = PxCreatePlane(*gPhysX.mPhysics,PxPlane(PxVec3(-1,0,0),L/2),*gPhysX.mDefaultMaterial);
	if (!anActor)
		ncc__error("Containment wall creation failed. Experiment aborted.\a");
	anActor->setName("~floodgate"); // ~name is a quick way to request invisibility
	gPhysX.mScene->addActor(*anActor);
	cliff::VIPs.floodGate = anActor->isRigidStatic();

	// For convenience, place the camera for a nice viewpoint
	gCamera.pos.y = H;
	gCamera.pos.z = 3*H;
}
void cliff::BuildCliff() // TODO: finish this
{

}

// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif