/////////////////////////////////////////////////////////////////////////////////
// Source file for project LabScaleImpacts. This project implements some low
// velocity impacts into regolith in uniform 1g gravity. Idea is of course to
// benchmark against lab tests.
//
// Author: Naor Movshovits (nmovshov at google dot com)
/////////////////////////////////////////////////////////////////////////////////

#include "ARSS.h" // all supporting global, non project-specific entities
#include "LabScaleImpacts.h" // global project-specific entities

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

// Experiment specific functions called from ARSS.cpp
void CreateExperiment()
{
	switch (labscale::eExperimentType)
	{
	case labscale::eHOLSAPPLE_1:
		labscale::CreateSpringerExperiment();
		break;
	case labscale::eBAD_EXPERIMENT_TYPE: // intentional fall through!
	default:
		ncc__error("Unknown experiment type. Experiment aborted.\a");
	}
}
void RebootExperiment()
{
	gSim.isRunning=false;
	DestroyPhysX();
	InitPhysX();
	CreateExperiment();
}
void LogExperiment()
{
	switch (labscale::eExperimentType)
	{
	case labscale::eHOLSAPPLE_1:
		labscale::LogSpringerExperiment();
		break;
	case labscale::eBAD_EXPERIMENT_TYPE:
		ncc__warning("Unknown experiment type. Nothing logged.");
		break;
	default:
		break;
	}
}
void ControlExperiment()
{
	switch (labscale::eExperimentType)
	{
	case labscale::eHOLSAPPLE_1:
		labscale::ControlSpringerExperiment();
		break;
	case labscale::eBAD_EXPERIMENT_TYPE:
		ncc__warning("Unknown experiment type.");
		break;
	default:
		break;
	}
}
void CustomizeScene(PxSceneDesc &sceneDesc)
{
	
}
void CustomizeRun(int argc, char** argv)
{
	// Read experiment-specific program options from file.
	if (!ConfigExperimentOptions())
		ncc__warning("Could not find ARSS config file. Attempting to continue with all default options.\a");
	
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
	if      (strcmp(buf,"holsapple_1")==0)
		labscale::eExperimentType=labscale::eHOLSAPPLE_1;
	else if (strcmp(buf,"springer")==0)
		labscale::eExperimentType=labscale::eHOLSAPPLE_1;
	else
		labscale::eExperimentType=labscale::eBAD_EXPERIMENT_TYPE;

	// Collider experiment parameters
	ncc::GetStrPropertyFromINIFile("experiment","spin_magnitude","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	//labscale::spinMag = atof(buf);

	// Physical parameters
	ncc::GetStrPropertyFromINIFile("units","little_g","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	labscale::units.littleG = atof(buf);
	
	return true;
}
void CustomizeGLUT()
{
	gCamera.pos.z=6;
}
void CustomizeHUD()
{
	labscale::hudMsgs.systemDiag1 = gHUD.hud.AddElement("",0.9,0.1);
	labscale::hudMsgs.systemDiag2 = gHUD.hud.AddElement("",0.9,0.2);
	labscale::hudMsgs.systemDiag3 = gHUD.hud.AddElement("",0.9,0.3);
	labscale::hudMsgs.actorDiag   = gHUD.hud.AddElement("",0.9,0.3);
}
void RefreshCustomHUDElements()
{
	char buf[MAX_CHARS_PER_NAME];
	int	  ch2px	    = 18; // hud uses 18pt font
	float px2width  = 1.0/glutGet(GLUT_WINDOW_WIDTH);
	float scrPos = glutGet(GLUT_WINDOW_WIDTH);

	switch (labscale::eExperimentType)
	{
	case labscale::eHOLSAPPLE_1:
		// Diagnostic stuff
		sprintf_s(buf,MAX_CHARS_PER_NAME,"Diagnostic 1 = %d",42);
		scrPos = 0.8;
		gHUD.hud.SetElement(labscale::hudMsgs.systemDiag1,buf,scrPos,0.1);
		break;
	}
}
void FireAction()
{
	//int x=0;
}
void PrintDebug()
{
	//int x=0;
}
void ApplyCustomInteractions()
{
	if (gSim.isRunning)
	{
		switch (labscale::eExperimentType)
		{
		case labscale::eHOLSAPPLE_1:
			break;
		}
	}
}
void RenderOtherStuff()
{
	switch (labscale::eExperimentType)
	{
	case labscale::eHOLSAPPLE_1:
		//DrawArrow(PxVec3(0),labscale::tumbler.L_now*10 ,0.12,4.0,ncc::rgb::rDarkRed);
		break;
	}
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

// Project namespace functions
void labscale::CreateSpringerExperiment()
/* Swing a mass on a spring. Check integration and scaling.*/
{
	// Put a box at the end of a pulled spring (sold separately)
	labscale::VIPs.ball1 = CreateRubbleGrain(PxVec3(10,0,0),eBOX_GRAIN,1,*gPhysX.mDefaultMaterial);

	// Move the camera to a better vantage point and turn on a grid
	gCamera.pos = PxVec3(0,0,14);
	gDebug.bXYGridOn = true;

	// Start a log
	if (gRun.outputFrequency)
	{
		ostringstream header;
		header << "# This is the run log of " << gRun.baseName << endl;
		header << "# Experiment type: SPRINGER (" << labscale::eExperimentType << ")" << endl;
		header << "# Time step used = " << gSim.timeStep << endl;
		header << "# Columns are (values in code units):" << endl;
		header << "# [t]    [x]    [v]" << endl;
		ofstream fbuf(gRun.outFile.c_str(),ios::trunc);
		if (!fbuf.is_open())
			ncc__error("Could not start a log. Experiment aborted.\a\n");
		fbuf << header.str() << endl;
	}

	// Start the action
	gSim.isRunning=true;
	gSim.bPause=false;
	gSim.codeTime = 0.0f;
	RefreshHUD();
}
void labscale::ControlSpringerExperiment()
{
	// Check for stop condition (target time interpreted as number of cycles)
	if (!gSim.targetTime && (gSim.codeTime >= 2*PxPi/PxSqrt(labscale::units.springK)))
		gSim.isRunning = false;
}
void labscale::LogSpringerExperiment()
{
	// Let's do this without worrying about minimizing access to disk - it will be so much easier!
	ofstream fbuf(gRun.outFile.c_str(),ios::app);
	if (!fbuf.is_open()) {
		ncc__warning("Could not open log. Nothing written!\a\n");
		return;
	}

	// Collect
	PxReal t = gSim.codeTime;
	PxReal x = labscale::VIPs.ball1->getGlobalPose().p.x;
	PxReal v = labscale::VIPs.ball1->getLinearVelocity().x;

	// Write
	fbuf.setf(ios::fixed);
	fbuf << setw(8) << t << "    " << setw(8) << x << "    " << setw(8) << v << endl;
	fbuf.close();
}

// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif