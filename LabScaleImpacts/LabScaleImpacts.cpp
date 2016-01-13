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
    // All experiments are in lab setup, give the lab gravity
    gPhysX.mScene->setGravity(PxVec3(0,-labscale::units.littleG,0));

    // Now dispatch to specific setup
    switch (labscale::eExperimentType)
    {
    case labscale::eHOLSAPPLE1:
        if (labscale::eExperimentSubtype==labscale::eFILL_BOX)
            labscale::CreateFillBoxExperiment();
        else if (labscale::eExperimentSubtype==labscale::ePENETRATOR)
            labscale::CreatePenetratorExperiment();
		else if (labscale::eExperimentSubtype==labscale::eTILT_BOX)
			labscale::CreateTiltBoxExperiment();
        else
            ncc__error("Unknown experiment type. Experiment aborted.\a");
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
    if (labscale::eExperimentType == labscale::eHOLSAPPLE1 && labscale::eExperimentSubtype == labscale::ePENETRATOR)
        labscale::LogPenetratorExperiment();
}
void ControlExperiment()
{
    switch (labscale::eExperimentType)
    {
    case labscale::eHOLSAPPLE1:
        if (labscale::eExperimentSubtype == labscale::eFILL_BOX)
            labscale::ControlFillBoxExperiment();
		else if (labscale::eExperimentSubtype == labscale::eTILT_BOX)
			labscale::ControlTiltBoxExperiment();
        else if (labscale::eExperimentSubtype == labscale::ePENETRATOR)
            labscale::ControlPenetratorExperiment();
        else ncc__error("Unknown experiment subtype.");
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

    // Experiment type and subtype
    ncc::GetStrPropertyFromINIFile("experiment","experiment_type","",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    if      (strcmp(buf,"holsapple1")==0)
        labscale::eExperimentType=labscale::eHOLSAPPLE1;
    else
        labscale::eExperimentType=labscale::eBAD_EXPERIMENT_TYPE;

    ncc::GetStrPropertyFromINIFile("experiment","experiment_subtype","",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    if      (strcmp(buf,"fill_container")==0)
        labscale::eExperimentSubtype=labscale::eFILL_BOX;
	else if (strcmp(buf,"tilt_container")==0)
		labscale::eExperimentSubtype=labscale::eTILT_BOX;
    else if (strcmp(buf,"penetrator")==0)
        labscale::eExperimentSubtype=labscale::ePENETRATOR;
    else
        labscale::eExperimentSubtype=labscale::eBAD_EXP_SUBTYPE;

    // Regolith container parameters
    ncc::GetStrPropertyFromINIFile("container","diameter","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    labscale::reg_box.diameter = atof(buf);
    ncc::GetStrPropertyFromINIFile("container","fill_height","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    labscale::reg_box.fillHeight = atof(buf);
	ncc::GetStrPropertyFromINIFile("container","pour_rate","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	labscale::reg_box.pourRate = atof(buf);
	ncc::GetStrPropertyFromINIFile("container","tilt_rate","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	labscale::reg_box.tiltRate = atof(buf);

    // Regolith parameters
    ncc::GetStrPropertyFromINIFile("regolith","diameter","0.0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    labscale::regolith.diameter = atof(buf);
    ncc::GetStrPropertyFromINIFile("regolith","material_density","0.0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    labscale::regolith.materialDensity = atof(buf);
    labscale::regolith.nbGrains = ncc::GetIntPropertyFromINIFile("regolith","nb_grains",1,gRun.iniFile.c_str());

    // Impactor parameters
    ncc::GetStrPropertyFromINIFile("impactor","diameter","0.0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    labscale::impactor.diameter = atof(buf);
    ncc::GetStrPropertyFromINIFile("impactor","material_density","0.0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    labscale::impactor.materialDensity = atof(buf);
    ncc::GetStrPropertyFromINIFile("impactor","speed","0.0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    labscale::impactor.speed = atof(buf);

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
    labscale::hudMsgs.systemDiag1 = gHUD.hud.AddElement("",0.8,0.04);
    labscale::hudMsgs.systemDiag2 = gHUD.hud.AddElement("",0.8,0.08);
    labscale::hudMsgs.systemDiag3 = gHUD.hud.AddElement("",0.8,0.12);
    labscale::hudMsgs.actorDiag   = gHUD.hud.AddElement("",0.8,0.16);
}
void RefreshCustomHUDElements()
{
    char buf[MAX_CHARS_PER_NAME];
    int	  ch2px	    = 18; // hud uses 18pt font
    float px2width  = 1.0/glutGet(GLUT_WINDOW_WIDTH);
    float scrPos = glutGet(GLUT_WINDOW_WIDTH);

    if (labscale::eExperimentType == labscale::eHOLSAPPLE1 && labscale::eExperimentSubtype == labscale::eFILL_BOX)
    {
		// Diagnostic 1: regolith particles
		PxU32 nbActors = gPhysX.mScene->getNbActors(gPhysX.roles.dynamics) - 1; // don't count the box
		PxU32 nbSleep = CountSleepers() - 1; // don't count the box
		sprintf_s(buf,MAX_CHARS_PER_NAME,"# particles (sleeping) = %u (%u)",nbActors,nbSleep);
		scrPos = 1.0 - strlen(buf)*ch2px*px2width*0.5;
		gHUD.hud.SetElement(labscale::hudMsgs.systemDiag1,buf,scrPos,0.04);

		// Diagnostic 2: surface level
        PxReal boxFloor = labscale::VIPs.container->getGlobalPose().p.y + labscale::reg_box.diameter/20; // width hard coded for now :(
        PxReal regSurf = labscale::getRegolithSurface() - boxFloor;
		sprintf_s(buf,MAX_CHARS_PER_NAME,"Regolith surface = %f cm (%d layers)",regSurf*100,int(regSurf/labscale::regolith.diameter));
		scrPos = 1.0 - strlen(buf)*ch2px*px2width*0.5;
		gHUD.hud.SetElement(labscale::hudMsgs.systemDiag2,buf,scrPos,0.08);
    }

    if (labscale::eExperimentType == labscale::eHOLSAPPLE1 && labscale::eExperimentSubtype == labscale::ePENETRATOR)
    {
        // Diagnostic 1: impactor speed
        PxReal v = labscale::VIPs.ball1->getLinearVelocity().magnitude();
        sprintf_s(buf,MAX_CHARS_PER_NAME,"|v| = %f",v);
        scrPos = 1.0 - strlen(buf)*ch2px*px2width*0.5;
        gHUD.hud.SetElement(labscale::hudMsgs.systemDiag1,buf,scrPos,0.04);

        // Diagnostic 2: penetration depth
        PxReal d = (labscale::VIPs.ball1->getGlobalPose().p.y - labscale::impactor.diameter/2.0 - labscale::impactor.iniSurface);
        sprintf_s(buf,MAX_CHARS_PER_NAME,"Penetration depth = %f cm",d*100);
        scrPos = 1.0 - strlen(buf)*ch2px*px2width*0.5;
        gHUD.hud.SetElement(labscale::hudMsgs.systemDiag2,buf,scrPos,0.08);
    }
    
}
void FireAction()
{
    if (labscale::eExperimentType==labscale::eHOLSAPPLE1 && labscale::eExperimentSubtype==labscale::ePENETRATOR)
    {
		labscale::VIPs.ball1->setRigidDynamicFlag(PxRigidDynamicFlag::eKINEMATIC, false);
		labscale::VIPs.ball1->setLinearVelocity(PxVec3(0,-labscale::impactor.speed,0));
    }

	if (labscale::eExperimentType==labscale::eHOLSAPPLE1 && labscale::eExperimentSubtype==labscale::eTILT_BOX)
	{
		labscale::reg_box.bTilt = !labscale::reg_box.bTilt;
	}
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
        case labscale::eHOLSAPPLE1:
            break;
        }
    }
}
void RenderOtherStuff()
{
    switch (labscale::eExperimentType)
    {
    case labscale::eHOLSAPPLE1:
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
void labscale::CreatePenetratorExperiment()
/* Throw a steel ball at glass regolith; measure penetration depth.*/
{
    // Load pre-settled regolith or create empty container
    if (gRun.loadSceneFromFile.empty())
        labscale::CreateRegolithContainer();
    else
    {
        if(!LoadSceneFromFile(gRun.loadSceneFromFile)) ncc__error("Failed to load partially filed box.\a\n");
        PxActor* theBox = FindNamedActor("the_box");
        if (theBox)
        {
            ColorActor(theBox, ncc::rgb::rRed);
            labscale::VIPs.container = theBox->isRigidDynamic();
        }
    }

    // Create impactor material (steel, for now)
    PxMaterial* steel = gPhysX.mPhysics->createMaterial(0.5, 0.5, 0.5);
    if (!steel)
        ncc__error("\aPxPhysics::createMaterial() failed!");
    PxReal rho = labscale::impactor.materialDensity;

    // Ready, aim impactor (will fire manually)
    PxReal radius = labscale::impactor.diameter/2;
	PxReal hLaunch = labscale::reg_box.diameter;
    labscale::VIPs.ball1 = CreateRubbleGrain(PxVec3(0,hLaunch,0),eSPHERE_GRAIN,radius,*steel,rho);
    labscale::VIPs.ball1->setName("impactor");
    labscale::VIPs.ball1->setRigidDynamicFlag(PxRigidDynamicFlag::eKINEMATIC, true);
    ColorActor(labscale::VIPs.ball1, ncc::rgb::gDarkOlive);
    labscale::impactor.iniSurface = getRegolithSurface(); // target surface in global frame

    // Adjust camera, grid, display
    gCamera.pos.x = 0.0;
    gCamera.pos.y = labscale::reg_box.fillHeight*1.4;
    gCamera.pos.z = labscale::reg_box.diameter*1.6;
    gDebug.bXYGridOn = true;

    // Start a log
    if (gRun.outputFrequency)
    {
        ostringstream header;
        header << "# This is the run log of " << gRun.baseName << endl;
        header << "# Experiment type: HOLSAPPLE1 (" << labscale::eExperimentType << ")" << endl;
        header << "# Time step used = " << gSim.timeStep << endl;
        header << "# Top of regolith surface = " << labscale::impactor.iniSurface << endl;
        header << "# Columns are (values in code units, d is penetration depth):" << endl;
        header << "# [t]    [d]    [x]    [y]    [z]" << endl;
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
void labscale::ControlPenetratorExperiment()
{
    // Finish when everyone is asleep, after impact!
	if (!(labscale::VIPs.ball1->getRigidDynamicFlags() & PxRigidDynamicFlag::eKINEMATIC))
	{
		if (CountSleepers() == gPhysX.mScene->getNbActors(gPhysX.roles.dynamics))
		{
		    gSim.isRunning = false;
		    SaveSceneToRepXDump();
		}
	}
}
void labscale::LogPenetratorExperiment()
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
    PxReal y = labscale::VIPs.ball1->getGlobalPose().p.y;
    PxReal z = labscale::VIPs.ball1->getGlobalPose().p.z;
    PxReal d = y - labscale::impactor.diameter/2.0 - labscale::impactor.iniSurface;

    // Write
    fbuf.setf(ios::fixed);
    fbuf << setw(8) << t << "    ";
    fbuf << setw(8) << showpos << d << "    ";
    fbuf << setw(8) << showpos << x << "    ";
    fbuf << setw(8) << showpos << y << "    ";
    fbuf << setw(8) << showpos << z << endl;
    fbuf.close();
}
void labscale::CreateFillBoxExperiment()
/* Put a box on the ground ready to be filled with regolith.*/
{
    // Put a box on the floor
	if (gRun.loadSceneFromFile.empty())
		CreateRegolithContainer();
	else
	{
		if(!LoadSceneFromFile(gRun.loadSceneFromFile)) ncc__error("Failed to load partially filed box.\a\n");
		PxActor* theBox = FindNamedActor("the_box");
		if (theBox)
		{
			ColorActor(theBox,ncc::rgb::rRed);
			labscale::VIPs.container = theBox->isRigidDynamic();
		}
	}

    // Adjust camera, grid, display
    gCamera.pos.x = 0.0;
    gCamera.pos.y = labscale::reg_box.fillHeight*1.4;
    gCamera.pos.z = labscale::reg_box.diameter*1.6;
    gDebug.bXYGridOn = true;
    
    // Start the action, regolith poured in runtime
    gSim.isRunning=true;
    gSim.bPause=false;
    gSim.codeTime = 0.0f;
    RefreshHUD();

}
void labscale::CreateRegolithContainer()
{
    // Geometry variables
    PxReal box_d = labscale::reg_box.diameter;
    PxReal box_h = labscale::reg_box.fillHeight*2;
    PxReal wall_dh = box_d/20;

    // We'll make the regolith container with a kinematic actor
    PxRigidDynamic* theBox = gPhysX.mPhysics->createRigidDynamic(PxTransform(PxVec3(0,0.05*wall_dh,0)));
    if (!theBox)
        ncc__error("Actor creation failed!");
    theBox->setRigidDynamicFlag(PxRigidDynamicFlag::eKINEMATIC, true);

    // Define sides
    PxBoxGeometry box_bottom(box_d/2,wall_dh/2,box_d/2);
    PxBoxGeometry box_side(wall_dh/2,box_h/2,box_d/2);
    PxMaterial* defmat=gPhysX.mDefaultMaterial;

    // Attach the sides, making front wall invisible
    theBox->createShape(box_bottom,*defmat,PxTransform(PxVec3(0,wall_dh/2,0))); // the bottom
    theBox->createShape(box_side,*defmat,PxTransform(PxVec3(-box_d/2,box_h/2,0),PxQuat(0,PxVec3(0,1,0)))); // left wall
    theBox->createShape(box_side,*defmat,PxTransform(PxVec3( box_d/2,box_h/2,0),PxQuat(0,PxVec3(0,1,0)))); // right wall
    theBox->createShape(box_side,*defmat,PxTransform(PxVec3(0,box_h/2,-box_d/2),PxQuat(PxPi/2,PxVec3(0,1,0)))); // back wall
    PxShape* fwall =  theBox->createShape(box_side,*defmat,PxTransform(PxVec3(0,box_h/2, box_d/2),PxQuat(PxPi/2,PxVec3(0,1,0)))); // front wall
    fwall->setName("~fwall");

    // Name, color, and register the container
    theBox->setName("the_box");
    ColorActor(theBox, ncc::rgb::rRed);
    gPhysX.mScene->addActor(*theBox);
    labscale::VIPs.container = theBox;

}
void labscale::ControlFillBoxExperiment()
{
    // Pour regolith, one by one every second
    int nbGrains = gPhysX.mScene->getNbActors(gPhysX.roles.dynamics) - 1;
    static PxReal poured_time = 0;
    if ((gSim.codeTime - poured_time) > (1/labscale::reg_box.pourRate) && nbGrains < labscale::regolith.nbGrains)
    {
        // Pour a grain
        PxRigidDynamic* grain = CreateRegolithGrain();
        RandLaunchActor(grain,2);
        PxVec3 v = grain->getLinearVelocity();
        grain->setLinearVelocity(PxVec3(v.x,0,v.z));

        // Reset timer
        poured_time = gSim.codeTime;
    }

    // When done save scene and stop
    if (nbGrains >= labscale::regolith.nbGrains && CountSleepers() == gPhysX.mScene->getNbActors(gPhysX.roles.dynamics))
    {
        gSim.isRunning = false;
        SaveSceneToRepXDump();
    }

    // Hack if rebooting (F10 was pressed) - not really important
    if ((gSim.codeTime - poured_time) < 0) poured_time = gSim.codeTime;
}
PxRigidDynamic * labscale::CreateRegolithGrain()
{
    PxMaterial* glass = gPhysX.mPhysics->createMaterial(0.5, 0.5, 0.5); // "glass"
    PxReal radius = labscale::regolith.diameter/2;
    PxReal rho = labscale::regolith.materialDensity;
    PxRigidDynamic* actor = CreateRubbleGrain(PxVec3(0,1.5*labscale::reg_box.fillHeight,0),eSPHERE_GRAIN,radius,*glass,rho);
    actor->setName("regolith");

    if (!actor)
        ncc__error("actor creations failed");
    return actor;
}
void labscale::CreateTiltBoxExperiment()
/* Tilt the box back and forth to smooth regolith layer*/
{
	// Put a box on the floor
	if (gRun.loadSceneFromFile.empty())
		CreateRegolithContainer();
	else
	{
		if(!LoadSceneFromFile(gRun.loadSceneFromFile)) ncc__error("Failed to load partially filed box.\a\n");
		PxActor* theBox = FindNamedActor("the_box");
		if (theBox)
		{
			ColorActor(theBox,ncc::rgb::rRed);
			labscale::VIPs.container = theBox->isRigidDynamic();
		}
	}

	// Adjust camera, grid, display
	gCamera.pos.x = 0.0;
	gCamera.pos.y = labscale::reg_box.fillHeight*1.4;
	gCamera.pos.z = labscale::reg_box.diameter*1.6;
	gDebug.bXYGridOn = true;

	// Start the action, tilt action handled in control function
	labscale::reg_box.bTilt = false;
	gSim.isRunning=true;
	gSim.bPause=false;
	gSim.codeTime = 0.0f;
	RefreshHUD();
}
void labscale::ControlTiltBoxExperiment()
{
	static PxReal dir = 1;
	if (labscale::reg_box.bTilt)
	{
		PxReal omega = (PxPi/4.0)*labscale::reg_box.tiltRate;
		PxReal dtau = omega*gSim.timeStep;
		PxTransform old = labscale::VIPs.container->getGlobalPose();
		PxReal cur = old.q.getAngle();
		bool inBounds = (cur <= PxPi/32.0) && (cur >= -PxPi/32.0); // reverse direction?
		if (!inBounds) dir*=-1;
		PxTransform rot(PxQuat(dir*dtau, PxVec3(0,0,1)));
		labscale::VIPs.container->setKinematicTarget(rot*old);
	}
}
PxReal labscale::getRegolithSurface()
{
    PxReal regSurf = -PX_MAX_REAL;
    PxU32 nbActors = gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,MAX_ACTORS_PER_SCENE);
    while (nbActors--)
    {
        PxRigidDynamic* actor = gPhysX.cast[nbActors]->isRigidDynamic();
        if (actor && strcmp(actor->getName(), "regolith") == 0)
        {
            PxReal y = actor->getGlobalPose().p.y;
            if (y>regSurf) regSurf = y;
        }
    }

    return (regSurf + labscale::regolith.diameter/2.0);
}

// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif