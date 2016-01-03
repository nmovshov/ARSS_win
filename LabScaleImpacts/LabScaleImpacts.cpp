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
    case labscale::eHOLSAPPLE1:
        if (labscale::eExperimentSubtype==labscale::eFILL_BOX)
            labscale::CreateFillBoxExperiment();
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
    switch (labscale::eExperimentType)
    {
    case labscale::eHOLSAPPLE1:
        labscale::LogHolsapple1Experiment();
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
    case labscale::eHOLSAPPLE1:
        labscale::ControlHolsapple1Experiment();
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
    if (strcmp(buf,"fill_container")==0)
        labscale::eExperimentSubtype=labscale::eFILL_BOX;
    else
        labscale::eExperimentSubtype=labscale::eBAD_EXP_SUBTYPE;

    // Regolith container parameters
    ncc::GetStrPropertyFromINIFile("container","diameter","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    labscale::reg_box.diameter = atof(buf);
    ncc::GetStrPropertyFromINIFile("container","fill_height","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    labscale::reg_box.fillHeight = atof(buf);

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
   // int	  ch2px	    = 18; // hud uses 18pt font
    //float px2width  = 1.0/glutGet(GLUT_WINDOW_WIDTH);
    float scrPos = glutGet(GLUT_WINDOW_WIDTH);

    switch (labscale::eExperimentType)
    {
    case labscale::eHOLSAPPLE1:
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
void labscale::CreateHolsapple1Experiment()
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
void labscale::ControlHolsapple1Experiment()
{
    // Check for stop condition (target time interpreted as number of cycles)
    if (false)
        gSim.isRunning = false;
}
void labscale::LogHolsapple1Experiment()
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
void labscale::CreateFillBoxExperiment()
/* Put a box on the ground ready to be filled with regolith.*/
{
    // Give the lab a floor, with gravity
    CreateGroundPlane();

    // Put a box on the floor
    CreateRegolithContainer();

    // Adjust camera, grid, display
    gCamera.pos.x = 0.0;
    gCamera.pos.y = 1.0;
    gCamera.pos.z = 2.0;
    gDebug.bXZGridOn = true;
    
    // Start the action
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
    PxRigidDynamic* theBox = gPhysX.mPhysics->createRigidDynamic(PxTransform(PxVec3(0,wall_dh,0)));
    if (!theBox)
        ncc__error("Actor creation failed!");
    theBox->setRigidDynamicFlag(PxRigidDynamicFlag::eKINEMATIC, true);

    // Define sides
    PxBoxGeometry box_bottom(box_d/2,wall_dh/2,box_d/2);
    PxBoxGeometry box_side(wall_dh/2,box_h/2,box_d/2);
    PxMaterial* defmat=gPhysX.mDefaultMaterial;

    // Attach the sides, making front wall invisible
    theBox->createShape(box_bottom,*defmat); // the bottom
    theBox->createShape(box_side,*defmat,PxTransform(PxVec3(-box_d/2,box_h/2,0),PxQuat(0,PxVec3(0,1,0)))); // left wall
    theBox->createShape(box_side,*defmat,PxTransform(PxVec3( box_d/2,box_h/2,0),PxQuat(0,PxVec3(0,1,0)))); // right wall
    theBox->createShape(box_side,*defmat,PxTransform(PxVec3(0,box_h/2,-box_d/2),PxQuat(PxPi/2,PxVec3(0,1,0)))); // back wall
    PxShape* fwall =  theBox->createShape(box_side,*defmat,PxTransform(PxVec3(0,box_h/2, box_d/2),PxQuat(PxPi/2,PxVec3(0,1,0)))); // front wall
    fwall->setName("~fwall");

    // Name, color, and register the container
    theBox->setName("the_box");
    gColors.colorBucket.push_back(vector<GLubyte>(3));
    gColors.colorBucket.back()[0] = ncc::rgb::rRed[0];
    gColors.colorBucket.back()[1] = ncc::rgb::rRed[1];
    gColors.colorBucket.back()[2] = ncc::rgb::rRed[2];
    theBox->userData = &(gColors.colorBucket.back()[0]);
    gPhysX.mScene->addActor(*theBox);
    labscale::VIPs.container = theBox;

}

// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif