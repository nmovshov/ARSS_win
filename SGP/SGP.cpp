/////////////////////////////////////////////////////////////////////////////////
// Header file for project SGP. Stands for Self-Gravitating Pile. This project is
// normally used to load a rubble pile from some source, impose some initial
// conditions, and then let it evolve.
//
// Author: Naor Movshovitz (nmovshov at google dot com)
/////////////////////////////////////////////////////////////////////////////////

#include "ARSS.h"
#include "SGP.h"
#ifdef HAVE_CUDA_TK
#include "cugrav.h"
#endif

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
    InitCUDA(10);
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
    if		(strcmp(buf,"make_sgp")==0)
        sgp::eExperimentType = sgp::eMAKE_SGP;
    else if (strcmp(buf,"orbit_sgp")==0)
        sgp::eExperimentType = sgp::eORBIT_SGP;
    else if (strcmp(buf,"load_sgp")==0)
        sgp::eExperimentType = sgp::eLOAD_SGP;
    else if (strcmp(buf,"test_scaling")==0)
        sgp::eExperimentType = sgp::eTEST_SCALING;
    else
        sgp::eExperimentType = sgp::eBAD_EXPERIMENT_TYPE;
    
    // The make_sgp subgroup includes parameters of the desired ellipsoid shape and rubble elements
    ncc::GetStrPropertyFromINIFile("experiment:make_sgp","sgp_mass","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::msgp.mass = atof(buf);

    ncc::GetStrPropertyFromINIFile("experiment:make_sgp","ellipsoid_long_axis","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::msgp.ellipsoid.longAxis = atof(buf);
    ncc::GetStrPropertyFromINIFile("experiment:make_sgp","ellipsoid_ab_ratio","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::msgp.ellipsoid.abAxesRatio = atof(buf);
    ncc::GetStrPropertyFromINIFile("experiment:make_sgp","ellipsoid_ac_ratio","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::msgp.ellipsoid.acAxesRatio = atof(buf);

    ncc::GetStrPropertyFromINIFile("experiment:make_sgp","gsd_size_scale","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::msgp.gsd.sizeScale = atof(buf);
    ncc::GetStrPropertyFromINIFile("experiment:make_sgp","gsd_type","uniform",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    if (strcmp(buf,"uniform")==0)
        sgp::msgp.gsd.type = sgp::msgp.gsd.eGSD_UNIFORM;
    else if (strcmp(buf,"identical")==0)
        sgp::msgp.gsd.type = sgp::msgp.gsd.eGSD_IDENTICAL;
    else
        sgp::msgp.gsd.type = sgp::msgp.gsd.eBAD_GSD_TYPE;
    
    ncc::GetStrPropertyFromINIFile("experiment:make_sgp","grain_density","1000",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::msgp.grain.density = atof(buf);
    ncc::GetStrPropertyFromINIFile("experiment:make_sgp","grain_shape","convex",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    if (strcmp(buf,"convex")==0)
        sgp::msgp.grain.shape = eCONVEX_GRAIN;
    else if (strcmp(buf,"sphere")==0)
        sgp::msgp.grain.shape = eSPHERE_GRAIN;
    else
        sgp::msgp.grain.shape = eBAD_RUBBLE_GRAIN_TYPE;

    // The load_sgp subgroup includes parameters for rescaling a saved sgp
    ncc::GetStrPropertyFromINIFile("experiment:load_sgp","sgp_remass","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::lsgp.remass = atof(buf);

    // The orbit_sgp subgroup includes parameters for orbit selection/generation
    ncc::GetStrPropertyFromINIFile("experiment:orbit_sgp","sgp_mass","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::orbit.sgpMass = atof(buf);
    ncc::GetStrPropertyFromINIFile("experiment:orbit_sgp","big_M","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::orbit.bigM = atof(buf);
    sgp::orbit.nbOrbits = ncc::GetIntPropertyFromINIFile("experiment:orbit_sgp","nb_orbits",1,gRun.iniFile.c_str());
    
    // Parameters of the grain size distribution (OBSOLETE REMOVE WHEN READY)
    /*OBSOLETE gsd group remove when ready*/
    ncc::GetStrPropertyFromINIFile("experiment","gsd_type","uniform",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    if		(strcmp(buf,"uniform")==0)
        sgp::gsd.type = sgp::gsd.eGSD_UNIFORM;
    else if	(strcmp(buf,"bimodal")==0)
        sgp::gsd.type = sgp::gsd.eGSD_BIMODAL;
    else
        sgp::gsd.type = sgp::gsd.eBAD_GSD_TYPE;
    ncc::GetStrPropertyFromINIFile("experiment","gsd_size1","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::gsd.size1 = atof(buf);
    ncc::GetStrPropertyFromINIFile("experiment","gsd_size2","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::gsd.size2 = atof(buf);
    sgp::gsd.totalNumber = ncc::GetIntPropertyFromINIFile("experiment","gsd_total_number",0,gRun.iniFile.c_str());
    sgp::gsd.numberRatio = ncc::GetIntPropertyFromINIFile("experiment","gsd_number_ratio",1,gRun.iniFile.c_str());
    ncc::GetStrPropertyFromINIFile("experiment","nucleus_radius", "0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::msgp.nucleusRadius = atof(buf);
    /*OBSOLETE gsd group remove when ready*/

    // The "scaled integration" tests are not really useful sgp experiments - more like debugging benchmarks
    ncc::GetStrPropertyFromINIFile("experiment:test_scaling","radius", "1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::sclTest.radius = atof(buf);
    ncc::GetStrPropertyFromINIFile("experiment:test_scaling","d_ini", "8",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::sclTest.dInitial = atof(buf);
    ncc::GetStrPropertyFromINIFile("experiment:test_scaling","density", "1000",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::sclTest.density = atof(buf);
    sgp::sclTest.nballs = ncc::GetIntPropertyFromINIFile("experiment:test_scaling","nballs",2,gRun.iniFile.c_str());

    // The diagnostics group keeps track of options used for real-time and/or post pre/post-processing
    ncc::GetStrPropertyFromINIFile("diagnostics","color_code","off",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    if (strcmp(buf,"off")==0)
        sgp::diag.eColorCodeType = sgp::diag.eNO_CCODE;
    else if (strcmp(buf,"gradient")==0)
        sgp::diag.eColorCodeType = sgp::diag.eGRADIENT;
    else if (strcmp(buf,"two_layers")==0)
        sgp::diag.eColorCodeType = sgp::diag.eTWO_LAYER;
    else if (strcmp(buf,"three_layers")==0)
        sgp::diag.eColorCodeType = sgp::diag.eTHREE_LAYER;
    else if (strcmp(buf,"surface")==0)
        sgp::diag.eColorCodeType = sgp::diag.eSURFACE;
    else if (strcmp(buf,"by_name")==0)
        sgp::diag.eColorCodeType = sgp::diag.eNAME;
    else
        sgp::diag.eColorCodeType = sgp::diag.eBAD_CCODE;

    ncc::GetStrPropertyFromINIFile("diagnostics","surface_thickness","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::diag.nbSurfaceThickness = atof(buf);

    // Code units and scaling
    ncc::GetStrPropertyFromINIFile("units","cu_length","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::cunits.length = atof(buf);
    ncc::GetStrPropertyFromINIFile("units","cu_mass","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::cunits.mass = atof(buf);
    ncc::GetStrPropertyFromINIFile("units","cu_time","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::cunits.time = atof(buf);
    ncc::GetStrPropertyFromINIFile("units","cu_big_g","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::cunits.bigG = atof(buf);

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
    sgp::hudMsgs.systemDiag1 = gHUD.hud.AddElement("",0.68,0.04);
    sgp::hudMsgs.systemDiag2 = gHUD.hud.AddElement("",0.68,0.08);
    sgp::hudMsgs.systemDiag3 = gHUD.hud.AddElement("",0.68,0.12);
    sgp::hudMsgs.systemDiag4 = gHUD.hud.AddElement("",0.68,0.16);
    sgp::hudMsgs.systemDiag5 = gHUD.hud.AddElement("",0.68,0.20);
    sgp::hudMsgs.systemDiag6 = gHUD.hud.AddElement("",0.68,0.24);
}
void RefreshCustomHUDElements()
{
    if (gHUD.verbosity < 0) return;
    switch (sgp::eExperimentType)
    {
    case sgp::eORBIT_SGP:
        sgp::RefreshOrbitSGPHUD();
        break;
    case sgp::eMAKE_SGP:
        sgp::RefreshMakeSGPHUD();
        break;
    case sgp::eLOAD_SGP:
        sgp::RefreshLoadSGPHUD();
        break;
    case sgp::eBAD_EXPERIMENT_TYPE: // intentional fall through
    default:
        break;
    }
}
void FireAction()
{
    if (sgp::eExperimentType==sgp::eMAKE_SGP)
    {
        cout << "Stopping short" << endl;
        DeadStop();
    }
    else if (sgp::eExperimentType==sgp::eLOAD_SGP)
    {
        cout << "Stopping short" << endl;
        DeadStop();
    }
    else if (sgp::eExperimentType==sgp::eORBIT_SGP)
    {
        sgp::orbit.bTrackingCamera = !sgp::orbit.bTrackingCamera;
        cout << "Tracking camera = " << sgp::orbit.bTrackingCamera << endl;
    }
}
void LogExperiment()
{
    switch (sgp::eExperimentType)
    {
    case sgp::eMAKE_SGP:
        sgp::LogMakeSGPExperiment();
        break;
    case sgp::eTEST_SCALING:
        sgp::LogTestScalingExperiment();
        break;
    case sgp::eORBIT_SGP:
        sgp::LogOrbitSGPExperiment();
        break;
    case sgp::eBAD_EXPERIMENT_TYPE:
        ncc__warning("Unkown experiment type. Nothing logged.");
        break;
    default:
        break;
    }
}
void PrintDebug()
{
    gCamera;
}
void ApplyCustomInteractions()
{
    sgp::GravitateSelf(true);
    if (sgp::eExperimentType == sgp::eORBIT_SGP) {
        sgp::ApplyTidingForce();
    }
}
void ControlExperiment()
{
    if (gSim.bPause) return;
    switch (sgp::eExperimentType)
    {
    case sgp::eTEST_SCALING:
        sgp::ControlTestScalingExperiment();
        break;
    case sgp::eMAKE_SGP:
        sgp::ControlMakeSGPExperiment();
        break;
    case sgp::eORBIT_SGP:
        sgp::ControlOrbitSGPExperiment();
        break;
    case sgp::eLOAD_SGP:
        sgp::ControlLoadSGPExperiment();
    default:
        break;
    }
}
void RenderOtherStuff()
{
    if (sgp::eExperimentType == sgp::eORBIT_SGP)
    {
        DrawLine(gExp.IOMs.systemCM, sgp::VIPs.gravitator->getGlobalPose().p, 1.0, ncc::rgb::oOrangeRed);
    }
}
void CreateExperiment()
{
    switch (sgp::eExperimentType)
    {
    case sgp::eMAKE_SGP:
        sgp::CreateMakeSGPExperiment();
        break;
    case sgp::eLOAD_SGP:
        sgp::CreateLoadSGPExperiment();
        break;
    case sgp::eORBIT_SGP:
        sgp::CreateOrbitSGPExperiment();
        break;
    case sgp::eTEST_SCALING:
        sgp::CreateTestScalingExperiment();
        break;
    case sgp::eBAD_EXPERIMENT_TYPE: // intentional fall through
    default:
        ncc__error("Unknown experiment type. Experiment aborted.\a");
    }
    
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
    PxVec3 d = gExp.IOMs.systemCM;
    RelocateScene(d);
    gCamera.pos += d;
}
void DownArrowAction()
{

}
void LeftArrowAction()
{
    RelocateScene(PxVec3(-1,0,0));
}
void RightArrowAction()
{
    RelocateScene(PxVec3(1,0,0));
}

// SGP namespace functions
void sgp::CreateMakeSGPExperiment()
/*
 * Make fresh pile with geometry specified by free parameters. Let it self-gravitate
 * to a steady state.
*/
{
    // Put rubble elements in initial, loose positions
    sgp::msgp.gsd.nbTotal = sgp::MakeLooseRubblePile(sgp::msgp.mass);

    // Color-code rubble
    if (sgp::diag.eColorCodeType)
        sgp::ColorCodeRubblePile();
    else
        ncc__warning("Unknown color-coding scheme - using default color");

    // Move the camera to a good location
    PxReal looseExtent = 0;
    FindExtremers();
    if (gExp.VIPs.extremers.rightmost)
        looseExtent = 2*gExp.VIPs.extremers.rightmost->getGlobalPose().p.x;
    gCamera.pos.z = looseExtent + 2*sgp::msgp.gsd.sizeScale;
    gCamera.zBufFar = 2*looseExtent;
    gCamera.zBufNear = 0.5*sgp::msgp.gsd.sizeScale;
    gCamera.speed = gCamera.speed*sgp::msgp.gsd.sizeScale;

    // Start a log
    if (gRun.outputFrequency)
    {
        time_t now = time(NULL);
        ostringstream header;
        header << "# This is the run log of " << gRun.baseName << " from " << ctime(&now); // ctime includes a newline
        header << "# Experiment type: MAKE_SGP (" << sgp::eExperimentType << ")" << endl;
        header << "# Rubble elements: number = " << sgp::msgp.gsd.nbTotal << ";  rho = " << sgp::msgp.grain.density << endl;
        header << "# Material: eps = " << gPhysX.mDefaultMaterial->getRestitution() << ";  mu = " << gPhysX.mDefaultMaterial->getDynamicFriction() << endl;
        header << "# Time step used = " << gSim.timeStep << " (cu)" << endl;
        header << "# Code units: 1 cu = [" << sgp::cunits.length << " m | " << sgp::cunits.mass << " kg | " << sgp::cunits.time << " s]" << endl;
        header << "# Scaled G = " << sgp::cunits.bigG << " (cu)" << endl;
        header << "# Columns are (values in code units):" << endl;
        header << "# [time]    [long axis]    [a/b axes ratio]    [a/c axes ratio]    [U]    [K]" << endl;
        ofstream fbuf(gRun.outFile.c_str(),ios::trunc);
        if (!fbuf.is_open())
            ncc__error("Could not start a log. Experiment aborted.\a\n");
        fbuf << header.str() << endl;
    }

    // Start the action
    gSim.isRunning=true;
    gSim.bPause=false;
    gSim.codeTime = 0.0;

}
void sgp::CreateLoadSGPExperiment()
{
    // Load a previously saved SGP
    if (!sgp::LoadSGP(gRun.loadSceneFromFile))
        ncc__error("Scene could not be loaded from repx file; experiment aborted.\a");

    // Optionally rescale SGP
    if (sgp::lsgp.remass > 0) sgp::ReMassSGP(sgp::lsgp.remass);

    // Move the camera to a good location
    SpyOnSGP();

    // Start a log
    if (gRun.outputFrequency)
    {
        time_t now = time(NULL);
        ostringstream header;
        header << "# This is the run log of " << gRun.baseName << " from " << ctime(&now); // ctime includes a newline
        header << "# Experiment type: LOAD_SGP (" << sgp::eExperimentType << ")" << endl;
        ofstream fbuf(gRun.outFile.c_str(),ios::trunc);
        if (!fbuf.is_open())
            ncc__error("Could not start a log. Experiment aborted.\a\n");
        fbuf << header.str() << endl;
    }

    // Start the action
    gSim.isRunning=true;
    gSim.bPause=true;
    gSim.codeTime = 0.0;
    
}
void sgp::CreateTestScalingExperiment()
{
    // Make two balls, unit radius unit density
    PxReal halfD = sgp::sclTest.dInitial*sgp::sclTest.radius/2;
    sgp::VIPs.lBall = CreateRubbleGrain(PxVec3(-halfD,0,0),eSPHERE_GRAIN,sgp::sclTest.radius,*gPhysX.mDefaultMaterial,sgp::sclTest.density);
    sgp::VIPs.rBall = CreateRubbleGrain(PxVec3(+halfD,0,0),eSPHERE_GRAIN,sgp::sclTest.radius,*gPhysX.mDefaultMaterial,sgp::sclTest.density);

    // Move the camera to a good location
    gCamera.pos.z = 1.25*sgp::sclTest.dInitial*sgp::sclTest.radius;
    gCamera.zBufFar = 2*gCamera.pos.magnitude();
    gCamera.speed *= sgp::sclTest.radius;

    // Start a log
    if (gRun.outputFrequency)
    {
        time_t now = time(NULL);
        ostringstream header;
        header << "# This is the run log of " << gRun.baseName << " from " << ctime(&now); // ctime includes a newline
        header << "# Experiment type: TEST_SCALING (" << sgp::eExperimentType << ")" << endl;
        header << "# Active force: point mass gravity" << endl;
        header << "# Actors: two spheres with r = " << sgp::sclTest.radius << " (cu) and rho = " << sgp::sclTest.density << " (cu)" << endl;
        header << "# Measured quantities: COM separation, COM relative velocity" << endl;
        header << "# Time step used = " << gSim.timeStep << " (cu)" << endl;
        header << "# Code units: 1 cu = [" << sgp::cunits.length << " m | " << sgp::cunits.mass << " kg | " << sgp::cunits.time << " s]" << endl;
        header << "# Scaled G = " << sgp::cunits.bigG << " (cu)" << endl;
        header << "# Columns are (values in code units):" << endl;
        header << "# [t]    [R]    [V]    [U]    [K]" << endl;
        ofstream fbuf(gRun.outFile.c_str(),ios::trunc);
        if (!fbuf.is_open())
            ncc__error("Could not start a log. Experiment aborted.\a\n");
        fbuf << header.str() << endl;
    }

    // Start the action
    gDebug.bXYGridOn = true;
    gSim.isRunning=true;
    gSim.bPause=false;
    gSim.codeTime = 0.0;
    gCUDA.cudaCapable=false; // TODO: remove when CUDA gravity is implemented

    return;
}
void sgp::GravitateSelf(bool bIgnoreKinematics/*=false*/)
{
    if (gCUDA.cudaCapable && gCUDA.enableGPU)
        sgp::GravitateOnGPU();
    else
        sgp::GravitateOnHost(bIgnoreKinematics);
}
void sgp::GravitateOnHost(bool bIgnoreKinematics/*=false*/)
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
        if (bIgnoreKinematics && (actor->getRigidDynamicFlags() & PxRigidDynamicFlag::eKINEMATIC)) {
            bodies[4*k+0] = 0;
        }
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

            float s = sgp::cunits.bigG * bodies[4*j+0] * bodies[4*k+0] * invDistCube;

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
        if (actor->getRigidDynamicFlags() & PxRigidDynamicFlag::eKINEMATIC) continue;
        PxVec3 F(forces[3*k+0],forces[3*k+1],forces[3*k+2]);
        actor->addForce(F);
    }

    // Clean up
    delete [] bodies;
    delete [] forces;

    return;
}
void sgp::GravitateOnGPU()
{
    // Prepare
    PxU32 nbActors = gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,MAX_ACTORS_PER_SCENE);
    if (nbActors<2) return;
    AllocateCUDAGlobals(nbActors);

    // Linearize actor positions
    for (PxU32 k=0; k<nbActors; k++)
    {
        PxRigidDynamic* actor = gPhysX.cast[k]->isRigidDynamic();
        PxTransform pose = actor->getGlobalPose().transform(actor->getCMassLocalPose());
        PxVec3 pos = pose.p;
        g_hostPositionsArr[k].x = pos.x;
        g_hostPositionsArr[k].y = pos.y;
        g_hostPositionsArr[k].z = pos.z;
        g_hostPositionsArr[k].w = actor->getMass();
    }

    // And now, the N^2 double loop, on the device
    GravitateOnDevice(nbActors);

    // Add accumulated forces to actors (don't forget big G)
    for (PxU32 k=0; k<nbActors; k++)
    {
        PxRigidDynamic* actor = gPhysX.cast[k]->isRigidDynamic();
        if (actor->getRigidDynamicFlags() & PxRigidDynamicFlag::eKINEMATIC) continue;
        PxVec3 F;
        F.x = g_hostForcesArr[k].x;
        F.y = g_hostForcesArr[k].y;
        F.z = g_hostForcesArr[k].z;
        F = sgp::cunits.bigG*F;
        actor->addForce(F);
    }

    // Clean up
    ReleaseCUDA();

    return;
}
PxU32 sgp::MakeLooseRubblePile(PxReal mass /*=0*/)
/* This function creates a loose rubble pile by placing rubble elements, called
 * grains, inside the volume of an imaginary ellipsoid. There will be ample space
 * between neighboring grains to ensure no initial overlap regardless of (possibly
 * randomized) shape and orientation. The total number of grains will be
 * calculated from the average grain size and desired size of the final rubble
 * pile. The initial loose rubble pile is usually allowed to settle under self
 * gravity and the dimensions of the final aggregate will be only roughly similar
 * to the specified desired ellipsoid.
 * 
 * In the first version all rubble elements are drawn from a uniform size
 * distribution, but their shapes are optionally randomized.
*/
{
    // 1. Begin by calculating expected number of grains
    PxReal a = sgp::msgp.ellipsoid.longAxis;
    PxReal b = a/sgp::msgp.ellipsoid.abAxesRatio;
    PxReal c = a/sgp::msgp.ellipsoid.acAxesRatio;
    PxReal ellipsoidVolume = 4.0/3.0*PxPi*a*b*c;
    PxReal rGrain = sgp::msgp.gsd.sizeScale;
    PxReal grainVolume = (rGrain)*(rGrain)*(rGrain);
    PxU32 nbGrains = ellipsoidVolume/grainVolume;
    if (nbGrains > MAX_ACTORS_PER_SCENE)
        ncc__error("Too many grains!");
    if (nbGrains < 12)
        ncc__warning("Too few grains!");

    // 2. Now calculate placement positions (this is the hard part)
    vector<PxVec3> positions(nbGrains);
    PxReal safeDL = 2.06*rGrain*PxMax(a/b, a/c);
    PxReal r = safeDL, teta = 0, phi = 0;
    PxU32 lastLayerOccupancy = 0;
    for (PxU32 k=0; k<positions.size(); k++)
    {
        PxReal dr = 0, dteta = 0, dphi=0;
        // try advancing on a slice
        if (teta) {
            dphi = safeDL/(r*sin(teta));
            phi += dphi;
            if (phi > (2*PxPi - dphi))
                phi = 0;
        }
        // if not possible, try a stack
        if (phi == 0) {
            dteta = safeDL/r;
            teta += dteta;
            if (teta > PxPi)
                teta = 0;
        }
        // if not possible, try a radius
        if (teta == 0) {
            dr = safeDL;
            r += dr;
            lastLayerOccupancy = 0;
        }
        // make Cartesian coordinates 
        PxReal x = r*sin(teta)*cos(phi);
        PxReal y = (b/a)*r*sin(teta)*sin(phi);
        PxReal z = (c/a)*r*cos(teta);
        positions[k] = PxVec3(x,y,z);
        lastLayerOccupancy++;
    }
    // and shave the last layer, for symmetry (symmetry eh...)
    while (lastLayerOccupancy--)
    {
        positions.pop_back();
    }

    // 3. Finally, create and place the grains

    // Make a random convex mesh
    vector<PxVec3> verts = MakeRandomVertexList();
    PxConvexMesh* theMesh = MakePxMeshFromVertexList(verts);
    PxReal grainMass = mass/positions.size();

    for (PxU32 k=0; k<positions.size(); k++)
    {
        // Create a convex geometry for the next grain
        PxMeshScale meshScale;
        meshScale.scale = PxVec3(rGrain);
        PxMat33 M = meshScale.toMat33();
        PxConvexMeshGeometry meshGeometry;

        if (sgp::msgp.gsd.type == sgp::msgp.gsd.eGSD_IDENTICAL) // (use a single shape)
        {
            meshGeometry = PxConvexMeshGeometry(theMesh,meshScale);
        } 
        else // (or make a new shape each time)
        {
            vector<PxVec3> newVerts = MakeRandomVertexList();
            PxConvexMesh* aNewMesh = MakePxMeshFromVertexList(newVerts);
            meshGeometry = PxConvexMeshGeometry(aNewMesh,meshScale);
        }

        // Create a convex actor from this geometry and add to the scene
        PxRigidDynamic* aGrain = PxCreateDynamic(*gPhysX.mPhysics,PxTransform(positions[k]),meshGeometry,*gPhysX.mDefaultMaterial,sgp::msgp.grain.density);
        if (aGrain)
        {
            // (manual override of some dynamic features, we do not normally use this)
            if (gPhysX.props.sleepThreshold > -1) aGrain->setSleepThreshold(0.5*gPhysX.props.sleepThreshold*gPhysX.props.sleepThreshold);
            if (gPhysX.props.angularDamping > -1) aGrain->setAngularDamping(gPhysX.props.angularDamping);
            if (gPhysX.props.linearDamping  > -1) aGrain->setLinearDamping(gPhysX.props.linearDamping);

            if (mass) {
                PxRigidBodyExt::setMassAndUpdateInertia(*aGrain, grainMass);
            }

            aGrain->setName("rubble");
            RandOrientActor(aGrain);
            RandSpinActor(aGrain,1);
            gPhysX.mScene->addActor(*aGrain);
        }
    }

    return positions.size();
}
bool sgp::MakeNewSGP()
/* **************OBSOLETE********************
 * This function creates a rubble pile by placing grains in a volume of an
 * imaginary spheroid, and lets them fall in. The shape of the spheroid is
 * adjusted based on the number and size scale of rubble grains. The individual
 * grain sizes will be drawn from a size distribution (sgp::gsd). An optional
 * special grain called the nucleus is placed at the center.
 * 
 * In the "uniform" variant, all grains share a convex mesh, just scaled
 * differently.
 * TODO: In the non-uniform variant, grains are generated individually, with a
 * size scale.
 * TODO: Implement more size distributions (currently bimodal)
 * TODOL Implement more nucleus options (currently capsule)
*/
{
    // Make a random convex mesh to be referenced by all grains
    vector<PxVec3> verts = MakeRandomVertexList();
    PxConvexMesh* theMesh = MakePxMeshFromVertexList(verts);

    // Calculate placement positions
    vector<PxVec3> positions(sgp::gsd.totalNumber);
    PxReal safeDL = PxMax(sgp::gsd.size1,sgp::gsd.size2) * 2.06;

    PxReal r = sgp::msgp.nucleusRadius, teta = 0, phi = 0; // the center is reserved for an optional kinematic nucleus
    PxU32 lastLayerOccupancy = 0;
    for (PxU32 k=0; k<positions.size(); k++)
    {
        PxReal dr = 0, dteta = 0, dphi=0;
        // try advancing on a slice
        if (teta) {
            dphi = safeDL / (r * sin(teta));
            phi += dphi;
            if (phi > 2*PxPi)
                phi = 0;
        }
        // if not possible, try a stack
        if (phi == 0) {
            dteta = safeDL / r;
            teta += dteta;
            if (teta > PxPi)
                teta = 0;
        }
        // if not possible, try a radius
        if (teta == 0) {
            dr = safeDL;
            r += dr;
            lastLayerOccupancy = 0;
        }

        PxReal x = r * sin(teta) * cos(phi);
        PxReal y = (1/sgp::msgp.ellipsoid.abAxesRatio) * r * sin(teta) * sin(phi);
        PxReal z = (1/sgp::msgp.ellipsoid.acAxesRatio) * r * cos(teta);
        positions[k] = PxVec3(x,y,z);
        lastLayerOccupancy++;
    }

    // shave the last layer
    while (lastLayerOccupancy--)
    {
        positions.pop_back();
    }

    // place the nucleus
    if (sgp::msgp.nucleusRadius)
        sgp::VIPs.nucleus = CreateRubbleGrain(PxVec3(0),eCAPSULE_GRAIN,sgp::msgp.nucleusRadius/2,*gPhysX.mDefaultMaterial);
    if (sgp::VIPs.nucleus)
    {
        ColorActor(sgp::VIPs.nucleus,ncc::rgb::rDarkRed);
        sgp::VIPs.nucleus->setRigidDynamicFlag(PxRigidDynamicFlag::eKINEMATIC,true);
    }

    // Place actors
    for (PxU32 k=0; k<positions.size(); k++)
    {
        // select a size for the next grain
        PxReal grainScale = sgp::gsd.size1;
        //bool isSize2 = (k % (sgp::gsd.numberRatio)) == 0;
        bool isSize2 = (sgp::gsd.type==sgp::gsd.eGSD_BIMODAL && k < sgp::gsd.numberRatio);
        if (sgp::gsd.type==sgp::gsd.eGSD_BIMODAL && isSize2)
            grainScale = sgp::gsd.size2;
        
        // create a convex geometry for the next grain
        PxMeshScale meshScale;
        meshScale.scale = PxVec3(grainScale);
        PxMat33 M = meshScale.toMat33();
        PxConvexMeshGeometry meshGeometry(theMesh,meshScale);

        // create a convex actor from this geometry
        PxRigidDynamic* aGrain = PxCreateDynamic(*gPhysX.mPhysics,PxTransform(positions[k]),meshGeometry,*gPhysX.mDefaultMaterial,gExp.defGrainDensity);
        if (aGrain)
        {
            if (gPhysX.props.sleepThreshold > -1) aGrain->setSleepThreshold(0.5*gPhysX.props.sleepThreshold*gPhysX.props.sleepThreshold);
            if (gPhysX.props.angularDamping > -1) aGrain->setAngularDamping(gPhysX.props.angularDamping);
            if (gPhysX.props.linearDamping  > -1) aGrain->setLinearDamping(gPhysX.props.linearDamping);

            if (isSize2)
            {
                aGrain->setName("size2");
                ColorActor(aGrain,ncc::rgb::oDarkOrange);
            }
            else
            {
                aGrain->setName("size1");
            }

            RandOrientActor(aGrain);

            gPhysX.mScene->addActor(*aGrain);
        }
    }
    return true;
}
void sgp::RefreshMakeSGPHUD()
{
    char buf[MAX_CHARS_PER_NAME];
    if (gPhysX.mScene->getNbActors(gPhysX.roles.dynamics)==0)
        return;

    // SGP info: element count, total mass, mean density
    UpdateIntegralsOfMotion(true);
    PxReal rhoBulk = sgp::SGPBulkDensity();
    sprintf(buf,"Rubble elements (\"grains\") = %u",gExp.rubbleCount - 1);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag1,buf);
    sprintf(buf,"M_tot = %0.2g; <rho> = %4.0f",gExp.IOMs.systemMass,rhoBulk);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag2,buf);

    // Dynamic shape info
    FindExtremers(true);
    PxVec3 rRight = gExp.VIPs.extremers.rightmost->getGlobalPose().transform(gExp.VIPs.extremers.rightmost->getCMassLocalPose()).p;
    PxVec3 rLeft  = gExp.VIPs.extremers.leftmost->getGlobalPose().transform(gExp.VIPs.extremers.leftmost->getCMassLocalPose()).p;
    PxVec3 rUp    = gExp.VIPs.extremers.upmost->getGlobalPose().transform(gExp.VIPs.extremers.upmost->getCMassLocalPose()).p;
    PxVec3 rDown  = gExp.VIPs.extremers.downmost->getGlobalPose().transform(gExp.VIPs.extremers.downmost->getCMassLocalPose()).p;
    PxVec3 rIn    = gExp.VIPs.extremers.inmost->getGlobalPose().transform(gExp.VIPs.extremers.inmost->getCMassLocalPose()).p;
    PxVec3 rOut   = gExp.VIPs.extremers.outmost->getGlobalPose().transform(gExp.VIPs.extremers.outmost->getCMassLocalPose()).p;
    PxReal a = rRight.x - rLeft.x;
    PxReal b = rOut.z - rIn.z;
    PxReal c = rUp.y - rDown.y;
    sprintf(buf,"Ellipsoid long (\"a\") axis = %0.2g (cu)",a);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag3,buf);
    sprintf(buf,"Ellipsoid a/b axes ratio = %-8.2f",a/b);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag4,buf);
    sprintf(buf,"Ellipsoid a/c axes ratio = %-8.2f",a/c);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag5,buf);
    PxReal meanR = PxPow(a*b*c/8,1.0/3.0);
    sprintf(buf,"Mean radius = %0.2f",meanR);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag6,buf);

    //// Surface g
    //PxReal gee = sgp::cunits.bigG*gExp.IOMs.systemMass/(a*a);
    //sprintf(buf,"Surface acceleration = %0.2g (cu)",gee);
    //gHUD.hud.SetElement(sgp::hudMsgs.systemDiag6,buf);
    
}
void sgp::LogMakeSGPExperiment()
{
    // Shape information
    FindExtremers();
    if (gExp.rubbleCount == 0) return;
    PxVec3 rRight = gExp.VIPs.extremers.rightmost->getGlobalPose().transform(gExp.VIPs.extremers.rightmost->getCMassLocalPose()).p;
    PxVec3 rLeft  = gExp.VIPs.extremers.leftmost->getGlobalPose().transform(gExp.VIPs.extremers.leftmost->getCMassLocalPose()).p;
    PxVec3 rUp    = gExp.VIPs.extremers.upmost->getGlobalPose().transform(gExp.VIPs.extremers.upmost->getCMassLocalPose()).p;
    PxVec3 rDown  = gExp.VIPs.extremers.downmost->getGlobalPose().transform(gExp.VIPs.extremers.downmost->getCMassLocalPose()).p;
    PxVec3 rIn    = gExp.VIPs.extremers.inmost->getGlobalPose().transform(gExp.VIPs.extremers.inmost->getCMassLocalPose()).p;
    PxVec3 rOut   = gExp.VIPs.extremers.outmost->getGlobalPose().transform(gExp.VIPs.extremers.outmost->getCMassLocalPose()).p;
    PxReal a = rRight.x - rLeft.x;
    PxReal b = rOut.z - rIn.z;
    PxReal c = rUp.y - rDown.y;

    // Potential energy information
    PxReal V = sgp::SystemPotentialEnergy();

    // Kinetic energy information
    UpdateIntegralsOfMotion();
    PxReal K = gExp.IOMs.systemKE*gExp.IOMs.systemMass;

    // Format and write it to log
    char buf[MAX_CHARS_PER_NAME];
    sprintf(buf,"%8f    %8g    %8f    %8f    %12.3g    %12.3g",gSim.codeTime,a,a/b,a/c,V,K);
    ncc::logEntry(gRun.outFile.c_str(),buf);
}
PxReal sgp::SystemPotentialEnergy()
{
    // Prepare
    PxReal V = 0.0;
    PxU32 nbActors = gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,MAX_ACTORS_PER_SCENE);
    if (nbActors<2) return V;
    float *bodies = new float[4*nbActors];
    if (bodies==NULL)
        ncc__error("Error allocating memory for body/forces arrays.");

    // Linearize actor positions/masses
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
        for (PxU32 k=0; k<j; k++)
        {
            PxReal x = bodies[4*k+1] - bodies[4*j+1];
            PxReal y = bodies[4*k+2] - bodies[4*j+2];
            PxReal z = bodies[4*k+3] - bodies[4*j+3];

            float distSqr = x*x + y*y + z*z;
            float invDist = 1.0f/sqrtf(distSqr);

            V -= bodies[4*j+0] * bodies[4*k+0] * invDist; // Multiply by G later...
        }
    }

    // Clean up
    delete [] bodies;

    return sgp::cunits.bigG*V;
}
void sgp::LogTestScalingExperiment()
{
    static vector<PxReal> t;
    static vector<PxReal> R;
    static vector<PxReal> V;
    static vector<PxReal> U;
    static vector<PxReal> K;

    if (gSim.isRunning) // While running collect measurements
    {
        UpdateIntegralsOfMotion();
        PxReal d = sgp::VIPs.rBall->getGlobalPose().p.x - sgp::VIPs.lBall->getGlobalPose().p.x;
        PxReal v = sgp::VIPs.rBall->getLinearVelocity().x - sgp::VIPs.lBall->getLinearVelocity().x;
        PxReal u = sgp::SystemPotentialEnergy();
        PxReal k = gExp.IOMs.systemKE*gExp.IOMs.systemMass;
        t.push_back(gSim.codeTime);
        R.push_back(d);
        V.push_back(v);
        U.push_back(u);
        K.push_back(k);
    } 
    else // When done, save collected measurements
    {
        ofstream fbuf(gRun.outFile.c_str(),ios::app);
        if (!fbuf.is_open()) {
            ncc__warning("Could not open log. Nothing written!\a\n");
            return;
        }
        fbuf.setf(ios::fixed);
        for (unsigned int k=0; k<t.size(); k++)
            fbuf << setw(8) << t[k] << "    " << setw(8) << R[k] << "    " << setw(8) << V[k] << "    " << U[k] << "    " << K[k] << "\n";
        fbuf.close();
    }

}
void sgp::ControlTestScalingExperiment()
{
    // Stop just before contact
    PxReal d = sgp::VIPs.rBall->getGlobalPose().p.x - sgp::VIPs.lBall->getGlobalPose().p.x;
    PxReal skin =  0.1*gPhysX.mPhysics->getTolerancesScale().length;
    if (((d < 2.0*sgp::sclTest.radius + skin) && !gSim.targetTime) || (gSim.targetTime && (gSim.codeTime >= gSim.targetTime - gSim.timeStep)))
        gSim.isRunning = false;
}
void sgp::ControlMakeSGPExperiment()
{
    // Finish when sgp settles down
    PxReal T = gExp.IOMs.systemKE*gExp.IOMs.systemMass;
    PxReal U = sgp::SystemPotentialEnergy();
    PxReal X = PxAbs(T/U);
    PxReal GT = 1.0/PxSqrt(sgp::cunits.bigG*sgp::msgp.grain.density);
    if (X < 1e-4 || (gSim.codeTime > 4*GT))
    {
        SaveSceneToRepXDump();
        gSim.isRunning = false;
    }
}
void sgp::ColorCodeRubblePile()
{
    // Short-circuit
    if (sgp::diag.eColorCodeType == sgp::diag.eNO_CCODE)
        return;

    // Choose paint for rust
    const GLubyte *color = ncc::rgb::oCoral;
    gColors.colorBucket.push_back(vector<GLubyte>(3));
    gColors.colorBucket.back()[0]=color[0];
    gColors.colorBucket.back()[1]=color[1];
    gColors.colorBucket.back()[2]=color[2];
    size_t rCIndex = gColors.colorBucket.size() - 1; // core color index

    // Choose paint for snow
    color = ncc::rgb::wBeige;
    gColors.colorBucket.push_back(vector<GLubyte>(3));
    gColors.colorBucket.back()[0]=color[0];
    gColors.colorBucket.back()[1]=color[1];
    gColors.colorBucket.back()[2]=color[2];
    size_t sCIndex = gColors.colorBucket.size() - 1; // mantle color index

    // Detect shape information
    UpdateIntegralsOfMotion();
    if (gExp.rubbleCount == 0) return;
    FindExtremers();
    PxVec3 X0 = gExp.IOMs.systemCM;
    PxVec3 rRight = gExp.VIPs.extremers.rightmost->getGlobalPose().transform(gExp.VIPs.extremers.rightmost->getCMassLocalPose()).p;
    PxVec3 rLeft  = gExp.VIPs.extremers.leftmost->getGlobalPose().transform(gExp.VIPs.extremers.leftmost->getCMassLocalPose()).p;
    PxVec3 rUp    = gExp.VIPs.extremers.upmost->getGlobalPose().transform(gExp.VIPs.extremers.upmost->getCMassLocalPose()).p;
    PxVec3 rDown  = gExp.VIPs.extremers.downmost->getGlobalPose().transform(gExp.VIPs.extremers.downmost->getCMassLocalPose()).p;
    PxVec3 rIn    = gExp.VIPs.extremers.inmost->getGlobalPose().transform(gExp.VIPs.extremers.inmost->getCMassLocalPose()).p;
    PxVec3 rOut   = gExp.VIPs.extremers.outmost->getGlobalPose().transform(gExp.VIPs.extremers.outmost->getCMassLocalPose()).p;
    PxReal a = (rRight.x - rLeft.x)/2;
    PxReal b = (rOut.z - rIn.z)/2;
    PxReal c = (rUp.y - rDown.y)/2;
    PxU32  nbRubble = gPhysX.mScene->getNbActors(gPhysX.roles.dynamics);

    // Paint actors by scheme (avoiding switch/case this time)
    if (sgp::diag.eColorCodeType == sgp::diag.eTWO_LAYER) // Half snow, half default
    {
        // Loop over rubble and color by position
        PxReal ca = a/2, cb = b/2, cc = c/2; // core semi-axes
        PxU32 nbActors = gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,MAX_ACTORS_PER_SCENE);
        while (nbActors--)
        {
            PxRigidDynamic* actor = gPhysX.cast[nbActors]->isRigidDynamic();
            if (actor->getRigidDynamicFlags() & PxRigidDynamicFlag::eKINEMATIC) continue;
            PxVec3 pos = actor->getGlobalPose().transform(actor->getCMassLocalPose()).p - X0;
            bool incore = ((pos.x*pos.x/ca/ca + pos.y*pos.y/cc/cc + pos.z*pos.z/cb/cb) < 1);
            if (incore)
            {
                actor->userData=&(gColors.colorBucket[sCIndex][0]); // yes I know I know :/
                actor->setName("rubble-core");
            }
            else
            {
                actor->setName("rubble-mantle");
            }
        }
    }
    if (sgp::diag.eColorCodeType == sgp::diag.eTHREE_LAYER) // Half cotton, half linen, half wool
    {
        // Loop over rubble and color by position
        PxReal ca = a/3, cb = b/3, cc = c/3; // core semi-axes
        PxReal ma = 2*a/3, mb = 2*b/3, mc = 2*c/3; // mantle semi-axes
        PxU32 nbActors = gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,MAX_ACTORS_PER_SCENE);
        while (nbActors--)
        {
            PxRigidDynamic* actor = gPhysX.cast[nbActors]->isRigidDynamic();
            if (actor->getRigidDynamicFlags() & PxRigidDynamicFlag::eKINEMATIC) continue;
            PxVec3 pos = actor->getGlobalPose().transform(actor->getCMassLocalPose()).p - X0;
            bool inmantle = ((pos.x*pos.x/ma/ma + pos.y*pos.y/mc/mc + pos.z*pos.z/mb/mb) < 1);
            bool incore = ((pos.x*pos.x/ca/ca + pos.y*pos.y/cc/cc + pos.z*pos.z/cb/cb) < 1);
            if (incore)
            {
                actor->userData=&(gColors.colorBucket[rCIndex][0]); // yes I know I know :/
                actor->setName("rubble-core");
            }
            else if (inmantle)
            {
                actor->userData=&(gColors.colorBucket[sCIndex][0]); // yes I know I know :/
                actor->setName("rubble-mantle");
            }
            else
            {
                actor->setName("rubble-crust");
            }
        }
    }
    else if (sgp::diag.eColorCodeType == sgp::diag.eSURFACE) // surface coat
    {
        // Estimate surface layer thickness (normalized)
        PxReal s = 1.0 - sgp::diag.nbSurfaceThickness*PxPow(nbRubble,-1.0/3.0);

        // Loop over rubble and color by position
        PxU32 nbActors = gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,MAX_ACTORS_PER_SCENE);
        while (nbActors--)
        {
            PxRigidDynamic* actor = gPhysX.cast[nbActors]->isRigidDynamic();
            if (actor->getRigidDynamicFlags() & PxRigidDynamicFlag::eKINEMATIC) continue;
            PxVec3 pos = actor->getGlobalPose().transform(actor->getCMassLocalPose()).p - X0;
            bool incore = ((pos.x*pos.x/a/a + pos.y*pos.y/c/c + pos.z*pos.z/b/b) < s);
            if (incore)
            {
                actor->userData=&(gColors.colorBucket[sCIndex][0]); // yes I know I know :/
                actor->setName("rubble-interior");
            }
            else
            {
                actor->setName("rubble-surface");
            }
        }
    }
    else if (sgp::diag.eColorCodeType == sgp::diag.eNAME)
    {
        // Loop over rubble and color by name
        PxU32 nbActors = gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,MAX_ACTORS_PER_SCENE);
        while (nbActors--)
        {
            PxRigidDynamic* actor = gPhysX.cast[nbActors]->isRigidDynamic();
            if (actor->getRigidDynamicFlags() & PxRigidDynamicFlag::eKINEMATIC) continue;
            const char* buf = actor->getName();
            if (strcmp(buf,"rubble-interior")==0)
                actor->userData=&(gColors.colorBucket[sCIndex][0]); // yes I know I know :/
            if (strcmp(buf,"rubble-mantle")==0)
                actor->userData=&(gColors.colorBucket[sCIndex][0]); // yes I know I know :/
            if (strcmp(buf,"rubble-core")==0)
                actor->userData=&(gColors.colorBucket[rCIndex][0]); // yes I know I know :/
        }
    }   
}
bool sgp::LoadSGP(string filename)
/*Load a scene and do some processing specific to SGPs*/
{
    // Load a previously serialized scene
    if (!LoadSceneFromFile(filename))
        return false;

    // Some dynamics parameters must be overridden actor-by-actor :/
    PxU32 nbActors = gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,MAX_ACTORS_PER_SCENE);
    while (nbActors--)
    {
        PxRigidDynamic* actor = gPhysX.cast[nbActors]->isRigidDynamic();
        if (gPhysX.props.angularDamping > 0)
            actor->setAngularDamping(gPhysX.props.angularDamping);
        if (gPhysX.props.linearDamping > 0)
            actor->setLinearDamping( gPhysX.props.linearDamping );
    }

    // Color-code rubble
    if (sgp::diag.eColorCodeType)
        sgp::ColorCodeRubblePile();
    else
        ncc__warning("Unknown color-coding scheme - using default color");

    // And update the IOMs just for good measure (and just for non-kinematic)
    UpdateIntegralsOfMotion(true);

    return true;
}
void sgp::CreateOrbitSGPExperiment()
{
    // Load a previously saved SGP
    if (!sgp::LoadSGP(gRun.loadSceneFromFile)) {
        //ncc__error("Could not load SGP from file; experiment aborted.\a");
        ncc__warning("Could not load SGP from file; experiment aborted.");
        CreateRubbleGrain(PxVec3(1,0,0),eSPHERE_GRAIN,0.5,*gPhysX.mDefaultMaterial,100);
        CreateRubbleGrain(PxVec3(-1,0,0),eSPHERE_GRAIN,0.5,*gPhysX.mDefaultMaterial,100);
    }
    DeadStop(); // stomp any residual velocities
    RecenterScene(); // put center-of-mass at origin
    if (sgp::orbit.sgpMass > 0) // set sgp total mass
        sgp::ReMassSGP(sgp::orbit.sgpMass);

    // Load orbit from run_base_name.orb
    sgp::orbit.orbFile = gRun.workingDirectory + "/" + gRun.baseName + ".orb";
    nr3::MatDoub raw;
    if (ncc::load(sgp::orbit.orbFile.c_str(), &raw, 5))
    {
        // Distribute to vectors (kind of pointless really)
        sgp::orbit.tvec.resize(raw.nrows());
        sgp::orbit.xvec.resize(raw.nrows());
        sgp::orbit.yvec.resize(raw.nrows());
        for (int k=0; k<raw.nrows(); k++) {
            sgp::orbit.tvec[k] = raw[k][0];
            sgp::orbit.xvec[k] = raw[k][1];
            sgp::orbit.yvec[k] = raw[k][2];
        }

        // Estimate simulation length
        double t0 = sgp::orbit.tvec[0];
        double tf = sgp::orbit.tvec[raw.nrows() - 1];
        cout << "Requested orbit integration time " << (tf-t0) << " (cu) in " << (int)((tf-t0)/gSim.timeStep) << " time steps." << endl;
        if ((tf-t0)/gSim.timeStep > 2e4)
            ncc__warning("Long orbit requested -- this might take a while\a");
        sgp::orbit.tStart = t0;
        sgp::orbit.tEnd = tf;

        // Save the initial location and orbit parameters
        sgp::orbit.X0.x = sgp::orbit.xvec[0];
        sgp::orbit.X0.y = sgp::orbit.yvec[0];
        double q2 = 1e200;
        for (int k=0; k<sgp::orbit.xvec.size(); k++)
        {
            double x = sgp::orbit.xvec[k];
            double y = sgp::orbit.yvec[k];
            double r2 = x*x + y*y;
            if (r2 < q2) q2 = r2;
        }
        sgp::orbit.periapse = sqrt(q2);
    } 
    else
    {
        ncc__error("Could not read orbit data, experiment aborted.")
    }
    
    // Create the gravitator (actor representing the primary) located at -X0
    PxRigidDynamic* center = CreateRubbleGrain(-sgp::orbit.X0,eSPHERE_GRAIN,1,*gPhysX.mDefaultMaterial);
    center->setRigidDynamicFlag(PxRigidDynamicFlag::eKINEMATIC, true);
    center->setMass(sgp::orbit.bigM); // we don't care about the inertia for this kinematic actor
    ColorActor(center, ncc::rgb::oOrange);
    sgp::VIPs.gravitator = center;

    // Move the camera to a good location
    SpyOnSGP();
    sgp::orbit.bTrackingCamera = false;

    // Start a log
    if (gRun.outputFrequency)
    {
        time_t now = time(NULL);
        ostringstream header;
        header << "# This is the run log of " << gRun.baseName << " from " << ctime(&now); // ctime includes a newline
        header << "# Experiment type: ORBIT_SGP (" << sgp::eExperimentType << ")" << endl;
        header << "# Central mass = " << sgp::orbit.bigM << " (cu)" << endl;
        header << "# Periapse = " << sgp::orbit.periapse << " (cu)" << endl;
        header << "# Rubble elements = " << gExp.rubbleCount << ";  Bulk density ~ " << SGPBulkDensity() << " (cu)" << endl;
        header << "# Time step used = " << gSim.timeStep << " (cu)" << endl;
        header << "# Code units: 1 cu = [" << sgp::cunits.length << " m | " << sgp::cunits.mass << " kg | " << sgp::cunits.time << " s]" << endl;
        header << "# Scaled G = " << sgp::cunits.bigG << " (cu)" << endl;
        header << "# Columns are (values in code units):" << endl;
        header << "# [time]    [CoM X]    [CoM Y]    [CoM Z]    [SGP bulk density]" << endl;
        ofstream fbuf(gRun.outFile.c_str(),ios::trunc);
        if (!fbuf.is_open())
            ncc__error("Could not start a log. Experiment aborted.\a\n");
        fbuf << header.str() << endl;
    }

    // Start the action
    gSim.isRunning=true;
    gSim.bPause=true;
    gSim.codeTime = 0.0;

}
void sgp::ControlOrbitSGPExperiment()
{
    // Remove CM drift
    PxVec3 d = sgp::FindSGPCenterOfMass();
    if (d.magnitude() > 0.01*sgp::orbit.periapse)
    {
        ncc__warning("CM drift detected (and corrected)");
    	RelocateScene(-d);
        PxVec3 V = gExp.IOMs.systemLM; // linear momentum per unit mass, aka CM velocity
        PxU32 nbActors = gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,MAX_ACTORS_PER_SCENE);
        while (nbActors--)
        {
            PxRigidDynamic* actor = gPhysX.cast[nbActors]->isRigidDynamic();
            if (actor)
            {
                if (actor->getRigidDynamicFlags() & PxRigidDynamicFlag::eKINEMATIC) continue;
                PxVec3 v = actor->getLinearVelocity();
                actor->setLinearVelocity(v - V);
            }
        }
    }

    // Put gravitator where it belongs
    static int orbStep = 0;
    PxReal orbTime = sgp::orbit.tStart + gSim.codeTime;
    while (orbStep < sgp::orbit.tvec.size() && sgp::orbit.tvec[orbStep] < orbTime) {
        orbStep++;
    }
    PxVec3 cforce(sgp::orbit.xvec[orbStep - 1], sgp::orbit.yvec[orbStep - 1], 0);
    PxTransform pose = sgp::VIPs.gravitator->getGlobalPose();
    pose.p = -cforce;
    sgp::VIPs.gravitator->setGlobalPose(pose);

    // Move the camera to a good location
    if (sgp::orbit.bTrackingCamera)
        SpyOnSGP();

    // Save and quit when done
    if (gSim.codeTime > (sgp::orbit.tEnd - sgp::orbit.tStart)) {
        SaveSceneToRepXDump();
        gSim.isRunning = false;
    }
}
void sgp::SpyOnSGP(PxReal f/*=1.0*/, bool bZoomOutOnly/*=true*/)
{
    // Move the camera to a good location
    static PxReal currentExtent = 0.0f;
    UpdateIntegralsOfMotion(true);
    FindExtremers(true);
    if (gExp.VIPs.extremers.rightmost)
    {
        PxVec3 X0 = gExp.IOMs.systemCM;
        PxVec3 rRight = gExp.VIPs.extremers.rightmost->getGlobalPose().transform(gExp.VIPs.extremers.rightmost->getCMassLocalPose()).p;
        PxVec3 rLeft  = gExp.VIPs.extremers.leftmost->getGlobalPose().transform(gExp.VIPs.extremers.leftmost->getCMassLocalPose()).p;
        PxVec3 rUp    = gExp.VIPs.extremers.upmost->getGlobalPose().transform(gExp.VIPs.extremers.rightmost->getCMassLocalPose()).p;
        PxVec3 rDown  = gExp.VIPs.extremers.downmost->getGlobalPose().transform(gExp.VIPs.extremers.leftmost->getCMassLocalPose()).p;
        PxReal pileExtent = PxMax((rUp.y - rDown.y),(rRight.x - rLeft.x));
        PxU32  nbRubble = gPhysX.mScene->getNbActors(gPhysX.roles.dynamics);
        PxReal gsize = PxPow(nbRubble,-1.0/3.0);

        gCamera.pos.x = X0.x;
        gCamera.pos.y = X0.y;
        if ((pileExtent > currentExtent) || (!bZoomOutOnly)) {
            currentExtent = pileExtent;
            gCamera.pos.z = X0.z + 1.4*pileExtent + 2*gsize;
        }
        gCamera.zBufFar = 20*pileExtent;
        gCamera.zBufNear = 0.5*gsize;
        gCamera.pos *= f;
    }
}
PxVec3 sgp::FindSGPCenterOfMass()
{
    // Right now it's just calling UpdateIntegralsOfMotion(true) so kind of redundant
    UpdateIntegralsOfMotion(true);
    return gExp.IOMs.systemCM;
}
PxReal sgp::SGPBulkDensity(bool bRoughGuess/*=true*/)
{
    PxReal rhoBulk = 1000;
    UpdateIntegralsOfMotion(true);
    if (gExp.rubbleCount < 3) return rhoBulk;
    FindExtremers(true);
    if (gExp.VIPs.extremers.rightmost)
    {
        PxVec3 rRight = gExp.VIPs.extremers.rightmost->getGlobalPose().transform(gExp.VIPs.extremers.rightmost->getCMassLocalPose()).p;
        PxVec3 rLeft  = gExp.VIPs.extremers.leftmost->getGlobalPose().transform(gExp.VIPs.extremers.leftmost->getCMassLocalPose()).p;
        PxVec3 rUp    = gExp.VIPs.extremers.upmost->getGlobalPose().transform(gExp.VIPs.extremers.upmost->getCMassLocalPose()).p;
        PxVec3 rDown  = gExp.VIPs.extremers.downmost->getGlobalPose().transform(gExp.VIPs.extremers.downmost->getCMassLocalPose()).p;
        PxVec3 rIn    = gExp.VIPs.extremers.inmost->getGlobalPose().transform(gExp.VIPs.extremers.inmost->getCMassLocalPose()).p;
        PxVec3 rOut   = gExp.VIPs.extremers.outmost->getGlobalPose().transform(gExp.VIPs.extremers.outmost->getCMassLocalPose()).p;
        PxReal a = (rRight.x - rLeft.x)/2;
        PxReal b = (rOut.z - rIn.z)/2;
        PxReal c = (rUp.y - rDown.y)/2;
        PxReal V = 4.0*PxPi/3.0*a*b*c;
        
        rhoBulk = gExp.IOMs.systemMass/V;
    }
    return rhoBulk;
}
void sgp::LogOrbitSGPExperiment()
{
    // Orbit information
    PxReal t = gSim.codeTime;
    PxVec3 X = sgp::FindSGPCenterOfMass() - sgp::VIPs.gravitator->getGlobalPose().transform(sgp::VIPs.gravitator->getCMassLocalPose()).p;
    
    // Pile information (put cluster count here in the future)
    PxReal rho = sgp::SGPBulkDensity();
    rho = 0; //DEBUG

    // Format and write it to log (yes each time, unbuffered)
    char buf[MAX_CHARS_PER_NAME];
    sprintf(buf,"%8f    %12.3g    %12.3g    %12.3g    %12.3g", t, X.x, X.y, X.z, rho);
    ncc::logEntry(gRun.outFile.c_str(),buf);
}
void sgp::RefreshOrbitSGPHUD()
{
    char buf[MAX_CHARS_PER_NAME];
    if (gPhysX.mScene->getNbActors(gPhysX.roles.dynamics)==0)
        return;

    // SGP info: element count, total mass, mean density
    UpdateIntegralsOfMotion(true);
    PxReal rhoBulk = sgp::SGPBulkDensity();
    sprintf(buf,"Rubble elements (\"grains\") = %u",gExp.rubbleCount - 1);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag1,buf);
    sprintf(buf,"M_tot = %0.2g; <rho> = %4.0f",gExp.IOMs.systemMass,rhoBulk);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag2,buf);

    // Dynamic shape info
    FindExtremers(true);
    PxVec3 rRight = gExp.VIPs.extremers.rightmost->getGlobalPose().transform(gExp.VIPs.extremers.rightmost->getCMassLocalPose()).p;
    PxVec3 rLeft  = gExp.VIPs.extremers.leftmost->getGlobalPose().transform(gExp.VIPs.extremers.leftmost->getCMassLocalPose()).p;
    PxVec3 rUp    = gExp.VIPs.extremers.upmost->getGlobalPose().transform(gExp.VIPs.extremers.upmost->getCMassLocalPose()).p;
    PxVec3 rDown  = gExp.VIPs.extremers.downmost->getGlobalPose().transform(gExp.VIPs.extremers.downmost->getCMassLocalPose()).p;
    PxVec3 rIn    = gExp.VIPs.extremers.inmost->getGlobalPose().transform(gExp.VIPs.extremers.inmost->getCMassLocalPose()).p;
    PxVec3 rOut   = gExp.VIPs.extremers.outmost->getGlobalPose().transform(gExp.VIPs.extremers.outmost->getCMassLocalPose()).p;
    PxReal a = (rRight.x - rLeft.x)/2;
    PxReal b = (rOut.z - rIn.z)/2;
    PxReal c = (rUp.y - rDown.y)/2;
    PxReal meanR = PxPow(a*b*c,1.0/3.0);
    sprintf(buf,"Mean radius = %0.2f",meanR);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag3,buf);

    // Orbit info
    PxVec3 X = sgp::FindSGPCenterOfMass() - sgp::VIPs.gravitator->getGlobalPose().transform(sgp::VIPs.gravitator->getCMassLocalPose()).p;
    PxReal r = X.magnitude();
    PxReal roche = 1.51*PxPow(sgp::orbit.bigM/rhoBulk,1.0/3.0);
    sprintf(buf,"Distance = %0.3g (cu)",r);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag4,buf);
    sprintf(buf,"Distance = %0.2f x q = %0.2f x roche",r/sgp::orbit.periapse, r/roche);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag5,buf);
}
bool sgp::GenerateOrbit()
{
    // placeholder to maybe generate orbits in-house
    bool success = false;
    return success;
}
void sgp::ReMassSGP(PxReal newMass)
{
    PxU32 nbGrains = gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,MAX_ACTORS_PER_SCENE);
    PxReal newMPerGrain = newMass/nbGrains;
    while (nbGrains--)
    {
        PxRigidDynamic* actor = gPhysX.cast[nbGrains]->isRigidDynamic();
        PxRigidBodyExt::setMassAndUpdateInertia(*actor,newMPerGrain);
    }
}
void sgp::ApplyTidingForce()
{
    PxReal G = sgp::cunits.bigG;
    PxReal M = sgp::orbit.bigM;
    PxVec3 D = sgp::VIPs.gravitator->getGlobalPose().p; // vector FROM origin TO planet
    PxReal d = D.magnitude();
    PxVec3 pseudoA = -G*M/(d*d*d)*D; // FROM planet TO origin
    PxU32 nbActors = gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,MAX_ACTORS_PER_SCENE);
    while (nbActors--)
    {
        PxRigidDynamic* actor = gPhysX.cast[nbActors]->isRigidDynamic();
        if (actor)
        {
            if (actor->getRigidDynamicFlags() & PxRigidDynamicFlag::eKINEMATIC) continue;
            PxReal m = actor->getMass();
            PxVec3 pseudoF = m*pseudoA;
            PxVec3 X = D - actor->getGlobalPose().p; // FROM actor TO planet
            PxReal x = X.magnitude();
            PxVec3 gravF = G*M*m/(x*x*x)*X;
            PxVec3 netF = pseudoF + gravF;
            actor->addForce(netF);
        }
    }
}
void sgp::ControlLoadSGPExperiment()
{
    
}
void sgp::RefreshLoadSGPHUD()
{
    char buf[MAX_CHARS_PER_NAME];
    if (gPhysX.mScene->getNbActors(gPhysX.roles.dynamics)==0)
        return;

    // SGP info: element count, total mass, mean density
    if (gHUD.verbosity > 0)
    {
	    UpdateIntegralsOfMotion(true);
	    sprintf(buf,"Rubble elements (\"grains\") = %u",gExp.rubbleCount - 1);
	    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag1,buf);
	    sprintf(buf,"M_tot = %0.2g",gExp.IOMs.systemMass);
        if (gHUD.verbosity > 1)
        {
            PxReal rhoBulk = sgp::SGPBulkDensity();
            sprintf(buf,"M_tot = %0.2g; <rho> = %4.0f",gExp.IOMs.systemMass,rhoBulk);
        }
	    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag2,buf);
    }

    // Dynamic shape info
    if (gHUD.verbosity > 1)
    {
	    FindExtremers(true);
	    PxVec3 rRight = gExp.VIPs.extremers.rightmost->getGlobalPose().transform(gExp.VIPs.extremers.rightmost->getCMassLocalPose()).p;
	    PxVec3 rLeft  = gExp.VIPs.extremers.leftmost->getGlobalPose().transform(gExp.VIPs.extremers.leftmost->getCMassLocalPose()).p;
	    PxVec3 rUp    = gExp.VIPs.extremers.upmost->getGlobalPose().transform(gExp.VIPs.extremers.upmost->getCMassLocalPose()).p;
	    PxVec3 rDown  = gExp.VIPs.extremers.downmost->getGlobalPose().transform(gExp.VIPs.extremers.downmost->getCMassLocalPose()).p;
	    PxVec3 rIn    = gExp.VIPs.extremers.inmost->getGlobalPose().transform(gExp.VIPs.extremers.inmost->getCMassLocalPose()).p;
	    PxVec3 rOut   = gExp.VIPs.extremers.outmost->getGlobalPose().transform(gExp.VIPs.extremers.outmost->getCMassLocalPose()).p;
	    PxReal a = rRight.x - rLeft.x;
	    PxReal b = rOut.z - rIn.z;
	    PxReal c = rUp.y - rDown.y;
	    sprintf(buf,"Ellipsoid long (\"a\") axis = %0.2g (cu)",a);
	    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag3,buf);
	    sprintf(buf,"Ellipsoid a/b axes ratio = %-8.2f",a/b);
	    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag4,buf);
	    sprintf(buf,"Ellipsoid a/c axes ratio = %-8.2f",a/c);
	    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag5,buf);
	    PxReal meanR = PxPow(a*b*c/8,1.0/3.0);
	    sprintf(buf,"Mean radius = %0.2f",meanR);
	    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag6,buf);
    }
}


// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif
