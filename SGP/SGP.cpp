/////////////////////////////////////////////////////////////////////////////////
// Header file for project SGP. Stands for Self-Gravitating Pile. This project is
// normally used to load a rubble pile from some source, impose some initial
// conditions, and then let it evolve.
//
// Author: Naor Movshovitz (nmovshov at google dot com)
/////////////////////////////////////////////////////////////////////////////////

#include "ARSS.h"
#include "SGP.h"

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
    if		(strcmp(buf,"make_sgp")==0)
        sgp::eExperimentType = sgp::eMAKE_SGP;
    else if (strcmp(buf,"load_sgp")==0)
        sgp::eExperimentType = sgp::eLOAD_SGP;
    else if (strcmp(buf,"test_scaling")==0)
        sgp::eExperimentType = sgp::eTEST_SCALING;
    else
        sgp::eExperimentType = sgp::eBAD_EXPERIMENT_TYPE;
    
    // The make_sgp subgroup includes parameters of the desired ellipsoid shape and rubble elements
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
    

    // Parameters of the grain size distribution (OBSOLETE REMOVE WHEN READY)
    /*OBSOLETE gsd group remove when ready*/
    ncc::GetStrPropertyFromINIFile("experiment","gsd_type","uniform",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    if		(strcmp(buf,"uniform")==0)
        sgp::gsd.type = sgp::gsd.eGSD_UNIFORM;
    else if	(strcmp(buf,"bimodal")==0)
        sgp::gsd.type = sgp::gsd.eGSD_BIMODAL;
    else
        sgp::gsd.type = sgp::gsd.eBAD_GSD_TYPE;
    /*OBSOLETE gsd group remove when ready*/

    ncc::GetStrPropertyFromINIFile("experiment","gsd_size1","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::gsd.size1 = atof(buf);
    ncc::GetStrPropertyFromINIFile("experiment","gsd_size2","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::gsd.size2 = atof(buf);

    sgp::gsd.totalNumber = ncc::GetIntPropertyFromINIFile("experiment","gsd_total_number",0,gRun.iniFile.c_str());
    sgp::gsd.numberRatio = ncc::GetIntPropertyFromINIFile("experiment","gsd_number_ratio",1,gRun.iniFile.c_str());

    ncc::GetStrPropertyFromINIFile("experiment","nucleus_radius", "0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::msgp.nucleusRadius = atof(buf);

    // The "scaled integration" tests are not really useful sgp experiments - more like debugging benchmarks
    ncc::GetStrPropertyFromINIFile("experiment","scl_test_radius", "1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::sclTest.radius = atof(buf);
    ncc::GetStrPropertyFromINIFile("experiment","scl_test_d_ini", "8",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
    sgp::sclTest.dInitial = atof(buf);

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
    sgp::hudMsgs.systemDiag1 = gHUD.hud.AddElement("",0.72,0.04);
    sgp::hudMsgs.systemDiag2 = gHUD.hud.AddElement("",0.72,0.08);
    sgp::hudMsgs.systemDiag3 = gHUD.hud.AddElement("",0.72,0.12);
    sgp::hudMsgs.systemDiag4 = gHUD.hud.AddElement("",0.72,0.16);
    sgp::hudMsgs.systemDiag5 = gHUD.hud.AddElement("",0.72,0.20);
}
void RefreshCustomHUDElements()
{
    switch (sgp::eExperimentType)
    {
    case sgp::eMAKE_SGP:
        sgp::RefreshMakeSGPHUD();
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
        DeadStop();
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
    sgp::GravitateSelf();
}
void ControlExperiment()
{
    switch (sgp::eExperimentType)
    {
    case sgp::eTEST_SCALING:
        sgp::ControlTestScalingExperiment();
        break;
    default:
        break;
    }
}
void RenderOtherStuff()
{

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

// SGP namespace functions
void sgp::CreateMakeSGPExperiment()
/*
 * Make fresh pile with geometry specified by free parameters. Let it self-gravitate
 * to a steady state.
*/
{
    // Put rubble elements in initial, loose positions
    sgp::MakeLooseRubblePile();

    // Move the camera to a good location
    PxReal looseExtent = 0;
    FindExtremers();
    if (gExp.VIPs.extremers.rightmost)
        looseExtent = 2*gExp.VIPs.extremers.rightmost->getGlobalPose().p.x;
    gCamera.pos.z = looseExtent + 2*sgp::msgp.gsd.sizeScale;
    gCamera.zBufFar = 2*looseExtent;
    gCamera.zBufNear = 0.5*sgp::msgp.gsd.sizeScale;

    // Start a log
    if (gRun.outputFrequency)
    {
        time_t now = time(NULL);
        ostringstream header;
        header << "# This is the run log of " << gRun.baseName << " from " << ctime(&now); // ctime includes a newline
        header << "# Experiment type: MAKE_SGP (" << sgp::eExperimentType << ")" << endl;
        header << "# Time step used = " << gSim.timeStep << " (cu)" << endl;
        header << "# Code units: 1 cu = [" << sgp::cunits.length << " m | " << sgp::cunits.mass << " kg | " << sgp::cunits.time << " s]" << endl;
        header << "# Scaled G = " << sgp::cunits.bigG << " (cu)" << endl;
        header << "# Columns are (values in code units):" << endl;
        header << "# [time]    [SGP a axis]    [SGP a/b axes ratio]    [SGP a/c axes ratio]    [system binding energy]" << endl;
        ofstream fbuf(gRun.outFile.c_str(),ios::trunc);
        if (!fbuf.is_open())
            ncc__error("Could not start a log. Experiment aborted.\a\n");
        fbuf << header.str() << endl;
    }

    // Start the action
    gSim.isRunning=true;
    gSim.bPause=true;
    gSim.codeTime = 0.0;
    gCUDA.cudaCapable=false; // TODO: remove when CUDA gravity is implemented

}
void sgp::CreateLoadSGPExperiment()
/* 
 * The strategy for loading is to choose between two branches. If a valid .repx file
 * is pointed to by the gRun parameter, the scene will be loaded from this file.
 * Then the physical properties of all dynamic actors will be updated based on current values.
 * Otherwise, we try to read initial conditions and shape information from ascii files.
*/
{
    if (gRun.loadSceneFromFile.empty())
    {
        //TODO: implement reading from ascii file
        ncc__error("Scene could not be loaded. Check that repx file exists. Experiment aborted.\a");
    } 
    else
    {
        if (!LoadSceneFromFile(gRun.loadSceneFromFile))
            ncc__error("Scene could not be loaded. Check that repx file exists. Experiment aborted.\a");
        PxU32 nbActors = gPhysX.mScene->getNbActors(gPhysX.roles.dynamics);
        gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,nbActors);
        while (nbActors--)
        {
            PxRigidDynamic* actor = gPhysX.cast[nbActors]->isRigidDynamic();
            actor->setAngularDamping(gPhysX.props.angularDamping);
            actor->setLinearDamping( gPhysX.props.linearDamping );
            const char* buf = actor->getName();
            if (buf && strcmp(buf,"size2")==0)
                ColorActor(actor,ncc::rgb::oDarkOrange);
            PxRigidDynamicFlags flags = actor->getRigidDynamicFlags();
            if (flags & PxRigidDynamicFlag::eKINEMATIC)
            {
                sgp::VIPs.nucleus = actor;
                ColorActor(actor,ncc::rgb::rDarkRed);
            }
        }
    }
}
void sgp::CreateTestScalingExperiment()
{
    // Make two balls, unit radius unit density
    PxReal halfD = sgp::sclTest.dInitial*sgp::sclTest.radius/2;
    sgp::VIPs.lBall = CreateRubbleGrain(PxVec3(-halfD,0,0),eSPHERE_GRAIN,sgp::sclTest.radius,*gPhysX.mDefaultMaterial,gExp.defGrainDensity);
    sgp::VIPs.rBall = CreateRubbleGrain(PxVec3(+halfD,0,0),eSPHERE_GRAIN,sgp::sclTest.radius,*gPhysX.mDefaultMaterial,gExp.defGrainDensity);

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
        header << "# Actors: two spheres with r = " << sgp::sclTest.radius << " (cu) and rho = " << gExp.defGrainDensity << " (cu)" << endl;
        header << "# Measured quantities: COM separation, COM relative velocity" << endl;
        header << "# Time step used = " << gSim.timeStep << " (cu)" << endl;
        header << "# Code units: 1 cu = [" << sgp::cunits.length << " m | " << sgp::cunits.mass << " kg | " << sgp::cunits.time << " s]" << endl;
        header << "# Scaled G = " << sgp::cunits.bigG << " (cu)" << endl;
        header << "# Columns are (values in code units):" << endl;
        header << "# [t]    [R]    [V]" << endl;
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
void sgp::GravitateSelf()
{
    if (gCUDA.cudaCapable)
        sgp::GravitateOnDevice();
    else
        sgp::GravitateOnHost();
}
void sgp::GravitateOnHost()
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
        PxVec3 F(forces[3*k+0],forces[3*k+1],forces[3*k+2]);
        actor->addForce(F);
    }

    // Clean up
    delete [] bodies;
    delete [] forces;

    return;
}
void sgp::GravitateOnDevice()
{

}
PxU32 sgp::MakeLooseRubblePile()
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

    // 2. Now calculate placement positions (this is the hard part)
    vector<PxVec3> positions(nbGrains);
    PxReal safeDL = 2.06*rGrain;
    PxReal r = safeDL, teta = 0, phi = 0;
    PxU32 lastLayerOccupancy = 0;
    for (PxU32 k=0; k<positions.size(); k++)
    {
        PxReal dr = 0, dteta = 0, dphi=0;
        // try advancing on a slice
        if (teta) {
            dphi = safeDL/(r*sin(teta));
            phi += dphi;
            if (phi > 2*PxPi)
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
        PxReal y = (a/b)*r*sin(teta)*sin(phi);
        PxReal z = (a/c)*r*cos(teta);
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

            aGrain->setName("rubble");
            RandOrientActor(aGrain);
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

    // Rubble element count
    sprintf(buf,"Element count = %u",gPhysX.mScene->getNbActors(gPhysX.roles.dynamics));
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag1,buf);
    if (gPhysX.mScene->getNbActors(gPhysX.roles.dynamics)==0)
        return;

    // Ellipsoid dimensions
    FindExtremers();
    PxVec3 rRight = gExp.VIPs.extremers.rightmost->getGlobalPose().transform(gExp.VIPs.extremers.rightmost->getCMassLocalPose()).p;
    PxVec3 rLeft  = gExp.VIPs.extremers.leftmost->getGlobalPose().transform(gExp.VIPs.extremers.leftmost->getCMassLocalPose()).p;
    PxVec3 rUp    = gExp.VIPs.extremers.upmost->getGlobalPose().transform(gExp.VIPs.extremers.upmost->getCMassLocalPose()).p;
    PxVec3 rDown  = gExp.VIPs.extremers.downmost->getGlobalPose().transform(gExp.VIPs.extremers.downmost->getCMassLocalPose()).p;
    PxVec3 rIn    = gExp.VIPs.extremers.inmost->getGlobalPose().transform(gExp.VIPs.extremers.inmost->getCMassLocalPose()).p;
    PxVec3 rOut   = gExp.VIPs.extremers.outmost->getGlobalPose().transform(gExp.VIPs.extremers.outmost->getCMassLocalPose()).p;
    PxReal a = rRight.x - rLeft.x;
    PxReal b = rOut.z - rIn.z;
    PxReal c = rUp.y - rDown.y;
    sprintf(buf,"Ellipsoid a axis = %g",a);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag2,buf);
    sprintf(buf,"Ellipsoid a/b axes ratio = %-8.2f",a/b);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag3,buf);
    sprintf(buf,"Ellipsoid a/c axes ratio = %-8.2f",a/c);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag4,buf);

    // Gravitational binding energy
    PxReal V = sgp::SystemPotentialEnergy();
    sprintf(buf,"PE_tot = %g",V);
    gHUD.hud.SetElement(sgp::hudMsgs.systemDiag5,buf);
}
void sgp::LogMakeSGPExperiment()
{
    // Shape information
    FindExtremers();
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

    // Format and write it to  log
    char buf[MAX_CHARS_PER_NAME];
    sprintf(buf,"%f    %g    %f    %f    %g",gSim.codeTime,a,a/b,a/c,V);
    ncc::logEntry(gRun.outFile.c_str(),buf);
}
physx::PxReal sgp::SystemPotentialEnergy()
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

            V -= bodies[4*j+0] * bodies[4*k+0] * invDist; // Multiply by 2*G later...
        }
    }

    // Clean up
    delete [] bodies;

    return 2*sgp::cunits.bigG*V;
}
void sgp::LogTestScalingExperiment()
{
    static vector<PxReal> t;
    static vector<PxReal> R;
    static vector<PxReal> V;

    if (gSim.isRunning) // While running collect measurements
    {
        PxReal d = sgp::VIPs.rBall->getGlobalPose().p.x - sgp::VIPs.lBall->getGlobalPose().p.x;
        PxReal v = sgp::VIPs.rBall->getLinearVelocity().x - sgp::VIPs.lBall->getLinearVelocity().x;
        t.push_back(gSim.codeTime);
        R.push_back(d);
        V.push_back(v);
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
            fbuf << setw(8) << t[k] << "    " << setw(8) << R[k] << "    " << setw(8) << V[k] << "\n";
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

// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif
