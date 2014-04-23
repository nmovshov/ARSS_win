///////////////////////////////////////////////////////////////////////////////
// This is a common implementation of global functions and variables for all projects
// in the ARSS solution. Also includes all necessary external libraries. Hook functions
// are declared that are implemented (customized) in individual projects.
// 
// Author: Me (Naor)
///////////////////////////////////////////////////////////////////////////////

#include "ARSS.h"

// Request lint level warnings with the LINT macro on Microsoft compilers
#ifdef LINT
#pragma warning(push,4)
#pragma warning(disable:4100) // unreferenced formal parameter
#endif

// GLUT globals (window manager, event manager, rendering manager)
arss_glut_camera gCamera;
arss_glut_colors gColors;
arss_glut_controls gControls;
arss_glut_hud gHUD;
// Run globals (file names, run name, frame information etc.)
arss_run gRun;
// DevIL globals (screen capture parameters)
arss_devil gDevil;
// PhysX globals (PhysX Foundation and Scene)
arss_physx gPhysX;
// CUDA globals
arss_cuda gCUDA;
// Simulation globals (like physical parameters and numerical choices)
arss_simulation gSim;
// Experiment globals common to all experiments (projects)
arss_experiment gExp;
// Debugging globals
arss_debug gDebug;
// Utility globals
arss_utils gUtils;

// Initialization functions
bool InitRun(int argc, char **argv)
/*
This function sets global run parameters based on the command line arguments, and
reads others from an input file. The parameters set here are those common to all
projects in the ARSS solution. Parameters specific to each project are read in CustomizeRun(),
an implementation of which must exists in each projects main file. Note that parameters
may be read from the input file elsewhere as well.
*/
{
	char buf[MAX_CHARS_PER_NAME];

	// Get the project's name and current working directory
	if (argc==1)	// Default for program invocation with no arguments
	{
		_getcwd(buf,MAX_CHARS_PER_NAME);
		gRun.workingDirectory=buf;
		gRun.baseName=argv[0];
		size_t idx;
		idx = gRun.baseName.find_last_of(".");
		if (idx!=std::string::npos)
			gRun.baseName.erase(idx);
		idx = gRun.baseName.find_last_of("/\\");
		if (idx!=std::string::npos)
			gRun.baseName.erase(0,idx+1);
	}
	else			// The first command line argument, if it exists, is always the name of the run
	{
		_getcwd(buf,MAX_CHARS_PER_NAME);
		gRun.workingDirectory=buf;
		gRun.workingDirectory+="/";
		gRun.workingDirectory+=argv[1];
		gRun.baseName=argv[1];
	}
	
	// Read common program options from a file named run_base_name.ini
	gRun.iniFile = gRun.workingDirectory + "/" + gRun.baseName + ".ini";				
	if (!ConfigARSSOptions())
		ncc__warning("Could not find ARSS config file. Attempting to continue with all default options.\a");

	// A file named run_base_name.log may be used to record run stats
	gRun.logFile = gRun.workingDirectory + "/" + gRun.baseName + ".log";

	// A file named run_base_name.out may be used to record per-project output
	gRun.outFile = gRun.workingDirectory + "/" + gRun.baseName + ".out";

	// Set project specific parameters (from command line arguments and run_base_name.ini)
	CustomizeRun(argc,argv);

	return true;
}
bool InitGlut(int argc, char **argv)
{
	// Create top window
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(960,720);
	glutCreateWindow(argv[0]);

	// Enable/Disable optional OpenGL features
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glEnable(GL_LIGHTING);
	glEnable(GL_COLOR_MATERIAL);
	glEnableClientState(GL_VERTEX_ARRAY);

	// Register callbacks
	glutDisplayFunc(RenderScene);
	glutReshapeFunc(ReshapeCallback);
	glutIdleFunc(IdleCallback);
	glutKeyboardFunc(ProcessAsciiKeys);
	glutKeyboardUpFunc(ProcessKeyRelease);
	glutSpecialFunc(ProcessSpecialKeys);
	glutMouseFunc(MouseCallback);
	glutMotionFunc(MotionCallback);
	glutPassiveMotionFunc(NULL);
	glutEntryFunc(NULL);
	glutVisibilityFunc(NULL);
	atexit(ExitCallback);

	// Disable generation of callbacks during key repeat
	glutIgnoreKeyRepeat(1);

	// Set default colors and shading and lighting
	glClearColor(gColors.background[0],gColors.background[1],gColors.background[2],1.0f);
	gColors.defaultActor[0]=gColors.defaultActor[1]=gColors.defaultActor[2]=ncc::rgb::grGray[0];
	glShadeModel(GL_SMOOTH);
	GLfloat AmbientColor[] = {0.1f,0.1f,0.1f,1.0f}; glLightfv(GL_LIGHT0,GL_AMBIENT, AmbientColor); // add some ambient light to the scene
	GLfloat DiffuseColor[] = {1.0f,1.0f,1.0f,1.0f}; glLightfv(GL_LIGHT0,GL_DIFFUSE, DiffuseColor); // the miner's headlight color
	glEnable(GL_LIGHT0);

	// Initialize camera values
	gCamera.pos.x = 0.0f;		gCamera.pos.y = 0.0f;		gCamera.pos.z = 0.0f;
	gCamera.forward.x = 0.0f;	gCamera.forward.y = 0.0f;	gCamera.forward.z = -1.0f;
	gCamera.up.x = 0.0f;		gCamera.up.y = 1.0f;		gCamera.up.z = 0.0f;
	gCamera.right.x = 1.0f;	gCamera.right.y = 0.0f;		gCamera.right.z = 0.0f;
	gCamera.speed = 1.0f;
	gCamera.yFOV = 60.0f;
	gCamera.aspectRatio = (float) glutGet(GLUT_WINDOW_WIDTH)/glutGet(GLUT_WINDOW_HEIGHT);

	// Initialize controls values
	for (int k=0; k<MAX_KEYBOARD_KEYS; k++) gControls.keys[k]=false;
	for (int k=0; k<3; k++) gControls.mbut[k]=false;
	gControls.mxco = 0.0f; gControls.myco = 0.0f;
	gControls.mrotspd = 0.008f;

	CustomizeGLUT(); // Allow experiment specific initialization
	return true;
}
bool InitDevIL()
{
	ilInit();
	iluInit();
	ilutRenderer(ILUT_OPENGL);
	ilGenImages(1,&gDevil.screenShot);
	ilBindImage(gDevil.screenShot);
	return true;
}
bool InitHUD()
{
	// Set up ARSS-wide HUD elements
	gHUD.FPS		= gHUD.hud.AddElement("FPS = ",0.04,0.04);
	gHUD.WALL_TIME	= gHUD.hud.AddElement("Wall time = ", 0.04,0.08);
	gHUD.DBG_MSG	= gHUD.hud.AddElement("",0.5,0.66);
	gHUD.PAUSED		= gHUD.hud.AddElement("",0.4,0.5);

	// Set up experiment specific HUD elements
	CustomizeHUD();
	return true;
}
bool InitPhysX()
{
	// First, create a PxFoundation object
	static PxDefaultAllocator		defAllocator;
	static PxDefaultErrorCallback	defErrCallback;
	PxFoundation* mFoundation = PxCreateFoundation(PX_PHYSICS_VERSION,defAllocator,defErrCallback);
	if (!mFoundation)
		ncc__error("\aPxCreateFoundation failed!");

	// Next, create the top level PxPhysics object
	PxTolerancesScale scale;
	if (gPhysX.scaling.length > 0) scale.length = gPhysX.scaling.length; // default = 1
	if (gPhysX.scaling.mass   > 0) scale.mass   = gPhysX.scaling.mass;   // default = 1000
	if (gPhysX.scaling.speed  > 0) scale.speed  = gPhysX.scaling.speed;  // default = 10
	PxPhysics* mPhysics = PxCreatePhysics(PX_PHYSICS_VERSION,*mFoundation,scale,false);
	if (!mPhysics)
		ncc__error("\aPxCreatePhysics failed!");

	// And the Cooking library
	PxCooking* mCooking = PxCreateCooking(PX_PHYSICS_VERSION,*mFoundation,PxCookingParams());
	if (!mCooking)
		ncc__error("\aPxCreateCooking failed!");

	// And the Extensions library
	if (!PxInitExtensions(*mPhysics))
		ncc__error("\aPxInitExtensions failed!");

	// Also, try to connect to PVD, if it's not there, no harm done
#ifdef _DEBUG
	if (InitPVD(mPhysics)==true)
		cout << "\nPVD connected" << endl;
	else
		cout << "\nPVD not connected" << endl;
#endif

	// Initialize miscellaneous flags and names
	gPhysX.roles.everyone	=	PxActorTypeSelectionFlag::eRIGID_DYNAMIC | PxActorTypeSelectionFlag::eRIGID_STATIC;
	gPhysX.roles.statics	=	PxActorTypeSelectionFlag::eRIGID_STATIC;
	gPhysX.roles.dynamics	=	PxActorTypeSelectionFlag::eRIGID_DYNAMIC;

	// Save the top level objects for later (just for kicks)
	gPhysX.mFoundation = mFoundation;
	gPhysX.mPhysics = mPhysics;
	gPhysX.mCooking = mCooking;

	// Create a default material
	PxReal mud, mus, eps;
	char buf[MAX_CHARS_PER_NAME];
	ncc::GetStrPropertyFromINIFile("materials","default_friction_dyn","0.5",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	mud=atof(buf);
	ncc::GetStrPropertyFromINIFile("materials","default_friction_sta","0.5",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	mus=atof(buf);
	ncc::GetStrPropertyFromINIFile("materials","default_restitution","0.8",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	eps=atof(buf);
	gPhysX.mDefaultMaterial = gPhysX.mPhysics->createMaterial(mus,mud,eps);
	if (!gPhysX.mDefaultMaterial)
		ncc__error("\aPxPhysics::createMaterial() failed!");

	// Next create the scene
	PxSceneDesc sceneDesc(gPhysX.mPhysics->getTolerancesScale());
	if (!sceneDesc.cpuDispatcher) // provide default single threaded dispatcher
	{
		PxDefaultCpuDispatcher* mDispatch = PxDefaultCpuDispatcherCreate(1);
		if (!mDispatch) ncc__error("PxDefaultCPUDispatcherCreate failed!");
		sceneDesc.cpuDispatcher = mDispatch;
	}
	if (!sceneDesc.filterShader) // provide default filter shader
		sceneDesc.filterShader = PxDefaultSimulationFilterShader;

	sceneDesc.gravity = PxVec3(0.0f,-9.81f,0.0f); // convenient default, but usually modified later
	sceneDesc.flags |= PxSceneFlag::eENABLE_ONE_DIRECTIONAL_FRICTION; // replace default, constraint-based friction, with coulomb friction
	if (gPhysX.props.bounceThreshold>0) sceneDesc.bounceThresholdVelocity = gPhysX.props.bounceThreshold; // zero bounce threshold is unsustainable

	CustomizeScene(sceneDesc); // allow each experiment (project) to modify the scene constants

	gPhysX.mScene = gPhysX.mPhysics->createScene(sceneDesc);
	if (!gPhysX.mScene) ncc__error("Scene creation failed!");

	return true;
}
bool InitPVD(PxPhysics* mPhysics)
{
	// Check if a connection manager exists
	if (mPhysics->getPvdConnectionManager() == NULL)
		return false;

	// Set connection parameters
	const char*     pvd_host_ip = "127.0.0.1";  // IP of the PC which is running PVD
	int             port        = 5425;         // TCP port to connect to, where PVD is listening
	unsigned int    timeout     = 100;          // timeout in milliseconds (make longer for remote connections)
	PxVisualDebuggerConnectionFlags connectionFlags = PxVisualDebuggerExt::getAllConnectionFlags();

	// And try to connect
	PVD::PvdConnection* theConnection = PxVisualDebuggerExt::createConnection(mPhysics->getPvdConnectionManager(),pvd_host_ip,port,timeout,connectionFlags);
	if (theConnection == NULL)
		return false;

	gPhysX.mPVDConnection=theConnection;
	return true;
}
bool InitCUDA(int devToUse)
{
	bool success=false;
#ifdef HAVE_CUDA_TK
	int nbDevs;
	cudaError_t err;
	err = cudaGetDeviceCount(&nbDevs);
	switch (err)
	{
	case cudaErrorInsufficientDriver:
		ncc__warning("\aA CUDA driver is not installed, or the installed driver is too old.\n\tFalling back on software mode.");
		break;
	case cudaErrorNoDevice:
		printf("\aNo CUDA capable device detected.\n\tFalling back on software mode.");
		break;
	case cudaSuccess:
		success=true;
		if (devToUse>=nbDevs)
		{
			devToUse=nbDevs-1;
			printf("\a\n\nYou only have %d GPUs, I'll just use the last one.\n",nbDevs);
		}
		cudaSetDevice(devToUse);
		cudaGetDeviceProperties(&gCUDA.gpuProp,devToUse);

		printf("\n\nAlo World, my name is %s and I'll be your GPU today. Here are my specs:\n",gCUDA.gpuProp.name);
		printf("\tCompute capability            =  %d.%d\n",gCUDA.gpuProp.major,gCUDA.gpuProp.minor);
		printf("\tTotal global memory           =  %d kB\n",gCUDA.gpuProp.totalGlobalMem/1024);
		printf("\tTotal shared memory per block =  %d kB\n",gCUDA.gpuProp.sharedMemPerBlock/1024);
		printf("\tMax Grid size                 =  %d x %d x %d\n",gCUDA.gpuProp.maxGridSize[0],gCUDA.gpuProp.maxGridSize[1],gCUDA.gpuProp.maxGridSize[2]);
		printf("\tMax threads per block         =  %d x %d x %d\n",gCUDA.gpuProp.maxThreadsDim[0],gCUDA.gpuProp.maxThreadsDim[1],gCUDA.gpuProp.maxThreadsDim[2]);
		printf("\tMax threads per block         =  %d\n",gCUDA.gpuProp.maxThreadsPerBlock);
		printf("\tMulti-processors count        =  %d\n",gCUDA.gpuProp.multiProcessorCount);
		printf("\tKernel execution may time out =  %d\n",gCUDA.gpuProp.kernelExecTimeoutEnabled);
		printf("\tCompute mode                  =  %d\n",gCUDA.gpuProp.computeMode);
		printf("\tECC enabled                   =  %d\n",gCUDA.gpuProp.ECCEnabled);
		
		if (gCUDA.gpuProp.major<2)
		{
			ncc__warning("Compute capability < 2 is not going to cut it. Falling back on software mode.");
			success=false;
		}

		break;
	}
#endif
	gCUDA.cudaCapable=success;
	return success;
}
bool InitExperiment()
{
	CreateExperiment();
	return true;
}
// Cleanup function
void ExitCallback()
{
	// TODO: log some run stats
	LogRunStats();

	// The plan was to call DestroyPhysX() here, but that causes an unexplained crash.
}
void DestroyPhysX()
{
	PxCloseExtensions();
	gPhysX.mPhysics->release();
	gPhysX.mCooking->release();
	gPhysX.mFoundation->release();
}

// GLUT callbacks
void RenderScene()
{
	// Clear the buffers
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Set up camera and lights
	ProcessCameraControls();
	GLfloat Position[] = {gCamera.pos.x,gCamera.pos.y,gCamera.pos.z,1.0f};	// miner's hat light model
	glLightfv(GL_LIGHT0,GL_POSITION,Position);								//
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(gCamera.yFOV,gCamera.aspectRatio,gCamera.zBufNear,gCamera.zBufFar);
	gluLookAt(	gCamera.pos.x,                     gCamera.pos.y,                     gCamera.pos.z,
				gCamera.pos.x + gCamera.forward.x, gCamera.pos.y + gCamera.forward.y, gCamera.pos.z + gCamera.forward.z,
				gCamera.up.x,                      gCamera.up.y,                      gCamera.up.z);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// Draw the scene
	RenderActors();
	RenderOtherStuff();

	// Render the HUD
	gHUD.hud.Render();

	// Render debug information
	if (gDebug.bBoxGridOn)	Draw3DBoxGrid();
	if (gDebug.bXZGridOn)	DrawXZGrid();
	if (gDebug.bXYGridOn)	DrawXYGrid();

	// Swap the buffers
	glutSwapBuffers();
}
void IdleCallback()
{
	static unsigned int hudTime = 0;
	static unsigned int hudFrame = 0;
	static unsigned int renderTime = 0;
	static unsigned int outputFrame = 0;
	static unsigned int captureFrame = 0;
	static unsigned int physicsTime = 0;

	// Accumulate wall clock time
	gSim.wallTime = glutGet(GLUT_ELAPSED_TIME);

	// Request a redisplay
	if (gSim.wallTime - renderTime > 1000*gRun.renderFrequency) {
		renderTime = gSim.wallTime;
		glutPostRedisplay();
	}

	// Refresh  the HUD
	if (gSim.wallTime - hudTime > 1000) {
		gSim.fps = (gSim.frame-hudFrame);
		hudTime = gSim.wallTime;
		hudFrame = gSim.frame;
		RefreshHUD();
	}

	// Output something
	if ((gRun.outputFrequency > 0) && (gSim.frame - outputFrame >= gRun.outputFrequency)) {
		outputFrame = gSim.frame;
		LogExperiment();
	}

	// Take a screen shot
	if ((gRun.captureFrequency > 0) && (gSim.frame - captureFrame >= gRun.captureFrequency)) {
		captureFrame = gSim.frame;
		CaptureScreen();
	}

	// Advance the simulation
	if (gSim.isRunning && !gSim.bPause && ((gRun.physicsFrequency == 0) || ((gSim.wallTime - physicsTime) > (1000/gRun.physicsFrequency) )  ) ) {
		physicsTime = gSim.wallTime;
		AdvanceSimulation(gSim.timeStep);
	}
}
void ReshapeCallback(int w, int h)
{
	if (h>0)
		gCamera.aspectRatio = (float)glutGet(GLUT_WINDOW_WIDTH)/glutGet(GLUT_WINDOW_HEIGHT);
	glViewport(0,0,w,h);
	glutPostRedisplay();
}
void ProcessAsciiKeys(unsigned char key, int x, int y)
{
	gControls.keys[key] = true;
	
	switch (key)
	{
	case 27: /*ESC*/
		DestroyPhysX();
		exit(0);
		break;
	case 32: /*SPACE*/
		FireAction();
		break;
	case 'x': /*x or ALT-x*/
		if (glutGetModifiers()==GLUT_ACTIVE_ALT)
			PrintDebug();
		break;
	case 'g': /*g or ALT-g*/
		if (glutGetModifiers()==GLUT_ACTIVE_ALT) gDebug.bXYGridOn = !gDebug.bXYGridOn;
		break;
	case '?': /*WakeUp*/
		Reveille();
		break;
	case 'p': /*Toggle Pause*/
		gSim.bPause=!gSim.bPause;
		RefreshHUD();
		break;
	//default: // uncomment to find out unknown ascii codes
	//	printf("%d\n",key);
	}
}
void ProcessKeyRelease(unsigned char key, int x, int y)
{
	gControls.keys[key] = false;
}
void ProcessSpecialKeys(int key, int x, int y)
{
	switch (key)
	{
	case GLUT_KEY_F10: // reboot experiment
		RebootExperiment();
		break;
	case GLUT_KEY_F9: // screen capture
		CaptureScreen();
		break;
	case GLUT_KEY_F7: // load scene from file
		LoadSceneFromFile(gRun.loadSceneFromFile);
		break;
	case GLUT_KEY_F5: // save scene to file
		SaveSceneToRepXDump();
		break;
	case GLUT_KEY_UP:
		UpArrowAction();
		break;
	case GLUT_KEY_DOWN:
		DownArrowAction();
		break;
	case GLUT_KEY_LEFT:
		LeftArrowAction();
		break;
	case GLUT_KEY_RIGHT:
		RightArrowAction();
		break;
	}
}
void MouseCallback(int but, int state, int x, int y)
{
	gControls.mxco=x;
	gControls.myco=y;
}
void MotionCallback(int x, int y)
{
	int dx = gControls.mxco - x;
	int dy = gControls.myco - y;

	gCamera.forward.normalize();
	gCamera.up.normalize();
	gCamera.right = gCamera.forward.cross(gCamera.up);
	gCamera.right.normalize(); // I know this seems unnecessary, but PxQuat expects a *very precisely unit* quaternion

	PxQuat qx(dx * 0.001f , gCamera.up);
	PxQuat qy(dy * 0.001f , gCamera.right);
	gCamera.forward = qx.rotate(gCamera.forward);
	gCamera.forward = qy.rotate(gCamera.forward);

	//gCamera.up = gCamera.right.cross(gCamera.forward); // more precise but more nausea inducing motion

	gControls.mxco = x;
	gControls.myco = y;
}
void ProcessCameraControls()
{
	float dt=0.08;
	for (int k=0; k<MAX_KEYBOARD_KEYS; k++)
	{
		if (!gControls.keys[k]) continue;
		switch (k)
		{
		case 'w':
			gCamera.pos += gCamera.forward * gCamera.speed * dt;
			break;
		case 's':
			gCamera.pos -= gCamera.forward * gCamera.speed * dt;
			break;
		case 'a':
			gCamera.pos -= gCamera.right * gCamera.speed * dt;
			break;
		case 'd':
			gCamera.pos += gCamera.right * gCamera.speed * dt;
			break;
		case 'q':
			gCamera.pos += gCamera.up * gCamera.speed * dt;
			break;
		case 'z':
			gCamera.pos -= gCamera.up * gCamera.speed * dt;
			break;
		}
	}
}

// All project functions
void AdvanceSimulation(PxReal dt)
/* This is the wrapper function for time stepping in all projects. It starts with an asynchronous
// call to PxScene::simulate(). This call will advance the world based on the current state.
// Next, a call to a generic, experiment specific function is made. Each experiment must provide
// a void ApplyCustomInteractions() function that may accumulate external forces to be used in the next
// time step. Here, external means anything that is not due to collisions between actors. This
// function may also manipulate actors, or the scene, directly. When this function returns,
// the results of the current state are fetched with a blocking call to PxScene::fetchResults,
// the time and frame count are advanced, and control is returned to the calling function.
*/
{
	if (dt>0.0f)
	{
		gPhysX.mScene->simulate(dt); // async call, returns immediately
		ApplyCustomInteractions(); // each project must provide its own version of this function
		gPhysX.mScene->fetchResults(true); // blocks until simulate returns
		gSim.frame++;
		gSim.codeTime+=dt;
	}
}
void RefreshHUD()
{
	char s[256];
	sprintf(s,"Wall time = %u",gSim.wallTime/1000U);
	gHUD.hud.SetElement(gHUD.WALL_TIME,s);
	sprintf(s,"FPS = %5.1f",gSim.fps);
	gHUD.hud.SetElement(gHUD.FPS,s);
	if (gSim.bPause)
		sprintf(s,"PAUSED - Hit `p' to unpause");
	else
		sprintf(s,"");
	gHUD.hud.SetElement(gHUD.PAUSED,s);

	RefreshCustomHUDElements(); // Project specific HUD elements, implemented in project source
}
void LogRunStats()
{

}
void CaptureScreen()
{
	char imgName[FILENAME_MAX];
	sprintf(imgName,"%s/%s.%d.jpg",gRun.workingDirectory.c_str(),gRun.baseName.c_str(),gSim.frame);
	ilutGLScreen();
	ilSaveImage(imgName);	
}
bool ConfigARSSOptions()
{
	// First check that an options file exists
	ifstream fp(gRun.iniFile.c_str());
	bool success = fp.good();
	fp.close();

	// Read in parameters by group, file.ini style
	char buf[MAX_CHARS_PER_NAME];

	// Program options group
	gRun.outputFrequency	= ncc::GetIntPropertyFromINIFile("program","output_frequency", 0,gRun.iniFile.c_str());
	gRun.captureFrequency	= ncc::GetIntPropertyFromINIFile("program","capture_frequency",0,gRun.iniFile.c_str());
	gRun.renderFrequency	= ncc::GetIntPropertyFromINIFile("program","render_frequency" ,0,gRun.iniFile.c_str());
	gRun.physicsFrequency	= ncc::GetIntPropertyFromINIFile("program","physics_frequency",0,gRun.iniFile.c_str());

	ncc::GetStrPropertyFromINIFile("program","scene_file","new",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if (strcmp(buf,"new")==0)
		gRun.loadSceneFromFile.clear();
	else
		gRun.loadSceneFromFile = gRun.workingDirectory + "/" + buf;

	// Simulation options group
	ncc::GetStrPropertyFromINIFile("simulation","time_step","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	gSim.timeStep = atof(buf);

	// GLUT options group
	ncc::GetStrPropertyFromINIFile("glut","z_buffer_far","10",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	gCamera.zBufFar = atof(buf);
	ncc::GetStrPropertyFromINIFile("glut","z_buffer_near","0.1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	gCamera.zBufNear = atof(buf);
	ncc::GetStrPropertyFromINIFile("glut","xy_grid","off",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if (strcmp(buf,"on")==0) gDebug.bXYGridOn = true;
	ncc::GetStrPropertyFromINIFile("glut","xz_grid","off",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if (strcmp(buf,"on")==0) gDebug.bXZGridOn = true;

	// PhysX options
	ncc::GetStrPropertyFromINIFile("PhysX","sleep_threshold",     "-1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str()); // default is scaled
	gPhysX.props.sleepThreshold = atof(buf);
	ncc::GetStrPropertyFromINIFile("PhysX","bounce_threshold",    "-1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str()); // default NOT scaled (=2)
	gPhysX.props.bounceThreshold = atof(buf);
	ncc::GetStrPropertyFromINIFile("PhysX","linear_damping",      "-1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str()); // default is zero
	gPhysX.props.linearDamping = atof(buf);
	ncc::GetStrPropertyFromINIFile("PhysX","angular_damping",     "-1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str()); // default is NOT scaled (=0.05)
	gPhysX.props.angularDamping = atof(buf);
	ncc::GetStrPropertyFromINIFile("PhysX","max_angular_velocity","-1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str()); // default is NOT scaled (=7)
	gPhysX.props.maxAngularVelocity = atof(buf);
	ncc::GetStrPropertyFromINIFile("PhysX","skin_width",          "-1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str()); // default is scaled
	gPhysX.props.skinWidth = atof(buf);

	// Scaling parameters
	ncc::GetStrPropertyFromINIFile("scaling","length","-1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	gPhysX.scaling.length = atof(buf);
	ncc::GetStrPropertyFromINIFile("scaling","mass",  "-1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	gPhysX.scaling.mass = atof(buf);
	ncc::GetStrPropertyFromINIFile("scaling","speed", "-1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	gPhysX.scaling.speed = atof(buf);

	// Experiment group options common to all ARSS projects - Grain type
	ncc::GetStrPropertyFromINIFile("experiment","grain_type","",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if		(strcmp(buf,"sphere")==0)
		gExp.defGrainType=eSPHERE_GRAIN;
	else if (strcmp(buf,"box")==0)
		gExp.defGrainType=eBOX_GRAIN;
	else if (strcmp(buf,"capsule")==0)
		gExp.defGrainType=eCAPSULE_GRAIN;
	else if (strcmp(buf,"convex")==0)
		gExp.defGrainType=eCONVEX_GRAIN;
	else if (strcmp(buf,"pyramid")==0)
		gExp.defGrainType=ePYRAMID_GRAIN;
	else if (strcmp(buf,"d12")==0)
		gExp.defGrainType=eD12_GRAIN;
	else if (strcmp(buf,"random_convex")==0)
		gExp.defGrainType=eRAND_CONVEX;
	else
		gExp.defGrainType=eBAD_RUBBLE_GRAIN_TYPE;

	// Experiment group options common to all ARSS projects - Grain size and scatter
	ncc::GetStrPropertyFromINIFile("experiment","grain_size","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	gExp.defGrainSize = atof(buf);
	ncc::GetStrPropertyFromINIFile("experiment","grain_size_scatter","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	gExp.defGrainSizeScatter = atof(buf);
	ncc::GetStrPropertyFromINIFile("experiment","uniform_rubble","false",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	if (strcmp(buf,"true")==0)
		gExp.bUniformRubble=true;
	else
		gExp.bUniformRubble=false;

	// Experiment group options common to all ARSS projects - Grain density
	ncc::GetStrPropertyFromINIFile("experiment","grain_density","1",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	gExp.defGrainDensity = atof(buf);


	return success;
}
void UpdateIntegralsOfMotion()
{
	// Reset global values
	gExp.IOMs.systemAM = PxVec3(0);
	gExp.IOMs.systemLM = PxVec3(0);
	gExp.IOMs.systemCM = PxVec3(0);
	gExp.IOMs.systemKE = 0.0f;
	gExp.IOMs.systemMass = 0.0f;
	PxReal M_system = 0;
	
	// Loop over dynamic actors
	PxU32 nbActors = gPhysX.mScene->getNbActors(gPhysX.roles.dynamics);
	gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,nbActors);
	for (PxU32 k=0; k<nbActors; k++)
	{
		PxRigidDynamic* actor = gPhysX.cast[k]->isRigidDynamic();
		if (actor)
		{
			PxTransform globalPose = actor->getGlobalPose().transform(actor->getCMassLocalPose());
			PxReal	M_actor = actor->getMass();
			PxVec3	v_actor = actor->getLinearVelocity(); 
			PxVec3	w_actor = actor->getAngularVelocity();
			PxVec3	r_actor = globalPose.p;
			PxMat33 I_actor = PxMat33::createDiagonal(actor->getMassSpaceInertiaTensor()); // in body frame!

			// system mass
			M_system += M_actor;

			// system center of mass
			gExp.IOMs.systemCM += M_actor * r_actor;

			// linear momentum
			gExp.IOMs.systemLM += M_actor * v_actor;

			// kinetic energy
			PxVec3 w_massFrame = globalPose.rotateInv(w_actor);
			gExp.IOMs.systemKE += 0.5 * M_actor * v_actor.magnitudeSquared();
			gExp.IOMs.systemKE += 0.5 * w_massFrame.dot(I_actor * w_massFrame);

			// angular momentum
			gExp.IOMs.systemAM += globalPose.rotate(I_actor * w_massFrame);
			gExp.IOMs.systemAM += M_actor * r_actor.cross(v_actor);

		}
	}
	if (M_system) {
		PxReal perUnitMass=(1.0f/M_system);
		gExp.IOMs.systemCM *= perUnitMass;
		gExp.IOMs.systemLM *= perUnitMass;
		gExp.IOMs.systemKE *= perUnitMass;
		gExp.IOMs.systemAM *= perUnitMass;
		gExp.IOMs.systemMass = M_system;
	}

}
void CreateGroundPlane()
{
	if (gPhysX.mPhysics && gPhysX.mDefaultMaterial)
	{
		gExp.VIPs.groundPlane = PxCreatePlane(*gPhysX.mPhysics,PxPlane(PxVec3(0,1,0),0),*gPhysX.mDefaultMaterial);
		if (!gExp.VIPs.groundPlane)
			ncc__error("ground plane creation failed!");
		gExp.VIPs.groundPlane->setName("ground");
		gPhysX.mScene->addActor(*gExp.VIPs.groundPlane);
	}
}
PxU32 CountSleepers()
{
	PxU32 sleepers = 0;
	PxU32 nbActors = gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,MAX_ACTORS_PER_SCENE);
	while (nbActors--)
	{
		PxRigidDynamic* actor = gPhysX.cast[nbActors]->isRigidDynamic();
		if (actor && actor->isSleeping())
			sleepers++;
	}
	return sleepers;
}
void Reveille()
{
	PxU32 nbActors = gPhysX.mScene->getNbActors(gPhysX.roles.dynamics);
	if (nbActors)
	{
		gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,nbActors);
		for (PxU32 k=0; k<nbActors; k++)
		{
			(gPhysX.cast[k]->isRigidDynamic())->wakeUp();
		}
	}
}
void DeadStop()
{
	PxU32 nbActors = gPhysX.mScene->getNbActors(gPhysX.roles.dynamics);
	gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,nbActors);
	while (nbActors--)
	{
		PxRigidDynamic* actor = gPhysX.cast[nbActors]->isRigidDynamic();
		if (actor->getRigidDynamicFlags() & PxRigidDynamicFlag::eKINEMATIC) continue;
		actor->setLinearVelocity(PxVec3(0));
		actor->setAngularVelocity(PxVec3(0));
	}
}
void ReCenterActors()
{
	UpdateIntegralsOfMotion();
	PxU32 nbActors = gPhysX.mScene->getNbActors(gPhysX.roles.dynamics);
	gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,nbActors);
	while (nbActors--)
	{
		PxRigidDynamic* actor = gPhysX.cast[nbActors]->isRigidDynamic();
		PxTransform pose = actor->getGlobalPose();
		pose.p -= gExp.IOMs.systemCM;
		actor->setGlobalPose(pose);
	}
}
void FindExtremers()
{
	gExp.VIPs.extremers.downmost=gExp.VIPs.extremers.upmost		= NULL;
	gExp.VIPs.extremers.leftmost=gExp.VIPs.extremers.rightmost	= NULL;
	gExp.VIPs.extremers.inmost=gExp.VIPs.extremers.outmost		= NULL;
	if (gPhysX.mScene)
	{
		PxReal right	= -PX_MAX_REAL,	up		= -PX_MAX_REAL,	out	= -PX_MAX_REAL;
		PxReal left		= +PX_MAX_REAL,	down	= +PX_MAX_REAL,	inn	= +PX_MAX_REAL;
		PxU32 nbActors = gPhysX.mScene->getNbActors(gPhysX.roles.dynamics);
		if (nbActors)
		{
			gPhysX.mScene->getActors(gPhysX.roles.dynamics,gPhysX.cast,nbActors);
			for (PxU32 k=0; k<nbActors; k++)
			{
				PxRigidDynamic* anActor = gPhysX.cast[k]->isRigidDynamic();
				if (anActor==NULL) continue; // shouldn't happen though
				PxVec3 pos = anActor->getGlobalPose().p;
				if (pos.x > right)	{right = pos.x;	gExp.VIPs.extremers.rightmost = anActor;}
				if (pos.x < left)	{left = pos.x;	gExp.VIPs.extremers.leftmost = anActor;}
				if (pos.y > up)		{up = pos.y;	gExp.VIPs.extremers.upmost = anActor;}
				if (pos.y < down)	{down = pos.y;	gExp.VIPs.extremers.downmost = anActor;}
				if (pos.z > out)	{out = pos.z;	gExp.VIPs.extremers.outmost = anActor;}
				if (pos.z < inn)	{inn = pos.z;	gExp.VIPs.extremers.inmost = anActor;}
			}
		}
	}
}
bool SaveSceneToRepXDump()
{
	// Create a default RepX collection
	repx::RepXCollection* theCollection = repx::createCollection(gPhysX.mPhysics->getTolerancesScale(),PxGetFoundation().getAllocatorCallback());
	repx::RepXIdToRepXObjectMap* theIdMap = repx::RepXIdToRepXObjectMap::create(PxGetFoundation().getAllocatorCallback());

	// Populate collection with all SKD and scene items
	repx::addItemsToRepX(*gPhysX.mPhysics,*gPhysX.mScene,*theIdMap,*theCollection);

	// Save collection to an output file
	char buf[MAX_CHARS_PER_NAME];
	string dumpFile = gRun.workingDirectory + "/" + gRun.baseName + ".";
	if (sprintf(buf,"%d",gSim.frame)>0) dumpFile += buf;
	dumpFile += ".repx";

	PxDefaultFileOutputStream buffer(dumpFile.c_str());
	theCollection->save(buffer);
	cout << "Scene dump saved to " << dumpFile << endl;

	// Clean up
	theCollection->destroy();
	theIdMap->destroy();

	return true;
}
bool LoadSceneFromFile(string sceneFile, bool dynamicsOnly/*=false*/)
{
	if (sceneFile.empty())
	{
		ncc__warning("Empty scene file, nothing loaded.");
		return false;
	}

	PxDefaultFileInputData buffer(sceneFile.c_str());
	if (!buffer.isValid())
	{
		ncc__warning("Scene File not found. Nothing loaded.");
		return false;
	}

	// Create a RepX collection from the scene file
	repx::RepXCollection* theCollection = repx::createCollection(buffer,PxGetFoundation().getAllocatorCallback());

	// Add actors from collection to scene. Note: a valid scene should be ready to receive them or a crash is likely
	PxStringTable* strTable = &PxStringTableExt::createStringTable(PxGetFoundation().getAllocatorCallback());
	if (dynamicsOnly)
	{
		// Deal with this later. Look in RepXUtility.h, struct RepXCoreItemAdder to see how.
	} 
	else
	{
		repx::addObjectsToScene(*theCollection,*gPhysX.mPhysics,*gPhysX.mCooking,*gPhysX.mScene,strTable,NULL);
	}
	
	// Clean up
	theCollection->destroy();

	return true;
}
PxRigidDynamic* CreateRubbleGrain(PxVec3 startPosition, RubbleGrainType grainType, PxReal sizeScale, PxMaterial& material, PxReal density/*=1*/, bool isUniform/*=true*/, PxReal sizeScatter/*=1.0*/, vector<PxVec3> *verts/*=NULL*/)
/*
This is the high level actor creator for standard rubble grains. When invoked, this
function will dispatch an appropriate PxActor creator based on grainType and then
validate the created object and add it to the scene. A handle to the created actor
is returned so that it may be assigned more properties (e.g. velocity, spin, etc.)
*/
{
	PxRigidDynamic* actor=NULL;

	// Dispatch a PxActor creator
	switch (grainType)
	{
	case eSPHERE_GRAIN:
		actor = CreateDynamicSphericalGrain(startPosition,sizeScale,material,density,isUniform,sizeScatter);
		break;
	case eBOX_GRAIN:
		actor = CreateDynamicBoxGrain(startPosition,sizeScale,material,density,isUniform,sizeScatter);
		break;
	case eCAPSULE_GRAIN:
		actor = CreateDynamicCapsuleGrain(startPosition,sizeScale,material,density,isUniform,sizeScatter);
		break;
	case eCONVEX_GRAIN:
		if (verts)
		{actor = CreateDynamicConvexGrain(startPosition,*verts,material,density);}
		else
		{ncc__warning("Missing vertex list for CONVEX grain type - nothing was created.\a");}
		break;
	case ePYRAMID_GRAIN:
		{
			vector<PxVec3> pVerts = MakePyramidVertexList(sizeScale,isUniform,sizeScatter);
			actor = CreateDynamicConvexGrain(startPosition,pVerts,material,density);
		}
		break;
	case eD12_GRAIN:
		{
			vector<PxVec3> dVerts = MakeD12VertexList(sizeScale);
			actor = CreateDynamicConvexGrain(startPosition,dVerts,material,density);
		}
		break;
	case eRAND_CONVEX:
		{
			vector<PxVec3> rVerts = MakeRandomVertexList(sizeScale);
			actor = CreateDynamicConvexGrain(startPosition,rVerts,material,density);
		}
		break;
	case eBAD_RUBBLE_GRAIN_TYPE: // intentional fall-through
	default:
		ncc__warning("Unknown grain type - nothing was created.\a");
		actor = NULL;
	}

	// Add a validated actor to the scene
	if (actor)
	{
		// Set important actor flags and properties. Most are read from the .ini file.
		if (gPhysX.props.sleepThreshold > -1) actor->setSleepThreshold(0.5*gPhysX.props.sleepThreshold*gPhysX.props.sleepThreshold);
		if (gPhysX.props.angularDamping > -1) actor->setAngularDamping(gPhysX.props.angularDamping);
		if (gPhysX.props.linearDamping  > -1) actor->setLinearDamping(gPhysX.props.linearDamping);
		if (gPhysX.props.skinWidth      > +0) // set shape contact offset
		{
			PxShape* shapes[MAX_SHAPES_PER_ACTOR];
			PxU32 nbShapes = actor->getShapes(shapes,MAX_SHAPES_PER_ACTOR);
			while (nbShapes--) shapes[nbShapes]->setContactOffset(gPhysX.props.skinWidth);
		}

		gPhysX.mScene->addActor(*actor);
	}

	return actor;
}
vector<PxVec3> MakePyramidVertexList(PxReal sizeScale/*=1.0f*/, bool isUniform/*=true*/, PxReal sideScatter/*=1.0f*/)
{
	vector<PxVec3> verts(5);
	PxReal p;
	p = PxAbs ( sizeScale * (1 + (!isUniform)*sideScatter*(-0.5+gUtils.uniformRNG->doub())) );
	verts[0] = PxVec3(0,p,0);
	p = PxAbs ( sizeScale * (1 + (!isUniform)*sideScatter*(-0.5+gUtils.uniformRNG->doub())) );
	verts[1] = PxVec3(p,0,0);
	p = PxAbs ( sizeScale * (1 + (!isUniform)*sideScatter*(-0.5+gUtils.uniformRNG->doub())) );
	verts[2] = PxVec3(-p,0,0);
	p = PxAbs ( sizeScale * (1 + (!isUniform)*sideScatter*(-0.5+gUtils.uniformRNG->doub())) );
	verts[3] = PxVec3(0,0,p);
	p = PxAbs ( sizeScale * (1 + (!isUniform)*sideScatter*(-0.5+gUtils.uniformRNG->doub())) );
	verts[4] = PxVec3(0,0,-p);
	return verts;
}
vector<PxVec3> MakeRandomVertexList( PxReal sizeScale/*=1.0f*/, PxU32 numVerts/*=12*/, PxReal gama/*=10*/ )
/*
Generate numVerts random vertices on a sphere, with no two closer the gama degrees.
*/
{
	gama*=PxPi/180.0;
	vector<PxVec3> verts(0);
	while (verts.size() < numVerts)
	{
		// generate candidate vector
		PxReal x = gUtils.normalRNG->dev();
		PxReal y = gUtils.normalRNG->dev();
		PxReal z = gUtils.normalRNG->dev();
		PxVec3 candVert(x,y,z); candVert.normalize(); candVert *= sizeScale;

		// reject the candidate is angular opening to an existing vector is less than gama
		bool tooClose=false;
		for (PxU32 k=0; k<verts.size(); k++)
		{
			PxReal dp = verts.at(k).dot(candVert);
			if (dp > cos(gama) * sizeScale *sizeScale) tooClose=true;
		}
		if (!tooClose) verts.push_back(candVert);
	}

	return verts;
}
vector<PxVec3> MakeD12VertexList( PxReal sizeScale/*=1.0f*/ )
{
	vector<PxVec3> verts(20);
	PxReal p=(1.0+sqrt(5.0))/2.0;
	PxReal q=1.0/p;
	PxReal fac=sizeScale*p/2.0; // sizeScale will be the edge length

	// Note: the actor frame and center-of-mass frame will have different origins
	verts[0]  = PxVec3(+1,+1,+1) * fac;
	verts[1]  = PxVec3(+1,+1,-1) * fac;
	verts[2]  = PxVec3(+1,-1,+1) * fac;
	verts[3]  = PxVec3(+1,-1,-1) * fac;
	verts[4]  = PxVec3(-1,+1,+1) * fac;
	verts[5]  = PxVec3(-1,-1,+1) * fac;
	verts[6]  = PxVec3(-1,+1,-1) * fac;
	verts[7]  = PxVec3(-1,-1,-1) * fac;
	verts[8]  = PxVec3(+0,+q,+p) * fac;
	verts[9]  = PxVec3(+0,-q,+p) * fac;
	verts[10] = PxVec3(+0,+q,-p) * fac;
	verts[11] = PxVec3(+0,-q,-p) * fac;
	verts[12] = PxVec3(+q,+p,+0) * fac;
	verts[13] = PxVec3(-q,+p,+0) * fac;
	verts[14] = PxVec3(+q,-p,+0) * fac;
	verts[15] = PxVec3(-q,-p,+0) * fac;
	verts[16] = PxVec3(+p,+0,+q) * fac;
	verts[17] = PxVec3(-p,+0,+q) * fac;
	verts[18] = PxVec3(+p,+0,-q) * fac;
	verts[19] = PxVec3(-p,+0,-q) * fac;

	return verts;
}
PxRigidDynamic* CreateDynamicSphericalGrain(PxVec3 position, PxReal radius, PxMaterial& material, PxReal density, bool isUniform, PxReal radiusScatter)
{
	PxReal rad = radius * (1 + (!isUniform)*radiusScatter*(-0.5+gUtils.uniformRNG->doub()));
	return PxCreateDynamic(*gPhysX.mPhysics,PxTransform(position),PxSphereGeometry(PxAbs(rad)),material,density);
}
PxRigidDynamic* CreateDynamicBoxGrain(PxVec3 position, PxReal side, PxMaterial& material, PxReal density, bool isUniform, PxReal sideScatter)
{
	PxReal hx, hy, hz;
	hx = PxAbs( side * (1 + (!isUniform)*sideScatter*(-0.5+gUtils.uniformRNG->doub())) );
	hy = PxAbs( side * (1 + (!isUniform)*sideScatter*(-0.5+gUtils.uniformRNG->doub())) );
	hz = PxAbs( side * (1 + (!isUniform)*sideScatter*(-0.5+gUtils.uniformRNG->doub())) );
	return PxCreateDynamic(*gPhysX.mPhysics,PxTransform(position),PxBoxGeometry(hx/2,hy/2,hz/2),material,density);
}
PxRigidDynamic* CreateDynamicCapsuleGrain(PxVec3 position, PxReal radius, PxMaterial& material, PxReal density, bool isUniform, PxReal radiusScatter)
{
	PxReal rad = radius * (1 + (!isUniform)*radiusScatter*(-0.5+gUtils.uniformRNG->doub()));
	PxReal halfHeight = 1.0f*rad;
	return PxCreateDynamic(*gPhysX.mPhysics,PxTransform(position),PxCapsuleGeometry(rad,halfHeight),material,density);
}
PxRigidDynamic* CreateDynamicConvexGrain(PxVec3 position, vector<PxVec3> verts, PxMaterial& material, PxReal density)
{
	if (verts.size()==0) return NULL;

	// Describe a convex mesh
	PxConvexMeshDesc		convDesc;
	convDesc.points.count	=	verts.size();
	convDesc.points.stride	=	sizeof(PxVec3);
	convDesc.points.data	=	verts.data();
	convDesc.flags			=	PxConvexFlag::eCOMPUTE_CONVEX;

	// Cook the descriptor into a mesh
	PxDefaultMemoryOutputStream buf;
	bool cooked = gPhysX.mCooking->cookConvexMesh(convDesc,buf);
	if (!cooked)
	{
		ncc__warning("Mesh cooking failed!\a");
		return NULL;
	}
	PxDefaultMemoryInputData input(buf.getData(),buf.getSize());
	PxConvexMesh* theMesh = gPhysX.mPhysics->createConvexMesh(input);

	// Finally, create a single shape actor referencing the mesh
	return PxCreateDynamic(*gPhysX.mPhysics,PxTransform(position),PxConvexMeshGeometry(theMesh),material,density);
}
PxRigidDynamic* CreateDynamicConvexGrain(PxVec3 position, PxConvexMesh* mesh, PxMaterial& material, PxReal density)
{
	if (mesh)
		return PxCreateDynamic(*gPhysX.mPhysics,PxTransform(position),PxConvexMeshGeometry(mesh),material,density);
	else
		return NULL;
}
PxConvexMesh* MakePxMeshFromVertexList(vector<PxVec3> verts)
{
	if (verts.size()==0) return NULL;

	// Describe a convex mesh
	PxConvexMeshDesc		convDesc;
	convDesc.points.count	=	verts.size();
	convDesc.points.stride	=	sizeof(PxVec3);
	convDesc.points.data	=	verts.data();
	convDesc.flags			=	PxConvexFlag::eCOMPUTE_CONVEX;

	// Cook the descriptor into a mesh
	PxDefaultMemoryOutputStream buf;
	bool cooked = gPhysX.mCooking->cookConvexMesh(convDesc,buf);
	if (!cooked)
	{
		ncc__warning("Mesh cooking failed!\a");
		return NULL;
	}
	PxDefaultMemoryInputData input(buf.getData(),buf.getSize());
	PxConvexMesh* theMesh = gPhysX.mPhysics->createConvexMesh(input);

	return theMesh;
}
void RandOrientActor(PxRigidActor* actor)
{
	if (!actor) return;
	
	PxReal x = gUtils.normalRNG->dev();
	PxReal y = gUtils.normalRNG->dev();
	PxReal z = gUtils.normalRNG->dev();
	PxReal t = gUtils.uniformRNG->doub() * PxPi;

	PxVec3 n(x,y,z); n.normalize();
	PxQuat q(t,n);
	PxVec3 p = actor->getGlobalPose().p;

	actor->setGlobalPose(PxTransform(p,q));
}
void RandLaunchActor(PxRigidDynamic* actor, PxReal kickMag)
{
	if (!actor) return;

	PxReal x = gUtils.normalRNG->dev();
	PxReal y = gUtils.normalRNG->dev();
	PxReal z = gUtils.normalRNG->dev();

	PxVec3 v(x,y,z); v.normalize();
	v *= kickMag;

	actor->setLinearVelocity(v);
}
void RandSpinActor(PxRigidDynamic* actor, PxReal spinMag)
{
	if (!actor) return;

	PxReal x = gUtils.normalRNG->dev();
	PxReal y = gUtils.normalRNG->dev();
	PxReal z = gUtils.normalRNG->dev();

	PxVec3 w(x,y,z); w.normalize();
	w *= spinMag;
	
	actor->setAngularVelocity(w);
}

// Rendering functions
void RenderActors()
{
	PxU32 nbActors = gPhysX.mScene->getNbActors(gPhysX.roles.everyone);
	if (nbActors)
	{
		gPhysX.mScene->getActors(gPhysX.roles.everyone,gPhysX.cast,nbActors);
		for (PxU32 k=0; k<nbActors; k++)
		{
			DrawActor(static_cast<PxRigidActor*>(gPhysX.cast[k]));
		}
	}
}
void DrawActor(PxRigidActor* actor)
/*
This is the high level drawing routine for rigid body objects in any ARSS project.
Most of the work here is selecting the rendering color to use based on actor properties.
Then the shapes associated with the actor are extracted and sent, one by one, to a shape
drawing routine, which may further modify color before dispatching the shape to a low-level
openGL geometry drawer.
*/
{
	// Actors' names may signify special treatment
	const char* name=actor->getName();
	if (name && name[0]==126) return; // a name that starts with a tilde (~) signifies invisibility
	bool wireframe = false;
	if (name && name[0]==35) wireframe = true; // a name that starts with a pound (#) signifies wireframe mode

	// Save current OpenGL attributes that may be messed up during drawing
	glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT);

	// Choose a display color - start with the default actor color
	glColor3ubv(gColors.defaultActor);

	// A sleeping actor is drawn in blue
	if (actor->isRigidDynamic()) {
		PxRigidDynamic* abody = static_cast<PxRigidDynamic*>(actor);
		if (abody->isSleeping()) glColor3ubv(ncc::rgb::bBlue);
	}

	// Sometimes we save a color in the actor's user data field, which should point to a vector<GLUbyte>
	if (actor->userData) {
		glColor3ubv(static_cast<GLubyte*>(actor->userData));
	}

	// Get and draw the actor's shapes (often, but not always, there is just one)
	PxShape* shapes[MAX_SHAPES_PER_ACTOR];
	PxU32 nbShapes = actor->getNbShapes();
	if (nbShapes)
	{
		actor->getShapes(shapes,MAX_SHAPES_PER_ACTOR);
		while (nbShapes--)
		{
			const char* name = shapes[nbShapes]->getName();
			if (name && name[0]==126) continue;
			bool wf = false;
			if (name && name[0]==35) wf = true;
			DrawShape(shapes[nbShapes], (wireframe || wf));
		}
	}

	// Restore OpenGL attributes to whatever they were (probably irrelevant but who knows)
	glPopAttrib();
}
void ColorActor(PxActor* actor, const GLubyte* color)
{
	gColors.colorBucket.push_back(vector<GLubyte>(3));
	gColors.colorBucket.back()[0]=color[0];
	gColors.colorBucket.back()[1]=color[1];
	gColors.colorBucket.back()[2]=color[2];
	actor->userData=&(gColors.colorBucket.back()[0]);
}

// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif