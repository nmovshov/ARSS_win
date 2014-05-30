///////////////////////////////////////////////////////////////////////////////
// This is a common header of global functions and variables for all projects in
// the ARSS solution. Also includes all necessary external libraries. Hook functions
// are declared that are implemented (customized) in individual projects.
// 
// Author: Me (Naor)
///////////////////////////////////////////////////////////////////////////////

#ifndef _ARSS_H_
#define _ARSS_H_

// Necessary headers
#ifdef __APPLE__
#include "GLUT/glut.h"
#else
#include "GL/glut.h"					// The GLUT window and event manager
#endif
#include "Drawers.h"					// OpenGL drawing routines
#include "ncclib.h"						// Naor's C++ collection, with some Numerical Recipes
#include "HUD.h"						// A quick-and-dirty HUD for GLUT
#include "IL/il.h"						// |
#include "IL/ilu.h"						// | The DevIL library used to capture screen shots to JPEG files
#include "IL/ilut.h"					// |
#include "PxPhysicsAPI.h"				// The entire PhysX API in a single header (including many functions we don't need but this is easiest).
using namespace physx;					//
#ifdef HAVE_CUDA_TK
#include "cuda.h"						// |
#include "cuda_runtime.h"				// | The CUDA API
#include "device_launch_parameters.h"	// |
#endif
// Named constants with ARSS scope
#define MAX_KEYBOARD_KEYS 256
#define MAX_ACTORS_PER_SCENE 16284
#define MAX_SHAPES_PER_ACTOR 64
#define MAX_CHARS_PER_NAME FILENAME_MAX
#define RNG_SEED 23
enum RubbleGrainType {eBAD_RUBBLE_GRAIN_TYPE,eBOX_GRAIN,eSPHERE_GRAIN,eCAPSULE_GRAIN,eCONVEX_GRAIN,ePYRAMID_GRAIN,eD12_GRAIN,eRAND_CONVEX};

// GLUT globals (window manager, event manager, rendering manager)
struct arss_glut_camera {
	PxVec3 pos, forward, right, up;
	float aspectRatio;
	float yFOV;
	float zBufNear;
	float zBufFar;
	float speed;
};
extern arss_glut_camera gCamera;
struct arss_glut_colors {
	GLfloat background[3];
	GLubyte defaultActor[3];
	vector< vector<GLubyte> > colorBucket;
};
extern arss_glut_colors gColors;
struct arss_glut_controls {
	bool keys[MAX_KEYBOARD_KEYS];
	bool mbut[3];
	float mxco, myco;
	float mrotspd;
};
extern arss_glut_controls gControls;
struct arss_glut_hud {
	HUD hud;
	unsigned int FPS;
	unsigned int WALL_TIME;
	unsigned int DBG_MSG;
	unsigned int PAUSED;
	unsigned int HELP;
};
extern arss_glut_hud gHUD;
// Run globals (file names, run name, misc information etc.)
struct arss_run {
	string baseName;
	string workingDirectory;
	string iniFile;
	string logFile;
	string outFile;
	string loadSceneFromFile;
	unsigned int outputFrequency;	// frames per output
	unsigned int captureFrequency;	// frames per screen grab (not yet implemented)
	unsigned int renderFrequency;	// seconds between render requests
	unsigned int physicsFrequency;	// maximum FPS, for debugging
};
extern arss_run gRun;
// DevIL globals
struct arss_devil {
	ILuint screenShot;
};
extern arss_devil gDevil;
// PhysX globals
struct arss_physx {
	struct {
		bool bAdaptiveForce;
		bool bEnergySleep;
		bool bGravity;
	} flags;
	struct {
		PxReal bounceThreshold;
		PxReal sleepThreshold;
		PxReal linearDamping;
		PxReal angularDamping;
		PxReal maxAngularVelocity;
		PxReal skinWidth;
	} props;
	struct {
		PxReal length;
		PxReal mass;
		PxReal speed;
	} scaling;
	struct {
		PxActorTypeSelectionFlags everyone;
		PxActorTypeSelectionFlags statics;
		PxActorTypeSelectionFlags dynamics;
	} roles;
	PxFoundation*	mFoundation;
	PxPhysics*		mPhysics;
	PxCooking*		mCooking;
	PxScene*		mScene;
	PxActor*		cast[MAX_ACTORS_PER_SCENE];
	PxMaterial*		mDefaultMaterial;
	PVD::PvdConnection*	mPVDConnection;
};
extern arss_physx gPhysX;
// CUDA globals
struct arss_cuda {
	bool			cudaCapable;
#ifdef HAVE_CUDA_TK
	cudaDeviceProp	gpuProp;
#endif
};
extern arss_cuda gCUDA;
// Simulation globals
struct arss_simulation {
	PxReal			timeStep;
	bool			isRunning;
	unsigned int	frame;
	PxReal			codeTime;
	unsigned int	wallTime;
	PxReal			fps;
	bool			bPause;
};
extern arss_simulation gSim;
// Experiment globals common to all experiments (projects)
struct arss_experiment {
	unsigned int		rubbleCount;
	bool				bUniformRubble;
	bool				bManualControl;
	vector<PxMaterial*>	materials;
	RubbleGrainType		defGrainType;
	PxReal				defGrainSize;
	PxReal				defGrainSizeScatter;
	PxReal				defGrainDensity;
	struct {
		PxActor*	groundPlane;
		PxActor*	selectedActor;
		struct {
			PxRigidDynamic* leftmost;
			PxRigidDynamic* rightmost;
			PxRigidDynamic* upmost;
			PxRigidDynamic* downmost;
			PxRigidDynamic* inmost;
			PxRigidDynamic* outmost;
		} extremers;
	} VIPs;
	struct {
		PxVec3	systemCM;
		PxVec3	systemAM;
		PxVec3	systemLM;
		PxReal	systemKE;
		PxReal	systemMass;
	} IOMs; // integrals of motion
};
extern arss_experiment gExp;
// Debugging globals
struct arss_debug {
	bool	bBoxGridOn;
	bool	bXZGridOn;
	bool	bXYGridOn;
};
extern arss_debug gDebug;
// Utility globals
struct arss_utils {
	nr3::Ran* uniformRNG;
	nr3::Normaldev* normalRNG;
	arss_utils() {uniformRNG = new nr3::Ran(RNG_SEED); normalRNG = new nr3::Normaldev(0.0f,1.0f,RNG_SEED);}
};
extern arss_utils gUtils;

// Initialization functions
bool InitRun(int,char**);
bool InitGlut(int,char**);
bool InitDevIL();
bool InitHUD();
bool InitPhysX();
bool InitPVD(PxPhysics*);
bool InitCUDA(int devToUse=0);
bool InitExperiment();

// Cleanup functions
void DestroyPhysX();
void ExitCallback();

// GLUT callbacks
void RenderScene();
void ReshapeCallback(int,int);
void IdleCallback();
void ProcessAsciiKeys(unsigned char,int,int);
void ProcessSpecialKeys(int,int,int);
void ProcessKeyRelease(unsigned char,int,int);
void MouseCallback(int,int,int,int);
void MotionCallback(int,int);
void ProcessCameraControls();

// All project functions
void AdvanceSimulation(PxReal);
void RefreshHUD();
void DisplayHelp();
void CaptureScreen();
bool ConfigARSSOptions();
void UpdateIntegralsOfMotion();
void CreateGroundPlane();
void LogRunStats();
PxU32 CountSleepers();
void Reveille();
void DeadStop();
void ReCenterActors();
void FindExtremers();
bool SaveSceneToRepXDump();
bool LoadSceneFromFile(string sceneFile, bool dynamicsOnly=false);
PxRigidDynamic* CreateRubbleGrain(PxVec3 startPosition, RubbleGrainType grainType, PxReal sizeScale, PxMaterial& material, PxReal density=1, bool isUniform=true, PxReal sizeScatter=1.0, vector<PxVec3> *verts=NULL);
PxRigidDynamic* CreateDynamicSphericalGrain(PxVec3 position, PxReal radius, PxMaterial& material, PxReal density, bool isUniform, PxReal radiusScatter);
PxRigidDynamic* CreateDynamicBoxGrain(PxVec3 position, PxReal side, PxMaterial& material, PxReal density, bool isUniform, PxReal sideScatter);
PxRigidDynamic* CreateDynamicCapsuleGrain(PxVec3 position, PxReal radius, PxMaterial& material, PxReal density, bool isUniform, PxReal radiusScatter);
PxRigidDynamic* CreateDynamicConvexGrain(PxVec3 position, vector<PxVec3> verts, PxMaterial& material, PxReal density);
PxRigidDynamic* CreateDynamicConvexGrain(PxVec3 position, PxConvexMesh* mesh, PxMaterial& material, PxReal density);
void RandOrientActor(PxRigidActor*);
void RandLaunchActor(PxRigidDynamic*,PxReal);
void RandSpinActor(PxRigidDynamic*,PxReal);
vector<PxVec3> MakePyramidVertexList(PxReal sizeScale=1.0f, bool isUniform=true, PxReal sideScatter=1.0f);
vector<PxVec3> MakeD12VertexList(PxReal sizeScale=1.0f);
vector<PxVec3> MakeRandomVertexList(PxReal sizeScale=1.0f, PxU32 numVerts=12, PxReal gama=10);
PxConvexMesh* MakePxMeshFromVertexList(vector<PxVec3> verts);

// Experiment specific functions (declared in ARSS.h but defined differently in each project's source)
void CustomizeRun(int,char**);
bool ConfigExperimentOptions();
void CustomizeScene(PxSceneDesc&);
void CustomizeGLUT();
void CustomizeHUD();
void RefreshCustomHUDElements();
void FireAction();
void PrintDebug();
void ApplyCustomInteractions();
void RenderOtherStuff();
void CreateExperiment();
void RebootExperiment();
void LoadExperiment();
void LogExperiment();
void UpArrowAction();
void DownArrowAction();
void LeftArrowAction();
void RightArrowAction();

// Rendering functions
void RenderActors();
void DrawActor(PxRigidActor* actor);
void ColorActor(PxActor*,const GLubyte*);

// Shape drawing routine defined in Drawers.cpp
//void DrawShape(PxShape*);

#endif // include guards