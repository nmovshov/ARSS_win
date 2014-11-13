//////////////////////////////////////////////////////////////////////////////////////////
// Source file for project Verify. This project implements a suite of tests of rigid body
// dynamics in PhysX. Several scenarios can be simulated, that have either an analytic
// solution or at least a conserved quantity that can be measured. The available tests
// are:
//
// * Slider - a box on an inclined plane, testing the friction model in PhysX
// 
// * Collider - a collision between shapes of varying complexity, testing the collision
//              detection and collision resolution in PhysX.
//
// * Ball on ground - test of internal gravity in PhysX
//
// * Ball on ball - test of integration of manually implemented gravity
//
// Author: Naor Movshovits (nmovshov at google dot com)
//////////////////////////////////////////////////////////////////////////////////////////

#include "ARSS.h" // all supporting global, non project-specific entities
#include "Verify.h" // global project-specific entities

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
	switch (verify::eExperimentType)
	{
	case verify::eCOLLIDER:
		verify::CreateColliderExperiment();	
		break;
	case verify::eTUMBLER:
		verify::CreateTumblerExperiment();
		break;
	case verify::eTUMBLERS:
		verify::CreateTumblersExperiment();
		break;
	case verify::eINCLINER:
		verify::CreateInclinerExperiment();
		break;
	case verify::eBALL_ON_GROUND:
		verify::CreateBallOnGroundExperiment();
		break;
	case verify::eBAD_EXPERIMENT_TYPE: // intentional fall through!
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
void LogExperiment()
{
	switch (verify::eExperimentType)
	{
	case verify::eTUMBLER:
	case verify::eTUMBLERS:
		verify::LogTumblerExperiment();
		break;
	case verify::eBALL_ON_GROUND:
		verify::LogBallOnGroundExperiment();
		break;
	case verify::eBAD_EXPERIMENT_TYPE:
		ncc__warning("Unknown experiment type. Nothing logged.");
		break;
	default:
		break;
	}
}
void ControlExperiment()
{
	switch (verify::eExperimentType)
	{
	case verify::eBALL_ON_GROUND:
		verify::ControlBallOnGroundExperiment();
		break;
	case verify::eBAD_EXPERIMENT_TYPE:
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
	if      (strcmp(buf,"collider")==0)
		verify::eExperimentType=verify::eCOLLIDER;
	else if (strcmp(buf,"incliner")==0)
		verify::eExperimentType=verify::eINCLINER;
	else if (strcmp(buf,"tumbler")==0)
		verify::eExperimentType=verify::eTUMBLER;
	else if (strcmp(buf,"tumblers")==0)
		verify::eExperimentType=verify::eTUMBLERS;
	else if (strcmp(buf,"ball_on_ground")==0)
		verify::eExperimentType=verify::eBALL_ON_GROUND;
	else if (strcmp(buf,"ball_on_ball")==0)
		verify::eExperimentType=verify::eBALL_ON_BALL;
	else
		verify::eExperimentType=verify::eBAD_EXPERIMENT_TYPE;

	// Collider experiment parameters
	ncc::GetStrPropertyFromINIFile("experiment","spin_magnitude","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	verify::spinMag = atof(buf);
	ncc::GetStrPropertyFromINIFile("experiment","kick_magnitude","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	verify::kickMag = atof(buf);

	// Physical parameters
	ncc::GetStrPropertyFromINIFile("units","little_g","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	verify::units.littleG = atof(buf);
	ncc::GetStrPropertyFromINIFile("units","big_g","0",buf,MAX_CHARS_PER_NAME,gRun.iniFile.c_str());
	verify::units.bigG = atof(buf);
	

	return true;
}
void CustomizeGLUT()
{
	gCamera.pos.z=6;
}
void CustomizeHUD()
{
	verify::hudMsgs.systemDiag1 = gHUD.hud.AddElement("",0.9,0.1);
	verify::hudMsgs.systemDiag2 = gHUD.hud.AddElement("",0.9,0.2);
	verify::hudMsgs.systemDiag3 = gHUD.hud.AddElement("",0.9,0.3);
	verify::hudMsgs.actorDiag   = gHUD.hud.AddElement("",0.9,0.3);
}
void RefreshCustomHUDElements()
{
	char buf[MAX_CHARS_PER_NAME];
	int	  ch2px	    = 18; // hud uses 18pt font
	float px2width  = 1.0/glutGet(GLUT_WINDOW_WIDTH);
	float scrPos = glutGet(GLUT_WINDOW_WIDTH);

	switch (verify::eExperimentType)
	{
	case verify::eCOLLIDER:
		UpdateIntegralsOfMotion();

		// System diagnostic 1 - system linear momentum
		sprintf_s(buf,MAX_CHARS_PER_NAME,"System LM = (%+04.6f,%+04.6f,%+04.6f)",gExp.IOMs.systemLM.x,gExp.IOMs.systemLM.y,gExp.IOMs.systemLM.z);
		scrPos = 1.0 - strlen(buf)*ch2px*px2width*0.6;
		gHUD.hud.SetElement(verify::hudMsgs.systemDiag1,buf,scrPos,0.1);

		// System diagnostic 2 - system angular momentum
		sprintf_s(buf,MAX_CHARS_PER_NAME,"System AM = (%+04.6f,%+04.6f,%+04.6f)",gExp.IOMs.systemAM.x,gExp.IOMs.systemAM.y,gExp.IOMs.systemAM.z);
		scrPos = 1.0 - strlen(buf)*ch2px*px2width*0.6;
		gHUD.hud.SetElement(verify::hudMsgs.systemDiag2,buf,scrPos,0.2);

		// System diagnostic 3 - system kinetic energy
		sprintf_s(buf,MAX_CHARS_PER_NAME,"System KE = (%+04.6f)",gExp.IOMs.systemKE);
		gHUD.hud.SetElement(verify::hudMsgs.systemDiag3,buf,scrPos,0.3);
	
		// Selected actor diagnostics
		break;

	case verify::eTUMBLER: // intentional fall-through
	case verify::eTUMBLERS:
		if (verify::VIPs.tumbler)
		{
			// System diagnostic 1 - the tumbler's angular momentum (manual)
			PxVec3 L = verify::tumbler.L_now;
			sprintf_s(buf,MAX_CHARS_PER_NAME,"L = (%+04.6f,%+04.6f,%+04.6f)",L.x,L.y,L.z);
			scrPos = 1.0 - strlen(buf)*ch2px*px2width*0.6;
			gHUD.hud.SetElement(verify::hudMsgs.systemDiag1,buf,scrPos,0.1);

			// System diagnostic 2 - the tumbler's angular momentum relative discrepancy
			PxVec3 L0 = verify::tumbler.L_true;
			PxReal L_err = (L-L0).magnitude() / L0.magnitude() * 100 * PxSign(L.dot(L0));
			if (L_err == HUGE_VAL) L_err = 0;
			sprintf_s(buf,MAX_CHARS_PER_NAME,"DL = %g%%",L_err);
			gHUD.hud.SetElement(verify::hudMsgs.systemDiag2,buf,scrPos,0.2);

			// More actor diagnostics

		}
		break;

	case verify::eINCLINER:
		// Current inclination angle
		sprintf_s(buf,MAX_CHARS_PER_NAME,"Inclination = %d deg",verify::incliner.inc);
		scrPos = 0.8;
		gHUD.hud.SetElement(verify::hudMsgs.systemDiag1,buf,scrPos,0.1);
		break;
	}
}
void FireAction()
{
	static int k=0;
	static int j=0;
	if (verify::eExperimentType == verify::eTUMBLER)
	{
		int idx[10] = {-10,-8,-6,-4,-2,2,4,6,8,10};
		if (j<10 && k<10)
		{
			PxVec3 pos(idx[j],idx[k],0);
			PxReal randSize = gExp.defGrainSize * (1 + (!gExp.bUniformRubble)*gExp.defGrainSizeScatter*(-0.5+gUtils.uniformRNG->doub()));
			PxRigidDynamic* tumbler = CreateRubbleGrain(pos,gExp.defGrainType,randSize,*gPhysX.mDefaultMaterial,gExp.defGrainDensity);
			if (tumbler)
			{
				tumbler->setName("tumbler");
				RandOrientActor(tumbler);
				RandSpinActor(tumbler,verify::spinMag*(1+gUtils.uniformRNG->doub()));
			}
			if (k<10) k++;
			if (k==10) {k=0; j++;}
		}
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
		switch (verify::eExperimentType)
		{
		case verify::eTUMBLER:
			verify::CalcTumblerDynamics();
			break;
		case verify::eTUMBLERS:
			UpdateIntegralsOfMotion();
			verify::tumbler.L_now = gExp.IOMs.systemAM;
			break;
		}
	}
}
void RenderOtherStuff()
{
	switch (verify::eExperimentType)
	{
	case verify::eTUMBLER:
		DrawArrow(PxVec3(0),verify::tumbler.L_now*10 ,0.12,4.0,ncc::rgb::rDarkRed);
		DrawArrow(PxVec3(0),verify::tumbler.L_true*10,0.12,4.0,ncc::rgb::gDarkGreen);
		break;
	}
}
void UpArrowAction()
{
	switch (verify::eExperimentType)
	{
	case verify::eINCLINER:
		verify::incliner.inc++;
		if (verify::incliner.inc>80) verify::incliner.inc=80;
		verify::InclineGravity(verify::incliner.inc);
		Reveille();
		break;
	}
}
void DownArrowAction()
{
	switch (verify::eExperimentType)
	{
	case verify::eINCLINER:
		verify::incliner.inc--;
		if (verify::incliner.inc<0) verify::incliner.inc=0;
		verify::InclineGravity(verify::incliner.inc);
		Reveille();
		break;
	}
}
void LeftArrowAction()
{

}
void RightArrowAction()
{

}

// Project namespace functions
void verify::CreateColliderExperiment()
/*
Create a spinning but stationary target and a moving projectile and observe the
collision. Check conservation of the integrals of motion.
*/
{
	PxRigidDynamic* target;
	PxRigidDynamic* bullet;
	gDebug.bXYGridOn=true;
	gPhysX.mScene->setGravity(PxVec3(0));
	
	bullet = CreateRubbleGrain(PxVec3(+2*gExp.defGrainSize,0,0),gExp.defGrainType,gExp.defGrainSize,*(gPhysX.mDefaultMaterial),gExp.defGrainDensity,gExp.bUniformRubble);
	target = CreateRubbleGrain(PxVec3(-2*gExp.defGrainSize,0,0),gExp.defGrainType,gExp.defGrainSize,*(gPhysX.mDefaultMaterial),gExp.defGrainDensity,gExp.bUniformRubble);

	if (bullet)
	{
		bullet->setName("bullet");
		bullet->setLinearVelocity(PxVec3(-verify::kickMag*gExp.defGrainSize,0,0));
		RandOrientActor(bullet);
	}

	if (target)
	{
		target->setName("target");
		RandOrientActor(target);
		RandSpinActor(target,verify::spinMag);
	}

	verify::VIPs.bullet = bullet;
	verify::VIPs.target = target;
}
void verify::CreateTumblerExperiment()
/*
Create a dynamic free tumbler and watch it tumble, hoping for a changing angular velocity
vector but stationary angular momentum vector. EDIT: this test proves that PhysX (<=3.2.0)
does not integrate Euler's equations of rigid body motion. Therefore the angular momentum
is not constant. Use this test to constrain the angular momentum error.
*/
{
	PxRigidDynamic* tumbler = NULL;
	gDebug.bXYGridOn=true;
	gPhysX.mScene->setGravity(PxVec3(0));

	PxReal randSize = gExp.defGrainSize * (1 + (!gExp.bUniformRubble)*gExp.defGrainSizeScatter*(-0.5+gUtils.uniformRNG->doub()));
	tumbler = CreateRubbleGrain(PxVec3(0),gExp.defGrainType,randSize,*(gPhysX.mDefaultMaterial),gExp.defGrainDensity);
	if (tumbler)
	{
		tumbler->setName("tumbler");
		RandOrientActor(tumbler);
		RandSpinActor(tumbler,verify::spinMag);
		//tumbler->setAngularVelocity(PxVec3(0,1,0));
	}

	verify::VIPs.tumbler = tumbler;
	verify::tumbler.handle = tumbler;
	CalcTumblerDynamics();
	verify::tumbler.L_true = verify::tumbler.L_now;
}
void verify::CreateTumblersExperiment()
/*
Create 100 dynamic free tumblers and record the system's total angular momentum.
*/
{
	gDebug.bXYGridOn=true;
	gPhysX.mScene->setGravity(PxVec3(0));
	int idx[10] = {-10,-8,-6,-4,-2,2,4,6,8,10};
	for (int j=0; j<10; j++) {
		for (int k=0; k<10; k++)
		{
			PxVec3 pos(idx[j],idx[k],0);
			PxReal randSize = gExp.defGrainSize * (1 + (!gExp.bUniformRubble)*gExp.defGrainSizeScatter*(-0.5+gUtils.uniformRNG->doub()));
			PxRigidDynamic* tumbler = CreateRubbleGrain(pos,gExp.defGrainType,randSize,*gPhysX.mDefaultMaterial,gExp.defGrainDensity);
			if (tumbler)
			{
				tumbler->setName("tumbler");
				RandOrientActor(tumbler);
				RandSpinActor(tumbler,verify::spinMag*(1+gUtils.uniformRNG->doub()));
			}

			UpdateIntegralsOfMotion();
			verify::tumbler.L_true = gExp.IOMs.systemAM;
			verify::VIPs.tumbler = tumbler;
		}
	}

	gCamera.pos.z = 20;
}
void verify::CreateInclinerExperiment()
/*
Create a box on a the ground plane. Tilt the gravity vector with the up/down arrows.
*/
{
	CreateGroundPlane();
	gCamera.pos.z=4.0;
	gCamera.pos.y=1.0;
	gDebug.bXZGridOn=true;
	gPhysX.mScene->setGravity(PxVec3(0,-10,0));

	PxRigidDynamic* incliner = NULL;
	incliner = CreateRubbleGrain(PxVec3(0,0.52*gExp.defGrainSize,0),eBOX_GRAIN,gExp.defGrainSize,*gPhysX.mDefaultMaterial);
	if (incliner)
	{
		incliner->setName("incliner");
	}
	verify::VIPs.incliner = incliner;

}
void verify::CalcTumblerDynamics()
{
	if (verify::VIPs.tumbler)
	{
		// Start by obtaining the tumbler's pose change during the time step
		static	PxTransform P0 = verify::VIPs.tumbler->getGlobalPose(); // first pose
				PxTransform P1 = verify::VIPs.tumbler->getGlobalPose(); // current pose

		// First, we calculate the angular velocity in global frame
		PxQuat Q0 = P0.q; // quaternion of initial pose
		PxQuat Q1 = P1.q; // quaternion of final pose
		PxQuat dQ = Q1 * Q0.getConjugate(); // the "difference" quaternion
		PxVec3 w;
		PxReal teta;
		dQ.toRadiansAndUnitAxis(teta,w);
		w = w * (teta/gSim.timeStep);
		// NOTE: uncomment to override with the SDK returned value for w
		w = verify::VIPs.tumbler->getAngularVelocity();
		verify::tumbler.w = w;

		// Next, we calculate the body frame inertia tensor resolved in global frame
		PxMat33 I = PxMat33::createDiagonal(verify::VIPs.tumbler->getMassSpaceInertiaTensor());
		PxMat33 S(verify::tumbler.handle->getGlobalPose().transform(verify::tumbler.handle->getCMassLocalPose()).q);
		PxMat33 Sm1 = S.getTranspose();
		I = S * I * Sm1;
		
		// And now we can calculate the global frame angular momentum and kinetic energy of rotation
		verify::tumbler.L_now = I * w;

		// Save current pose for next time
		P0 = P1;
	}
}
void verify::InclineGravity( PxReal deg )
{
	deg *= PxPi/180.f;
	PxVec3 g(10*sin(deg),-10*cos(deg),0);
	gPhysX.mScene->setGravity(g);
}
void verify::CreateBallOnGroundExperiment()
/*
Drop a ball on the ground. Check integration and proper scaling.
*/
{
	// Make the ground and gravity
	CreateGroundPlane();
	gPhysX.mScene->setGravity(PxVec3(0,verify::units.littleG,0));

	// Put a ball in the air
	verify::VIPs.ball1 = CreateRubbleGrain(PxVec3(0,6,0),eSPHERE_GRAIN,1,*gPhysX.mDefaultMaterial);

	// Move the camera to a better vantage point and turn on a grid
	gCamera.pos = PxVec3(0,2,6);
	gDebug.bXZGridOn = true;

	// Start a log
	if (gRun.outputFrequency)
	{
		ostringstream header;
		header << "# This is the run log of " << gRun.baseName << endl;
		header << "# Time step used = " << gSim.timeStep << endl;
		header << "# Columns are (values in code units):" << endl;
		header << "# [t]    [y]    [v]" << endl;
		ofstream fbuf(gRun.outFile.c_str(),ios::trunc);
		if (!fbuf.is_open())
			ncc__error("Could not start a log. Experiment aborted.\a\n");
		fbuf << header.str() << endl;
	}

	// Start the action
	gSim.isRunning = true;
	gSim.bPause = false;

}
void verify::LogTumblerExperiment()
{
	static ofstream fp(gRun.outFile.c_str(),ios::trunc);
	if (fp.good())
	{
		fp.width(5);
		fp << gSim.frame << "\t";
		if (verify::eExperimentType == verify::eTUMBLER || verify::eExperimentType == verify::eTUMBLERS)
		{
			UpdateIntegralsOfMotion();
			fp << setiosflags(ios::fixed);
			fp << gExp.IOMs.systemAM.x << "\t" << gExp.IOMs.systemAM.y << "\t" << gExp.IOMs.systemAM.z << "\t" << gExp.IOMs.systemMass << endl;
		}
	}
	else
	{
		ncc__error("Could not access output file. Experiment aborted.\a");
	}
}
void verify::LogBallOnGroundExperiment()
{
	// Let's do this without worrying about minimizing access to disk - it will be so much easier!
	ofstream fbuf(gRun.outFile.c_str(),ios::app);
	if (!fbuf.is_open()) {
		ncc__warning("Could not open log. Nothing written!\a\n");
		return;
	}

	// Collect
	PxReal t = gSim.codeTime;
	PxReal y = verify::VIPs.ball1->getGlobalPose().p.y;
	PxReal v = verify::VIPs.ball1->getLinearVelocity().y;

	// Write
	fbuf.setf(ios::fixed);
	fbuf << setw(8) << t << "    " << setw(8) << y << "    " << setw(8) << v << endl;
	fbuf.close();
}
void verify::ControlBallOnGroundExperiment()
{
	// Check for stop condition
	PxReal y = verify::VIPs.ball1->getGlobalPose().p.y;
	PxReal d = gPhysX.scaling.length*0.04;
	if ((!gSim.targetTime && ((y - 1) <= d)) || (gSim.targetTime && (gSim.codeTime >= gSim.targetTime - gSim.timeStep)))
		gSim.isRunning = false;
}

// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif