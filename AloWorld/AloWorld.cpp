///////////////////////////////////////////////////////////////////////////////
// This project exists to test the proper "installation" of the ARSS package,
// testing compilation and linking with the external libraries. If everything
// is installed correctly, you should get a window showing ????????????????????
// and a HUD displaying wall time and Frames-per-second counter (which should be
// very close to your monitor's refresh rate if v-sync is enabled for your video
// card. In the terminal window, a message is displayed showing some of your system's
// environment settings, and a separate message enumerating your GPU's capabilities.
///////////////////////////////////////////////////////////////////////////////

#include "ARSS.h" // all supporting global, non project-specific entities

// Request lint level warnings with the LINT macro on Microsoft compilers
#ifdef LINT
#pragma warning(push,4)
#pragma warning(disable:4100) // unreferenced formal parameter
#endif

PxRigidDynamic* gBall;

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

// Experiment specific functions called from ARSS.cpp
void CreateExperiment()
{
	gDebug.bXZGridOn=true;

	gExp.VIPs.groundPlane = PxCreatePlane(*gPhysX.mPhysics,PxPlane(PxVec3(0,1,0),0),*gPhysX.mDefaultMaterial);
	if (!gExp.VIPs.groundPlane)
		ncc__error("ground plane creation failed!");
	gExp.VIPs.groundPlane->setName("ground");
	gPhysX.mScene->addActor(*gExp.VIPs.groundPlane);

	PxRigidStatic* basket = gPhysX.mPhysics->createRigidStatic(PxTransform(PxVec3(0)));
	if (!basket)
		ncc__error("basket fail!");
	PxMaterial* defmat=gPhysX.mDefaultMaterial;
	PxBoxGeometry base(0.5,0.1,0.5);
	PxBoxGeometry pole(0.05,1,0.05);
	PxBoxGeometry board(0.5,0.5,0.01);
	PxBoxGeometry hoopel(0.01,0.01,0.15);
	basket->createShape(base,*defmat,PxTransform(PxVec3(0,0.1,0)));
	basket->createShape(pole,*defmat,PxTransform(PxVec3(0,1,0)));
	PxShape* sboard = basket->createShape(board,*defmat,PxTransform(PxVec3(0,2,0.05)));
	PxShape* shoopel1 = basket->createShape(hoopel,*defmat,PxTransform(PxVec3(-0.15,2,0.15+0.06)));
	PxShape* shoopel2 = basket->createShape(hoopel,*defmat,PxTransform(PxVec3(+0.15,2,0.15+0.06)));
	PxShape* shoopel3 = basket->createShape(hoopel,*defmat,PxTransform(PxVec3(+0.00,2,0.30+0.05),PxQuat(PxPi/2,PxVec3(0,1,0))));

	gPhysX.mScene->addActor(*basket);
	
	gColors.colorBucket.push_back(vector<GLubyte>(3));
	gColors.colorBucket.back()[0]=ncc::rgb::yLightYellow[0];
	gColors.colorBucket.back()[1]=ncc::rgb::yLightYellow[1];
	gColors.colorBucket.back()[2]=ncc::rgb::yLightYellow[2];
	sboard->userData=&(gColors.colorBucket.back()[0]);

	gColors.colorBucket.push_back(vector<GLubyte>(3));
	gColors.colorBucket.back()[0]=ncc::rgb::grBlack[0];
	gColors.colorBucket.back()[1]=ncc::rgb::grBlack[1];
	gColors.colorBucket.back()[2]=ncc::rgb::grBlack[2];
	shoopel1->userData=&(gColors.colorBucket.back()[0]);
	shoopel2->userData=&(gColors.colorBucket.back()[0]);
	shoopel3->userData=&(gColors.colorBucket.back()[0]);

	gSim.isRunning = true;
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

}
void CustomizeScene(PxSceneDesc &sceneDesc)
{
	sceneDesc.gravity=PxVec3(0.0,-2.8,0.0);
}
void CustomizeRun(int argc, char** argv)
{
	gSim.isRunning=true;
	gSim.timeStep=1.0f/60.0f;
}
void CustomizeGLUT()
{
	gCamera.pos.z=4.0;
	gCamera.pos.y=1.0;
}
void CustomizeHUD()
{
	
}
void RefreshCustomHUDElements()
{

}
void FireAction()
{
	PxRigidDynamic* ball = CreateRubbleGrain(gCamera.pos+PxVec3(0,0,-0.3),gExp.defGrainType,gExp.defGrainSize,*gPhysX.mDefaultMaterial,gExp.defGrainDensity);
	if (ball)
	{
		ball->setLinearVelocity(gCamera.forward*2.5);
		gColors.colorBucket.push_back(vector<GLubyte>(3));
		gColors.colorBucket.back()[0]=ncc::rgb::oDarkOrange[0];
		gColors.colorBucket.back()[1]=ncc::rgb::oDarkOrange[1];
		gColors.colorBucket.back()[2]=ncc::rgb::oDarkOrange[2];
		ball->userData=&(gColors.colorBucket.back()[0]);
	}
	gBall=ball;

}
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
void PrintDebug()
{
	
}
void ApplyCustomInteractions()
{

}
void RenderOtherStuff()
{
	//if (gBall) DrawActorAxes(gBall,gExp.defGrainSize*2);
}

// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif