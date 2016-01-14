// ===============================================================================
//  a set of experiments with tidal encounters, using GPU kernels for gravity
// ===============================================================================

#include "nshared.h"
#include "cuShared.h"
#include "cuTidal.h"

// Experiment globals
enum ExperimentType		gExperimentType;
enum GravityType		gRavityType;
NxReal gT_tot, gT_tot0; // total scene kinetic energy
NxVec3 gL_tot, gL_tot0; // total scene angular momentum
NxVec3 gP_tot;			// total scene linear momentum
NxVec3 gSystemCM, gSystemCM0;
NxReal gTidingMass;
NxReal gGravSoftFactor;
NxVec3 gTidingCenter(0.0);
NxReal *gOrbit=NULL;
unsigned int gOrbitLength; // b/c I'm too lazy to change gOrbit into a vector

void CreateExperiment()
{
	switch (gExperimentType)
	{
	case TIDAL_ENCOUNTER:
		{
			gScene->releaseActor(*groundPlane); // get rid of ground plane
			CreateInitialPile(); // set initial conditions
#ifdef USE_CUDA_GRAV
			AllocateCUDAGlobals(RUBBLE_SIZE);
#endif
			ReadOrbit();

			// initialize debugging data
			UpdateConservables();
			gL_tot0=gL_tot;
			gSystemCM=LocateSystemCM();
			gSystemCM0=gSystemCM;
			
			isRunning=true;
			if (!isManualControl) {bDebugWireframeMode=false; bPause=false;}
			break;
		}
	case FREE_PILE:
		{
			gScene->releaseActor(*groundPlane); // get rid of ground plane
			CreateInitialPile();
#ifdef USE_CUDA_GRAV
			AllocateCUDAGlobals(RUBBLE_SIZE);
#endif

			//initialize debugging data
			UpdateConservables();
			gL_tot0=gL_tot;
			gSystemCM=LocateSystemCM();
			gSystemCM0=gSystemCM;

			isRunning=true;
			if (!isManualControl) {bDebugWireframeMode=false; bPause=false;}
			break;
		}

	default:				
		{printf("Error: unimplemented experiment type\a\n"); isRunning=false;	break;}
	}
	if (gColorScatter) ActorColorsScatter(RUBBLE_SIZE);
}
bool ReadOrbit()
{
	// open orbit file
	char orbfile[255];
	sprintf(orbfile,"%s\\%s\\%s.orb",_getcwd(NULL,0),gRunBaseName,gRunBaseName);
	ifstream fp(orbfile);
	if (fp.fail()){cout<<"ERROR:could not open orbit file\a"<<endl; return false;}
	
	// peek into the .orb file to find header lines.
	char line[255];
	unsigned int headlines=0;
	while (fp.good()){
		fp.getline(line,255);
		headlines++;
		if (strcmp(line,"---END HEADER---")==0) break;
	}
	if (fp.eof()) headlines=0; // explicit header/data separator not found. assume no header.
	fp.close();

	vector<vector<float>*>* temporb=RSReadMatrixFromFile(orbfile,headlines);
	// verify orbit time vector vs. gDeltaTime
	for (int k=0; k<temporb->size()-1; k++) {
		if (NxMath::abs(temporb->at(k+1)->at(0)-temporb->at(k)->at(0)-gDeltaTime) > 1e-2) nerror("delta tee mismatch between .ini (specified) and .orb (implied)"); //omg single precis is so bad!
	}
	// we have a temporary orbit in a <vector> of <vector>s, but i'd rather have a simple linear float array.
	gOrbit = new (nothrow) NxReal[temporb->size()*3];
	gOrbitLength=temporb->size();
	if (gOrbit==NULL) nerror("failed to allocate memory\a");
	// gOrbit will be filled with t,x,y,t,x,y,t,x,y...
	for (int k=0; k<temporb->size(); k++) {
		gOrbit[3*k+0]=temporb->at(k)->at(0);
		gOrbit[3*k+1]=temporb->at(k)->at(1);
		gOrbit[3*k+2]=temporb->at(k)->at(2);
	}
	delete temporb;

	return true;
}
void InitializeHUD()
{
	bHardwareScene = (gScene->getSimType() == NX_SIMULATION_HW);

	hud.Clear();

	char s[255];

	// Add hardware/software to HUD - string id 0
	if (bVerboseHUD)
	{
		if (bHardwareScene)
			hud.AddDisplayString("Hardware Scene", 0.84f, 0.96f);
		else
			hud.AddDisplayString("Software Scene", 0.84f, 0.96f);
		if (bRealTime)
			hud.SetDisplayColor(0,1,0,0,1);
	} 
	else	hud.AddDisplayString("",0.84f,0.96f);

	// Add pause to HUD - string id 1
	if (bPause)  
		hud.AddDisplayString("Paused - Hit \"p\" to Unpause", 0.3f, 0.55f);
	else
		hud.AddDisplayString("", 0.0f, 0.0f);

	// Add experiment title - string id 2
	if (bVerboseHUD)	sprintf(s,"Scene: %s",gRunBaseName);
	else				sprintf(s,"");
	hud.AddDisplayString(s,0.01,0.96);
	hud.SetDisplayColor(2,0,0,1,1);

	// Add NbActors - string id 3
	if (bVerboseHUD)	sprintf(s,"%s%d","Actors: ",gScene->getNbDynamicShapes());
	else				sprintf(s,"");
	hud.AddDisplayString(s,0.01,0.87);

	// Add rubble type - string id 4
	if (bVerboseHUD)		sprintf(s,"Rubble type: %d",gRubbleType);
	else					sprintf(s,"");
	hud.AddDisplayString(s,0.01,0.93);

	// Add Experiment type - string id 5
	if (bVerboseHUD)			sprintf(s,"Experiment Type: %d",gExperimentType);
	else						sprintf(s,"");
	hud.AddDisplayString(s,0.01,0.90);

	// Add time - string id 6
	if (bVerboseHUD)
	{
		sprintf(s,"Wall time: %d Code time: %d",0,0);
		hud.AddDisplayString(s,0.01,0.84);
	} 
	else
	{
		sprintf(s,"t = %d",0);
		hud.AddDisplayString(s,0.01,0.96);
		hud.SetDisplayColor(6,0.0f,0.0f,0.0f,1.0);
	}

	// Add FPS - string id 7
	if (bVerboseHUD)	sprintf(s,"FPS: %0.1f(%0.1f)",0.0,0.0);
	else				sprintf(s,"");
	hud.AddDisplayString(s,0.01,0.81);

	// Add selected actor name - string id 8
	if (bVerboseHUD)
	{
		if (gSelectedActor)
			sprintf(s,"Selected actor: %s",gSelectedActor->getName());
		else
			sprintf(s,"Selected actor: NULL");
	} 
	else	sprintf(s,"");
	hud.AddDisplayString(s,0.98-strlen(s)/100.0f,0.93);

	// Add selected actor specs - string id 9
	sprintf(s,"");
	hud.AddDisplayString(s,0.98-strlen(s)/100.0f,0.90);

	// Add scene diagnostics - string id 10
	sprintf(s,"");
	hud.AddDisplayString(s,0.98,0.87);
	if (!bVerboseHUD)
		hud.SetDisplayColor(10,0.0f,0.0f,0.0f,1.0);

}
void UpdateHUD()
{
	char s[255];

	// update scene diagnostics
	UpdateConservables();
	gSystemCM=LocateSystemCM();
	unsigned int step=NxMath::max(gFrameCounter-1,0); // -1 b/c this is usually called from main loop after frame counter is advanced
	if (gExperimentType==TIDAL_ENCOUNTER){
	float whereami=0,whereamiR=0;
	if (step>=0 && step<gOrbitLength) {
		whereami=NxMath::sqrt(pow(gOrbit[3*step+1],2)+pow(gOrbit[3*step+2],2));
		whereamiR=whereami/(1.51*NxMath::pow(gTidingMass/gDefaultDensity,(1.f/3.f))); // in Roche limit units
		sprintf(s,"R = %0.3g km (%0.3g Roche limits)  ",whereami/1000,whereamiR);
		if (bVerboseHUD)	hud.SetDisplayString(10,s,0.98-strlen(s)/100.0f,0.84);
		else				hud.SetDisplayString(10,s,0.98-strlen(s)/100.0f,0.96);
	}
	}

	// update selected actor specs
	if (gSelectedActor && gSelectedActor->isDynamic())
	{
		NxVec3 W=gSelectedActor->getAngularVelocity();
		NxVec3 L=gSelectedActor->getAngularMomentum();//resolved in body space
		gSelectedActor->getGlobalOrientation().multiplyByTranspose(L,L);//resolved in world space
		NxMat33 Imin1;
		gSelectedActor->getGlobalInertiaTensor().getInverse(Imin1);
		Imin1.multiply(L,W);
		sprintf(s,"L=(%0.3g,%0.3g,%0.3g)    ",L.x,L.y,L.z);
	}
	else
		sprintf(s,"");
	hud.SetDisplayString(9,s,0.98-strlen(s)/100.0f,0.90);

	// update selected actor name
	if (bVerboseHUD)
	{
		if (gSelectedActor)
			sprintf(s,"Selected actor: %s",gSelectedActor->getName());
		else
			sprintf(s,"Selected actor: NULL");
	} 
	else	sprintf(s,"");
	hud.SetDisplayString(8,s,0.98-strlen(s)/100.0f,0.93);

	// update hardware/software
	if (bVerboseHUD)
	{
		if (bHardwareScene)
			hud.SetDisplayString(0,"Hardware Scene", 0.84f, 0.96f);
		else
			hud.SetDisplayString(0,"Software Scene", 0.84f, 0.96f);
	}
	
	// update pause/go
	if (bPause)
		hud.SetDisplayString(1, "Paused - Hit \"p\" to Unpause", 0.3f, 0.55f);
	else
		hud.SetDisplayString(1, "", 0.0f, 0.0f);
	
	if (bVerboseHUD)
	{
		// update NbActors
		int ignoredones=0;
		sprintf(s,"%s%d(%d)","Actors: ",gScene->getNbDynamicShapes()-ignoredones,NxMath::max(NxI32(0),NxI32(CountSleepers())));
		hud.SetDisplayString(3,s,0.01,0.87);

		// update rubble type
		if (isUniformRubble)
			sprintf(s,"Rubble type: %d",gRubbleType);
		else
			sprintf(s,"Rubble type: %d",gRubbleType);
		hud.SetDisplayString(4,s,0.01,0.93);

		// update experiment type
		if (isManualControl)
			sprintf(s,"Experiment Type: %d (manual)",gExperimentType);
		else
			sprintf(s,"Experiment Type: %d (auto)",gExperimentType);
		hud.SetDisplayString(5,s,0.01,0.90);
	}

	// update time
	if (bVerboseHUD)
	{
		if (gExperimentType==TIDAL_ENCOUNTER)
			sprintf(s,"Wall time: %d Code time: %0.2g (%0.2g)",NxMath::trunc(gWallTime),gOrbit[3*gFrameCounter],gDeltaTime);
		else
			sprintf(s,"Wall time: %d Code time: %0.2g (%0.2g)",NxMath::trunc(gWallTime),gFrameCounter*gDeltaTime,gDeltaTime);
		hud.SetDisplayString(6,s,0.01,0.84);
	} 
	else
	{
		if (gExperimentType==TIDAL_ENCOUNTER)
			sprintf(s,"time = %0.2g hr",gOrbit[3*gFrameCounter]/3600);
		else
			sprintf(s,"time = %0.2g hr",gFrameCounter*gDeltaTime/3600);
		hud.SetDisplayString(6,s,0.01,0.96);
	}

	// update FPS
	if (bVerboseHUD)
	{
		sprintf(s,"FPS: %0.1f(%0.1f)",gFPS,gFrameCounter/(NxMath::max(gWallTime,1.f)));
		hud.SetDisplayString(7,s,0.01,0.81);

		// update experiment title
		if (isRunning)
			hud.SetDisplayColor(2,1,0,0,1);
		else
			hud.SetDisplayColor(2,0,1,0,1);
	}

    return;
}
void AddSpecialForces()
{
	// but not *just* forces. this is a guaranteed point of entry from the main loop (RenderCallback)
	switch (gExperimentType)
	{
	case TIDAL_ENCOUNTER: {
		if (gFrameCounter>gOrbitLength)
			StopRun("orbit completed");
		else
			Gravitate();
		break;
	}
	case FREE_PILE:
		{
			Gravitate();
			break;
		}
	default: break;
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
void Gravitate()
{
	switch (gRavityType)
	{
	case ALL_PAIRS: {GravitateByPairs(); break;} //Note: this requires N^2 looping
	default: break;
	}
}
void GravitateByPairs()
{
	NxU32 nbActors = gScene->getNbActors();
	NxActor** actors = gScene->getActors();
	NxVec3 rj, r21, force;
	NxReal mj, mk;
	NxReal softsqr = pow(gGravSoftFactor * gRainSize, 2);
#ifdef USE_CUDA_GRAV
	// For memory coherency we'll extract a linear array of body descriptors to send over the bus
	for (int j=0; j<nbActors; j++){
		rj=actors[j]->getCMassGlobalPosition();
		mj=actors[j]->getMass();
		if (actors[j]->readBodyFlag(NX_BF_KINEMATIC)) mj=0;
		g_hostPositionsArr[j].x=rj.x;
		g_hostPositionsArr[j].y=rj.y;
		g_hostPositionsArr[j].z=rj.z;
		g_hostPositionsArr[j].w=mj;
	}

	// 1. self-gravity of the pile term
	// calculate self-gravity on the device here NOTE: returns immediately, get forces back later
	GravitateOnDevice(nbActors,1U);
	
	if (gExperimentType==TIDAL_ENCOUNTER) {
		// The tide force and fictitious force can be copied from the CPU only code
		float *bodyDescripts = new float[4*nbActors]; // the small performance gain of making these global is not worth the obscurity
		float *forces        = new float[3*nbActors]; //
		for (int j=0; j<nbActors; j++){
			bodyDescripts[4*j+0]=actors[j]->getCMassGlobalPosition().x;
			bodyDescripts[4*j+1]=actors[j]->getCMassGlobalPosition().y;
			bodyDescripts[4*j+2]=actors[j]->getCMassGlobalPosition().z;
			bodyDescripts[4*j+3]=actors[j]->getMass();
			if (actors[j]->readBodyFlag(NX_BF_KINEMATIC)) bodyDescripts[4*j+3] = 0;

			forces[3*j+0]=0; forces[3*j+1]=0; forces[3*j+2]=0;

			// 2. pull from tiding body term
			unsigned int step=gFrameCounter-1; // -1 b/c Gravitate is called from main loop after frame counter is advanced
			gTidingCenter.x=-gOrbit[3*step+1]; gTidingCenter.y=-gOrbit[3*step+2]; // gTidingCenter now points from (approx) CM to the planet
			rj.x=gTidingCenter.x-bodyDescripts[4*j+0];//
			rj.y=gTidingCenter.y-bodyDescripts[4*j+1];// rj now points from body j to the planet
			rj.z=gTidingCenter.z-bodyDescripts[4*j+2];//

			rj.setMagnitude(gUniversalBigG*gTidingMass*bodyDescripts[4*j+3]/(rj.magnitudeSquared()+4*gRainSize*gRainSize));

			forces[3*j+0]+=rj.x; forces[3*j+1]+=rj.y; forces[3*j+2]+=rj.z;

			// 3. fictitious force term on the accelerated frame
			rj=-gTidingCenter; // rj now points from planet to (approx) CM
			rj.setMagnitude(gUniversalBigG*gTidingMass*bodyDescripts[4*j+3]/(rj.magnitudeSquared()+4*gRainSize*gRainSize));
			forces[3*j+0]+=rj.x; forces[3*j+1]+=rj.y; forces[3*j+2]+=rj.z;

			force.x=forces[3*j+0]; force.y=forces[3*j+1]; force.z=forces[3*j+2];
			actors[j]->addForce(force);
		}

	delete [] bodyDescripts;
	delete [] forces;
	}

	GravitateOnDevice(nbActors,2U); // get forces back from device (don't forget to add universal G)
	for (int j=0; j<nbActors; j++){
		if (actors[j]->readBodyFlag(NX_BF_KINEMATIC)) continue;
		force.x=g_hostForcesArr[j].x;
		force.y=g_hostForcesArr[j].y;
		force.z=g_hostForcesArr[j].z;
		force=force*gUniversalBigG;
		actors[j]->addForce(force);
	}

#else
	float *bodyDescripts = new float[4*nbActors]; // the small performance gain of making these global is not worth the obscurity
	float *forces        = new float[3*nbActors]; //
	for (int j=0; j<nbActors; j++){
		bodyDescripts[4*j+0]=actors[j]->getCMassGlobalPosition().x;
		bodyDescripts[4*j+1]=actors[j]->getCMassGlobalPosition().y;
		bodyDescripts[4*j+2]=actors[j]->getCMassGlobalPosition().z;
		bodyDescripts[4*j+3]=actors[j]->getMass();
		if (actors[j]->readBodyFlag(NX_BF_KINEMATIC)) bodyDescripts[4*j+3] = 0;
	}

	for (int j=0; j<nbActors; j++){
		forces[3*j+0]=0; forces[3*j+1]=0; forces[3*j+2]=0;
		// 1. self-gravity of the pile term
		for (int k=0; k<j; k++){

			r21.x = bodyDescripts[4*k+0]-bodyDescripts[4*j+0];
			r21.y = bodyDescripts[4*k+1]-bodyDescripts[4*j+1];
			r21.z = bodyDescripts[4*k+2]-bodyDescripts[4*j+2];
			
			float distSqr = r21.x*r21.x + r21.y*r21.y + r21.z*r21.z + softsqr;
			float distSix = distSqr * distSqr * distSqr;
			float invDistCube = 1.0f/sqrtf(distSix);
			
			float s = gUniversalBigG*bodyDescripts[4*j+3]*bodyDescripts[4*k+3]*invDistCube;
			
			forces[3*j+0]+=r21.x*s; forces[3*j+1]+=r21.y*s; forces[3*j+2]+=r21.z*s;
			forces[3*k+0]-=r21.x*s; forces[3*k+1]-=r21.y*s; forces[3*k+2]-=r21.z*s;
		}

		if (gExperimentType==TIDAL_ENCOUNTER) {
		// 2. pull from tiding body term
		unsigned int step=gFrameCounter-1; // -1 b/c Gravitate is called from main loop after frame counter is advanced
		gTidingCenter.x=-gOrbit[3*step+1]; gTidingCenter.y=-gOrbit[3*step+2]; // gTidingCenter now points from (approx) CM to the planet
		rj.x=gTidingCenter.x-bodyDescripts[4*j+0];//
		rj.y=gTidingCenter.y-bodyDescripts[4*j+1];// rj now points from body j to the planet
		rj.z=gTidingCenter.z-bodyDescripts[4*j+2];//

		rj.setMagnitude(gUniversalBigG*gTidingMass*bodyDescripts[4*j+3]/(rj.magnitudeSquared()+softsqr));

		forces[3*j+0]+=rj.x; forces[3*j+1]+=rj.y; forces[3*j+2]+=rj.z;

		// 3. fictitious force term on the accelerated frame
		rj=-gTidingCenter; // rj now points from planet to (approx) CM
		rj.setMagnitude(gUniversalBigG*gTidingMass*bodyDescripts[4*j+3]/(rj.magnitudeSquared()+softsqr));
		forces[3*j+0]+=rj.x; forces[3*j+1]+=rj.y; forces[3*j+2]+=rj.z;
		}
	}

	for (int j=0; j<nbActors; j++){
		force.x=forces[3*j+0]; force.y=forces[3*j+1]; force.z=forces[3*j+2];
		actors[j]->addForce(force);
	}	
	
	delete [] bodyDescripts;
	delete [] forces;
#endif
	return;
}
void RenderSpecialStuff()
{
	if (bTrackingCamera && !bPause) ResetCamera();
	if (gExperimentType==TIDAL_ENCOUNTER) {
		DrawArrow(gSystemCM,gTidingCenter,NxVec3(1.,0.,0.));
		//MarkClusters();
	}
}
void MarkClusters()
{
	/*Find clusters
	% I use eclazz to find the equivalence classes of the dataset of point
	% coordinates with the equivalence relation inClust such that particle i is
	% inClust of particle j if the distance between their centers is less than
	% the sum of their radii plus a small "skin" correction. This is clearly
	% reflexive and symmetric, and eclazz assumes it is transitive so it doesn't
	% have to be defined as such explicitly.*/
	VecInt clustVec(gScene->getNbActors());
	eclazz(clustVec,inClust);

	/*clustVec now contains cluster labels. next we need to identify clusters with
	more than 2 grains, find the size of an inscribing sphere around the clust CM, and draw it*/
	vector<int> clusts(1); // make a unique(clustVec)
	clusts[0]=clustVec[0];


	/*glPushMatrix();
	SetupGLMatrix(pos,NxMat33(NX_IDENTITY_MATRIX));
	glScalef(r,r,r);
	glutWireSphere(1.0,12,12);
	glPopMatrix();*/
}
bool inClust(const Int a1, const Int a2)
{
	NxActor** actors = gScene->getActors();
	NxActor*  act1=actors[a1-1];
	NxActor*  act2=actors[a2-1];
	if (!act1->isDynamic() || !act2->isDynamic()) return false;
	NxReal dist=(act1->getCMassGlobalPosition()-act2->getCMassGlobalPosition()).magnitude();
	NxReal vol1=act1->getMass()/gDefaultDensity;
	NxReal vol2=act2->getMass()/gDefaultDensity;
	NxReal ssize=NxMath::pow(vol1*3.f/4.f/3.14f,1.0f/3.0f)+NxMath::pow(vol2*3.f/4.f/3.14f,1.0f/3.0f);;
	if (dist<ssize*1) return true;
	else return false;
}
void eclazz(VecInt_O &nf, Bool equiv(const Int, const Int))
{
	Int kk,jj,n=nf.size();
	nf[0]=0;
	for (jj=1;jj<n;jj++) {
		nf[jj]=jj;
		for (kk=0;kk<jj;kk++) {
			nf[kk]=nf[nf[kk]];
			if (equiv(jj+1,kk+1)) nf[nf[nf[kk]]]=jj;
		}
	}
	for (jj=0;jj<n;jj++) nf[jj]=nf[nf[jj]];
}
void UpdateConservables()
{
	NxU32 nbActors = gScene->getNbActors();
	NxActor** actors = gScene->getActors();
	gT_tot=0.0f; gL_tot.zero(); gP_tot.zero();
	while (nbActors--)
	{
		NxActor* actor = *actors++;
		if (!actor->isDynamic()) continue;
		switch (gExperimentType)
		{
		default:
			{
				gP_tot+=actor->getLinearMomentum();
				gT_tot+=actor->computeKineticEnergy();
				NxVec3 L=actor->getAngularMomentum();
				actor->getGlobalOrientation().multiplyByTranspose(L,L);
				gL_tot+=L;
				gL_tot+=actor->getCMassGlobalPosition().cross(actor->getLinearMomentum());
				break;
			}
		
		}
	}
}
bool CreateInitialPile()
{
	/* Initial conditions are required for most experiments. The initial conditions file should
		contain a rectangular matrix consisting of columns for x,y,z,vx,vy,vz and one row per 
		actor. The number of rows in the ic file will supersede the rubble_size line in the ini
		file and RUBBLE_SIZE will be replaced. The initial conditions file should be in the run
		directory and has the run base name with an extension .ic.
	*/

	// peek into the ic file to find header lines.
	char icfile[255];
	sprintf(icfile,"%s\\%s\\%s.ic",_getcwd(NULL,0),gRunBaseName,gRunBaseName);
	ifstream fp(icfile);
	if (fp.fail()){cout<<"ERROR:could not open initial conditions file\a"<<endl; return false;}
	char line[255];
	unsigned int headlines=0;
	while (fp.good()){
		fp.getline(line,255);
		headlines++;
		if (strcmp(line,"---END HEADER---")==0) break;
	}
	if (fp.eof()) headlines=0; // explicit header/data separator not found. assume no header.
	fp.close();

	// read initial conditions into a temporary holder. expected size is nb_grains-by-6.
	vector<vector<float>*>* ics=RSReadMatrixFromFile(icfile,headlines);
	if (ics==NULL){cout<<"ERROR:failed to read initial conditions"<<endl; return false;}
	if ((ics->at(0)->size())!=6){cout<<"ERROR:wrong number of initial conditions per actor"<<endl; return false;}
	RUBBLE_SIZE=ics->size(); // replacing number possibly specified in ini file.

	// place an actor with appropriate initial conditions
	for (int k=0; k<RUBBLE_SIZE; k++){
		NxVec3 pos(ics->at(k)->at(0),ics->at(k)->at(1),ics->at(k)->at(2));
		NxVec3 vel(ics->at(k)->at(3),ics->at(k)->at(4),ics->at(k)->at(5));
		bool success=PlaceActor(pos,vel,k);
		if (!success){cout<<"ERROR:actor placement failed, possibly bad IC\a"<<endl; return false;}
	}

	// all actors placed. get rid of temp holder and return.
	delete ics;
	return true;
}
bool PlaceActor(NxVec3 position, NxVec3 velocity, unsigned int actind)
/*	This function is an intermediate step between CreateInitialPile that defines that initial positions
	and velocities of each actor in the pile, and the actor creation routine, which may be simply the
	familiar CreateRubbleGrain if actors should have a generic (or random) shape, or the more specialized
	CreateSpecificActor, that reads the shape information from a shapes file. The branching is decided
	simply by looking for the existence of a shapes file. EDIT: a shape file for non-uniform spheres is now
	implemented, as is a mixed polyhedra-non-uniform-sphere file. Branching determined by the section
	separator: SPHERE or SHAPE.
*/
{
	bool success=false;
	static bool skipSHFile=false;

	// attempt to open a shapes file
	char shfile[255];
	sprintf(shfile,"%s\\%s\\%s.sh",_getcwd(NULL,0),gRunBaseName,gRunBaseName);
	static ifstream fp(shfile);
	if (fp.fail() || skipSHFile) {// the easy default path. creating actors based on grain_type form ini file.
		NxActor* actor=CreateRubbleGrain(position);
		if (actor) actor->setLinearVelocity(velocity);
		if (actor) success=true;
	}
	else // this is the more complex path. read shapes from the shapes file.
	{
		// sigh. ok, here we go. the .sh file remains open so we just need to skip the header once:
		static char line[255];
		if (actind==0){
			while (fp.good()){
				fp.getline(line,255);
				if (strcmp(line,"SHAPE")==0 || strcmp(line,"SPHERE")==0) break;
			}
			if (fp.eof()) {cout<<"ERROR:shape segments not found in shape file."<<endl; fp.close(); skipSHFile=true; return false;}
		}

		// otherwise, we should be in position to start reading vertex info.
		if (strcmp(line,"SHAPE")==0)
		{
			float xx,yy,zz;
			NxVec3 vertex;
			vector<NxVec3> verts;
			while (fp.good())
			{
				fp.getline(line,255);
				if (strcmp(line,"\0")==0) continue; //skip empty lines. (but lines with nothing but ws are a problem.)
				if (strcmp(line,"SHAPE")==0 || strcmp(line,"SPHERE")==0 || fp.eof()) break;
				if (sscanf_s(line,"%f%f%f",&xx,&yy,&zz)<3) {cout<<"ERROR:bad format in shapes file."<<endl; fp.close(); return false;}
				vertex.x=xx; vertex.y=yy; vertex.z=zz;
				verts.push_back(vertex);
			}
			if (fp.eof() && actind<(RUBBLE_SIZE-1)){
				cout<<"WARNING:not enough shapes in shape file. additional shapes determined by grain_type"<<endl;
				fp.close();
				skipSHFile=true;
			}

			// well well. i do believe we have extracted the vertex information for shape actind. now what?
			NxActor* actor=CreateSpecificConvexShape(position,gDefaultDensity,verts); // someone else's problem now.
			if (actor){
				actor->setLinearVelocity(velocity);
				success=true;
			}
		}
		else if (strcmp(line,"SPHERE")==0)
		{
			float rr;
			fp.getline(line,255);
			if (sscanf_s(line,"%f",&rr)<1) {cout<<"ERROR:bad format in shapes file."<<endl; fp.close(); return false;}
			while (fp.good())
			{
				fp.getline(line,255);
				if (strcmp(line,"\0")==0) continue; //skip empty lines. (but lines with nothing but ws are a problem.)
				if (strcmp(line,"SHAPE")==0 || strcmp(line,"SPHERE")==0 || fp.eof()) break;
				cout<<"ERROR:bad format in shapes file."<<endl; fp.close(); return false;
			}
			if (fp.eof() && actind<(RUBBLE_SIZE-1)){
				cout<<"WARNING:not enough shapes in shape file. additional shapes determined by grain_type"<<endl;
				fp.close();
				skipSHFile=true;
			}
			NxActor* actor=CreateSphericalGrain(position,gDefaultDensity,rr*2);
			if (actor){
				actor->setLinearVelocity(velocity);
				success=true;
			}
		}
		else
		{
			cout<<"WARNING:unknown *shape* in shapes file. additional shapes determined by grain_type"<<endl;
			fp.close();
			skipSHFile=true;
		}
	}

	// done.
	return success;
}
void ReCreateExperiment()
{
	ReadOrbit();
	UpdateConservables();
	gSystemCM=LocateSystemCM();
	isRunning=true;
	if (!isManualControl) {bDebugWireframeMode=false; bPause=false;}
}
void FireAction()
{
	gScene->simulate(0);
	gScene->checkResults(NX_ALL_FINISHED,true); // make sure time step is done
	Reveille(); // wake up everyone
}
void AutoExperimentControl()
{
	gScene->simulate(0);
	gScene->checkResults(NX_ALL_FINISHED,true); // make sure time step is done
	return;
}
void StopRun(char reason[255])
{
	isRunning=false;
	DumpToFile();
	//CaptureScene(); don't know why this doesn't work
	printf("Run Stopped. (Reason: %s)\n",reason);
	char logFile[255];
	sprintf(logFile,"%s\\%s\\%s.log",_getcwd(NULL,0),gRunBaseName,gRunBaseName);
	LogRunStats(logFile);
	LogExperimentStats(logFile);
	WritePrivateProfileString("experiment","stop_reason",reason,logFile);
	SetThreadExecutionState(ES_CONTINUOUS); // save electricity for the uni
	if (bQuitWhenDone) exit(0);
	UpdateHUD();
}
void LogExperimentStats(char logFile[255])
{
	char s[255];

	// experiment options
	switch (gExperimentType)
	{
	default:				{sprintf(s,"%s","?"); break;}
	}
	WritePrivateProfileString("experiment","experiment_type",s,logFile);
	sprintf(s,"%d",gRubbleType);
	WritePrivateProfileString("experiment","rubble_type",s,logFile);
	sprintf(s,"%d",RUBBLE_SIZE);
	WritePrivateProfileString("experiment","pile_size",s,logFile);
	sprintf(s,"%f",gRainSize);
	WritePrivateProfileString("experiment","monomer_size",s,logFile);
	sprintf(s,"%f",gDefaultDensity);
	WritePrivateProfileString("experiment","monomer_density",s,logFile);
	sprintf(s,"%d",isUniformRubble);
	WritePrivateProfileString("experiment","uniform_pile",s,logFile);
	sprintf(s,"%f",gDynamicFriction);
	WritePrivateProfileString("materials","default_dfriction",s,logFile);
	sprintf(s,"%f",gStaticFriction);
	WritePrivateProfileString("materials","default_sfriction",s,logFile);
	sprintf(s,"%f",gRestitution);
	WritePrivateProfileString("materials","default_restitution",s,logFile);
	
	return;
}
void ResetCamera()
{
	// Reset camera position. Move camera to where it can watch the whole scene.
	NxVec3 camPos;
	float viewAngle=NxMath::degToRad(30.0f);
	switch (gExperimentType)
	{
	case TIDAL_ENCOUNTER:
		{
			gSystemCM=LocateSystemCM();
			NxActor* extremers[6]; extremers[0]=NULL;
			GetExtremeActors(extremers);
			if (extremers[0]==NULL) return;
			float left=extremers[0]->getCMassGlobalPosition().x+gRainSize;
			float right=extremers[1]->getCMassGlobalPosition().x-gRainSize;
			float up=extremers[2]->getCMassGlobalPosition().y+gRainSize;
			float down=extremers[3]->getCMassGlobalPosition().y-gRainSize;
			float in=extremers[4]->getCMassGlobalPosition().z;
			float out=extremers[5]->getCMassGlobalPosition().z;

			if (bStalkingCamera)
			{
				camPos=gSystemCM;
				gCameraForward=(gTidingCenter-camPos);
				gCameraForward.normalize();
				gCameraRight.cross(gCameraForward,NxVec3(0,0,1));
				gCameraForward.normalize();

				camPos-=30*gRainSize*gCameraForward;
				camPos-=8*gRainSize*gCameraRight;
				
				float dist=(camPos-gSystemCM).magnitude();
				while ((up-down)/2/dist>NxMath::tan(viewAngle) || (in-out)/2/dist>NxMath::tan(viewAngle))
				{
					camPos*=1.01;
					dist=(camPos-gSystemCM).magnitude();
				}

				gCameraPos=camPos;
				gCameraAspectRatio=1.0f;
				gCameraSpeed=gRainSize;
			} 
			else
			{
				camPos.x=gSystemCM.x;
				camPos.y=gSystemCM.y;
				camPos.z=out-20*gRainSize;
				camPos.z=NxMath::min(camPos.z,out-((left-right)/2/NxMath::tan(viewAngle)));
				camPos.z=NxMath::min(camPos.z,out-((up-down)/2/NxMath::tan(viewAngle)));
				camPos.z=NxMath::max(camPos.z,-0.98f*gCameraZbufFar);
				gCameraPos=camPos;
				gCameraForward=NxVec3(0,0,1);
				gCameraRight=NxVec3(-1,0,0);
				gCameraAspectRatio=1.0f;
				gCameraSpeed=gRainSize;
			}
			break;
		}
	default:
		{
			NxActor* extremers[6]; extremers[0]=NULL;
			GetExtremeActors(extremers);
			if (extremers[0]==NULL) return;
			float left=extremers[0]->getGlobalPosition().x;
			float right=extremers[1]->getGlobalPosition().x;
			float up=extremers[2]->getGlobalPosition().y;
			float down=extremers[3]->getGlobalPosition().y;
			float in=extremers[4]->getGlobalPosition().z;
			float out=extremers[5]->getGlobalPosition().z;
			camPos.x=right+(left-right)/2;
			camPos.y=down+(up-down)/2;
			camPos.z=out-(in-out)-2;
			camPos.z=NxMath::min(camPos.z,-((left-right)/2/NxMath::tan(viewAngle)));
			camPos.z=NxMath::min(camPos.z,-((up-down)/2/NxMath::tan(viewAngle)));
			camPos.z=NxMath::min(camPos.z,-(NxMath::pow(RUBBLE_SIZE,1.0f/3.0f)+2*gRainSize));

			gCameraPos=camPos;
			gCameraForward=NxVec3(0,0,1);
			gCameraRight=NxVec3(-1,0,0);
			gCameraAspectRatio=1.0f;
			break;
		}
	}
}
void PrinterDebug()
{
	if (gSelectedActor==NULL) return;
	NxVec3 v=gSelectedActor->getLinearVelocity();
	return;
}
void ConfigExperimentFromFile(char iniFile[255])
{
	FILE  *chk;
	if (fopen_s(&chk,iniFile,"r")!=0) {printf("Error: Could not find ini file.\a\n"); exit(1);}

	char s[255];

	// configure physical properties
	{
		GetPrivateProfileString("physical","little_g","9.8",s,255,iniFile);
		gDefaultGravity=NxVec3(0,-atof(s),0);
		GetPrivateProfileString("physical","big_G","0",s,255,iniFile);
		gUniversalBigG=atof(s);
		GetPrivateProfileString("physical","gravity_softening_factor","0",s,255,iniFile);
		gGravSoftFactor=atof(s);
	}

	// read experiment-specific parameters
	{
		GetPrivateProfileString("experiment","experiment_type","NULL",s,255,iniFile);
		if		(strcmp(s,"free_pile")==0)
			gExperimentType=FREE_PILE;
		else if (strcmp(s,"tidal_encounter")==0)
			gExperimentType=TIDAL_ENCOUNTER;
		else
			gExperimentType=BAD_EXPERIMENT_TYPE;

		GetPrivateProfileString("experiment","tiding_mass","0",s,255,iniFile);
		gTidingMass=atof(s);

		GetPrivateProfileString("physical","gravity_type","NULL",s,255,iniFile);
		if		(strcmp(s,"all_pairs")==0)
			{gRavityType=ALL_PAIRS; bGravity=false;}
		else
			gRavityType=BAD_GRAVITY_TYPE;

		GetPrivateProfileString("experiment","rubble_type","",s,255,iniFile);
		if		(strcmp(s,"boxes")==0)
			gRubbleType=BOXES;
		else if (strcmp(s,"spheres")==0)
			gRubbleType=SPHERES;
		else if (strcmp(s,"icosahedra")==0)
			gRubbleType=ICOSAHEDRA;
		else if (strcmp(s,"tetrahedra")==0)
			gRubbleType=TETRAHEDRA;
		else if (strcmp(s,"mixed")==0)
			gRubbleType=MIXED;
		else if (strcmp(s,"convex")==0)
			gRubbleType=CONVEX_POLYHEDRA;
		else if (strcmp(s,"d12")==0)
			gRubbleType=D12;
		else if (strcmp(s,"capsules")==0)
			gRubbleType=CAPSULES;
		else if (strcmp(s,"mixed")==0)
			gRubbleType=MIXED;
		else
			gRubbleType=BAD_RUBBLE_TYPE;

		RUBBLE_SIZE=GetPrivateProfileInt("experiment","rubble_size",0,iniFile);

		GetPrivateProfileString("experiment","grain_size","0",s,255,iniFile);
		gRainSize=atof(s);

		GetPrivateProfileString("experiment","grain_density","1000",s,255,iniFile);
		gDefaultDensity=atof(s);

		gDefaultNumVerts=GetPrivateProfileInt("experiment","num_verts",12,iniFile);

		GetPrivateProfileString("experiment","vertex_spacing","10",s,255,iniFile);
		gDefaultGam=atof(s);

		GetPrivateProfileString("experiment","uniform_rubble","false",s,255,iniFile);
		if (strcmp(s,"true")==0)
			isUniformRubble=true;
		else
			isUniformRubble=false;

		GetPrivateProfileString("experiment","control_mode","manual",s,255,iniFile);
		if (strcmp(s,"manual")==0)
			isManualControl=true;
		else
			isManualControl=false;
	}
}
void InitRun(int argc, char** argv)
{
	char iniFile[255];
	if (argc==1)
		sprintf(iniFile,"%s\\tidal.ini",_getcwd(NULL,0));
	else
		sprintf(iniFile,"%s\\%s\\%s.ini",_getcwd(NULL,0),argv[1],argv[1]);
	ConfigSDKFromFile(iniFile);
	ConfigExperimentFromFile(iniFile);
	if (argc>1 && strcmp(gRunBaseName,argv[1])!=0)
		//cout<<"WARNING:run-base-name mismatch"<<endl;
		nerror("run-base-name mismatch\a");
	_mkdir(gRunBaseName); // ok if dir exists, nothing bad will happen.
	if (gFrameCounter==0)
	{
		bFreshRun=true;
	}
	else
	{
		char s[255];
		sprintf(s,"%s\\dump.%u",gRunBaseName,gFrameCounter);
		bFreshRun=!LoadFromFile(s,false);
		if (bFreshRun)
			printf("Failed to load initial conditions. Starting a fresh run.\n");
	}
}
int main(int argc, char** argv)
{
	InitWin7();
	InitDevIL();
	InitRun(argc,argv);
    InitGlut(argc, argv);
	if (bFreshRun)
		InitNx();
	else
		ReInitNx();
    glutMainLoop();
	ReleaseNx();
	return 0;
}
