// ===============================================================================
//  a set of experiments with asteroid instrument pods
// ===============================================================================

#include "nshared.h"
#include "pods.h"

// Experiment globals
enum ExperimentType		gExperimentType;
enum LanderType			gLanderType;
NxActor*	gContainerCell[4]={NULL,NULL,NULL,NULL}; // NxPlanes will be used to contain a finite cell
NxReal		gCellSize;
NxReal		gCellMatParam[2];
NxReal		gLanderUpMatParam[2];
NxReal		gLanderDownMatParam[2];
NxMaterial* gCellMaterial=NULL;
NxMaterial* gLanderUpMaterial=NULL;
NxMaterial* gLanderDownMaterial=NULL;
NxVec3		gLanderSize; // for a shaped charge this is (width,depth,height), for a roly-poly this is (radius,offset,r_bot/r_top)
NxReal		gLanderMass;
NxReal		gRPWeight;
NxReal		gLanderMassRatio;
NxActor*	gLander=NULL;
NxReal		gLandingRoughness;
NxReal		gRavitatorThresholdMass;
bool		bVibrator;
NxReal		gVibrateAmp;
NxReal		gVibrateFreq;
Ran			gRNG(17);
vector<NxActor*> gRavitators; // hold unknown number of big boulders
float		gsd[3]; // holds the parameters of the grain size distribution


void CreateExperiment()
{
	switch (gExperimentType)
	{
	case LANDER:
		  {
			  gScene->releaseActor(*groundPlane); groundPlane=NULL; // will be recreated with the container cell
			  CreateExperimentMaterials();
			  CreateContainerCell(gCellSize);
			  gPhysicsSDK->setParameter(NX_VISUALIZE_ACTOR_AXES,0);
			  gPhysicsSDK->setParameter(NX_VISUALIZE_COLLISION_FNORMALS,0);
			  gPhysicsSDK->setParameter(NX_VISUALIZE_COLLISION_SHAPES,1);
			  gPhysicsSDK->setParameter(NX_VISUALIZE_WORLD_AXES,0);
			  //gLander=CreateLander(NxVec3(0.0f,0.0f,0.0f));
			  gLander=CreateLander(NxVec3(0.0f,gCellSize/2.0f,0.0f));
			  IdentifyGravitators();
			  if (gLander) {
			  ReOrientActor(gLander);
			  SpinActor(gLander,gLandingRoughness*sqrt(gCellSize*gDefaultGravity.magnitude()/(gLander->getMassSpaceInertiaTensor().magnitude())));
			  LaunchActor(gLander,gLandingRoughness*sqrt(gCellSize*gDefaultGravity.magnitude()));
			  gSelectedActor=gLander;
			  gPhysicsSDK->setParameter(NX_VISUALIZE_ACTOR_AXES,gLanderSize.magnitude());
			  isRunning=true;
			  }
			  else {
				  printf("Error: Unable to create lander\a\n");
				  isRunning=false;
			  }
			  break;
		  }
	case LAY_SUBSTRATE:
		{
			gScene->releaseActor(*groundPlane); groundPlane=NULL; // will be recreated with the container cell
			CreateExperimentMaterials();
			CreateContainerCell(gCellSize);
			gPhysicsSDK->setParameter(NX_VISUALIZE_ACTOR_AXES,0);
			gPhysicsSDK->setParameter(NX_VISUALIZE_COLLISION_FNORMALS,0);
			gPhysicsSDK->setParameter(NX_VISUALIZE_COLLISION_SHAPES,1);
			gPhysicsSDK->setParameter(NX_VISUALIZE_WORLD_AXES,0);
			//ThrowStone();
			isRunning=true;
			break;
		}
	default:				
		{printf("Error: unimplemented experiment type\a\n"); isRunning=false; break;}
	}
}
NxActor* CreateLander(NxVec3 position)
{
	NxActor* actor=NULL;
	switch (gLanderType)
	{
		case CAPSULE_LANDER:
			actor = CreateCapsuleGrain(position,gLanderMass/pow(gLanderSize.x,3.0f),gLanderSize.x);
			break;
		case SPHERE_LANDER:
			actor = CreateSphericalGrain(position,gLanderMass/pow(gLanderSize.x,3.0f),gLanderSize.x);
			break;
		case ROLY_POLY:
			{// Two slightly offset spheres, one heavy, one slightly smaller, and hollow
				NxActorDesc actorDesc;
				NxBodyDesc bodyDesc;
				NxReal radius=gLanderSize.x;

				// The bottom "roly"
				NxSphereShapeDesc rolyDesc;
				rolyDesc.radius=radius*gLanderSize.z; // relative to top hemisphere
				rolyDesc.localPose.t=NxVec3(0.0f,radius,0.0f);
				rolyDesc.materialIndex=gLanderDownMaterial->getMaterialIndex(); // note: both hemispheres use "down" material
				rolyDesc.mass=1; // relative to top hemisphere
				rolyDesc.skinWidth=radius/60;
				rolyDesc.name="downside";
				actorDesc.shapes.pushBack(&rolyDesc);

				// The top "poly"
				NxSphereShapeDesc polyDesc;
				polyDesc.radius=radius;
				polyDesc.localPose.t=NxVec3(0.0f,radius+gLanderSize.y,0.0f);
				polyDesc.materialIndex=gLanderDownMaterial->getMaterialIndex(); // note: both hemispheres use "down" material
				polyDesc.mass=1;
				polyDesc.skinWidth=radius/60;
				polyDesc.name="upside";
				actorDesc.shapes.pushBack(&polyDesc);

				// The weight
				NxSphereShapeDesc wDesc;
				wDesc.radius=radius/10;
				wDesc.localPose.t=NxVec3(0.0f,wDesc.radius,0.0f);
				wDesc.materialIndex=gLanderDownMaterial->getMaterialIndex(); // note: this is irrelevant
				wDesc.mass=gRPWeight; // relative to top hemisphere
				wDesc.skinWidth=radius/120;
				wDesc.name="weight";
				actorDesc.shapes.pushBack(&wDesc);

				// the lander
				bodyDesc.mass=gLanderMass;
				actorDesc.body=&bodyDesc;
				actorDesc.globalPose.t=position;
				assert(actorDesc.isValid());
				actor=gScene->createActor(actorDesc);
				assert(actor);
				actor->setName("~lander");
				//actor->raiseBodyFlag(NX_BF_KINEMATIC);

				break;
			}
		case SHAPED_CHARGE:
			{// Two flat plates, like a two-layer ruler, one heavy and absorbing, one hollow and elastic
				NxActorDesc actorDesc;
				NxBodyDesc bodyDesc;
				NxVec3 boxDim=gLanderSize/2; // remember NxShapeDesc.dimension is the "radius"
				boxDim.y/=2; // remember we're going to stack two of these

				// The down side
				NxBoxShapeDesc downsideDesc;
				downsideDesc.dimensions.set(boxDim);
				downsideDesc.localPose.t=NxVec3(0.0f,boxDim.y,0.0f);
				downsideDesc.materialIndex=gLanderDownMaterial->getMaterialIndex();
				downsideDesc.density=gLanderMassRatio;
				downsideDesc.skinWidth=boxDim.y/60;
				downsideDesc.name="downside";
				actorDesc.shapes.pushBack(&downsideDesc);

				// The up side
				NxBoxShapeDesc upsideDesc;
				upsideDesc.dimensions.set(boxDim);
				upsideDesc.localPose.t=NxVec3(0.0f,3.0*boxDim.y,0.0f);
				upsideDesc.materialIndex=gLanderUpMaterial->getMaterialIndex();
				upsideDesc.density=1;
				upsideDesc.skinWidth=boxDim.y/60;
				upsideDesc.name="upside";
				actorDesc.shapes.pushBack(&upsideDesc);

				// The Lander
				bodyDesc.mass=gLanderMass;
				actorDesc.body = &bodyDesc;
				actorDesc.globalPose.t	= position;
				assert(actorDesc.isValid());
				actor = gScene->createActor(actorDesc);	
				assert(actor);
				actor->setName("~lander");

				break;
			}
		case CUBESAT:
			{// JPL's new fad, comes in units of 10x10x10 cm
				NxActorDesc actorDesc;
				NxBodyDesc bodyDesc;
				NxVec3 boxDim(0.05); // remember NxShapeDesc.dimension is the "radius"

				// U1 - the dead weight
				NxBoxShapeDesc U1Desc;
				U1Desc.dimensions.set(boxDim);
				U1Desc.localPose.t=NxVec3(0.0f,boxDim.y,0.0f);
				U1Desc.materialIndex=gLanderDownMaterial->getMaterialIndex();
				U1Desc.mass=10*gLanderMass;
				U1Desc.skinWidth=boxDim.y/60;
				U1Desc.name="U1";
				actorDesc.shapes.pushBack(&U1Desc);

				// U2
				NxBoxShapeDesc U2Desc;
				U2Desc.dimensions.set(boxDim);
				U2Desc.localPose.t=NxVec3(0.0f,3.0*boxDim.y,0.0f);
				U2Desc.materialIndex=gLanderUpMaterial->getMaterialIndex();
				U2Desc.mass=gLanderMass;
				U2Desc.skinWidth=boxDim.y/60;
				U2Desc.name="U2";
				actorDesc.shapes.pushBack(&U2Desc);

				// U3
				NxBoxShapeDesc U3Desc;
				U3Desc.dimensions.set(boxDim);
				U3Desc.localPose.t=NxVec3(0.0f,5.0*boxDim.y,0.0f);
				U3Desc.materialIndex=gLanderUpMaterial->getMaterialIndex();
				U3Desc.mass=gLanderMass;
				U3Desc.skinWidth=boxDim.y/60;
				U3Desc.name="U3";
				actorDesc.shapes.pushBack(&U3Desc);

				// The Cubesat
				bodyDesc.mass=gLanderMass;
				actorDesc.body = &bodyDesc;
				actorDesc.globalPose.t	= position;
				assert(actorDesc.isValid());
				actor = gScene->createActor(actorDesc);	
				assert(actor);
				actor->setName("~lander");
				break;
			}
		default:
			printf("Error: unknown lander type\n");
			break;
	}

	return actor;
}
bool ThrowStone()
{
	bool success=false;
	// figure out grain size. draw a size from a distribution dN = C * D^-alpha dD (Asphaug 2009 AnnRev)
	NxReal drawnSize = DrawFromBoundedPowerLaw(gsd[0],gsd[1],gsd[2]); // CHOOSE GOOD VALUES HERE
	
	// make a random grain, use default material
	NxActor* actor=NULL;
	actor=CreateConvexPolyGrain(NxVec3(0,2,0),gDefaultDensity,drawnSize,false,(8 + gRNG.int32()%8));
	if (actor) {
		actor->setName("rubble");
		if (gRavitatorThresholdMass>0 && actor->getMass()>gRavitatorThresholdMass)
			actor->setName("boulder");
	}
	// randomize and let fly
	if (actor) {
		ReOrientActor(actor);
		SpinActor(actor,1);
		LaunchActor(actor,1);
	}
	success=true;
	return success;
}
NxReal DrawFromBoundedPowerLaw( float a0, float a1, float alfa )
{
	// using the inverse transformation law for probabilities. see for example, NR3 p.362, although the explanation sucks.
	// for the pdf p(a)=C*a^-alfa, between a0 and a1.
	float t=gRNG.doub();
	float onemalfa=1-alfa;
	float a = NxMath::pow(t*(NxMath::pow(a1,onemalfa)-NxMath::pow(a0,onemalfa)) + NxMath::pow(a0,onemalfa) , 1.0f/(onemalfa));
	return a;
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
	separator: SPHERE or SHAPE. EDIT: this function should not be called except by CreateInitialPile.
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
void CreateExperimentMaterials()
{
	NxMaterialDesc matDesc;
	
	// cell material, smooth and bouncy
	matDesc.restitution=gCellMatParam[0];
	matDesc.dynamicFriction=gCellMatParam[1];
	matDesc.staticFriction=gCellMatParam[1];
	matDesc.frictionCombineMode=NX_CM_MIN;
	matDesc.restitutionCombineMode=NX_CM_MAX;
	gCellMaterial=gScene->createMaterial(matDesc);

	// lander up side material, smooth and bouncy
	matDesc.restitution=gLanderUpMatParam[0];
	matDesc.dynamicFriction=gLanderUpMatParam[1];
	matDesc.staticFriction=gLanderUpMatParam[1];
	matDesc.frictionCombineMode=NX_CM_AVERAGE;
	matDesc.restitutionCombineMode=NX_CM_MAX;
	gLanderUpMaterial=gScene->createMaterial(matDesc);

	// lander down side material, rough and crunchy
	matDesc.restitution=gLanderDownMatParam[0];
	matDesc.dynamicFriction=gLanderDownMatParam[1];
	matDesc.staticFriction=gLanderDownMatParam[1];
	matDesc.frictionCombineMode=NX_CM_MAX;
	matDesc.restitutionCombineMode=NX_CM_MIN;
	gLanderDownMaterial=gScene->createMaterial(matDesc);

	return;
}
void CreateContainerCell(NxReal width)
{
	// We need to create 4 vertical, static planes (facing inward, why not?) plus a ground plane
	NxPlaneShapeDesc planeDesc;
	NxActorDesc actorDesc;

	// ground plane, facing up, same material as grains
	planeDesc.normal=NxVec3(0,1,0);
	planeDesc.d=0;
	planeDesc.materialIndex=0; // use default material for ground plane
	actorDesc.shapes.push_back(&planeDesc);
	groundPlane=gScene->createActor(actorDesc);
	groundPlane->setName("ground");
	actorDesc.shapes.pop_back();

	// left wall, facing right
	planeDesc.normal=NxVec3(-1,0,0);
	planeDesc.d=-0.5*width;
	planeDesc.materialIndex=gCellMaterial->getMaterialIndex();
	actorDesc.shapes.push_back(&planeDesc);
	gContainerCell[0]=gScene->createActor(actorDesc);
	gContainerCell[0]->setName("left wall");
	actorDesc.shapes.pop_back();

	// right wall, facing left
	planeDesc.normal=NxVec3(1,0,0);
	planeDesc.d=-0.5*width;
	planeDesc.materialIndex=gCellMaterial->getMaterialIndex();
	actorDesc.shapes.push_back(&planeDesc);
	gContainerCell[1]=gScene->createActor(actorDesc);
	gContainerCell[1]->setName("right wall");
	actorDesc.shapes.pop_back();

	// in wall, facing out
	planeDesc.normal=NxVec3(0,0,1);
	planeDesc.d=-0.5*width;
	planeDesc.materialIndex=gCellMaterial->getMaterialIndex();
	actorDesc.shapes.push_back(&planeDesc);
	gContainerCell[2]=gScene->createActor(actorDesc);
	gContainerCell[2]->setName("in wall");
	actorDesc.shapes.pop_back();

	// out wall, facing in
	planeDesc.normal=NxVec3(0,0,-1);
	planeDesc.d=-0.5*width;
	planeDesc.materialIndex=gCellMaterial->getMaterialIndex();
	actorDesc.shapes.push_back(&planeDesc);
	gContainerCell[3]=gScene->createActor(actorDesc);
	gContainerCell[3]->setName("out wall");
	actorDesc.shapes.pop_back();

	return;
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
} 
else
{
	hud.AddDisplayString("",0.84f,0.96f);
}
	if (bRealTime)
		hud.SetDisplayColor(0,1,0,0,1);

	// Add pause to HUD - string id 1
	if (bPause)  
		hud.AddDisplayString("Paused - Hit \"p\" to Unpause", 0.3f, 0.55f);
	else
		hud.AddDisplayString("", 0.0f, 0.0f);

	// Add experiment title - string id 2
	sprintf(s,"POD Experiment: %s",gRunBaseName);
	hud.AddDisplayString(s,0.01,0.96);
	hud.SetDisplayColor(2,0,0,1,1);

	// Add NbActors - string id 3
	sprintf(s,"%s%d","Actors: ",gScene->getNbDynamicShapes());
	hud.AddDisplayString(s,0.01,0.87);

	// Add rubble type - string id 4
	if (bVerboseHUD)
		sprintf(s,"Rubble type: %d",gRubbleType);
	else
		sprintf(s,"");
	hud.AddDisplayString(s,0.01,0.93);

	// Add Experiment type - string id 5
	if (bVerboseHUD)
		sprintf(s,"Experiment Type: %d",gExperimentType);
	else
		sprintf(s,"");
	hud.AddDisplayString(s,0.01,0.90);

	// Add time - string id 6
	sprintf(s,"Wall time: %d s, Code time: %d s",0,0);
	hud.AddDisplayString(s,0.01,0.84);

	// Add FPS - string id 7
	sprintf(s,"FPS: %0.1f(%0.1f)",0.0,0.0);
	hud.AddDisplayString(s,0.01,0.81);

	// Add selected actor name - string id 8
	if (bVerboseHUD)
{
	if (gSelectedActor)
			sprintf(s,"Selected actor: %s",gSelectedActor->getName());
		else
			sprintf(s,"Selected actor: NULL");
} 
else
{
	sprintf(s,"");
}
	hud.AddDisplayString(s,0.98-strlen(s)/100.0f,0.93);

	// Add selected actor specs - string id 9
	sprintf(s,"");
	hud.AddDisplayString(s,0.98-strlen(s)/100.0f,0.90);

	// Add scene diagnostics - string id 10
	sprintf(s,"");
	hud.AddDisplayString(s,0.98,0.87);

}
void UpdateHUD()
{
	char s[255];

	// update scene diagnostics

	// update selected actor specs
	if (bVerboseHUD)
	{
		if (gSelectedActor && gSelectedActor->isDynamic())
		{
			NxVec3 W=gSelectedActor->getAngularVelocity();
			NxVec3 L=gSelectedActor->getAngularMomentum();//resolved in body space
			gSelectedActor->getGlobalOrientation().multiplyByTranspose(L,L);//resolved in world space
			NxMat33 Imin1;
			gSelectedActor->getGlobalInertiaTensor().getInverse(Imin1);
			Imin1.multiply(L,W);
			sprintf(s,"W=(%0.3g,%0.3g,%0.3g)    ",W.x,W.y,W.z);
		}
		else
			sprintf(s,"");
		hud.SetDisplayString(9,s,0.98-strlen(s)/100.0f,0.90);
	}

	// update selected actor name
	if (bVerboseHUD)
	{
		if (gSelectedActor)
			sprintf(s,"Selected actor: %s",gSelectedActor->getName());
		else
			sprintf(s,"Selected actor: NULL");
		hud.SetDisplayString(8,s,0.98-strlen(s)/100.0f,0.93);
	}

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

	// update NbActors
	int ignoredones=5;
	sprintf(s,"%s%d(%d sleeping)","Actors: ",gScene->getNbActors()-ignoredones,NxMath::max(NxI32(0),NxI32(CountSleepers())));
	hud.SetDisplayString(3,s,0.01,0.87);

	// update rubble type
	if (bVerboseHUD)
		sprintf(s,"Rubble type: %d",gRubbleType);
	else
		sprintf(s,"");
	hud.SetDisplayString(4,s,0.01,0.93);

	// update experiment type
	if (bVerboseHUD)
		sprintf(s,"Experiment Type: %d",gExperimentType);
	else
		sprintf(s,"");
	hud.SetDisplayString(5,s,0.01,0.90);

	// update time
	if (bVerboseHUD)
		sprintf(s,"Wall time: %d s, Code time: %0.2g s (dt=%0.2g s)",NxMath::trunc(gWallTime),gCodeTime,gDeltaTime);
	else
		sprintf(s,"");
	hud.SetDisplayString(6,s,0.01,0.84);

	// update FPS
	if (bVerboseHUD)
		sprintf(s,"FPS: %0.1f(%0.1f)",gFPS,gFrameCounter/(NxMath::max(gWallTime,1.f)));
	else
		sprintf(s,"FPS: %0.1f",gFPS);
	hud.SetDisplayString(7,s,0.01,0.81);

	// update experiment title
	if (isRunning)
		hud.SetDisplayColor(2,1,0,0,1);
	else
		hud.SetDisplayColor(2,0,1,0,1);

	return;
}
void AddSpecialForces()
{
	// but not *just* forces. this is a guaranteed point of entry from the main loop (RenderCallback)
	switch (gExperimentType)
	{
	case LANDER:
		if (bVibrator)
		{
			NxReal t=gFrameCounter*gDeltaTime;
			NxReal fx=gVibrateAmp*NxMath::cos(gVibrateFreq*t);
			NxReal fy=gVibrateAmp*NxMath::sin(gVibrateFreq*t);
			NxVec3 F(fx,fy,0);
			gLander->addForce(F);
		}
		GravitateToGravitators();
		break;
	default: break;
	}
}
void UpArrowAction()
{
	bVibrator=true;
}
void DownArrowAction()
{
	bVibrator=false;
}
void LeftArrowAction()
{
}
void RightArrowAction()
{
}
void RenderSpecialStuff()
{
	if (bTrackingCamera && !bPause) ResetCamera();
	DrawLanders();
}
void DrawLanders()
{
	// draw the landers in 2 colors
	NxU32 nbActors = gScene->getNbActors();
	NxActor** actors = gScene->getActors();
	while (nbActors--)
	{
		NxActor* actor = *actors++;
		const char *name=actor->getName();
		if (strcmp(name,"~lander")==0) {
			NxShape*const* shapes;
			shapes = actor->getShapes();
			NxU32 nShapes = actor->getNbShapes();
			//glDisable(GL_LIGHTING);
			char sname[256]="";
			while (nShapes--)
			{
				sprintf(sname,"%s",shapes[nShapes]->getName());
				if      (strcmp(sname,"upside")==0)
					glColor4f(0.0f,1.0f,0.0f,1.0f);
				else if (strcmp(sname,"downside")==0)
					glColor4f(1.0f,0.0f,0.0f,1.0f);
				else if (strcmp(sname,"U1")==0)
					glColor4f(1.0f,0.0f,0.0f,1.0f);
				else if (strcmp(sname,"U2")==0)
					glColor4f(0.0f,1.0f,0.0f,1.0f);
				else if (strcmp(sname,"U3")==0)
					glColor4f(0.0f,0.0f,1.0f,1.0f);
				DrawShape(shapes[nShapes], false);
			}
			//glEnable(GL_LIGHTING);
		}
	}
}
void ReCreateExperiment()
{
	gScene->simulate(0);
	gScene->checkResults(NX_ALL_FINISHED,true); // make sure time step is done
	
	// Re-set the crucial SDK parameters
	if (gPhysicsSDK)
	{
		gPhysicsSDK->setParameter(NX_VISUALIZE_ACTOR_AXES,0);
		gPhysicsSDK->setParameter(NX_VISUALIZE_WORLD_AXES,0);
		if (gBounceEps>0)
			gPhysicsSDK->setParameter(NX_BOUNCE_THRESHOLD,-gBounceEps);
		if (gSkinWidth>0)
			gPhysicsSDK->setParameter(NX_SKIN_WIDTH, gSkinWidth);
	}
	if (gScene)
	{
		if (bVarTimeStep)
			gScene->setTiming(gDeltaTime/4.0,4,NX_TIMESTEP_VARIABLE);
		else
			gScene->setTiming(gDeltaTime/4.0,4,NX_TIMESTEP_FIXED);
		gScene->setGravity(gDefaultGravity*bGravity);
	}

	// Next get to cleaning up the scene that was maybe loaded from file.
	if (gScene)
	{
		// Get rid of defined materials (who knows where they've been...)
		NxMaterial* mats[20]; // I hardly think we'll have more than 20
		NxU32 iter=0; // don't worry about it
		int nmats = gScene->getMaterialArray(mats,20,iter);
		while (nmats--)
		{
			gScene->releaseMaterial(*(mats[nmats]));
		}
		gCellMaterial=NULL;
		gLanderUpMaterial=NULL;
		gLanderDownMaterial=NULL;

		CreateExperimentMaterials(); // make shiny new materials

		// Now, look over all the actors, who knows where they've been. Find out then.
		NxU32 nbActors = gScene->getNbActors();
		NxActor** actors = gScene->getActors();
		while (nbActors--)
		{
			NxActor* actor = *actors++;
			const char *name=actor->getName();
			NxShape*const* shapes;
			shapes = actor->getShapes();
			NxU32 nShapes = actor->getNbShapes();

			// If it's a rubble grain, make sure it uses the default material
			if      (strcmp(name,"rubble")==0 || strcmp(name,"boulder")==0)
				while (nShapes--) shapes[nShapes]->setMaterial(0);

			// If it's a cell container, make sure it uses the cell material, and remember its pointer
			else if (strcmp(name,"left wall")==0) {
				shapes[0]->setMaterial(gCellMaterial->getMaterialIndex());
				gContainerCell[0]=actor;
			}
			else if (strcmp(name,"right wall")==0) {
				shapes[0]->setMaterial(gCellMaterial->getMaterialIndex());
				gContainerCell[1]=actor;
			}
			else if (strcmp(name,"in wall")==0) {
				shapes[0]->setMaterial(gCellMaterial->getMaterialIndex());
				gContainerCell[2]=actor;
			}
			else if (strcmp(name,"out wall")==0) {
				shapes[0]->setMaterial(gCellMaterial->getMaterialIndex());
				gContainerCell[3]=actor;
			}

			// If it's the ground plane it also uses default material
			else if (strcmp(name,"ground")==0) {
				shapes[0]->setMaterial(0);
				groundPlane=actor;
			}

			// If it's the lander, it has two shapes made of different material, but the lander will be recreated elsewhere...

			// Anything else, get rid of
			else
				gScene->releaseActor(*actor);
		
		}

		// If the container is damaged, recreate it
		for (int k=0; k<4; k++)
		{
			if (gContainerCell[k]==NULL) {
				for (int j=0; j<4; j++)
				{
					if (gContainerCell[j]) {
						gScene->releaseActor(*gContainerCell[j]);
						gContainerCell[j]=NULL;
					}
				}
				CreateContainerCell(gCellSize);
				break;
			}
		}
	}

	// OK, all clean, now create the experiment as needed.
	switch (gExperimentType) {
	case LANDER:
		gLander=CreateLander(NxVec3(0.0f,gCellSize/2,0.0f));
		IdentifyGravitators();
		if (gLander) {
			ReOrientActor(gLander);
			SpinActor(gLander,gLandingRoughness*sqrt(gCellSize*gDefaultGravity.magnitude()/(gLander->getMassSpaceInertiaTensor().magnitude())));
			LaunchActor(gLander,gLandingRoughness*sqrt(gCellSize*gDefaultGravity.magnitude()));
			gSelectedActor=gLander;
			//gPhysicsSDK->setParameter(NX_VISUALIZE_ACTOR_AXES,gLanderSize.magnitude()); // not sure why this is a problem when reading from previous xml
			isRunning=true;
		}
		else {
			printf("Error: Unable to create lander\a\n");
			isRunning=false;
		}
		break;
	case LAY_SUBSTRATE:
		{
			isRunning=true;
			break;
		}
	default:
		printf("Error: Unknown experiment type\a\n");
		isRunning=false;
		break;
	}
	if (!isManualControl){bPause=false;}
}
void IdentifyGravitators()
{
	gRavitators.clear();
	NxU32 nbActors = gScene->getNbActors();
	NxActor** actors = gScene->getActors();
	while (nbActors--)
	{
		NxActor* actor = *actors++;
		if (actor->getMass()>gRavitatorThresholdMass)
			gRavitators.push_back(actor);
	}
}
void GravitateToGravitators()
{
	NxU32 nbActors = gScene->getNbActors();
	NxActor** actors = gScene->getActors();
	NxReal m1,m2,G=gUniversalBigG,R;
	NxVec3 r,F;
	while (nbActors--)
	{
		NxActor* actor = *actors++;
		const char* name=actor->getName();
		if (strcmp(name,"~lander")==0)
		{
			F.zero();
			m1=actor->getMass();
			for (int k=0; k<gRavitators.size(); k++)
			{
				r = (gRavitators[k]->getCMassGlobalPosition() - actor->getCMassGlobalPosition());
				R = r.magnitudeSquared();
				m2= gRavitators[k]->getMass();
				F = r;
				F.setMagnitude(G*m1*m2/R);
				actor->addForce(F);
			}
		}
	}
}
void FireAction()
{
	gScene->simulate(0);
	gScene->checkResults(NX_ALL_FINISHED,true); // make sure time step is done
	if (gExperimentType==LANDER)
	{
		Reveille();
		gLander=CreateLander(NxVec3(0.0f));
		gLander->setGlobalPosition(NxVec3(0.0f,gCellSize/2,0.0f));
		ReOrientActor(gLander);
		SpinActor(gLander,gLandingRoughness*sqrt(gCellSize*gDefaultGravity.magnitude()/(gLander->getMassSpaceInertiaTensor().magnitude())));
		LaunchActor(gLander,gLandingRoughness*sqrt(gCellSize*gDefaultGravity.magnitude()));
	}
	else
	{
		Reveille();
	}
}
void AutoExperimentControl()
{
	gScene->simulate(0);
	gScene->checkResults(NX_ALL_FINISHED,true); // make sure time step is done
	switch (gExperimentType)
	{
	case LAY_SUBSTRATE:
		if (gScene->getNbDynamicShapes()<RUBBLE_SIZE) ThrowStone();
		else
			if (CountSleepers()==RUBBLE_SIZE) StopRun("everyone's asleep");
		break;
	}
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
	case LANDER:
		{
			camPos.x=0;
			camPos.y=gCellSize*0.48+gLanderSize.x;
			camPos.z=-gCellSize*0.8-6.0*gLanderSize.x;
			break;
		}
	case LAY_SUBSTRATE:
		{
			camPos.x=0;
			camPos.y=gCellSize*0.48+gRainSize;
			camPos.z=-gCellSize-2*gRainSize;
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
			break;
		}
	}
	gCameraPos=camPos;
	gCameraForward=NxVec3(0,-0.20,1);
	gCameraRight=NxVec3(-1,0,0);
	gCameraAspectRatio=1.0f;
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
		bGravity=true;

		GetPrivateProfileString("physical","big_G","1",s,255,iniFile);
		gUniversalBigG=atof(s);
	}

	// configure special materials
	{
		GetPrivateProfileString("material","cell_restitution","0.5",s,255,iniFile);
		gCellMatParam[0]=atof(s);
		GetPrivateProfileString("material","cell_friction","0.5",s,255,iniFile);
		gCellMatParam[1]=atof(s);

		GetPrivateProfileString("material","lander_upside_restitution","0.5",s,255,iniFile);
		gLanderUpMatParam[0]=atof(s);
		GetPrivateProfileString("material","lander_upside_friction","0.5",s,255,iniFile);
		gLanderUpMatParam[1]=atof(s);
		GetPrivateProfileString("material","lander_downside_restitution","0.5",s,255,iniFile);
		gLanderDownMatParam[0]=atof(s);
		GetPrivateProfileString("material","lander_downside_friction","0.5",s,255,iniFile);
		gLanderDownMatParam[1]=atof(s);
	}

	// read experiment-specific parameters
	{
		GetPrivateProfileString("experiment","experiment_type","NULL",s,255,iniFile);
		if		(strcmp(s,"lander")==0)
			gExperimentType=LANDER;
		else if (strcmp(s,"lay_substrate")==0)
			gExperimentType=LAY_SUBSTRATE;
		else
			gExperimentType=BAD_EXPERIMENT_TYPE;

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

		// Lander parameters
		GetPrivateProfileString("experiment","lander_type","",s,255,iniFile);
		if		(strcmp(s,"capsule")==0)
			gLanderType=CAPSULE_LANDER;
		else if (strcmp(s,"sphere")==0)
			gLanderType=SPHERE_LANDER;
		else if (strcmp(s,"shaped_charge")==0)
			gLanderType=SHAPED_CHARGE;
		else if (strcmp(s,"roly_poly")==0)
			gLanderType=ROLY_POLY;
		else if (strcmp(s,"cubesat")==0)
			gLanderType=CUBESAT;
		else
			gLanderType=BAD_LANDER_TYPE;

		GetPrivateProfileString("experiment","lander_size.w","0",s,255,iniFile);
		gLanderSize.x=atof(s);
		GetPrivateProfileString("experiment","lander_size.h","0",s,255,iniFile);
		gLanderSize.y=atof(s);
		GetPrivateProfileString("experiment","lander_size.d","0",s,255,iniFile);
		gLanderSize.z=atof(s);
		if (gLanderType==ROLY_POLY) {
		GetPrivateProfileString("experiment","roly_poly_radius","0",s,255,iniFile);
		gLanderSize.x=atof(s);
		GetPrivateProfileString("experiment","roly_poly_offset","0",s,255,iniFile);
		gLanderSize.y=atof(s);
		GetPrivateProfileString("experiment","roly_poly_hemisphere_ratio","1",s,255,iniFile);
		gLanderSize.z=atof(s);
		GetPrivateProfileString("experiment","roly_poly_weight_ratio","0",s,255,iniFile);
		gRPWeight=atof(s);
		}

		GetPrivateProfileString("experiment","lander_mass","0",s,255,iniFile);
		gLanderMass=atof(s);
		GetPrivateProfileString("experiment","lander_mass_ratio","1",s,255,iniFile);
		gLanderMassRatio=atof(s);

		GetPrivateProfileString("experiment","landing_roughness","0",s,255,iniFile);
		gLandingRoughness=atof(s);

		GetPrivateProfileString("experiment","vibrator_amplitude","0",s,255,iniFile);
		gVibrateAmp=atof(s);
		GetPrivateProfileString("experiment","vibrator_frequency","0",s,255,iniFile);
		gVibrateFreq=atof(s);

		// rubble parameters
		RUBBLE_SIZE=GetPrivateProfileInt("experiment","rubble_size",0,iniFile);

		GetPrivateProfileString("experiment","grain_size","0",s,255,iniFile);
		gRainSize=atof(s);

		GetPrivateProfileString("experiment","gsd.amin","0.1",s,255,iniFile);
		gsd[0]=atof(s);
		GetPrivateProfileString("experiment","gsd.amax","1.0",s,255,iniFile);
		gsd[1]=atof(s);
		GetPrivateProfileString("experiment","gsd.alpha","4.0",s,255,iniFile);
		gsd[2]=atof(s);

		GetPrivateProfileString("experiment","gravitator_threshold_mass","0",s,255,iniFile);
		gRavitatorThresholdMass=atof(s);

		GetPrivateProfileString("experiment","cell_size","1",s,255,iniFile);
		gCellSize=atof(s);

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
		sprintf(iniFile,"%s\\pods.ini",_getcwd(NULL,0));
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
void PrinterDebug()
{
	if (gLander==NULL) return;
	NxReal m=gLander->getMass();
	NxVec3 I=gLander->getMassSpaceInertiaTensor();
	NxMat33 II=gLander->getGlobalInertiaTensor();
	
	return;
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