//////////////////////////////////////////////////////////////////////////////////////////
// Header file for project Verify. This project implements a suite of tests of rigid body
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

#ifndef VERIFY_H
#define VERIFY_H
namespace verify
{
enum {eBAD_EXPERIMENT_TYPE, eCOLLIDER, eINCLINER, eTUMBLER, eTUMBLERS, eBALL_ON_GROUND,
      eSPRINGER_EXPERIMENT} eExperimentType;

PxReal spinMag; // Spin magnitude for target in collider experiment
PxReal kickMag; // Kick magnitude for bullet in collider experiment

PxTransform R1, R0; // throwaway transform holders with static duration
PxVec3 V;           // throwaway vector with static duration

struct {
	unsigned int actorDiag;
	unsigned int systemDiag1;
	unsigned int systemDiag2;
	unsigned int systemDiag3;
} hudMsgs;

struct {
	PxRigidDynamic* bullet;
	PxRigidDynamic* target;
	PxRigidDynamic* tumbler;
	PxRigidDynamic* ball1;
	PxRigidDynamic* ball2;
	vector <PxRigidDynamic*> Tumblers;
	PxRigidDynamic* incliner;
} VIPs;

struct {
	PxVec3 w;
	PxVec3 L_now, L_true;
	PxReal K_now, K_true;
	PxRigidDynamic* handle;
} tumbler;

struct {
	int inc;
	PxRigidDynamic* handle;
} incliner;

struct {
	PxReal littleG;
	PxReal bigG;
} units;

void CreateColliderExperiment();
void CreateTumblerExperiment();
void CreateTumblersExperiment();
void CreateInclinerExperiment();
void CreateBallOnGroundExperiment();
void CreateSpringerExperiment();
void CalcTumblerDynamics();
void InclineGravity(PxReal deg);
void LogTumblerExperiment();
void LogBallOnGroundExperiment();
void LogSpringerExperiment();
void ControlBallOnGroundExperiment();
void ControlSpringerExperiment();
};

#endif