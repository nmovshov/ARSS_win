///////////////////////////////////////////////////////////////////////////////
// Header file for project Verify. This project implements a suite of tests of
// rigid body dynamics in PhysX. Several scenarios can be simulated, that have
// either an analytic solution or at least a conserved quantity that can be measured.
// The available tests are:
// * Slider - a box on an inclined plane, testing the friction model in PhysX
// * Collider - a collision between shapes of varying complexity, testing the collision
// detection and collision resolution in PhysX.
// 
// Usage: An options file specifies the type of experiment to run and the experiment
// parameters.
// 
// Author: Me (Naor)
///////////////////////////////////////////////////////////////////////////////

#ifndef VERIFY_H
#define VERIFY_H
namespace verify
{
enum {eBAD_EXPERIMENT_TYPE, eCOLLIDER, eINCLINER, eTUMBLER, eTUMBLERS} eExperimentType;
PxReal spinMag; // Spin magnitude for target in collider experiment
PxReal kickMag; // Kick magnitude for bullet in collider experiment
PxTransform R1, R0; // throwaway transform holders with static duration
PxVec3 V; // throwaway vector with static duration
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
void CreateColliderExperiment();
void CreateTumblerExperiment();
void CreateTumblersExperiment();
void CreateInclinerExperiment();
void CalcTumblerDynamics();
void InclineGravity(PxReal deg);
};

#endif