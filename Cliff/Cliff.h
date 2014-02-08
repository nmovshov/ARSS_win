///////////////////////////////////////////////////////////////////////////////
// Header file for project Cliff. This project implements modeling of several variants
// of the cliff collapse problem. A review of the problem description and the expected
// relevant model parameters can be found in Holsapple (2012) and references therein.
//
// Author: Me (Naor)
///////////////////////////////////////////////////////////////////////////////

#ifndef CLIFF_H
#define CLIFF_H
namespace cliff
{
	enum {eBAD_EXPERIMENT_TYPE, eRECT_SMOOTH_BASE_COLLAPSE, eRECT_ROUGH_BASE_COLLAPSE, eRECT_FILL} eExperimentType;

	struct {
		unsigned int systemDiag;
	} hudMsgs;
	
	struct {
		PxRigidStatic* floodGate;
	} VIPs;

	struct {
		PxU32 H0;
		PxU32 L0;
		PxU32 NMax;
	} params;

	void CreateContainment();
	void BuildCliff();
};
#endif