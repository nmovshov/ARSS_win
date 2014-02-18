///////////////////////////////////////////////////////////////////////////////
// Header file for project BoxOfRox. 
//
// Authors:  Viranga and Naor
///////////////////////////////////////////////////////////////////////////////

#ifndef BOXOFROX_H
#define BOXOFROX_H
namespace rox
{
	enum {eBAD_EXPERIMENT_TYPE, eFILL_BOX, eSHAKE_BOX} eExperimentType;

	struct {
		unsigned int systemDiag;
	} hudMsgs;

	struct {
		PxRigidDynamic* theBox;
		PxRigidDynamic* nucleus;
	} VIPs;

	struct {
		PxReal boxSize;
		PxReal shakeMagnitude;
		PxReal nucleusRadius;
	} params;

	struct {
		enum {eBAD_GRAIN_TYPE, eGRAIN_UNIFORM, eGRAIN_BIMODAL} type;
		PxReal size2;
		PxReal size1;
		PxU32 numberRatio;
		PxU32 totalNumber;
	} grain; // grain size distribution

	struct {
		PxReal length;
		PxReal mass;
		PxReal velocity;
		PxReal bigG;
	} units;

	void CreateContainment();
	void CreateFillBoxExperiment();
	void GravitateSelf();
	void GravitateOnHost();
	void GravitateOnDevice();

};
#endif