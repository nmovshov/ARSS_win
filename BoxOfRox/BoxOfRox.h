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
	enum {eBAD_BOX_DESIGN, eAOSAT1} eBoxDesign;

	struct {
		unsigned int systemDiag;
	} hudMsgs;

	struct {
		PxRigidDynamic* theBox;
	} VIPs;

	struct {
		PxReal boxSize;
		PxReal shakeMagnitude;
	} params;

	struct {
		enum {eBAD_GRAIN_TYPE, eGRAIN_UNIFORM, eGRAIN_BIMODAL} type;
		PxReal size2;
		PxReal size1;
		PxU32 numberRatio;
		PxU32 totalNumber;
	} regolith;

	struct {
		PxReal length;
		PxReal mass;
		PxReal velocity;
		PxReal bigG;
	} units;

	void CreateTheBox();
	void CreateAOSAT1();
	void CreateFillBoxExperiment();
	void GravitateSelf();
	void GravitateOnHost();
	void GravitateOnDevice();

};
#endif