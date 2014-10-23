///////////////////////////////////////////////////////////////////////////////
// Header file for project SGP. Stands for Self-Gravitating Pile. This project
// is normally used to load a rubble pile from some source, with some initial conditions,
// and then let it evolve. There are variants for things like impacting the pile
// with something.
//
// Author: Me (Naor)
///////////////////////////////////////////////////////////////////////////////

#ifndef SGP_H
#define SGP_H
namespace sgp
{
	enum {eBAD_EXPERIMENT_TYPE, eMAKE_SGP, eLOAD_SGP, eSHAKE_SGP, eKICK_SGP} eExperimentType;

	struct {
		unsigned int systemDiag;
	} hudMsgs;

	struct {
		PxReal shakeMagnitude;
		PxReal kickMagnitude;
		PxReal nucleusRadius;
	} params;

	struct {
		enum {eBAD_GSD_TYPE, eGSD_UNIFORM, eGSD_BIMODAL} type;
		PxReal size2;
		PxReal size1;
		PxU32  numberRatio;
		PxU32  totalNumber;
	} gsd; // grain size distribution

	struct {
		PxRigidDynamic* nucleus;
		PxRigidDynamic* kicker;
	} VIPs;
	struct {
		PxReal length;
		PxReal mass;
		PxReal velocity;
		PxReal bigG;
	} units;

	void CreateMakeSGPExperiment();
	void CreateLoadSGPExperiment();
	bool MakeNewSGP();
	void GravitateSelf();
	void GravitateOnHost();
	void GravitateOnDevice();
};
#endif