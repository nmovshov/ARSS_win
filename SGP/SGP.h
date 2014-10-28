///////////////////////////////////////////////////////////////////////////////
// Header file for project SGP. Stands for Self-Gravitating Pile. This project
// is normally used to load a rubble pile from some source, impose some initial
// conditions, and then let it evolve.
//
// Author: Me (Naor)
///////////////////////////////////////////////////////////////////////////////

#ifndef SGP_H
#define SGP_H
namespace sgp
{
	enum {eBAD_EXPERIMENT_TYPE, eMAKE_SGP, eLOAD_SGP} eExperimentType;

	struct {
		unsigned int systemDiag;
		unsigned int nbActors;
	} hudMsgs;

	struct {
		PxReal nucleusRadius;
		struct {
			PxReal abAxesRatio;
			PxReal acAxesRatio;
		} ellipsoid;
	} params;

	struct {
		enum {eBAD_GSD_TYPE, eGSD_UNIFORM, eGSD_BIMODAL} type;
		PxReal size2;
		PxReal size1;
		PxU32  numberRatio; // used for bimodal
		PxU32  totalNumber;
	} gsd; // grain size distribution

	struct {
		PxRigidDynamic* nucleus;
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