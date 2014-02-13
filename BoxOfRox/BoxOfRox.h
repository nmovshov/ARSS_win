///////////////////////////////////////////////////////////////////////////////
// Header file for project BoxOfRox. 
//
// Authors:  Viranga and Naor
///////////////////////////////////////////////////////////////////////////////

#ifndef BOXOFROX_H
#define BOXOFROX_H
namespace rox
{
	enum {eBAD_EXPERIMENT_TYPE, eFILL_BOX} eExperimentType;

	struct {
		unsigned int systemDiag;
	} hudMsgs;

	struct {
		PxRigidDynamic* theBox;
	} VIPs;

	struct {
		PxReal boxSize;
	} params;

	struct {
		enum {eBAD_GSD_TYPE, eGSD_UNIFORM} type;
	} gsd; // grain size distribution

	void CreateContainment();
	void CreateFillBoxExperiment();
};
#endif