///////////////////////////////////////////////////////////////////////////////
// Header file for project BoxOfRox. 
//
// Author:  Viranga
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

	} params;

	struct {
		enum {eBAD_GSD_TYPE, eGSD_UNIFORM} type;
	} gsd; // grain size distribution

	
};
#endif