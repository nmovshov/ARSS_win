/////////////////////////////////////////////////////////////////////////////////
// Header file for project SGP. Stands for Self-Gravitating Pile. This project is
// normally used to load a rubble pile from some source, impose some initial
// conditions, and then let it evolve.
//
// Author: Naor Movshovitz (nmovshov at google dot com)
/////////////////////////////////////////////////////////////////////////////////

#ifndef SGP_H
#define SGP_H
namespace sgp
{
    // Implement different scenarios with different experiment labels. The dispatcher
    // function CreateExperiment will call different routines, basically totally
    // different programs that share mostly just the dispatching and logging mechanism.
    enum {eBAD_EXPERIMENT_TYPE, eMAKE_SGP, eLOAD_SGP, eTEST_SCALING} eExperimentType;

    // Some throwaway objects with static duration
    PxTransform R1, R0; // throwaway transform holders
    PxVec3 V;           // throwaway vector

    // Really dumb way to implement hud messages; but it works
    struct {
        unsigned int systemDiag1;
        unsigned int systemDiag2;
        unsigned int systemDiag3;
        unsigned int systemDiag4;
        unsigned int systemDiag5;
    } hudMsgs;

    // Named actors
    struct {
        PxRigidDynamic* nucleus;
        PxRigidDynamic* lBall;
        PxRigidDynamic* rBall;
    } VIPs;

    // Misc params
    struct {
        PxReal nucleusRadius;
        struct {
            PxReal abAxesRatio;
            PxReal acAxesRatio;
        } ellipsoid;
        PxReal defaultDensity;
    } params;

    // Grain Size Distribution
    struct {
        enum {eBAD_GSD_TYPE, eGSD_UNIFORM, eGSD_BIMODAL} type;
        PxReal size2;
        PxReal size1;
        PxU32  numberRatio; // used for bimodal
        PxU32  totalNumber;
    } gsd;

    // Code units
    struct {
        PxReal length; // cu length in meters
        PxReal mass; // cu mass in kg
        PxReal time; // cu time in seconds
        PxReal bigG;
    } cunits;

    // Namespace functions
    void CreateMakeSGPExperiment();
    void CreateLoadSGPExperiment();
    void CreateTestScalingExperiment();
    bool MakeNewSGP();
    void GravitateSelf();
    void GravitateOnHost();
    void GravitateOnDevice();
    void RefreshMakeSGPHUD();
    void LogMakeSGPExperiment();
    void LogTestScalingExperiment();
    void ControlTestScalingExperiment();
    PxReal SystemPotentialEnergy();
};
#endif