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
    enum {eBAD_EXPERIMENT_TYPE, eMAKE_SGP, eLOAD_SGP, eTEST_SCALING, eORBIT_SGP} eExperimentType;

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
        unsigned int systemDiag6;
    } hudMsgs;

    // Named actors
    struct {
        PxRigidDynamic* nucleus;
        PxRigidDynamic* lBall;
        PxRigidDynamic* rBall;
        PxRigidDynamic* gravitator;
    } VIPs;

    // Parameters for make_sgp experiment
    struct {
        PxReal nucleusRadius; // obsolete
        PxReal mass;
        struct {
            PxReal longAxis;    // the longest semi-principle-axis, called a
            PxReal abAxesRatio; // a/b semi-principle-axes ratio
            PxReal acAxesRatio; // a/c semi-principle-axes ratio
        } ellipsoid;
        struct {
            enum {eBAD_GSD_TYPE, eGSD_UNIFORM, eGSD_IDENTICAL} type;
            PxReal sizeScale;
            PxU32 nbTotal;
        } gsd;
        struct {
            RubbleGrainType shape;
            PxReal density;
        } grain;
    } msgp;

    // Parameters for load_sgp experiment
    struct {
        PxReal remass;
        PxReal rescale;
    } lsgp;

    // Parameters for orbit_sgp experiment
    struct {
        PxReal bigM; // primary mass
        PxReal sgpMass; // overrides loaded sgp
        PxReal sgpRadius; // overrides loaded sgp
        string orbFile;
        nr3::VecDoub tvec;
        nr3::VecDoub xvec;
        nr3::VecDoub yvec;
        PxReal periapse;
        bool bTrackingCamera;
        PxReal tStart;
        PxReal tEnd;
        PxVec3 X0;
        PxVec3 V0;
        int nbOrbits; // elliptical orbits only
    } orbit;

    // Parameters for "scaled integration" test
    struct {
        PxReal radius;
        PxU32 dInitial; // initial separation in radii
        PxReal density;
        int nballs;
    } sclTest;

    // Grain Size Distribution (obsolete now, keeping to not break old code right away)
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

    // Diagnostics
    struct {
        enum {eBAD_CCODE, eNO_CCODE, eGRADIENT, eTWO_LAYER, eTHREE_LAYER, eSURFACE, eNAME} eColorCodeType;
        PxReal nbSurfaceThickness;
    } diag;

    // Namespace functions
    void CreateMakeSGPExperiment();
    void CreateLoadSGPExperiment();
    void CreateOrbitSGPExperiment();
    void CreateTestScalingExperiment();
    void ControlMakeSGPExperiment();
    void ControlLoadSGPExperiment();
    void ControlOrbitSGPExperiment();
    void ControlTestScalingExperiment();
    void LogMakeSGPExperiment();
    void LogOrbitSGPExperiment();
    void LogTestScalingExperiment();
    void RefreshMakeSGPHUD();
    void RefreshLoadSGPHUD();
    void RefreshOrbitSGPHUD();

    bool MakeNewSGP();
    PxU32 MakeLooseRubblePile(PxReal mass=0);
    void GravitateSelf(bool bIgnoreKinematics=false);
    void GravitateOnHost(bool bIgnoreKinematics=false);
    void GravitateOnGPU(); 
    void ApplyTidingForce();
    PxReal SystemPotentialEnergy();
    void ColorCodeRubblePile();
    bool LoadSGP(string filename);
    void AsciizeSGP(string filename);
    void SpyOnSGP(PxReal f=1.0, bool bZoomOutOnly=true);
    PxVec3 FindSGPCenterOfMass();
    PxReal SGPBulkDensity(bool bRoughGuess=true);
    bool GenerateOrbit();
    bool LoadOrbitFromFile();
    void ReMassSGP(PxReal newMass);
    bool RescaleSGP(PxReal factor);
};
#endif