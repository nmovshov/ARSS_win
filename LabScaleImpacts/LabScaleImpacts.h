/////////////////////////////////////////////////////////////////////////////////
// Header file for project LabScaleImpacts. This project implements some low
// velocity impacts into regolith in uniform 1g gravity. Idea is of course to
// benchmark against lab tests.
//
// Author: Naor Movshovits (nmovshov at google dot com)
/////////////////////////////////////////////////////////////////////////////////

#ifndef LABSCALE_H
#define LABSCALE_H
namespace labscale
{
// Implement different scenarios with different experiment labels. The dispatcher
// function CreateExperiment will call different routines, basically totally
// different programs that share mostly just the dispatching and logging mechanism.
enum {eBAD_EXPERIMENT_TYPE, eHOLSAPPLE1} eExperimentType;
enum {eBAD_EXP_SUBTYPE, eFILL_BOX} eExperimentSubtype;

// Some throwaway objects with static duration
PxTransform R1, R0; // throwaway transform holders
PxVec3 V;           // throwaway vector

// Really dumb way to implement hud messages; but it works
struct {
    unsigned int actorDiag;
    unsigned int systemDiag1;
    unsigned int systemDiag2;
    unsigned int systemDiag3;
} hudMsgs;

// Named actors
struct {
    PxRigidDynamic* bullet;
    PxRigidDynamic* container;
    PxRigidDynamic* ball1;
} VIPs;

// A regolith container
struct {
    PxReal diameter;
    PxReal fillHeight;
} reg_box;

// Measured and diagnostic quantities

// Units
struct {
    PxReal littleG;
} units;

void CreateHolsapple1Experiment();
void CreateFillBoxExperiment();
void LogHolsapple1Experiment();
void ControlHolsapple1Experiment();
void CreateRegolithContainer();
};

#endif