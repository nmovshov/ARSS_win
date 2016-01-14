#ifndef CUTIDAL_H
#define CUTIDAL_H

// experiment specific functions
void CreateExperiment();
void ReCreateExperiment();
void ResetCamera();
void InitRun(int, char**);
void StopRun(char reason[255]="unknown");
void ConfigExperimentFromFile(char iniFile[255]);
void AutoExperimentControl();
void LogExperimentStats(char logFile[255]);
void AddSpecialForces();
void Gravitate();
void GravitateByPairs();
void MarkClusters();
void eclazz(VecInt_O &nf, Bool equiv(const Int, const Int));
bool inClust(const Int, const Int);
void UpdateConservables();
bool ReadOrbit();
bool CreateInitialPile();
bool PlaceActor(NxVec3 position, NxVec3 velocity, unsigned int actind);

// some types
enum ExperimentType		{BAD_EXPERIMENT_TYPE, FREE_PILE, TIDAL_ENCOUNTER};
enum GravityType		{BAD_GRAVITY_TYPE, ALL_PAIRS};

int main(int argc, char** argv);

#endif