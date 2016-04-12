///////////////////////////////////////////////////////////////////////////////
// Global variables and functions required for GPU operation
///////////////////////////////////////////////////////////////////////////////
#ifndef CUGRAV_H
#define CUGRAV_H
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// global device variables
extern float4 *g_devicePositionsArr;
extern float4 *g_hostPositionsArr;
extern float3 *g_deviceForcesArr;
extern float3 *g_hostForcesArr;

// global device functions functions
bool AllocateCUDAGlobals(int nbThings);
void ReleaseCUDA();
void GravitateOnDevice(unsigned int, unsigned int mode=3);

#endif
