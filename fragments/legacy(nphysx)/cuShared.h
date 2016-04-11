#ifndef CUSHARED_H
#define CUSHARED_H
// ===============================================================================
// Global variables and functions required for GPU operation
// ===============================================================================
#include <stdio.h>
#include "cuda.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// global device and device related variables
#define MAX_EXPECTED_ACTORS 10000
extern cudaDeviceProp gpuProps;

extern float4 *g_devicePositionsArr;
extern float4 *g_hostPositionsArr;
extern float3 *g_deviceForcesArr;
extern float3 *g_hostForcesArr;


// functions
bool AllocateCUDAGlobals(int nbThings);
void ReleaseCUDA();
void GravitateOnDevice(unsigned int, unsigned int mode=3);

#endif
