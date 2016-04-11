#include "cuShared.h"
#include <math.h>

#define TILE_WIDTH 64
#define EPS2 0.00125f

// global device variables and functions
float4 *g_devicePositionsArr=NULL;
float4 *g_hostPositionsArr=NULL;
float3 *g_deviceForcesArr=NULL;
float3 *g_hostForcesArr=NULL;

__global__ void GravityKernel(unsigned int bodies, float4 *bodyDescripts, float3 *forcesToUpodate);
__global__ void calculate_forces(float4 *bodyDescripts, float3 *forcesToUpdate, unsigned int N);
__device__ float3 bodyBodyInteraction(float4 bi, float4 bj, float3 ai);

bool AllocateCUDAGlobals(int nbThings)
{
	// Not just allocate. This is a convenient place to create CUDA debugging information and look into capabilities of physical devices
	cudaDeviceProp gpuProps;
	int devToUse=1;
	cudaSetDevice(devToUse);
	cudaGetDeviceProperties(&gpuProps,devToUse);

	printf("Device name=%s\n",gpuProps.name);
	printf("Compute capability=%d.%d\n",gpuProps.major,gpuProps.minor);
	printf("Total global memory=%d kB\n",gpuProps.totalGlobalMem/1024);
	printf("Total shared memory per block=%d kB\n",gpuProps.sharedMemPerBlock/1024);
	printf("Max Grid size=%d x %d x %d\n",gpuProps.maxGridSize[0],gpuProps.maxGridSize[1],gpuProps.maxGridSize[2]);
	printf("Max threads per block=%d x %d x %d\n",gpuProps.maxThreadsDim[0],gpuProps.maxThreadsDim[1],gpuProps.maxThreadsDim[2]);
	printf("Max threads per block=%d\n",gpuProps.maxThreadsPerBlock);
	printf("Multi-processors count=%d\n",gpuProps.multiProcessorCount);
	printf("Kernel execution may time out=%d\n",gpuProps.kernelExecTimeoutEnabled);
	printf("Compute mode=%d\n",gpuProps.computeMode);
	printf("ECC enabled=%d\n",gpuProps.ECCEnabled);

	// allocate memory for device and host arrays of actor positions and forces
	cudaError_t err;
	size_t nbytes4=nbThings * sizeof(float) * 4;
	size_t nbytes3=nbThings * sizeof(float) * 3;
	err = cudaMalloc((void **) &g_devicePositionsArr,nbytes4);
	if (err != cudaSuccess) {
		printf("Error in cudaMalloc()...%s\a\n",cudaGetErrorString(err));
		exit(1);
	}
	err = cudaMalloc((void **) &g_deviceForcesArr,nbytes3);
	if (err != cudaSuccess) {
		printf("Error in cudaMalloc()...%s\a\n",cudaGetErrorString(err));
		exit(1);
	}

	g_hostPositionsArr = (float4 *) malloc(nbytes4);
	g_hostForcesArr = (float3 *) malloc(nbytes3);

	//err = cudaFuncSetCacheConfig(GravityKernel,cudaFuncCachePreferL1);
	//err = cudaFuncSetCacheConfig(calculate_forces,cudaFuncCachePreferShared);

	return true;
}

void ReleaseCUDA()
{
	cudaError_t err;
	err=cudaFree(g_devicePositionsArr);
	err=cudaFree(g_deviceForcesArr);
	free(g_hostPositionsArr);
	free(g_hostForcesArr);
}

void GravitateOnDevice(unsigned int bodies, unsigned int mode)
	// mode controls the async behavior. mode=1 will set up and launch the gravity kernel and then return. mode=2 will complete the gravity
	// kernel and retrieve the forces. mode=3 will do everything.
{
	cudaError err;
	size_t nbytes3 = bodies * sizeof(float) * 3;
	size_t nbytes4 = bodies * sizeof(float) * 4;

	if (mode==1 || mode==3)
	{
		// copy the body descriptions to device global memory
		err = cudaMemcpy(g_devicePositionsArr,g_hostPositionsArr,nbytes4,cudaMemcpyHostToDevice);

		if (err!=cudaSuccess){
			printf("cuda memcpy failed!\n");
			exit(1);
		}

		// reset device forces
		err = cudaMemset(g_deviceForcesArr,0,nbytes3);

		// launch a device force calculation
		int nt=(TILE_WIDTH<bodies)?TILE_WIDTH:bodies;
		int nb = (int) ceil((float)bodies/(float)nt); // recommended to use a good multiple of 32 for bodies
		dim3 grid(nb,1,1);
		dim3 block(nt,1,1);


		calculate_forces<<<grid,block>>>(g_devicePositionsArr,g_deviceForcesArr,bodies);	
		//GravityKernel<<<grid,block>>>(bodies, g_devicePositionsArr, g_deviceForcesArr);
		err=cudaGetLastError();
		if (err!=cudaSuccess)
		{
			printf("kernel launch failed!\a\n");
			exit(1);
		}
	}

	if (mode==2 || mode==3)
	{
		// copy updated forces from device to host
		err = cudaDeviceSynchronize();
		err = cudaMemcpy(g_hostForcesArr,g_deviceForcesArr,nbytes3,cudaMemcpyDeviceToHost);
		if (err!=cudaSuccess){
			printf("cuda memcpy failed!\n");
			exit(1);
		}
	}

	return;
}

__global__ void GravityKernel(unsigned int bodies, float4 *bodyDescripts, float3 *forcesToUpdate)
{
	int tid = (blockIdx.x * blockDim.x) + threadIdx.x;
	if (tid<bodies)
	{
		float4 rp = bodyDescripts[tid];
		float4 r;
		float distSqr, invDist, invDistCube, s;

		for (int k=0; k<bodies; k++)
		{
			if (k==tid) continue;
			r = bodyDescripts[k];
			r.x -= rp.x;
			r.y -= rp.y;
			r.z -= rp.z;

			distSqr = r.x*r.x + r.y*r.y + r.z*r.z;
			invDist = 1.0f/sqrtf(distSqr);
			invDistCube = invDist*invDist*invDist;

			s = r.w * rp.w * invDistCube; // universal big G added outside

			forcesToUpdate[tid].x += r.x * s;
			forcesToUpdate[tid].y += r.y * s;
			forcesToUpdate[tid].z += r.z * s;
		}
	}
	return;
}

__device__ float3 bodyBodyInteraction(float4 bi, float4 bj, float3 ai)
{
	float3 r;
	// r_ij [3 flops]
	r.x = bj.x - bi.x;
	r.y = bj.y - bi.y;
	r.z = bj.z - bi.z;

	// distSqr = dot(r_ij, r_ij) + EPS^2  [6 FLOPS]
	float distSqr = r.x*r.x + r.y*r.y + r.z*r.z + EPS2;

	// invDistCube =1/distSqr^(3/2)  [4 FLOPS (2 mul, 1 sqrt, 1 inv)]
	float distSixth = distSqr * distSqr * distSqr;
	float invDistCube = 1.0f/ sqrtf (distSixth);
	
	// s = m_j * invDistCube [1 FLOP]
	float s = bi.w * bj.w * invDistCube; // universal big G added outside
	
	// a_i =  a_i + s * r_ij [6 FLOPS]
	ai.x += r.x * s;
	ai.y += r.y * s;
	ai.z += r.z * s;
	
	return ai;
}

__global__ void calculate_forces(float4 *bodyDescripts, float3 *forcesToUpdate, unsigned int N)
{
	__shared__ float4 shPosition[TILE_WIDTH];

	float4 myPosition;
	int j, tile;
	float3 acc = {0.0f, 0.0f, 0.0f};
	int gtid = blockIdx.x * blockDim.x + threadIdx.x;

	myPosition=bodyDescripts[gtid];

	for (j = 0, tile = 0; j < N; j+=TILE_WIDTH, tile++) {
		int idx = tile * blockDim.x + threadIdx.x;
		shPosition[threadIdx.x] = bodyDescripts[idx];
		__syncthreads();
		for (int k=0; k<TILE_WIDTH; k++)
		{
			acc=bodyBodyInteraction(myPosition,shPosition[k],acc);
		}
		__syncthreads();
	}
	// Save the result in global memory for the integration step.
	//float4 acc4 = {acc.x, acc.y, acc.z, 0.0f};
	//globalA[gtid] = acc4;
	forcesToUpdate[gtid]=acc;
}