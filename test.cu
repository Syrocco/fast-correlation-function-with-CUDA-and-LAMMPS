#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

__global__ void hello_cuda() {
        printf("hello from GPU\n");
}

int main() {
        printf("hello from CPU\n");
        hello_cuda <<<1, 1>>> ();
        gpuErrchk( cudaPeekAtLastError() );
        gpuErrchk( cudaDeviceSynchronize() );

        gpuErrchk(cudaDeviceReset());
        printf("bye bye from CPU\n");
        return 0;
}