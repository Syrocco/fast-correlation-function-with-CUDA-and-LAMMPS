
extern "C" { 
#include "parser.h"
}
#include <stdlib.h>
#include <math.h>

#define floatingType double

typedef struct {   
    float r;
    float g;
} corr;

//__global__ void histogramKernel(float* x, float* y, float* qa, floatingType* g6Temp, int* count, int N, float halfL, int n) {
__global__ void histogramKernel(float* x, float* y, float* qa, floatingType* g6Temp, int* count, int N, float halfL, int n) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (tid < N) {
        printf("%d\n", tid);
        
        float xi = x[tid];
        float yi = y[tid];
        float qi = qa[tid];

        
        for (int j = 0; j < tid; j++) {
            
        
            float dx = xi - x[j];
            float dy = yi - y[j];
            if (dx >= halfL)
                dx -= 2*halfL;
            else if (dx < -halfL)
                dx += 2*halfL;
            if (dy >= halfL)
                dy -= 2*halfL;
            else if (dy < -halfL)
                dy += 2*halfL;
            float r = sqrt(dx*dx + dy*dy);

    
            if (r < halfL) {
                int idx = (int)(n * (r / halfL));
                
                float qj = qa[j];

                atomicAdd(&g6Temp[idx], (floatingType)2*cos(qj - qi));
                atomicAdd(&count[idx], 2);
            }
        }
    }
}

corr* correlation(Dump* dump, int n, int start, int end, int step) {
    char hboxx;
    float L = get_boxx(&hboxx, 1, dump);
    float halfL = 0.5 * L;
    int N = get_natoms(dump);

    float* x = (float*)calloc(N, sizeof(float));
    float* y = (float*)calloc(N, sizeof(float));
    float* qa = (float*)calloc(N, sizeof(float));
    corr* g6 = (corr*)calloc(n, sizeof(corr));

    int loopN = 0;
    for (int frame = start; frame < end; frame += step) {
        printf("%d/%d\r", frame, end);
        fflush(stdout);
        loopN++;
        jump_to_frame(frame, dump);

        get_floatatomprop("x", x, N, dump);
        get_floatatomprop("y", y, N, dump);
        get_floatatomprop("qA6", qa, N, dump);

        floatingType* g6Temp = (floatingType*)calloc(n, sizeof(floatingType));
        int* count = (int*)calloc(n, sizeof(int));

        float* device_x;
        float* device_y;
        float* device_qa;
        floatingType* device_g6Temp;
        int* device_count;

        // Allocate memory on the GPU
        cudaMalloc((void**)&device_x, N * sizeof(float));
        cudaMalloc((void**)&device_y, N * sizeof(float));
        cudaMalloc((void**)&device_qa, N * sizeof(float));
        cudaMalloc((void**)&device_g6Temp, n * sizeof(floatingType));
        cudaMalloc((void**)&device_count, n * sizeof(int));

        // Copy data from host to device
        cudaMemcpy(device_x, x, N * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(device_y, y, N * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(device_qa, qa, N * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemset(device_g6Temp, 0, n * sizeof(floatingType));
        cudaMemset(device_count, 0, n * sizeof(int));

        // Define the number of threads per block and the number of blocks
        int threadsPerBlock = 256;
        int numBlocks = (N + threadsPerBlock - 1) / threadsPerBlock;

        // Launch the CUDA kernel to calculate the histogram
        histogramKernel<<<numBlocks, threadsPerBlock>>>(device_x, device_y, device_qa, device_g6Temp, device_count, N, halfL, n);

        // Copy the histogram result back from the device to the host
        cudaMemcpy(g6Temp, device_g6Temp, n * sizeof(floatingType), cudaMemcpyDeviceToHost);
        cudaMemcpy(count, device_count, n * sizeof(int), cudaMemcpyDeviceToHost);

        for (int i = 0; i < n; i++){
            if (count[i] > 1)
                g6[i].g += g6Temp[i]/count[i];
            g6Temp[i] = 0;
            count[i] = 0;
        }
        // Clean up device memory
        cudaFree(device_x);
        cudaFree(device_y);
        cudaFree(device_qa);
        cudaFree(device_g6Temp);
        cudaFree(device_count);
    }

    for (int i = 0; i < n; i++){
        g6[i].r = i*(halfL/n);
        g6[i].g /= loopN;
    }


    return g6;
}

void print(corr* g, int n){
    FILE* file;

    file = fopen ("data2.txt", "w+");
    for (int i = 0; i < n; i++)
        fprintf(file, "%f %f\n", g[i].r, g[i].g);
}

int main(){
    char* path_in = "/home/syrocco/Documents/correlation function/N_99960dtnoise_0.300res_0.950gamma_0.300T_0.010phi_0.765000rat_1.000vo_3.500ao_1.500delta_0.030Lx_640.704Ly_640.704q_0.000v_1.dump";
    Dump* dump = dump_open(path_in, 'r');
    int nframes = dump->nframes;
    int n = 1000;
    corr* g6 = correlation(dump, n, 0, nframes, 1);
    print(g6, n);
}

