extern "C" { 
#include "parser.h"
}
#include <stdlib.h>
#include <math.h>

typedef struct {   
    float r;
    float g;
} corr;

int a, b, c;
void pars(char* data){
    char input[100];
    char *token;
    strcpy(input, data);

    
    input[strcspn(input, "\n")] = '\0';

  
    token = strtok(input, ":");

    if (token != NULL) {
        a = atoi(token);
        token = strtok(NULL, ":");
        if (token != NULL) {
            b = atoi(token);
            token = strtok(NULL, ":");
            if (token != NULL) {
                c = atoi(token);
            }
            else {
                printf("Invalid input: Not enough values.\n");
            }
        } 
        else {
            printf("Invalid input: Not enough values.\n");
        }
    } 
    else {
        printf("Invalid input: Not enough values.\n");
    }
}

#define floatingType double
//__global__ void histogramKernel(float* x, float* y, float* qa, floatingType* g6Temp, int* count, int N, float halfL, int n) {
__global__ void histogramKernel(float* x, float* y, float* qa, floatingType* g6Temp, int* count, int N, float halfLx, float halfLy, int n) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (tid < N) {
        printf("%d\n", tid);
        
        float xi = x[tid];
        float yi = y[tid];
        float qi = qa[tid];

        
        for (int j = 0; j < tid; j++) {
            
        
            float dx = xi - x[j];
            float dy = yi - y[j];
            if (dx >= halfLx)
                dx -= 2*halfLx;
            else if (dx < -halfLx)
                dx += 2*halfLx;
            if (dy >= halfLy)
                dy -= 2*halfLy;
            else if (dy < -halfLy)
                dy += 2*halfLy;
            float r = sqrt(dx*dx + dy*dy);

    
            if (r < min(halfLx, halfLy)) {
                int idx = (int)(n * (r / min(halfLx, halfLy)));
                
                float qj = qa[j];

                atomicAdd(&g6Temp[idx], (floatingType)2*cos(qj - qi));
                atomicAdd(&count[idx], 2);
            }
        }
    }
}

corr* correlation(Dump* dump, int n) {
    char hboxx;
    char hboxy;
    float Lx = get_boxx(&hboxx, 1, dump);
    float halfLx = 0.5 * Lx;
    float Ly = get_boxy(&hboxy, 1, dump);
    float halfLy = 0.5 * Ly;
    int N = get_natoms(dump);

    float* x = (float*)calloc(N, sizeof(float));
    float* y = (float*)calloc(N, sizeof(float));
    float* qa = (float*)calloc(N, sizeof(float));
    corr* g6 = (corr*)calloc(n, sizeof(corr));

    int loopN = 0;
     if (a < 0){
        a = 0;
    }
    if (b < 1){
        b = dump->nframes;
    }
    if (c < 1){
        c = 1;
    }

    for (int frame = a; frame < b; frame += c) {
        printf("%d/%d\r", frame, b);
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
        histogramKernel<<<numBlocks, threadsPerBlock>>>(device_x, device_y, device_qa, device_g6Temp, device_count, N, halfLx, halfLy, n);

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
        g6[i].r = i*(min(halfLx, halfLy)/n);
        g6[i].g /= loopN;
    }


    return g6;
}


int main(int argc, char** argv){

    float factor = atof(argv[1]);
    char* path_in = argv[2];
    char* path_out = argv[3];
    pars(argv[4]);
    size_t length = strlen(path_in);
    Dump* dump = dump_open(path_in, 'r');
    double sig = 2;
    if ('L' == path_in[length - 1]){
        sig = 0.0025;
    }
    char hboxx;
    int n = get_boxx(&hboxx, 1, dump)/(sig*2*factor);
    corr* g = correlation(dump, n);
    FILE* file;

    file = fopen (path_out, "w+");
    for (int i = 0; i < n; i++)
        fprintf(file, "%f %f\n", g[i].r, g[i].g);
    fclose(file);
}

