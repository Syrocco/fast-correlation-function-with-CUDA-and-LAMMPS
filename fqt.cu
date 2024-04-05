extern "C" { 
#include "parser.h"
}
#include <stdlib.h>
#include <math.h>
FILE* fileout;

int index(int i, int j, int k, int ntime, int nq){
    return i + k*nq*ntime + j*ntime;
}
__device__ int index2(int i, int j, int k, int ntime, int nq){
    return i + k*nq*ntime + j*ntime;
}
int a, b, c;
void pars(char* data){
    char input[100];
    char *token;
    strcpy(input, data);

    // Remove the newline character from the input
    input[strcspn(input, "\n")] = '\0';

    // Use strtok to split the string
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

__global__ void compute_structf_kernel(float* xpos, float* ypos, int nscat, float* qx, int nqx, float* qy, int nqy, float* structfR, float* structfI)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < nqx && j < nqy)
    {   
        
        
        float im = 0;
        float re = 0;

        for (int k = 0; k < nscat; ++k)
        {
            float qr = qx[i] * xpos[k] + qy[j] * ypos[k];

            re += cos(qr);
            im += sin(qr);
        }

        structfR[i * nqy + j] = re;
        structfI[i * nqy + j] = im;
    }
}

__global__ void calculateFQT(float *FQT, float *Srt, float *Sit, int timeSize, int nq) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < timeSize && j < timeSize && j >= i) {
        printf("computing FQT: %d / %d  and %d / %d\r", i, timeSize, j, timeSize);
        for (int xx = 0; xx < nq; xx++) {
            for (int yy = 0; yy < nq; yy++) {
                
                atomicAdd(&FQT[index2(j - i, xx, yy, timeSize, nq)],
                          Srt[index2(i, xx, yy, timeSize, nq)] * Srt[index2(j, xx, yy, timeSize, nq)] +
                          Sit[index2(i, xx, yy, timeSize, nq)] * Sit[index2(j, xx, yy, timeSize, nq)]);
            }
        }
    }
}

void compute(Dump* dump, float qmax) {
    char hboxx;
    float L = get_boxx(&hboxx, 1, dump);
    float xx = 2*M_PI/L;
    int nq = 2*qmax/xx + 1;
    int N = get_natoms(dump);
    
    float* x = (float*)calloc(N, sizeof(float));
    float* y = (float*)calloc(N, sizeof(float));
    float* q = (float*)calloc(nq, sizeof(float));
    float* Si = (float*)calloc(nq*nq, sizeof(float));
    float* Sr = (float*)calloc(nq*nq, sizeof(float));

    for (int i = 0; i < nq; ++i){
		q[i] = xx * (i - (nq - 1)/2);
    }
    

    if (a < 0){
        a = 0;
    }
    if (b < 1){
        b = dump->nframes;
    }
    if (c < 1){
        c = 1;
    }
    int timeSize = (int)((b - a) / c);
    float* Sit = (float*)calloc(timeSize*nq*nq, sizeof(float));
    float* Srt = (float*)calloc(timeSize*nq*nq, sizeof(float));
    float* FQT = (float*)calloc(timeSize*nq*nq, sizeof(float));
    int* count = (int*)calloc(timeSize, sizeof(int));

    fprintf(fileout, "%d\n", timeSize);
    for (int i = 0; i < nq; ++i)
        fprintf(fileout, "%g ", q[i]);
    fprintf(fileout, "\n");

    jump_to_frame(0, dump);
    double tinit = get_timestep(dump);
    for (int j = 0; j < dump->nframes; ++j){
        jump_to_frame(j, dump);
        fprintf(fileout, "%lf ", get_timestep(dump) - tinit);
    }
    fprintf(fileout, "\n");

    float* device_x;
    float* device_y;
    float* device_q;
    float* device_Si;
    float* device_Sr;
    // Allocate memory on the GPU
    cudaMalloc((void**)&device_x, N * sizeof(float));
    cudaMalloc((void**)&device_y, N * sizeof(float));
    cudaMalloc((void**)&device_q, nq * sizeof(float));
    cudaMalloc((void**)&device_Sr, nq*nq * sizeof(float));
    cudaMalloc((void**)&device_Si, nq*nq * sizeof(float));
    cudaMemcpy(device_q, q, nq * sizeof(float), cudaMemcpyHostToDevice);
    

    // Define the number of threads per block and the number of blocks
    dim3 blockDim(32, 32); // Adjust block dimensions as needed
    dim3 gridDim((nq + blockDim.x - 1) / blockDim.x, (nq + blockDim.y - 1) / blockDim.y);
    
    int cc = 0;
    for (int i = a; i < b; i = i + c){
        printf("computing structure factor: %d / %d \r", i, b);
        jump_to_frame(i, dump);

        get_floatatomprop("x", x, N, dump);
        get_floatatomprop("y", y, N, dump);
        

        

        

        // Copy data from host to device
        cudaMemset(device_Sr, 0, nq * nq * sizeof(float));
        cudaMemset(device_Si, 0, nq * nq * sizeof(float));
        cudaMemcpy(device_x, x, N * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(device_y, y, N * sizeof(float), cudaMemcpyHostToDevice);


        compute_structf_kernel<<<gridDim, blockDim>>>(device_x, device_y, N, device_q, nq, device_q, nq, device_Sr, device_Si);

        // Copy the result back to the host
        cudaMemcpy(Sr, device_Sr, nq * nq * sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(Si, device_Si, nq * nq * sizeof(float), cudaMemcpyDeviceToHost);
        
        for (int l = 0; l < nq; ++l)
        {
            for (int j = 0; j < nq; ++j){
                Sit[index(cc, l, j, timeSize, nq)] = Si[l*nq + j];
                Srt[index(cc, l, j, timeSize, nq)] = Sr[l*nq + j];
            }   
    	}
        cc++;
    }
    printf("\nDONE!\n");

    float *d_FQT, *d_Srt, *d_Sit;


    // Allocate device memory
    cudaMalloc((void **)&d_FQT, timeSize * nq * nq * sizeof(float));
    cudaMalloc((void **)&d_Srt, timeSize * nq * nq * sizeof(float));
    cudaMalloc((void **)&d_Sit, timeSize * nq * nq * sizeof(float));


    cudaMemset(d_FQT, 0, timeSize * nq * nq * sizeof(float));


    cudaMemcpy(d_Srt, Srt, timeSize * nq * nq* sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Sit, Sit, timeSize * nq * nq* sizeof(float), cudaMemcpyHostToDevice);

    dim3 blockDimN(32, 32); // Adjust block dimensions as needed
    dim3 gridDimN((timeSize + blockDimN.x - 1) / blockDimN.x, (timeSize + blockDimN.y - 1) / blockDimN.y);
    calculateFQT<<<gridDimN, blockDimN>>>(d_FQT, d_Srt, d_Sit, timeSize, nq);

    // Copy results back to host
    cudaMemcpy(FQT, d_FQT, timeSize* nq * nq * sizeof(float), cudaMemcpyDeviceToHost);


    for (int t = 0; t < timeSize; t++){
        for (int i = 0; i < nq; ++i){
            for (int j = 0; j < nq; ++j){
                
                fprintf(fileout, "%g ", FQT[index(t, i, j, timeSize, nq)]/(N*(timeSize - t)));
            }
            fprintf(fileout, "\n");
        }
    }
    fclose(fileout);
    printf("\nDONE!\n");
}


int main(int argc, char** argv){
    float qmax = atof(argv[1]);
    char* path_in = argv[2];
    char* path_out = argv[3];
    pars(argv[4]);
    fileout = fopen(path_out, "w");
    Dump* dump = dump_open(path_in, 'r');
    compute(dump, qmax);
}