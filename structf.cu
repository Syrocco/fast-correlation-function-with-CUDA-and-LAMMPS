extern "C" { 
#include "parser.h"
}
#include <stdlib.h>
#include <math.h>

FILE* fileout;

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

__global__ void compute_structf_kernel(float* xpos, float* ypos, int nscat, float* qx, int nqx, float* qy, int nqy, float* structf)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i < nqx && j < nqy)
    {   
        
        printf("%d / %d \r", j, nqx);
        float im = 0;
        float re = 0;

        for (int k = 0; k < nscat; ++k)
        {
            float qr = qx[i] * xpos[k] + qy[j] * ypos[k];

            re += cos(qr);
            im += sin(qr);
        }

        structf[i * nqy + j] = (re * re + im * im) / nscat;
    }
}

void compute(Dump* dump, float qmax) {
    char hboxx;
    float Lx = get_boxx(&hboxx, 1, dump);
    char hboxy;
    float Ly = get_boxy(&hboxy, 1, dump);
    float xx = 2*M_PI/Lx;
    float yy = 2*M_PI/Ly;
    int nqx = 2*qmax/xx + 1;
    int nqy = 2*qmax/yy + 1;
    int N = get_natoms(dump);

    float* x = (float*)calloc(N, sizeof(float));
    float* y = (float*)calloc(N, sizeof(float));
    float* qx = (float*)calloc(nqx, sizeof(float));
    float* qy = (float*)calloc(nqy, sizeof(float));
    float* S = (float*)calloc(nqy*nqx, sizeof(float));

    for (int i = 0; i < nqx; ++i){
        qx[i] = xx * ( i - (nqx - 1) / 2);
    }
    for (int i = 0; i < nqy; ++i){
        qy[i] = yy * ( i - (nqy - 1) / 2);
    }
    

    if (a < 0){
        a = 0;
    }
    if (b < 1){
        b = dump->nframes + b;
    }
    if (c < 1){
        c = 1;
    }
    fprintf(fileout, "%d\n", (int)((b - a) / c));
    for (int i = 0; i < nqx; ++i)
        fprintf(fileout, "%g ", qx[i]);
    fprintf(fileout, "\n");

    for (int j = 0; j < nqy; ++j)
        fprintf(fileout, "%g ", qy[j]);
    fprintf(fileout, "\n");

    float* device_x;
    float* device_y;
    float* device_qx;
    float* device_qy;
    float* device_S;
    // Allocate memory on the GPU
    cudaMalloc((void**)&device_x, N * sizeof(float));
    cudaMalloc((void**)&device_y, N * sizeof(float));
    cudaMalloc((void**)&device_qx, nqx * sizeof(float));
    cudaMalloc((void**)&device_qy, nqy * sizeof(float));
    cudaMalloc((void**)&device_S, nqx*nqy * sizeof(float));
    cudaMemcpy(device_qy, qy, nqy * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(device_qx, qx, nqx * sizeof(float), cudaMemcpyHostToDevice);
    

    // Define the number of threads per block and the number of blocks
    dim3 blockDim(32, 32); // Adjust block dimensions as needed
    dim3 gridDim((nqx + blockDim.x - 1) / blockDim.x, (nqy + blockDim.y - 1) / blockDim.y);
    
    
    for (int i = a; i < b; i = i + c){
        
        jump_to_frame(i, dump);

        get_floatatomprop("x", x, N, dump);
        get_floatatomprop("y", y, N, dump);
        

        

        

        // Copy data from host to device
        cudaMemset(device_S, 0, nqy * nqx * sizeof(float));
        cudaMemcpy(device_x, x, N * sizeof(float), cudaMemcpyHostToDevice);
        
        cudaMemcpy(device_y, y, N * sizeof(float), cudaMemcpyHostToDevice);


        compute_structf_kernel<<<gridDim, blockDim>>>(device_x, device_y, N, device_qx, nqx, device_qy, nqy, device_S);

        // Copy the result back to the host
        cudaMemcpy(S, device_S, nqx * nqy * sizeof(float), cudaMemcpyDeviceToHost);

        

        /* Write the structure factor as a nqx * nqy matrix. */
        for (int i = 0; i < nqx; ++i)
        {
            for (int j = 0; j < nqy; ++j)
                fprintf(fileout, "%g ", S[i*nqy + j]);

            fprintf(fileout, "\n");
    	}
    }
    fclose(fileout);
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
