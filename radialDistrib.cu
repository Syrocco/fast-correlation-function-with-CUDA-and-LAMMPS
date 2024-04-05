
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

#define floatingType double


__global__ void histogramKernel(float* x, float* y, floatingType* count, int N, float halfL, int n) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    float qx = 0;
    //float qy = 3.33503;
    float qy = 2691.98;
    if (tid < N) {
        printf("%d / %d \r", tid, N);
        float xi = x[tid];
        float yi = y[tid];

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
            float r = sqrt(dx * dx + dy * dy);

            if (r < halfL) {
                int idx = (int)(n * (r / halfL));
                atomicAdd(&count[idx], (floatingType)2*cos(qx*dx + qy*dy));
            }
        }
    }
}

corr* correlation(Dump* dump, int n) {
    char hboxx;
    float L = get_boxx(&hboxx, 1, dump);
    float halfL = 0.5 * L;
    int N = get_natoms(dump);

    float* x = (float*)calloc(N, sizeof(float));
    float* y = (float*)calloc(N, sizeof(float));
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

    float* device_x;
    float* device_y;
    floatingType* device_count;
    // Allocate memory on the GPU
    cudaMalloc((void**)&device_x, N * sizeof(float));
    cudaMalloc((void**)&device_y, N * sizeof(float));
    cudaMalloc((void**)&device_count, n * sizeof(floatingType));
    floatingType* count = (floatingType*)calloc(n, sizeof(floatingType));
    for (int frame = a; frame < b; frame += c) {
        printf("%d/%d\r", a, b);
        fflush(stdout);
        loopN++;
        jump_to_frame(frame, dump);

        get_floatatomprop("x", x, N, dump);
        get_floatatomprop("y", y, N, dump);

        

        // Copy data from host to device
        cudaMemcpy(device_x, x, N * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(device_y, y, N * sizeof(float), cudaMemcpyHostToDevice);
        cudaMemset(device_count, 0, n * sizeof(floatingType));

        // Define the number of threads per block and the number of blocks
        int threadsPerBlock = 256;
        int numBlocks = (N + threadsPerBlock - 1)/threadsPerBlock;

        // Launch the CUDA kernel to calculate the histogram
        histogramKernel<<<numBlocks, threadsPerBlock>>>(device_x, device_y, device_count, N, halfL, n);

        // Copy the histogram result back from the device to the host
        cudaMemcpy(count, device_count, n * sizeof(floatingType), cudaMemcpyDeviceToHost);

        for (int i = 0; i < n; i++){
            
                //printf("%lf ", g6Temp[i]);
            g6[i].g += count[i];
            count[i] = 0;
        }


    }
    // Clean up device memory
    cudaFree(device_x);
    cudaFree(device_y);
    cudaFree(device_count);
    free(count);
    for (int i = 0; i < n; i++){
        g6[i].r = (i + 0.5)*(halfL/n);
        g6[i].g /= loopN;
        g6[i].g = L*L*g6[i].g/(2*M_PI*N*(N - 1)*g6[i].r*(g6[1].r - g6[0].r));
        printf("%lf ", g6[i].g);
    }


    return g6;
}

void print(corr* g, int n){

}

int main(int argc, char** argv){
    float factor = atof(argv[1]);
    char* path_in = argv[2];
    char* path_out = argv[3];
    pars(argv[4]);
    size_t length = strlen(path_in);
    Dump* dump = dump_open(path_in, 'r');
    int nframes = dump->nframes;
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





