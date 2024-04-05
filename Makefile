
NVCC = nvcc

NVCC_FLAGS = -O3 -arch=sm_61

CUDA_FILES = vStructf.cu structf.cu radialDistrib.cu fqt.cu bondOrientational.cu

all: $(CUDA_FILES:.cu=)

%: %.cu
	$(NVCC) parser.c $(NVCC_FLAGS) $< -o $@



