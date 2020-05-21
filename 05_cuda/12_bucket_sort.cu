#include <cstdio>
#include <cstdlib>
#include <vector>
using namespace std;
__global__ void BucketSort(int * bucket,int *key,int N,int range){
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  bucket[i%5] = 0;
  __syncthreads();
  if(key[i] == 0){
    atomicAdd(&bucket[0],1);
  }else if(key[i] == 1){
    atomicAdd(&bucket[1],1);
  }else if(key[i] == 2){
    atomicAdd(&bucket[2],1);
  }else if(key[i] == 3){
    atomicAdd(&bucket[3],1);
  }else{
    atomicAdd(&bucket[4],1);
  }

  if(i < bucket[0]){
    key[i] = 0;
  }else if(i < bucket[0]+bucket[1]){
    key[i] = 1;
  }else if(i < bucket[0]+bucket[1]+bucket[2]){
    key[i] = 2;
  }else if(i < bucket[0]+bucket[1]+bucket[2]+bucket[3]){
    key[i] = 3;
  }else if(i < N){
    key[i] = 4;
  }

}

int main() {
  int n = 50;
  int range = 5;
  int *key;
  cudaMallocManaged(&key,n*sizeof(int));
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  int *bucket;
  cudaMallocManaged(&bucket,range*sizeof(int));
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
  BucketSort<<<1,n>>>(bucket,key,n,range);
  cudaDeviceSynchronize();
  cudaFree(bucket);

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
  cudaFree(key);
}
