#include <cstdio>
#include <cstdlib>
#include <vector>
#include <omp.h>

int main() {
  omp_set_num_threads(8);
  int n = 50;
  int range = 5;
  std::vector<int> key(n);
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  std::vector<int> bucket(range),vec(range); 
#pragma omp parallel for shared(bucket)
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
  for (int i=0; i<n; i++) {
#pragma omp atomic update
    bucket[key[i]]++;
  }
/*
  for (int i=0, j=0; i<range; i++) {
    for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
    }
  }
*/


#pragma omp parallel
  for(int j=0;j<range;j++){
#pragma omp for 
    for(int i=0;i<j;i++){
#pragma omp atomic update
      vec[j] += bucket[i];
    }
#pragma omp for
    for(int i=vec[j];i<vec[j]+bucket[j];i++){
      key[i] = j;
    }
  }

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}
