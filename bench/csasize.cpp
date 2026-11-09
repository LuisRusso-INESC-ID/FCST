#include <sdsl/suffix_arrays.hpp>
#include <iostream>
using namespace sdsl;
int main(int, char** argv){
  csa_wt<> a; construct(a, argv[1], 1);
  csa_wt<wt_huff<rrr_vector<63>>,32,64> b; construct(b, argv[1], 1);
  csa_sada<> c; construct(c, argv[1], 1);
  double n=a.size();
  printf("%s csa_wt(huff,32/64)=%.2f csa_wt(huff-rrr,32/64)=%.2f csa_sada=%.2f bits/char\n", argv[1], size_in_bytes(a)*8/n, size_in_bytes(b)*8/n, size_in_bytes(c)*8/n);
}
