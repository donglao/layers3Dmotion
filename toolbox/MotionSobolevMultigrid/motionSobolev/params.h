#ifndef PARAMS_H
#define PARAMS_H

// #define WRITE_PRINT_STATEMENTS

#define EPS 1e-6
//#define TRANSLATION_EPS 1e-3
#define TRANSLATION_EPS 0.5

#define CG_ERR_TOL 1e-4
//#define CG_ERR_TOL 5e-2
#define MAX_ENERGY_DIFF_MOTION 1e-4
// number of blocks for parallelization
#define Nblocks 16

#define BLOCKS(block) int block=0; block<Nblocks; block++
#define PIXELS(p) int p=start; p<end; p++ 
#define BLOCK_ENDS(N)                                         \
  int start = block*(N)/Nblocks;                              \
  int end   = block==Nblocks-1 ? (N) : (block+1)*(N)/Nblocks

#endif
