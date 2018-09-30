#ifndef DEBUG_C
#define DEBUG_C

#include "multigrid.h"
#include "mex.h"

template<class T>
void multigrid<T>::plotImage(double *I, int level)
{
  mxArray *prhs[1];
  double *image_mx;

  prhs[0]=mxCreateNumericMatrix(sizes[level][0], sizes[level][1], 
				(mxClassID)mxDOUBLE_CLASS, mxREAL);
  
  int GridSize=sizes[level][0]*sizes[level][1];
  
  image_mx=(double*)mxGetData(prhs[0]);

  for (int p=0; p<GridSize; p++) image_mx[p]=I[p];
  
  mexCallMATLAB(0, 0, 1, prhs, "plotImage");

  mxDestroyArray(prhs[0]);

  return;
}

template<class T>
void multigrid<T>::plotImage(vector2D<double> *I, int level, const char *str)
{
  mxArray *prhs[3];
  double *image_x_mx, *image_y_mx;

  prhs[0]=mxCreateNumericMatrix(sizes[level][0], sizes[level][1], 
        (mxClassID)mxDOUBLE_CLASS, mxREAL);
  prhs[1]=mxCreateNumericMatrix(sizes[level][0], sizes[level][1], 
      (mxClassID)mxDOUBLE_CLASS, mxREAL);
  prhs[2]=mxCreateString( str );

  int GridSize=sizes[level][0]*sizes[level][1];
  
  image_x_mx=(double*)mxGetData(prhs[0]);
  image_y_mx=(double*)mxGetData(prhs[1]);

  for (int p=0; p<GridSize; p++) {
    image_x_mx[p]=I[p].x;
    image_y_mx[p]=I[p].y;
  }

  mexCallMATLAB(0, 0, 3, prhs, "plotImage2");

  mxDestroyArray(prhs[0]);

  return;
}

#endif
