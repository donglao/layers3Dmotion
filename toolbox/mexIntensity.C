#include "mex.h"
#include "matrix.h"
#include "stdio.h"
#include "string.h"
#include "stdlib.h"

#define min(a,b)(((a)<(b))?(a):(b))
#define max(a,b)(((a)>(b))?(a):(b))


typedef double ARRAYTYPE;

// [intensity_likelihood] = mexIntensity(image, region, edge, threshold, size_a, size_b)
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    ARRAYTYPE *image, *edge;
    bool *region;
    int XSize, YSize, size_a, size_b, edge_count;
    ARRAYTYPE threshold;
    const mwSize *size  = mxGetDimensions( prhs[0] );
    XSize    = size[0];
    YSize    = size[1];
    
    image     = (ARRAYTYPE*)mxGetData(prhs[0]);
    region  = mxGetLogicals(prhs[1]);
    edge = (ARRAYTYPE*)mxGetData(prhs[2]);
    const mwSize *temp_size2  = mxGetDimensions( prhs[2] );
    edge_count = temp_size2[0];
    threshold = (double)mxGetScalar(prhs[3]);
    size_a               =    (int)mxGetScalar(prhs[4]);
    size_b              =    (int)mxGetScalar(prhs[5]);
    
    
    
    
    mwSize dims[2];
    dims[0]=XSize; dims[1]=YSize;
    plhs[0] = mxCreateNumericArray( (mwSize)2, dims, mxDOUBLE_CLASS, mxREAL );
    
    
    double *intensity_seg   = (double*)mxGetData(plhs[0]);

    
    int a,b,c,d, similar, total;
    double difference;
    double temp_f;
    int ind, ind_org, i, j, m;
    
   
    for (int m = 0; m < edge_count; m++)
    {   
        i = (int)edge[m]-1;
        j = (int)edge[m + edge_count]-1;

                similar = 0;
                total = 0;
                a = max(i-size_a,0);
                b = min(i+size_a,XSize-1);
                c = max(j-size_b,0);
                d = min(j+size_b,YSize-1);
                for ( int k = a; k < b+1; k++)
                    for (int l = c; l < d+1; l++)
                    {
                        if ( region[k + l*XSize] )
                        {
                            total += 1;
                            difference = 0;
                            for (int tempi = 0; tempi < 3; tempi++)
                            {
                                ind = k + l*XSize + tempi*XSize*YSize;
                                ind_org = i + j*XSize + tempi*XSize*YSize;
                                
                                temp_f = (image[ind] - image[ind_org]);
//                                 printf("%f\n", temp_f);
                                difference += temp_f*temp_f;
                            }
                            if ( difference < threshold)
                                similar += 1;
                            
                        }
                    }
                
                intensity_seg[i + j*XSize] = -(double)(similar)/(double)(total+20);
    }
    
    return;
}
