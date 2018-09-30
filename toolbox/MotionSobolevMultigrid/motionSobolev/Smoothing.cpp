#include "Smoothing.h"
#include "math.h"
#include <stdio.h>
#include <string.h>

#define T( I, r, c, d, nRow, nCol )     I[ (r) + (c) * nRow + (d) * nRow * nCol ]
#define ROW( nIndex, nRow, nCol )       (nIndex % nRow)
#define COL( nIndex, nRow, nCol )       ( ((nIndex - (nIndex % nRow)) / nRow) % nCol )

#ifndef MAX
    #define MAX( a, b ) ( ( (a) > (b) ) ? (a) : (b) )
#endif

#ifndef MIN
    #define MIN( a, b ) ( ( (a) < (b) ) ? (a) : (b) )
#endif


CSmoothing::CSmoothing()
{

}

CSmoothing::~CSmoothing()
{

}

void    CSmoothing::Blur( double *pImageBlur, const double *pImage, int nRow, int nCol, int nVec, double dSigma )
{
    if( dSigma == 0 ) {

        memcpy( pImageBlur, pImage, sizeof(double) * nRow * nCol * nVec );
        return;
    }

    int     nSizeRow = floor( dSigma * 3 );
    int     nSizeCol = floor( dSigma * 3 );
    int     hr, hc;
    double  hg, sum_hg = 0;

    int     nRowKernel  = 2 * nSizeRow + 1;
    int     nColKernel  = 2 * nSizeCol + 1;
    double  *pKernel    = new double [ nRowKernel * nColKernel ];

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hr = r - nSizeRow;
            hc = c - nSizeCol;

            hg = exp( - ((hr*hr) + (hc*hc)) / (2*dSigma*dSigma) );

            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg;
            sum_hg += hg;
        }
    }

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hg = T( pKernel, r, c, 0, nRowKernel, nColKernel );
            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg / sum_hg;
        }
    }

    double  dI, dK;
    double  dSumKernel = 0, dSumBlur = 0;
    int     nTop, nBottom, nLeft, nRight;
    int     kr, kc;

    for( int v = 0; v < nVec; v++ ) {
        for( int r = 0; r < nRow; r++ ) {
            for( int c = 0; c < nCol; c++ ) {

                dSumBlur    = 0;
                dSumKernel  = 0;

                nTop        = MAX( 0, r-nSizeRow );
                nBottom     = MIN( nRow-1, r+nSizeRow );
                nLeft       = MAX( 0, c-nSizeCol );
                nRight      = MIN( nCol-1, c+nSizeCol );

                for( int rr = nTop; rr <= nBottom; rr++ ) {
                    for( int cc = nLeft; cc <= nRight; cc++ ) {

                        kr = rr - r + nSizeRow;
                        kc = cc - c + nSizeCol;

                        dI = T( pImage,rr,cc,v,nRow,nCol );
                        dK = T( pKernel,kr,kc,0,nRowKernel,nColKernel );

                        dSumBlur += dI * dK;
                        dSumKernel  += dK;
                    }
                }

                T( pImageBlur,r,c,v,nRow,nCol ) = dSumBlur / dSumKernel;
            }
        }
    }
}

void    CSmoothing::Blur( float *pImageBlur, const float *pImage, int nRow, int nCol, int nVec, float dSigma )
{
    if( dSigma == 0 ) {

        memcpy( pImageBlur, pImage, sizeof(float) * nRow * nCol * nVec );
        return;
    }

    int     nSizeRow = floor( dSigma * 3 );
    int     nSizeCol = floor( dSigma * 3 );
    int     hr, hc;
    float  hg, sum_hg = 0;

    int     nRowKernel  = 2 * nSizeRow + 1;
    int     nColKernel  = 2 * nSizeCol + 1;
    float  *pKernel    = new float [ nRowKernel * nColKernel ];

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hr = r - nSizeRow;
            hc = c - nSizeCol;

            hg = exp( - ((hr*hr) + (hc*hc)) / (2*dSigma*dSigma) );

            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg;
            sum_hg += hg;
        }
    }

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hg = T( pKernel, r, c, 0, nRowKernel, nColKernel );
            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg / sum_hg;
        }
    }

    float  dI, dK;
    float  dSumKernel = 0, dSumBlur = 0;
    int     nTop, nBottom, nLeft, nRight;
    int     kr, kc;

    for( int v = 0; v < nVec; v++ ) {
        for( int r = 0; r < nRow; r++ ) {
            for( int c = 0; c < nCol; c++ ) {

                dSumBlur    = 0;
                dSumKernel  = 0;

                nTop        = MAX( 0, r-nSizeRow );
                nBottom     = MIN( nRow-1, r+nSizeRow );
                nLeft       = MAX( 0, c-nSizeCol );
                nRight      = MIN( nCol-1, c+nSizeCol );

                for( int rr = nTop; rr <= nBottom; rr++ ) {
                    for( int cc = nLeft; cc <= nRight; cc++ ) {

                        kr = rr - r + nSizeRow;
                        kc = cc - c + nSizeCol;

                        dI = T( pImage,rr,cc,v,nRow,nCol );
                        dK = T( pKernel,kr,kc,0,nRowKernel,nColKernel );

                        dSumBlur += dI * dK;
                        dSumKernel  += dK;
                    }
                }

                T( pImageBlur,r,c,v,nRow,nCol ) = dSumBlur / dSumKernel;
            }
        }
    }
}

void    CSmoothing::BlurMask( float *pImageBlur, const float *pImage, const bool *pMask, int nRow, int nCol, int nVec, float dSigma )
{
    if( dSigma == 0 ) {

        for( int v = 0; v < nVec; v++ ) {
            for( int r = 0; r < nRow; r++ ) {
                for( int c = 0; c < nCol; c++ ) {

                    if( T( pMask,r,c,0,nRow,nCol ) ) {

                        T( pImageBlur,r,c,v,nRow,nCol ) = T( pImage,r,c,v,nRow,nCol );
                    }
                }
            }
        }

        return;
    }

    int     nSizeRow = floor( dSigma * 3 );
    int     nSizeCol = floor( dSigma * 3 );
    int     hr, hc;
    float  hg, sum_hg = 0;

    int     nRowKernel  = 2 * nSizeRow + 1;
    int     nColKernel  = 2 * nSizeCol + 1;
    float  *pKernel    = new float [ nRowKernel * nColKernel ];

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hr = r - nSizeRow;
            hc = c - nSizeCol;

            hg = exp( - ((hr*hr) + (hc*hc)) / (2*dSigma*dSigma) );

            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg;
            sum_hg += hg;
        }
    }

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hg = T( pKernel, r, c, 0, nRowKernel, nColKernel );
            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg / sum_hg;
        }
    }

    float  dI, dK;
    float  dSumKernel = 0, dSumBlur = 0;
    int     nTop, nBottom, nLeft, nRight;
    int     kr, kc;

    for( int v = 0; v < nVec; v++ ) {
        for( int r = 0; r < nRow; r++ ) {
            for( int c = 0; c < nCol; c++ ) {

                if( T( pMask,r,c,0,nRow,nCol ) ) {

                    dSumBlur    = 0;
                    dSumKernel  = 0;

                    nTop        = MAX( 0, r-nSizeRow );
                    nBottom     = MIN( nRow-1, r+nSizeRow );
                    nLeft       = MAX( 0, c-nSizeCol );
                    nRight      = MIN( nCol-1, c+nSizeCol );

                    for( int rr = nTop; rr <= nBottom; rr++ ) {
                        for( int cc = nLeft; cc <= nRight; cc++ ) {

                            if( T( pMask,rr,cc,0,nRow,nCol ) ) {

                                kr = rr - r + nSizeRow;
                                kc = cc - c + nSizeCol;

                                dI = T( pImage,rr,cc,v,nRow,nCol );
                                dK = T( pKernel,kr,kc,0,nRowKernel,nColKernel );

                                dSumBlur += dI * dK;
                                dSumKernel  += dK;
                            }
                        }
                    }

                    T( pImageBlur,r,c,v,nRow,nCol ) = dSumBlur / dSumKernel;
                }
            }
        }
    }
}


void    CSmoothing::BlurLabel( float *pImageBlur, const float *pImage, const int *pLabel, int nLabel, int nRow, int nCol, int nVec, float dSigma )
{
    if( dSigma == 0 ) {

        for( int v = 0; v < nVec; v++ ) {
            for( int r = 0; r < nRow; r++ ) {
                for( int c = 0; c < nCol; c++ ) {

                    if( T( pLabel,r,c,0,nRow,nCol ) == nLabel ) {

                        T( pImageBlur,r,c,v,nRow,nCol ) = T( pImage,r,c,v,nRow,nCol );
                    }
                }
            }
        }

        return;
    }

    int     nSizeRow = floor( dSigma * 3 );
    int     nSizeCol = floor( dSigma * 3 );
    int     hr, hc;
    float   hg, sum_hg = 0;

    int     nRowKernel  = 2 * nSizeRow + 1;
    int     nColKernel  = 2 * nSizeCol + 1;
    float   *pKernel    = new float [ nRowKernel * nColKernel ];

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hr = r - nSizeRow;
            hc = c - nSizeCol;

            hg = exp( - ((hr*hr) + (hc*hc)) / (2*dSigma*dSigma) );

            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg;
            sum_hg += hg;
        }
    }

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hg = T( pKernel, r, c, 0, nRowKernel, nColKernel );
            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg / sum_hg;
        }
    }

    float   dI, dK;
    float   dSumKernel = 0, dSumBlur = 0;
    int     nTop, nBottom, nLeft, nRight;
    int     kr, kc;

    for( int v = 0; v < nVec; v++ ) {
        for( int r = 0; r < nRow; r++ ) {
            for( int c = 0; c < nCol; c++ ) {

                if( T( pLabel,r,c,0,nRow,nCol ) == nLabel ) {

                    dSumBlur    = 0;
                    dSumKernel  = 0;

                    nTop        = MAX( 0, r-nSizeRow );
                    nBottom     = MIN( nRow-1, r+nSizeRow );
                    nLeft       = MAX( 0, c-nSizeCol );
                    nRight      = MIN( nCol-1, c+nSizeCol );

                    for( int rr = nTop; rr <= nBottom; rr++ ) {
                        for( int cc = nLeft; cc <= nRight; cc++ ) {

                            if( T( pLabel,rr,cc,0,nRow,nCol ) == nLabel ) {

                                kr = rr - r + nSizeRow;
                                kc = cc - c + nSizeCol;

                                dI = T( pImage,rr,cc,v,nRow,nCol );
                                dK = T( pKernel,kr,kc,0,nRowKernel,nColKernel );

                                dSumBlur += dI * dK;
                                dSumKernel  += dK;
                            }
                        }
                    }

                    T( pImageBlur,r,c,v,nRow,nCol ) = dSumBlur / dSumKernel;
                }
            }
        }
    }
}


void    CSmoothing::BlurLabelRegion( float *pImageBlur, const float *pImage, const int *pLabel, int nLabel, const int *pIndexRegion, int nNumIndex, int nRow, int nCol, int nVec, float dSigma )
{
    int nIndex, r, c;

    if( dSigma == 0 ) {

        for( int v = 0; v < nVec; v++ ) {

            for( int n = 0; n < nNumIndex; n++ ) {

                //nIndex  = pIndexRegion[n];
                //r       = nIndex % nRow;
                //c       = ( (nIndex - r) / nRow ) % nCol;

                nIndex  = pIndexRegion[n];
                r       = ROW(nIndex, nRow, nCol);
                c       = COL(nIndex, nRow, nCol);

                if( T( pLabel,r,c,0,nRow,nCol ) == nLabel ) {

                    T( pImageBlur,r,c,v,nRow,nCol ) = T( pImage,r,c,v,nRow,nCol );
                }
            }
        }

        return;
    }

    int     nSizeRow = floor( dSigma * 3 );
    int     nSizeCol = floor( dSigma * 3 );
    int     hr, hc;
    float   hg, sum_hg = 0;

    int     nRowKernel  = 2 * nSizeRow + 1;
    int     nColKernel  = 2 * nSizeCol + 1;
    float   *pKernel    = new float [ nRowKernel * nColKernel ];

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hr = r - nSizeRow;
            hc = c - nSizeCol;

            hg = exp( - ((hr*hr) + (hc*hc)) / (2*dSigma*dSigma) );

            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg;
            sum_hg += hg;
        }
    }

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hg = T( pKernel, r, c, 0, nRowKernel, nColKernel );
            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg / sum_hg;
        }
    }

    float   dI, dK;
    float   dSumKernel = 0, dSumBlur = 0;
    int     nTop, nBottom, nLeft, nRight;
    int     kr, kc;

    for( int v = 0; v < nVec; v++ ) {

        for( int n = 0; n < nNumIndex; n++ ) {

            //nIndex  = pIndexRegion[n];
            //r       = nIndex % nRow;
            //c       = ( (nIndex - r) / nRow ) % nCol;

            nIndex  = pIndexRegion[n];
            r       = ROW(nIndex, nRow, nCol);
            c       = COL(nIndex, nRow, nCol);

            if( T( pLabel,r,c,0,nRow,nCol ) == nLabel ) {

                dSumBlur    = 0;
                dSumKernel  = 0;

                nTop        = MAX( 0, r-nSizeRow );
                nBottom     = MIN( nRow-1, r+nSizeRow );
                nLeft       = MAX( 0, c-nSizeCol );
                nRight      = MIN( nCol-1, c+nSizeCol );

                for( int rr = nTop; rr <= nBottom; rr++ ) {
                    for( int cc = nLeft; cc <= nRight; cc++ ) {

                        if( T( pLabel,rr,cc,0,nRow,nCol ) == nLabel ) {

                            kr = rr - r + nSizeRow;
                            kc = cc - c + nSizeCol;

                            dI = T( pImage,rr,cc,v,nRow,nCol );
                            dK = T( pKernel,kr,kc,0,nRowKernel,nColKernel );

                            dSumBlur += dI * dK;
                            dSumKernel  += dK;
                        }
                    }
                }

                T( pImageBlur,r,c,v,nRow,nCol ) = dSumBlur / dSumKernel;
            }
        }
    }
}

void    CSmoothing::BlurMaskRegion( float *pImageBlur, const float *pImage, const bool *pMask, const int *pIndexRegion, int nNumIndex, int nRow, int nCol, int nVec, float dSigma )
{
    int nIndex, r, c;

    if( dSigma == 0 ) {

        for( int v = 0; v < nVec; v++ ) {

            for( int n = 0; n < nNumIndex; n++ ) {

                //nIndex  = pIndexRegion[n];
                //r       = nIndex % nRow;
                //c       = ( (nIndex - r) / nRow ) % nCol;

                nIndex  = pIndexRegion[n];
                r       = ROW(nIndex, nRow, nCol);
                c       = COL(nIndex, nRow, nCol);

                if( T( pMask,r,c,0,nRow,nCol ) ) {

                    T( pImageBlur,r,c,v,nRow,nCol ) = T( pImage,r,c,v,nRow,nCol );
                }
            }
        }

        return;
    }

    int     nSizeRow = floor( dSigma * 3 );
    int     nSizeCol = floor( dSigma * 3 );
    int     hr, hc;
    float   hg, sum_hg = 0;

    int     nRowKernel  = 2 * nSizeRow + 1;
    int     nColKernel  = 2 * nSizeCol + 1;
    float   *pKernel    = new float [ nRowKernel * nColKernel ];

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hr = r - nSizeRow;
            hc = c - nSizeCol;

            hg = exp( - ((hr*hr) + (hc*hc)) / (2*dSigma*dSigma) );

            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg;
            sum_hg += hg;
        }
    }

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hg = T( pKernel, r, c, 0, nRowKernel, nColKernel );
            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg / sum_hg;
        }
    }

    float   dI, dK;
    float   dSumKernel = 0, dSumBlur = 0;
    int     nTop, nBottom, nLeft, nRight;
    int     kr, kc;

    for( int v = 0; v < nVec; v++ ) {

        for( int n = 0; n < nNumIndex; n++ ) {

            //nIndex  = pIndexRegion[n];
            //r       = nIndex % nRow;
            //c       = ( (nIndex - r) / nRow ) % nCol;

            nIndex  = pIndexRegion[n];
            r       = ROW(nIndex, nRow, nCol);
            c       = COL(nIndex, nRow, nCol);

            if( T( pMask,r,c,0,nRow,nCol ) ) {

                dSumBlur    = 0;
                dSumKernel  = 0;

                nTop        = MAX( 0, r-nSizeRow );
                nBottom     = MIN( nRow-1, r+nSizeRow );
                nLeft       = MAX( 0, c-nSizeCol );
                nRight      = MIN( nCol-1, c+nSizeCol );

                for( int rr = nTop; rr <= nBottom; rr++ ) {
                    for( int cc = nLeft; cc <= nRight; cc++ ) {

                        if( T( pMask,rr,cc,0,nRow,nCol ) ) {

                            kr = rr - r + nSizeRow;
                            kc = cc - c + nSizeCol;

                            dI = T( pImage,rr,cc,v,nRow,nCol );
                            dK = T( pKernel,kr,kc,0,nRowKernel,nColKernel );

                            dSumBlur += dI * dK;
                            dSumKernel  += dK;
                        }
                    }
                }

                T( pImageBlur,r,c,v,nRow,nCol ) = dSumBlur / dSumKernel;
            }
        }
    }
}

void    CSmoothing::BlurLabelAll( float *pImageBlur, const float *pImage, const int *pLabel, int nRow, int nCol, int nVec, float dSigma )
{
    int nIndex, r, c;

    if( dSigma == 0 ) {

        memcpy(pImageBlur, pImage, sizeof(float) * nRow * nCol * nVec);
        return;
    }

    int     nSizeRow = floor( dSigma * 3 );
    int     nSizeCol = floor( dSigma * 3 );
    int     hr, hc;
    float   hg, sum_hg = 0;

    int     nRowKernel  = 2 * nSizeRow + 1;
    int     nColKernel  = 2 * nSizeCol + 1;
    float   *pKernel    = new float [ nRowKernel * nColKernel ];

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hr = r - nSizeRow;
            hc = c - nSizeCol;

            hg = exp( - ((hr*hr) + (hc*hc)) / (2*dSigma*dSigma) );

            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg;
            sum_hg += hg;
        }
    }

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hg = T( pKernel, r, c, 0, nRowKernel, nColKernel );
            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg / sum_hg;
        }
    }

    float   dI, dK;
    float   dSumKernel = 0, dSumBlur = 0;
    int     nTop, nBottom, nLeft, nRight;
    int     nLabel, kr, kc;

    for( int v = 0; v < nVec; v++ ) {
        for( int r = 0; r < nRow; r++ ) {
            for( int c = 0; c < nCol; c++ ) {

                nLabel = T( pLabel,r,c,0,nRow,nCol );

                dSumBlur    = 0;
                dSumKernel  = 0;

                nTop        = MAX( 0, r-nSizeRow );
                nBottom     = MIN( nRow-1, r+nSizeRow );
                nLeft       = MAX( 0, c-nSizeCol );
                nRight      = MIN( nCol-1, c+nSizeCol );

                for( int rr = nTop; rr <= nBottom; rr++ ) {
                    for( int cc = nLeft; cc <= nRight; cc++ ) {

                        if( T( pLabel,rr,cc,0,nRow,nCol ) == nLabel ) {

                            kr = rr - r + nSizeRow;
                            kc = cc - c + nSizeCol;

                            dI = T( pImage,rr,cc,v,nRow,nCol );
                            dK = T( pKernel,kr,kc,0,nRowKernel,nColKernel );

                            dSumBlur += dI * dK;
                            dSumKernel  += dK;
                        }
                    }
                }

                T( pImageBlur,r,c,v,nRow,nCol ) = dSumBlur / dSumKernel;
            }
        }
    }
}

void    CSmoothing::BlurRegionBand( float *pImageBlur, float *pMeanRegion, const float *pImage, const bool *pMaskRegion, const int *pIndexRegion, const int &nNumIndexRegion, const bool *pMaskBand, const int *pIndexBand, const int &nNumIndexBand,  const int &nRow, const int &nCol, const int &nVec, const float &dSigma )
{
    int     nIndex, r, c;
    float   *pSumRegion     = new float [nVec];

    memset(pSumRegion, 0, sizeof(float) * nVec);

    if( dSigma == 0 ) {

        for( int v = 0; v < nVec; v++ ) {

            for( int n = 0; n < nNumIndexRegion; n++ ) {

                nIndex  = pIndexRegion[n];
                r       = ROW(nIndex, nRow, nCol);
                c       = COL(nIndex, nRow, nCol);

                T( pImageBlur,r,c,v,nRow,nCol ) = T( pImage,r,c,v,nRow,nCol );
                pSumRegion[v] += T( pImageBlur,r,c,v,nRow,nCol );
            }

            for( int m = 0; m < nNumIndexBand; m++ ) {

                nIndex  = pIndexBand[m];
                r       = ROW(nIndex, nRow, nCol);
                c       = COL(nIndex, nRow, nCol);

                T( pImageBlur,r,c,v,nRow,nCol ) = T( pImage,r,c,v,nRow,nCol );
            }

            pMeanRegion[v] = pSumRegion[v] / nNumIndexRegion;
        }

        return;
    }

    int     nSizeRow = floor( dSigma * 3 );
    int     nSizeCol = floor( dSigma * 3 );
    int     hr, hc;
    float   hg, sum_hg = 0;

    int     nRowKernel  = 2 * nSizeRow + 1;
    int     nColKernel  = 2 * nSizeCol + 1;
    float   *pKernel    = new float [ nRowKernel * nColKernel ];

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hr = r - nSizeRow;
            hc = c - nSizeCol;

            hg = exp( - ((hr*hr) + (hc*hc)) / (2*dSigma*dSigma) );

            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg;
            sum_hg += hg;
        }
    }

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hg = T( pKernel, r, c, 0, nRowKernel, nColKernel );
            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg / sum_hg;
        }
    }

    float   dI, dK;
    float   dSumKernel = 0, dSumBlur = 0;
    int     nTop, nBottom, nLeft, nRight;
    int     kr, kc;

    for( int v = 0; v < nVec; v++ ) {

        for( int n = 0; n < nNumIndexRegion; n++ ) {

            nIndex  = pIndexRegion[n];
            r       = ROW(nIndex, nRow, nCol);
            c       = COL(nIndex, nRow, nCol);

            dSumBlur    = 0;
            dSumKernel  = 0;

            nTop        = MAX( 0, r-nSizeRow );
            nBottom     = MIN( nRow-1, r+nSizeRow );
            nLeft       = MAX( 0, c-nSizeCol );
            nRight      = MIN( nCol-1, c+nSizeCol );

            for( int rr = nTop; rr <= nBottom; rr++ ) {
                for( int cc = nLeft; cc <= nRight; cc++ ) {

                    if( T( pMaskRegion,rr,cc,0,nRow,nCol ) ) {

                        kr = rr - r + nSizeRow;
                        kc = cc - c + nSizeCol;

                        dI = T( pImage,rr,cc,v,nRow,nCol );
                        dK = T( pKernel,kr,kc,0,nRowKernel,nColKernel );

                        dSumBlur += dI * dK;
                        dSumKernel  += dK;
                    }
                }
            }

            T( pImageBlur,r,c,v,nRow,nCol ) = dSumBlur / dSumKernel;
            pSumRegion[v] += T( pImageBlur,r,c,v,nRow,nCol );
        }

        for( int m = 0; m < nNumIndexBand; m++ ) {

            nIndex  = pIndexBand[m];
            r       = ROW(nIndex, nRow, nCol);
            c       = COL(nIndex, nRow, nCol);

            dSumBlur    = 0;
            dSumKernel  = 0;

            nTop        = MAX( 0, r-nSizeRow );
            nBottom     = MIN( nRow-1, r+nSizeRow );
            nLeft       = MAX( 0, c-nSizeCol );
            nRight      = MIN( nCol-1, c+nSizeCol );

            for( int rr = nTop; rr <= nBottom; rr++ ) {
                for( int cc = nLeft; cc <= nRight; cc++ ) {

                    if( T( pMaskRegion,rr,cc,0,nRow,nCol ) || T( pMaskBand,rr,cc,0,nRow,nCol ) ) {

                        kr = rr - r + nSizeRow;
                        kc = cc - c + nSizeCol;

                        dI = T( pImage,rr,cc,v,nRow,nCol );
                        dK = T( pKernel,kr,kc,0,nRowKernel,nColKernel );

                        dSumBlur += dI * dK;
                        dSumKernel  += dK;
                    }
                }
            }

            T( pImageBlur,r,c,v,nRow,nCol ) = dSumBlur / dSumKernel;
        }

        pMeanRegion[v] = pSumRegion[v] / nNumIndexRegion;
    }

    delete [] pSumRegion;
}

void    CSmoothing::BlurRegion( float *pImageBlur, float *pMeanRegion, const float *pImage, const bool *pMaskRegion, const int *pIndexRegion, const int &nNumIndexRegion, const int &nRow, const int &nCol, const int &nVec, const float &dSigma )
{
    int     nIndex, r, c;
    float   *pSumRegion     = new float [nVec];

    memset(pSumRegion, 0, sizeof(float) * nVec);

    if( dSigma == 0 ) {

        for( int v = 0; v < nVec; v++ ) {

            for( int n = 0; n < nNumIndexRegion; n++ ) {

                nIndex  = pIndexRegion[n];
                r       = ROW(nIndex, nRow, nCol);
                c       = COL(nIndex, nRow, nCol);

                T( pImageBlur,r,c,v,nRow,nCol ) = T( pImage,r,c,v,nRow,nCol );
                pSumRegion[v] += T( pImageBlur,r,c,v,nRow,nCol );
            }

            pMeanRegion[v] = pSumRegion[v] / nNumIndexRegion;
        }

        return;
    }

    int     nSizeRow = floor( dSigma * 3 );
    int     nSizeCol = floor( dSigma * 3 );
    int     hr, hc;
    float   hg, sum_hg = 0;

    int     nRowKernel  = 2 * nSizeRow + 1;
    int     nColKernel  = 2 * nSizeCol + 1;
    float   *pKernel    = new float [ nRowKernel * nColKernel ];

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hr = r - nSizeRow;
            hc = c - nSizeCol;

            hg = exp( - ((hr*hr) + (hc*hc)) / (2*dSigma*dSigma) );

            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg;
            sum_hg += hg;
        }
    }

    for( int r = 0; r < nRowKernel; r++ ) {
        for( int c = 0; c < nColKernel; c++ ) {

            hg = T( pKernel, r, c, 0, nRowKernel, nColKernel );
            T( pKernel, r, c, 0, nRowKernel, nColKernel ) = hg / sum_hg;
        }
    }

    float   dI, dK;
    float   dSumKernel = 0, dSumBlur = 0;
    int     nTop, nBottom, nLeft, nRight;
    int     kr, kc;

    for( int v = 0; v < nVec; v++ ) {

        for( int n = 0; n < nNumIndexRegion; n++ ) {

            nIndex  = pIndexRegion[n];
            r       = ROW(nIndex, nRow, nCol);
            c       = COL(nIndex, nRow, nCol);

            dSumBlur    = 0;
            dSumKernel  = 0;

            nTop        = MAX( 0, r-nSizeRow );
            nBottom     = MIN( nRow-1, r+nSizeRow );
            nLeft       = MAX( 0, c-nSizeCol );
            nRight      = MIN( nCol-1, c+nSizeCol );

            for( int rr = nTop; rr <= nBottom; rr++ ) {
                for( int cc = nLeft; cc <= nRight; cc++ ) {

                    if( T( pMaskRegion,rr,cc,0,nRow,nCol ) ) {

                        kr = rr - r + nSizeRow;
                        kc = cc - c + nSizeCol;

                        dI = T( pImage,rr,cc,v,nRow,nCol );
                        dK = T( pKernel,kr,kc,0,nRowKernel,nColKernel );

                        dSumBlur += dI * dK;
                        dSumKernel  += dK;
                    }
                }
            }

            T( pImageBlur,r,c,v,nRow,nCol ) = dSumBlur / dSumKernel;
            pSumRegion[v] += T( pImageBlur,r,c,v,nRow,nCol );
        }

        pMeanRegion[v] = pSumRegion[v] / nNumIndexRegion;
    }

    delete [] pSumRegion;
}
