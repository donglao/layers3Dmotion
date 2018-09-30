#if !defined(_SMOOTHING_H_)
#define _SMOOTHING_H_

class CSmoothing
{

// constructor and destructor
public:
    CSmoothing();
    virtual ~CSmoothing();

// member function
public:
    static  void    Blur( float *pImageBlur, const float *pImage, int nRow, int nCol, int nVec, float dSigma );
    static  void    Blur( double *pImageBlur, const double *pImage, int nRow, int nCol, int nVec, double dSigma );

    static  void    BlurLabel( float *pImageBlur, const float *pImage, const int *pLabel, int nLabel, int nRow, int nCol, int nVec, float dSigma );
    static  void    BlurMask( float *pImageBlur, const float *pImage, const bool *pMask, int nRow, int nCol, int nVec, float dSigma );

    static  void    BlurLabelRegion( float *pImageBlur, const float *pImage, const int *pLabel, int nLabel, const int *pIndexRegion, int nNumIndex, int nRow, int nCol, int nVec, float dSigma );
    static  void    BlurMaskRegion( float *pImageBlur, const float *pImage, const bool *pMask, const int *pIndexRegion, int nNumIndex, int nRow, int nCol, int nVec, float dSigma );

    static  void    BlurRegionBand( float *pImageBlur, float *pMeanRegion, const float *pImage, const bool *pMaskRegion, const int *pIndexRegion, const int &nNumIndexRegion, const bool *pMaskBand, const int *pIndexBand, const int &nNumIndexBand,  const int &nRow, const int &nCol, const int &nVec, const float &dSigma );
    static  void    BlurRegion( float *pImageBlur, float *pMeanRegion, const float *pImage, const bool *pMaskRegion, const int *pIndexRegion, const int &nNumIndexRegion, const int &nRow, const int &nCol, const int &nVec, const float &dSigma );

    static  void    BlurLabelAll( float *pImageBlur, const float *pImage, const int *pLabel, int nRow, int nCol, int nVec, float dSigma );

protected:

// member variable
protected:

};

#endif // !defined(_SMOOTHING_H_)
