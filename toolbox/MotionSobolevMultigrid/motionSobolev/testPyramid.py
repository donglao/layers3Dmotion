from numpy import *
import matplotlib.pyplot as plt
from skimage import morphology
import pdb

dData = "/Users/sundarg/work/Datasets_Vision/moseg_dataset/marple7/";
strI0 = "marple7_001.jpg";
strI1 = "marple7_002.jpg";
strGT = "marple7_001.pgm";

I0     = plt.imread( "%s%s" %(dData, strI0) );
I1     = plt.imread( "%s%s" %(dData, strI1) );
region = plt.imread( "%s%s" %(dData, strGT) );

J0     = float64(I0) / 255;
J1     = float64(I0) / 255;
region = int8( region==255 );

strel  = morphology.square( 30 );
region = morphology.dilation( region, strel );

plt.subplot(1,3,1);
plt.imshow(I0); 
plt.title('I0');
plt.subplot(1,3,2);
plt.imshow(region);
plt.title('region');
plt.subplot(1,3,3);
plt.imshow(I1);
plt.title('I1');