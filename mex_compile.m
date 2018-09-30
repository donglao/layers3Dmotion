cd ./toolbox/MotionSobolevMultigrid/motionSobolev

parallelize=1;

mex -O mexComputePropagatedSegmentation.C Smoothing.cpp

if 1,
if parallelize,
	mex -I./multigrid2/ mexMotionEstimation.C Motion_Sobolev.C PoissonMultigrid.C ConnectedComponentDetection.C CXXFLAGS="\$CXXFLAGS -O3 -fopenmp -ffast-math" LDFLAGS="\$LDFLAGS -fopenmp" -DPARALLELIZE_CODE %-DDEBUG
%	mex mexMotionEstimation2.C Motion_Sobolev.C conjugateGradient.C ConnectedComponentDetection.C CXXFLAGS="\$CXXFLAGS -O3 -fopenmp -ffast-math" LDFLAGS="\$LDFLAGS -fopenmp" -DPARALLELIZE_CODE

else
	mex mexMotionEstimation.C Motion_Sobolev.C conjugateGradient.C ConnectedComponentDetection.C CXXFLAGS="\$CXXFLAGS -O3 -ffast-math"
	mex mexMotionEstimation2.C Motion_Sobolev.C conjugateGradient.C ConnectedComponentDetection.C CXXFLAGS="\$CXXFLAGS -O3 -ffast-math"	
end
end

if 0,
if parallelize,
	mex mexMotionEstimation.C Motion_Sobolev.C conjugateGradient_no_virtual.C ConnectedComponentDetection.C CXXFLAGS="\$CXXFLAGS -O3 -fopenmp -ffast-math" LDFLAGS="\$LDFLAGS -fopenmp" -DPARALLELIZE_CODE -v
else
	mex mexMotionEstimation.C Motion_Sobolev.C conjugateGradient_no_virtual.C ConnectedComponentDetection.C CXXFLAGS="\$CXXFLAGS -O3 -ffast-math"
end
end

cd ../../

mex mexIntensity.C

cd ..