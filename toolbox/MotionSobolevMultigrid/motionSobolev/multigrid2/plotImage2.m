function plotImage2(Ix,Iy,str)
  
  figure;
  subplot(121); imagesc(Ix); colormap gray; axis image; colorbar;
  subplot(122); imagesc(Iy); colormap gray; axis image; colorbar;
%  subplot(133); imagesc(sqrt(Ix.*Ix+Iy.*Iy)); colormap gray; axis image; colorbar;
  title(str);
  drawnow;
  pause;