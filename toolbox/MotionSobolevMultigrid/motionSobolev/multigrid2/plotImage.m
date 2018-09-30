function plotImage(I)
  
  figure; imagesc(I>=0.5); colormap gray; axis image; colorbar;
  drawnow;
  pause;