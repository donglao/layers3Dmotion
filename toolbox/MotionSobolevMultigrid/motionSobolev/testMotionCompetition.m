dWork = '/Users/hong/work/FastLabelScaleSpace/';
dData = [dWork 'data/videos/moseg_dataset/marple7/'];


nIters = 1000;

I0=imread( strcat( dData, 'marple7_001.jpg') );
I1=imread( strcat( dData, 'marple7_003.jpg') );

region = imread( strcat( dData,'marple7_001.pgm' ) );

I0 = imresize( I0, 0.5 );
I1 = imresize( I1, 0.5 );
region = imresize( region, 0.5 );

label = region == 255;

if 0
se =  strel( 'square', 15 );
label = int8( imdilate( label, se ) );
end

%label = int8(region == 255);

J0 = double( I0 )/255;
J1 = double( I1 )/255;

figure;
MotionCompetition( J0, J1, label, nIters );
