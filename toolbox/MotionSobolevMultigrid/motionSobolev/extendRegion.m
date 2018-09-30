function [region_ext,band,ind]=extendRegion(label,labelNo,band_length)

	se=strel('disk',band_length,0);
	se=strel('square',3);
	region=double(label==labelNo);
	region_ext=imdilate(region,se);
	band=region_ext-region;
	region_ext=int32(region_ext);
	ind=find(band==1);