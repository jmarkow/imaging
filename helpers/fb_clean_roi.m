function [ROI INCLUDE] = fb_clean_roi(ROI)
%
%
%

nrois=length(ROI.coordinates);
counter=1;

nrois=length(ROI.coordinates);
TO_DEL=[];

for i=1:nrois
	for j=setdiff(1:nrois,i)
		if all(ROI.stats(j).Centroid==ROI.stats(i).Centroid)
			TO_DEL=[TO_DEL j];
		end
	end
end

ROI.stats(TO_DEL)=[];
ROI.coordinates(TO_DEL)=[];

INCLUDE=setdiff(1:nrois,TO_DEL);
