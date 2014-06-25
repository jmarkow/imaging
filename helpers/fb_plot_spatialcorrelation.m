function [roi_dist,dff_dist]=fb_plot_spatialcorrelation(ROI,ROI_AVE,varargin)
%
%
%
%
%
%



% select file to load

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})	
		case 'fig_num'
			fig_num=varargin{i+1};
			
	end
end

% get the pairwise distance between all ROIs

centroid=cat(1,ROI.stats(:).Centroid);

% for now, take euclidean distance between the centroid of each ROI

roi_dist=pdist(centroid,'euclidean');

% average dff traces across trials

dff=mean(ROI_AVE.interp_dff,3);
dff=zscore(dff')';

% activity difference, correlation

dff_dist=pdist(dff,'seuclidean');

fig_num=figure();
scatter(roi_dist(:),dff_dist(:));


