function [roi_dist,dff_dist]=fb_plot_spatialcorrelation(ROI,ROI_AVE,varargin)
%
%
%
%
%
%



% select file to load

nparams=length(varargin);
exclude=[];
use_com=0;

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})	
		case 'fig_num'
			fig_num=varargin{i+1};
		case 'exclude'
			exclude=varargin{i+1};	
		case 'use_com'
			use_com=varargin{i+1};
	end	
end

% get the pairwise distance between all ROIs

centroid=cat(1,ROI.stats(:).Centroid);

% for now, take euclidean distance between the centroid of each ROI


% average dff traces across trials

dff=mean(ROI_AVE.interp_dff,3);
dff_max=max(dff,[],2);

if ~isempty(exclude)

	to_del=find(dff_max<exclude);
	centroid(to_del,:)=[];
	dff(to_del,:)=[];

end

%dff=zscore(dff')';

% activity difference, correlation

roi_dist=pdist(centroid,'euclidean');

% take the center of mass for each dff trace

if use_com
	[nrois,nsamples]=size(dff)

	ind=1:nsamples;
	com=zeros(1,nrois);

	for i=1:nrois
		com(i)=sum(ind.*dff(i,:))./sum(dff(i,:));
	end

	dff=com(:);
end

dff_dist=pdist(dff,'euclidean');

fig_num=figure();
scatter(roi_dist(:),dff_dist(:));


