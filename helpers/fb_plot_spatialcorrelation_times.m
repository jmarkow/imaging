function [dist_t,dist_s,data_y]=fb_plot_spatialcorrelation(ROI,PEAKS,varargin)
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
range=[-inf inf];


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
		case 'range'
			range=varargin{i+1};
	end	
end



% get the pairwise distance between all ROIs

centroid=cat(1,ROI.stats(:).Centroid);
nrois=length(ROI.stats);
% for now, take euclidean distance between the centroid of each ROI

for i=1:nrois
	PEAKS{i}(PEAKS{i}<range(1)|PEAKS{i}>range(2))=[];
end

if ~isempty(exclude)

	to_del=find(dff_max<exclude);
	centroid(to_del,:)=[];
	dff(to_del,:)=[];

end

% activity difference, correlation

%roi_dist=pdist(centroid,'euclidean');

% take the center of mass for each dff trace

combos=nchoosek([1:nrois],2)

dist_t=[];
dist_s=[];

for i=1:size(combos,1)

	% take minimum distance

	m1=combos(i,1);
	m2=combos(i,2);

	pk_times=[PEAKS{m1} PEAKS{m2}]'
		
	if length(PEAKS{m1})<1|length(PEAKS{m2})<1
		continue;
	end

	total_min=inf;

	tmp=[];
	for j=1:length(PEAKS{m1})

		for k=1:length(PEAKS{m2})
		
			tmp(end+1)=abs(PEAKS{m1}(j)-PEAKS{m2}(k));

			if tmp<total_min
				total_min=tmp(end);
			end
		end

	end

	%dist_t(end+1)=total_min;
	dist_t(end+1)=min(tmp);
	dist_s(end+1)=sqrt(sum((centroid(m1,:)-centroid(m2,:)).^2));

end

%fig_num=figure();
%scatter(dist_s(:),dist_t(:));

% bin by dist_s

bins=unique([0:100:1000]);

[n,bins_x]=histc(dist_s,bins);

for i=1:length(bins)
	data_y{i}=dist_t(bins_x==i);
end

[r,p]=corrcoef(dist_t,dist_s)


