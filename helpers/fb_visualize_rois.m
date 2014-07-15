function fig_num=fb_visualize_rois(ROI,varargin)
%
%
%
%
%
%



nparams=length(varargin);

roi_map=[1 0 1];
bg_map='gray';

label_fontsize=25;
label_color=[1 1 0];
clims=[0 1];
label=0;
fig_num=[];
filled=0;
weights=[];
weights_map='winter';
ref_image=1;
ncolors=[];

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})	
		case 'roi_map'
			roi_map=varargin{i+1};
		case 'label_fontsize'
			label_fontsize=varargin{i+1};
		case 'label_color'
			label_color=varargin{i+1};
		case 'label'
			label=varargin{i+1};
		case 'resize'
			resize=varargin{i+1};
		case 'clims'
			clims=varargin{i+1};
		case 'scale_bar'
			scale_bar=varargin{i+1};
		case 'fig_num'
			fig_num=varargin{i+1};
		case 'filled'
			filled=varargin{i+1};
		case 'weights'
			weights=varargin{i+1};
		case 'weights_map'
			weights_map=varargin{i+1};
		case 'ref_image'
			ref_image=varargin{i+1};
		case 'ncolors'
			ncolors=varargin{i+1};
			
	end
end


nrois=length(ROI.coordinates);

if isempty(fig_num)
	fig_num=figure();
end

if ref_image
	imagesc(ROI.reference_image);
	colormap(gray);
	axis off;
	hold on;
end

% scale bar?

% scale weights, 64 colors

if ~iscell(weights)
	nweights=length(weights);
else
	nweights=sum(cellfun(@length,weights));
end

if ~isempty(weights)
	if ~isempty(ncolors)
		weights_map=eval([ weights_map '(' num2str(ncolors) ')' ]);
	else
		weights_map=eval([ weights_map '(' num2str(nweights) ')' ]);
		ncolors=nweights;
	end

end

weights_map
% map weights to colors

if ~iscell(weights)
	weights=(weights-min(weights))./(max(weights)-min(weights));
	weights=ceil(weights.*(ncolors-1)+1);
else


	tmp=cat(2,weights{:});
	weights_min=min(tmp)
	weights_max=max(tmp)

	for i=1:length(weights)
		weights{i}=(weights{i}-weights_min)./(weights_max-weights_min);
		weights{i}=ceil(weights{i}.*(ncolors-1)+1);
	end
end


if ~isfield(ROI.stats,'ConvexHull')
	for i=1:nrois
		k=convhull(ROI.coordinates{i}(:,1),ROI.coordinates{i}(:,2));
		ROI.stats(i).ConvexHull=ROI.coordinates{i}(k,:);
	end
end


counter=1;	

for i=1:nrois

	tmp=ROI.stats(i).ConvexHull;
	tmp_c=ROI.stats(i).Centroid;
	npoints=size(tmp,1);

	if filled
		if ~iscell(weights)
			if ~isempty(weights)
				fill(tmp(:,1),tmp(:,2),weights_map(weights(i),:))
			else
				fill(tmp(:,1),tmp(:,2),roi_map(counter,:))
			end
		else
			if ~isempty(weights{i})

				split=length(weights{i});
				split_points=floor(npoints/split);
				split_segments=floor(linspace(1,npoints,split+1));

				for j=2:split+1
					slice=split_segments(j-1):split_segments(j);
					fill([tmp_c(1);tmp(slice,1)],[tmp_c(2);tmp(slice,2)],...
						weights_map(weights{i}(j-1),:),'edgecolor','none');
				end
			else
				plot(tmp(:,1),tmp(:,2),'-','linewidth',1,'color',roi_map(counter,:));
			end
		end

	else	
		plot(tmp(:,1),tmp(:,2),'-','linewidth',1,'color',roi_map(counter,:));
	end

	hold on;

	if label

		x=mean(ROI.coordinates{i}(:,1));
		y=mean(ROI.coordinates{i}(:,2));

		text(x,y,sprintf('%i',i),...
			'color',label_color,'fontsize',label_fontsize,'fontweight','bold');
	end
	
	if counter<size(roi_map,1)
		counter=counter+1;
	else
		counter=1;
	end

end

if ~ref_image
	set(gca,'ydir','rev');
end
