function fig_num=fb_visualize_rois(ROI,varargin)
%
%
%
%
%
%



% select file to load

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
			
	end
end


nrois=length(ROI.coordinates);

if isempty(fig_num)
	fig_num=figure();
end

if ref_image
	imagesc(ROI.reference_image);
	colormap(gray);
	caxis([clims]);
	axis off;
	hold on;
end

% scale bar?

% scale weights, 64 colors

length(weights)

if ~isempty(weights)
	weights_map=eval([ weights_map '(' num2str(length(weights)) ')' ]);
end

% map weights to colors

weights=(weights-min(weights))./(max(weights)-min(weights));
weights=ceil(weights.*(length(weights)-1)+1);

if strcmp(lower(ROI.type(1)),'m') | strcmp(lower(ROI.type(1)),'a')

	counter=1;	
	
	for i=1:nrois

		tmp=ROI.stats(i).ConvexHull;

		if filled
			if ~isempty(weights)
				fill(tmp(:,1),tmp(:,2),weights_map(weights(i),:))
			else
				fill(tmp(:,1),tmp(:,2),roi_map(counter,:))
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


end

