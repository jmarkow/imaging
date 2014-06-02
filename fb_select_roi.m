function [EXTRACTED_ROI STATS]=fb_select_roi(varargin)
%fb_select_roi selects an arbitrary number of roi's for plotting
%
%
%
%


% select file to load

nparams=length(varargin);
baseline=2; % 0 for mean, 1 for median, 2 for trimmed mean
filt_rad=12; % gauss filter radius
filt_alpha=4; % gauss filter alpha
lims=2; % contrast prctile limits
roi_color=[1 1 0];
save_dir='roi';
per=2; % baseline percentile (0 for min)
ave_scale=40; % for adaptive threshold, scale for averaging filter
resize_correct=1; % correction of parameters for resized movies
activity_colormap='gray'; % colormap for activity
mode='dff';
resize=.5;
fig_resize=0;
label_color=[1 .5 0];

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'activity_colormap'
			activity_colormap=varargin{i+1};
		case 'baseline'
			baseline=varargin{i+1};
		case 'filt_rad'
			filt_rad=varargin{i+1};
		case 'trim_per'
			trim_per=varargin{i+1};
		case 'filt_alpha'
			filt_alpha=varargin{i+1};
		case 'save_dir'
			save_dir=varargin{i+1};
		case 'per'
			per=varargin{i+1};
		case 'lims'
			lims=varargin{i+1};
		case 'mode'
			mode=varargin{i+1};
		case 'resize'
			resize=varargin{i+1};
		case 'roi_color'
			roi_color=varargin{i+1};
		case 'fig_resize'
			fig_resize=varargin{i+1};
		case 'label_color'
			label_color=varargin{i+1};
	end
end

disp('Loading data...');

resize_skip=0;
im_resize=1; % if im_resize does not exist as a variable, the data has not been resized!

[filename,pathname]=uigetfile({'*.mat';'*.tif'},'Pick a mat file to extract the image data from',pwd);
[path,file,ext]=fileparts(filename);

if strcmp(ext,'.mat')
	load(fullfile(pathname,filename),'mov_data','im_resize');
end

if ~exist('mov_data','var')

	disp('Retrieving tiff data...');
	
	% assume we're in the mat directory, now drop back to retrieve_mov

	if resize~=1
		[mov_data,frame_idx]=fb_retrieve_mov(fullfile(pathname,filename),'im_resize',resize);
		resize_skip=1;
		im_resize=resize;
	end
end

[rows,columns,frames]=size(mov_data);

% resize if we want

if resize~=1 & ~resize_skip

	disp(['Resizing movie data by factor of ' num2str(resize)]);
	frameone=imresize(mov_data(:,:,1),resize);
	[new_rows,new_columns]=size(frameone);

	new_mov=zeros(new_rows,new_columns,frames);

	for i=1:frames		
		new_mov(:,:,i)=imresize(mov_data(:,:,i),resize);
	end
	
	im_resize=im_resize.*resize;
	mov_data=new_mov;
	clear new_mov;
end

if resize_correct & im_resize~=1

	disp('Correcting parameters since file has been downsampled...');
	filt_rad=round(filt_rad.*im_resize);
	filt_alpha=filt_alpha.*im_resize;

end

[rows,columns,frames]=size(mov_data);
mkdir(save_dir);

% maximum projection
% convert mov_data to df/f

disp('Filtering images, this may take a minute...');

h=fspecial('gaussian',filt_rad,filt_alpha);

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for j=1:frames
	fprintf(1,formatstring,round((j/frames)*100));	
	mov_data(:,:,j)=imfilter(mov_data(:,:,j),h,'circular');
end

fprintf(1,'\n');

raw_proj=max(mov_data,[],3);

% downsample to get the color limits

clims=prctile(raw_proj(:),[lims 100-lims]);

baseline=repmat(prctile(mov_data,per,3),[1 1 frames]);
dff=((mov_data-baseline)./baseline).*100;
dff=single(dff); % convert to single before flattening to preserve memory

dff_clims=prctile(dff(:),[lims 100-lims]);

clear tmp;

switch lower(mode(1))
	case 'd'
		
	case 'r'
		dff=mov_data;
		dff_clims=clims;
end

% convert to uint8, faster loading for animation

dff=min(dff,dff_clims(2)); % clip to max
dff=max(dff-dff_clims(1),0); % clip min
dff=dff./(dff_clims(2)-dff_clims(1)); % normalize to [0,1]
dff=uint8(dff.*255);

slider_fig=figure();
h_dff=image(dff(:,:,1));caxis([0 255]);
set(gca,'xlimmode','manual','ylimmode','manual',...
	'zlimmode','manual','climmode','manual','alimmode','manual');
set(slider_fig,'DoubleBuffer','off');
colormap('gray(255)');

if ~fig_resize
	set(slider_fig,'resize','off','position',[50 50 columns rows+75]);
	set(gca,'units','pixels','pos',[ 0 75 columns rows ]);
end

hold on;
axis off;

hsl = uicontrol(slider_fig,'Style','slider','Min',1,'Max',frames,...
          'SliderStep',[1/frames 1/frames],'Value',1,...
             'Units','Normalized','Position',[.1 .05 .8 .05]);
set(hsl,'Callback',{@slider_callback,slider_fig,dff,h_dff})

% return a cell array with the ROI

EXTRACTED_ROI={}; % indices for the ROI
centroid=[]; % keep the centroids for deleting

[xi,yi]=meshgrid(1:columns,1:rows); % collect all coordinates into xi and yi
exit_flag=0;
counter=1;

h_ellipse=imellipse(get(slider_fig,'currentaxes'));

xlimits=xlim();
ylimits=ylim();

xlimits(1)=xlimits(1)-10;
xlimits(2)=xlimits(2)+10;

ylimits(1)=ylimits(1)-10;
ylimits(2)=ylimits(2)+10;

fcn = makeConstrainToRectFcn('imellipse',xlimits,ylimits);
setPositionConstraintFcn(h_ellipse,fcn);

ellipses={};
pl_ellipses=[];
pl_centroids=[];
diameter=[];

counter=1;
del_flag=0;

while 1>0

	h=wait(h_ellipse);

	% convert the ROI into all pixels inside the polygon
	
	if isempty(h), break; end

	% if any element of h is outside the plotting window, enter delete mode

	if del_flag

		% loop through idx, any matches result in deleting that roi
		% is a centroid inside the ellipse?

		if length(centroid)>0	
		
			del=inpolygon(centroid(:,2),centroid(:,1),h(:,1),h(:,2));
			idx=find(del);

			% clean up

			EXTRACTED_ROI(idx)=[];
			centroid(idx,:)=[];
			ellipses(idx)=[];
			diameter(idx)=[];

			delete(pl_ellipses(idx));
			pl_ellipses(idx)=[];

		end
	end

	if any(h(:,1)>columns) | any(h(:,2)>rows) | any(h(:)<0)
	
		% placing the ellipse on the border changes delete mode

		if del_flag
			disp('Exiting delete mode');			
		else
			disp('Entering delete mode');
		end

		del_flag=~del_flag; % flip del_flag
		continue;
	end

	% don't plot anything if we're in delete mode

	if del_flag
		continue;
	end

	roi=inpolygon(xi,yi,h(:,1),h(:,2));
	[idx]=find(roi);

	% xi=columns yi=rows
	% collect the roi

	EXTRACTED_ROI{end+1}=[ yi(idx) xi(idx) ];
	centroid(end+1,:)=[ mean(yi(idx)) mean(xi(idx)) ];

	dist=pdist(h,'euclidean');
	diameter(end+1)=max(dist);

	set(0,'CurrentFigure',slider_fig);
	
	ellipses{end+1}=h;
	pl_ellipses(end+1)=plot(h(:,1),h(:,2),'-','linewidth',1.5,'color',roi_color);

	% what's inside of the ROI?  this could also be used to normalize fluorescene per 
	% Svoboda et al. (take an annulus surrounding the ROI)

	counter=counter+1; % increment the colormap

end

save_fig=figure();
h_rawproj=imagesc(raw_proj);caxis([clims]);
hold on;
colormap(activity_colormap);
axis off;

for i=1:length(ellipses)
	h=ellipses{i};
	plot(h(:,1),h(:,2),'-','linewidth',1.5,'color',roi_color);
	text(centroid(i,2),centroid(i,1),sprintf('%i',i),...
		'color',label_color,'fontsize',30*im_resize);
end

if resize~=1
	
	disp('Putting ROI coordinates back into the original movie frame...');

	for i=1:length(EXTRACTED_ROI)
		EXTRACTED_ROI{i}=EXTRACTED_ROI{i}.*(1/resize);
		centroid(i,:)=centroid(i,:).*(1/resize);
	end
end

% draw rois onto max_proj and be done!

STATS=struct('centroid',centroid,'diameter',diameter);

save(fullfile(save_dir,'roi_data_manual.mat'),'EXTRACTED_ROI','STATS');
fb_multi_fig_save(save_fig,save_dir,'roi_map_manual','tiff','res',100);
close([save_fig]);

end

function slider_callback(hObject,eventdata,fig,dff,h_dff)

val=get(hObject,'Value');
set(0,'CurrentFigure',fig)
set(h_dff,'cdata',dff(:,:,round(val)));
setappdata(fig,'framenumber',val);

end
