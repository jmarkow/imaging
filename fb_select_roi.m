function EXTRACTED_ROI=fb_select_roi(DIR,varargin)
%fb_select_roi selects an arbitrary number of roi's for plotting
%
%
%
%


% select file to load

nparams=length(varargin);
baseline=2; % 0 for mean, 1 for median, 2 for trimmed mean
filt_rad=30; % gauss filter radius
filt_alpha=2; % gauss filter alpha
lims=2; % contrast prctile limits
roi_map=lines(40);
save_dir='roi';
per=2; % baseline percentile (0 for min)
ave_scale=40; % for adaptive threshold, scale for averaging filter
resize_correct=1; % correction of parameters for resized movies
activity_colormap='gray'; % colormap for activity
mode='dff';

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
		case 'w'
			w=varargin{i+1};
		case 'erode_scale'
			erode_scale=varargin{i+1};
		case 'dilate_scale'
			dilate_scale=varargin{i+1};
		case 'ave_scale'
			ave_scale=varargin{i+1};
		case 'per'
			per=varargin{i+1};
		case 'lims'
			lims=varargin{i+1};
		case 'mode'
			mode=varargin{i+1};
	end
end


if nargin<1 | isempty(DIR), DIR=pwd; end

im_resize=1; % if im_resize does not exist as a variable, the data has not been resized!
[filename,pathname]=uigetfile({'*.mat'},'Pick a mat file to extract the image data from',fullfile(DIR,'..'));
load(fullfile(pathname,filename),'mov_data','im_resize');

%mov_data=mov_data(:,:,1:5);

if resize_correct & im_resize~=1

	disp('Correcting parameters since file has been downsampled...');
	ave_scale=round(ave_scale.*im_resize);
	filt_rad=round(filt_rad.*im_resize);
	filt_alpha=filt_alpha.*im_resize;

end


mkdir(save_dir);

% maximum projection
% convert mov_data to df/f

disp('Filtering images, this may take a minute...');

[rows,columns,frames]=size(mov_data);
h=fspecial('gaussian',filt_rad,filt_alpha);

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for j=1:frames
	fprintf(1,formatstring,round((j/frames)*100));	
	mov_data(:,:,j)=imfilter(mov_data(:,:,j),h,'circular');
end

fprintf(1,'\n');

raw_proj=max(mov_data,[],3);
clims(1)=prctile(raw_proj(:),lims);
clims(2)=prctile(raw_proj(:),100-lims);

baseline=repmat(prctile(mov_data,per,3),[1 1 frames]);
dff=((mov_data-baseline)./baseline).*100;

dff_clims(1)=prctile(dff(:),lims);
dff_clims(2)=prctile(dff(:),100-lims);

switch lower(mode(1))
	case 'd'
		
	case 'r'
		dff=mov_data;
		dff_clims=clims;
end

save_fig=figure();
h_rawproj=imagesc(raw_proj);caxis([clims]);
hold on;
colormap(activity_colormap);
axis off;

slider_fig=figure();
h_dff=imagesc(dff(:,:,1));caxis([dff_clims]);
hold on;
colormap(activity_colormap);
axis off;

hsl = uicontrol(slider_fig,'Style','slider','Min',1,'Max',frames,...
          'SliderStep',[1/frames 1/frames],'Value',1,...
             'Units','Normalized','Position',[.2 .05 .6 .05]);
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

pl_ellipses=[];
pl_centroids=[];
pl_ellipses_raw=[];
txt_centroid=[];

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
			delete(pl_centroids(idx));
			delete(pl_ellipses(idx));
			delete(txt_centroid(:));
			delete(pl_ellipses_raw(idx));

			txt_centroid=[];
			pl_centroids(idx)=[];
			pl_ellipses(idx)=[];
			pl_ellipses_raw(idx)=[];

			% renumber so we're consistent

			for i=1:size(centroid)
				set(0,'CurrentFigure',save_fig);
				txt_centroid(end+1)=text(centroid(i,2),centroid(i,1),sprintf('%i',i),...
					'color',[1 0 0],'fontsize',30*im_resize),
			end
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

	set(0,'CurrentFigure',slider_fig);
	pl_ellipses(end+1)=plot(h(:,1),h(:,2),'-','linewidth',1.5,'color',roi_map(counter,:));

	set(0,'CurrentFigure',save_fig);
	pl_ellipses_raw(end+1)=plot(h(:,1),h(:,2),'-','linewidth',1.5,'color',roi_map(counter,:));
	txt_centroid(end+1)=text(centroid(end,2),centroid(end,1),sprintf('%i',size(centroid,1)),...
		'color',[1 0 0],'fontsize',30*im_resize),

	set(0,'CurrentFigure',slider_fig);

	% what's inside of the ROI?  this could also be used to normalize fluorescene per 
	% Svoboda et al. (take an annulus surrounding the ROI)
	set(0,'CurrentFigure',slider_fig);
	pl_centroids(end+1)=...
		plot(centroid(end,2),centroid(end,1),'k.','color',roi_map(counter,:),'markersize',45.*im_resize);

	counter=counter+1; % increment the colormap

	if counter>size(roi_map,1)
		counter=counter-size(roi_map,1);
	end

end

% draw rois onto max_proj and be done!

save(fullfile(save_dir,'ave_roi_manual.mat'),'EXTRACTED_ROI','centroid');
fb_multi_fig_save(save_fig,save_dir,'roi_map_manual','tiff','res',100);
close([save_fig]);

end

function slider_callback(hObject,eventdata,fig,dff,h_dff)

val=get(hObject,'Value');
set(0,'CurrentFigure',fig)
set(h_dff,'cdata',dff(:,:,round(val)));
setappdata(fig,'framenumber',val);

end
