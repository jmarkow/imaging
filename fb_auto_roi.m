function [ROI]=fb_auto_roi(DIR,varargin)
%fb_select_roi selects an arbitrary number of roi's for plotting
%
%
%
% LOGIC:
% 1) spatially smooth each from with a small Gaussian filter (2 micron seems standard)
% 2) compute df/f
%
%


% select file to load

nparams=length(varargin);
baseline=2; % 0 for mean, 1 for median, 2 for trimmed mean
filt_rad=25; % gauss filter radius
filt_alpha=20; % gauss filter alpha
lims=.25; % contrast prctile limits
roi_map='lines';
save_dir='auto_roi';
per=3; % baseline percentile (0 for min)
ave_scale=80; % for adaptive threshold, scale for averaging filter
med_scale=20; % for removing speckle noise from maximum projection
resize_correct=1; % correction of parameters for resized movies
use_xcorr=0; % doesn't seem to improve matters, but leaving in for completeness
roi_ver='.01a';
label_fontsize=50;
label_color=[1 1 0];

% parameters used for morphological opening

erode_scale=2; % scale (in pxs) for erosion
dilate_scale=2; % scale (in pxs) for dilation

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})	
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
	end
end


if nargin<1 | isempty(DIR), DIR=pwd; end

im_resize=1; % if im_resize does not exist as a variable, the data has not been resized!

%TODO: allow for multiple file selection, concatenate movies

[filename,pathname]=uigetfile({'*.mat'},'Pick a mat file to extract the image data from',fullfile(DIR,'..'));
load(fullfile(pathname,filename),'mov_data','im_resize');

if resize_correct & im_resize~=1

	disp('Correcting parameters since file has been downsampled...');
	ave_scale=round(ave_scale.*im_resize);
	med_scale=round(med_scale.*im_resize);
	filt_rad=round(filt_rad.*im_resize);
	filt_alpha=filt_alpha.*im_resize;
	label_fontsize=round(label_fontsize.*im_resize);

end

mkdir(save_dir);

disp('Filtering images, this may take a minute...');

[rows,columns,frames]=size(mov_data);
%h=fspecial('gaussian',filt_rad,filt_alpha);
h=fspecial('disk',filt_rad);

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

raw_proj=max(mov_data,[],3);

for j=1:frames
	fprintf(1,formatstring,round((j/frames)*100));	
	mov_data(:,:,j)=imfilter(mov_data(:,:,j),h,'circular');
end

fprintf(1,'\n');

clims=prctile(raw_proj(:),[lims 100-lims]);

baseline=repmat(prctile(mov_data,per,3),[1 1 frames]);
dff=((mov_data-baseline)./baseline).*100;

dff_mu=zeros(size(dff));

h=fspecial('average',ave_scale);

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for j=1:frames
	fprintf(1,formatstring,round((j/frames)*100));	
	dff_mu(:,:,j)=imfilter(dff(:,:,j),h,'circular');
end

fprintf(1,'\n');

pad_pxs=floor(ave_scale/4); % don't take anything too close to the edge
dff=(dff-dff_mu);
dff=smooth3(dff,'box',[1 1 9]);

if use_xcorr
	max_proj=CrossCorrImage(dff); % from LabRigger, appears to work pretty well!
					   % also could be used to clean up initial pass
else
	max_proj=max(dff,[],3);
end

clear dff;
clear dff_mu;

max_proj=max(max_proj./max(max_proj(:)),0); % convert to [0,1]
max_proj=medfilt2(max_proj,[med_scale med_scale]); 
max_proj=max_proj(pad_pxs:end-pad_pxs,pad_pxs:end-pad_pxs);

% pad

pad_y=zeros(pad_pxs,size(max_proj,2));

max_proj=[ pad_y;max_proj;pad_y ];

pad_x=zeros(size(max_proj,1),pad_pxs);

max_proj=[ pad_x max_proj pad_x ];

% scale the raw projection from [0,1]

raw_proj=min(raw_proj,clims(2)); % clip to max
raw_proj=max(raw_proj-clims(1),0); % clip min
raw_proj=raw_proj./(clims(2)-clims(1)); % normalize to [0,1]

%raw_proj=min(raw_proj,1);


[rows,columns]=size(max_proj);

slmin=min(max_proj(:));
slmax=max(max_proj(:));

% display the anatomical "max projection"

slider_fig=figure('Name',['Auto ROI ver' roi_ver],'NumberTitle','off');
cm1=gray;

% convert the indexed images to rgb

im1_rgb=ind2rgb(round(raw_proj.*size(cm1,1)),cm1);
im2_rgb=ind2rgb(max_proj>0,[0 0 0;1 0 0]);
image(im1_rgb);

% now add the mask, change the transparency using the slider

hold on;
h=image(im2_rgb);
set(h,'AlphaData',(max_proj>0));
set(gca,'position',[0 .2 1 .75],'units','normalized')
axis off;

hed1=uicontrol(slider_fig,'Style','edit','units','normalized','string',num2str(erode_scale),'position',[.45 .01 .1 .1]);
hed1_label=uicontrol(slider_fig,'Style','text','units','normalized','string','Erode scale','position',[.45 .12 .1 .025]);
hed2=uicontrol(slider_fig,'Style','edit','units','normalized','string',num2str(dilate_scale),'position',[.3 .01 .1 .1]);
hed2_label=uicontrol(slider_fig,'Style','text','units','normalized','string','Dilate scale','position',[.3 .12 .1 .025]);

hsl = uicontrol(slider_fig,'Style','slider','units','normalized','Min',slmin,'Max',slmax,...
          'SliderStep',[1e-3 1e-3],'Value',slmin,...
             'Position',[.05 .01 .2 .1]);
set(hsl,'Callback',{@slider_callback,slider_fig,max_proj,h,hed1,hed2})

done = uicontrol(slider_fig,'units','normalized','Position',[.6 .01 .15 .1],'String','Done',...
              'Callback','uiresume(gcbf)');

% TODO: replace with uiwait!

uiwait(gcf);
threshold=get(hsl,'value');
close(slider_fig);

% convert threshold to [0,1)

threshold=threshold;
new_image=im2bw(max_proj,threshold);

% basic morphological operators before conncomp 

new_image=bwmorph(new_image,'clean');

se=strel('disk',erode_scale);
new_image=imerode(new_image,se);

se=strel('disk',dilate_scale);
new_image=imdilate(new_image,se);

% collect some basic stats if we want to exclude later

conn_comp=bwconncomp(new_image);
roi_map=eval([ roi_map '(' num2str(length(conn_comp.PixelIdxList)) ')' ]);

STATS=regionprops(conn_comp,'eccentricity','majoraxislength',...
	'minoraxislength','convexhull','centroid','boundingbox');

% get coordinates and draw ROIs on maximum projection

save_fig=figure();
image(uint8(raw_proj.*255));
colormap(gray(256));
axis off;
hold on;

% structure for further analysis

ROI.coordinates={};
ROI.stats=STATS;
ROI.type='auto'; % auto
ROI.reference_image=raw_proj; % for drawing later

clearvars STATS type reference_image;

for i=1:length(conn_comp.PixelIdxList);

	[yi,xi]=ind2sub(size(max_proj),conn_comp.PixelIdxList{i});
	ROI.coordinates{i}=[xi(:) yi(:)]; % get the coordinates
	tmp=ROI.stats(i).ConvexHull;
	plot(tmp(:,1),tmp(:,2),'-','linewidth',2,'color',roi_map(i,:));
	tmp=ROI.stats(i).BoundingBox;

	x=tmp(1)+tmp(3)/2;
	y=tmp(2)+tmp(4)/2;

	text(x,y,sprintf('%i',i),...
		'color',label_color,'fontsize',label_fontsize,'fontweight','bold');

end

pause();
fb_multi_fig_save(save_fig,save_dir,'roi_map_auto','tiff','res','100');
close(save_fig);
save(fullfile(save_dir,'roi_data_auto.mat'),'ROI');

% add cleanup feature

end

function slider_callback(hObject,eventdata,fig,max_proj,h,erode_box,dilate_box)

alpha=.4;

val=get(hObject,'Value');

erode_scale=round(str2num(get(erode_box,'string')));
set(erode_box,'string',num2str(erode_scale));

dilate_scale=round(str2num(get(dilate_box,'string')));
set(dilate_box,'string',num2str(dilate_scale));

new_image=max_proj>val;
new_image=bwmorph(new_image,'clean');

se=strel('disk',erode_scale);
new_image=imerode(new_image,se);

se=strel('disk',dilate_scale);
new_image=imdilate(new_image,se);

set(0,'CurrentFigure',fig)
set(h,'AlphaData',new_image.*alpha);

title(val);
setappdata(fig,'threshold',val);

end
