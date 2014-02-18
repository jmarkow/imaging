function [EXTRACTED_ROI,STATS]=fb_auto_roi(DIR,varargin)
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
filt_rad=30; % gauss filter radius
filt_alpha=10; % gauss filter alpha
lims=2; % contrast prctile limits
roi_map=colormap('lines');
save_dir='roi';
per=2; % baseline percentile (0 for min)
ave_scale=40; % for adaptive threshold, scale for averaging filter
med_scale=22; % for removing speckle noise from maximum projection
resize_correct=1; % correction of parameters for resized movies

% parameters used for morphological opening

erode_scale=3; % scale (in pxs) for erosion
dilate_scale=3; % scale (in pxs) for dilation

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

[filename,pathname]=uigetfile({'*.mat'},'Pick a mat file to extract the image data from',fullfile(DIR,'..'));
load(fullfile(pathname,filename),'mov_data','im_resize');

if resize_correct & im_resize~=1

	disp('Correcting parameters since file has been downsampled...');
	ave_scale=round(ave_scale.*im_resize);
	med_scale=round(med_scale.*im_resize);
	filt_rad=round(filt_rad.*im_resize);
	filt_alpha=filt_alpha.*im_resize;

end

%if exist(save_dir,'dir') rmdir(save_dir,'s'); end

mkdir(save_dir);

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

clear mov_data;

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
max_proj=max(dff,[],3);

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

EXTRACTED_ROI={};

[rows,columns]=size(max_proj);

slmin=min(max_proj(:));
slmax=max(max_proj(:));

cm1=colormap('gray');

% convert the indexed images to rgb

im1_rgb=ind2rgb(round(raw_proj.*size(cm1,1)),cm1);
im2_rgb=ind2rgb(max_proj>0,[0 0 0;1 0 0]);

% display the anatomical "max projection"

slider_fig=figure();
image(im1_rgb);

% now add the mask, change the transparency using the slider

hold on;
h=image(im2_rgb);
set(h,'AlphaData',(max_proj>0));

hsl = uicontrol(slider_fig,'Style','slider','Min',slmin,'Max',slmax,...
          'SliderStep',[1e-3 1e-3],'Value',slmin,...
             'Position',[20 10 200 20]);
set(hsl,'Callback',{@slider_callback,slider_fig,max_proj,h,erode_scale,dilate_scale})

disp('Press any key in the command window to use current threshold');
pause();
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

STATS=regionprops(conn_comp,'eccentricity','majoraxislength',...
	'minoraxislength','convexhull','centroid','boundingbox');

% get coordinates and draw ROIs on maximum projection

save_fig=figure('Visible','off');
imshow(raw_proj);
colormap(gray);
hold on;

for i=1:length(conn_comp.PixelIdxList);

	[xi,yi]=ind2sub(size(max_proj),conn_comp.PixelIdxList{i});
	EXTRACTED_ROI{i}=[xi(:) yi(:)]; % get the coordinates
	tmp=STATS(i).ConvexHull;
	plot(tmp(:,1),tmp(:,2),'-','linewidth',2,'color',roi_map(i,:));
	tmp=STATS(i).BoundingBox;

	x=tmp(1)+tmp(3)/2;
	y=tmp(2)+tmp(4)/2;

	text(x,y,[num2str(i)],'FontSize',12,'FontName','Helvetica','color','r','FontWeight','bold');
end

fb_multi_fig_save(save_fig,save_dir,'roi_map','tiff','res','100');
close(save_fig);
save(fullfile(save_dir,'roi_data.mat'),'EXTRACTED_ROI','STATS');

end

function slider_callback(hObject,eventdata,fig,max_proj,h,erode_scale,dilate_scale)

alpha=.4;

val=get(hObject,'Value');

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
