function [EXTRACTED_ROI,STATS]=fb_auto_roi(DIR,varargin)
%fb_select_roi selects an arbitrary number of roi's for plotting
%
%
%
%


% select file to load

nparams=length(varargin);
baseline=2; % 0 for mean, 1 for median, 2 for trimmed mean
filt_rad=30; % disk filter radius
filt_alpha=10;
lims=1;
save_dir='proc';
roi_map=colormap('lines');
save_dir='roi';

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
	end
end


if nargin<1 | isempty(DIR), DIR=pwd; end

[filename,pathname]=uigetfile({'*.mat'},'Pick a mat file to extract the image data from',fullfile(DIR,'..'));
load(fullfile(pathname,filename),'mov_data');

mkdir(save_dir);

% maximum projection
% convert mov_data to df/f

%mov_data=imresize(mov_data,.25);
%ccimage=xcorr_image(mov_data); %WAY TOO SLOW IF IMAGE IS NOT SUBSAMPLED

% convert to df/f then take maximum projection

[rows,columns,frames]=size(mov_data);

h=fspecial('gaussian',filt_rad,filt_alpha);

mov_data=imfilter(mov_data,h);
baseline=repmat(prctile(mov_data,8,3),[1 1 frames]);

dff=((mov_data-baseline)./baseline).*100;
max_proj=max(dff,[],3);
max_proj=max(max_proj./max(max_proj(:)),0); % convert to [0,1]
max_proj=medfilt2(max_proj,[20 20]); 

raw_proj=max(mov_data,[],3);
clims(1)=prctile(raw_proj(:),lims);
clims(2)=prctile(raw_proj(:),100-lims);

EXTRACTED_ROI={};

[rows,columns]=size(max_proj);

slmin=min(max_proj(:));
slmax=max(max_proj(:));

slider_fig=figure();
imagesc(max_proj>0);

hsl = uicontrol(slider_fig,'Style','slider','Min',slmin,'Max',slmax,...
          'SliderStep',[1e-3 1e-3],'Value',slmin,...
             'Position',[20 10 200 20]);
set(hsl,'Callback',{@slider_callback,slider_fig,max_proj})

disp('Press any key to use current threshold');
pause();
threshold=get(hsl,'value');
close(slider_fig);

% convert threshold to [0,1)

threshold=threshold;
new_image=im2bw(max_proj,threshold);

% basic morphological operators before conncomp 

new_image=bwmorph(new_image,'fill');
new_image=bwmorph(new_image,'clean');
new_image=bwmorph(new_image,'dilate');

% collect some basic stats if we want to exclude later

conn_comp=bwconncomp(new_image);
STATS=regionprops(conn_comp,'eccentricity','majoraxislength','minoraxislength','convexhull');

% get coordinates and draw ROIs on maximum projection

save_fig=figure('Visible','off');
imshow(raw_proj);caxis(clims);
colormap(gray);
hold on;

for i=1:length(conn_comp.PixelIdxList);
	[xi,yi]=ind2sub(size(max_proj),conn_comp.PixelIdxList{i});
	EXTRACTED_ROI{i}=[xi(:) yi(:)]; % get the coordinates
	tmp=STATS(i).ConvexHull;
	plot(tmp(:,1),tmp(:,2),'-','linewidth',2,'color',roi_map(i,:));
end

fb_multi_fig_save(save_fig,save_dir,'roi_map','tiff','res','100');
close(save_fig);
save(fullfile(save_dir,'roi_data.mat'),'EXTRACTED_ROI','STATS');

end

function slider_callback(hObject,eventdata,fig,max_proj)

val=get(hObject,'Value');
set(0,'CurrentFigure',fig)
imagesc(max_proj>val)
title(val);
setappdata(fig,'threshold',val);

end
