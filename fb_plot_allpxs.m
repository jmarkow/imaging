function fb_plot_allpxs(MOV_DATA,MIC_DATA,FRAME_IDX,varargin)
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
filt_rad=60; % gauss filter radius
filt_alpha=20; % gauss filter alpha
lims=5; % contrast prctile limits
cmap=colormap('jet');
save_dir='roi';
per=0; % baseline percentile (0 for min)
fs=24.414e3;

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
		case 'per'
			per=varargin{i+1};
		case 'lims'
			lims=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'cmap'
			cmap=varargin{i+1};
	end
end


[rows,columns,frames]=size(MOV_DATA);

h=fspecial('gaussian',filt_rad,filt_alpha);
MOV_DATA=imfilter(MOV_DATA,h,'circular');

baseline=repmat(prctile(MOV_DATA,per,3),[1 1 frames]);
dff=((MOV_DATA-baseline)./baseline).*100;

% convert frames to times

frame_idx=FRAME_IDX./fs;

% take the center of mass across dim 3 (time) for each point in space

com_idx=zeros(1,1,frames);

for i=1:frames
	com_idx(:,:,i)=i;
end

com_idx=repmat(com_idx,[rows columns 1]);

mass=sum(dff,3);
com_dff=sum((dff.*com_idx),3)./mass;

max_proj=max(dff,[],3);

clims(1)=prctile(max_proj(:),lims);
clims(2)=prctile(max_proj(:),100-lims);

norm_max_proj=min(max_proj,clims(2));
norm_max_proj=max(norm_max_proj-clims(1),0);
norm_max_proj=norm_max_proj./(clims(2)-clims(1));

% map to [0,1] for ind2rgb

cm1=cmap;
%cm1=colormap('jet');

clims(1)=min(com_dff(:));
clims(2)=max(com_dff(:));

norm_dff=min(com_dff,clims(2)); % clip to max
norm_dff=max(norm_dff-clims(1),0); % clip min
norm_dff=norm_dff./(clims(2)-clims(1)); % normalize to [0,1]

im1_rgb=ind2rgb(round(norm_dff.*size(cm1,1)),cm1);

figure();
h=image(im1_rgb);
set(h,'AlphaData',norm_max_proj);
set(gca,'color',[0 0 0]);

