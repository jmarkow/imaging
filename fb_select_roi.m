function EXTRACTED_ROI=fb_select_roi(DIR,varargin)
%fb_select_roi selects an arbitrary number of roi's for plotting
%
%
%
%


% select file to load

nparams=length(varargin);
baseline=2; % 0 for mean, 1 for median, 2 for trimmed mean
filt_rad=15; % disk filter radius
filt_alpha=10;
lims=1;
trim_per=20;
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
		case 'roi_map'
			roi_map=varargin{i+1};
	end
end


if nargin<1 | isempty(DIR), DIR=pwd; end

[filename,pathname]=uigetfile({'*.mat'},'Pick a mat file to extract the image data from',fullfile(DIR,'..'));
load(fullfile(pathname,filename),'mov_data');

mkdir(save_dir);

% maximum projection
% convert mov_data to df/f

h = fspecial('gaussian',filt_rad, filt_alpha); % filter for df/f
mov_data=imfilter(mov_data,h);

if baseline==0
	norm_fact=mean(mov_data,3);
elseif baseline==1
	norm_fact=median(mov_data,3);
else
	norm_fact=trimmean(mov_data,trim_per,'round',3);
end

norm_fact=repmat(mean(mov_data,3),[1 1 size(mov_data,3)]); % now convert to df/f
maxproj=max((mov_data-norm_fact)./norm_fact,[],3); % take maximum projection

% get climits

climits(1)=min(maxproj(:));
climits(2)=max(maxproj(:));

% normalize for initial plot

maxproj=max((maxproj-climits(1))./climits(2),0);
overview_fig=figure('Toolbar','none','Menubar','none');
overview_img=imshow(maxproj);
hold on;

overview_scroll=imscrollpanel(overview_fig,overview_img);
imoverview(overview_img);
contrast_handle=imcontrast(overview_img);

uiwait(contrast_handle);

clims=caxis();

ellipse_handle=imellipse(get(overview_fig,'CurrentAxes'));

% return a cell array with the ROI

EXTRACTED_ROI={};

[rows,columns]=size(maxproj);

[xi,yi]=meshgrid(1:columns,1:rows); % collect all coordinates into xi and yi
exit_flag=0;
counter=1;

save_fig=figure('Visible','off');
imshow(maxproj);caxis([clims]);
hold on;

while 1>0

	h=wait(ellipse_handle);

	% convert the ROI into all pixels inside the polygon
	
	if isempty(h), break; end

	set(0,'CurrentFigure',overview_fig);
	plot(h(:,1),h(:,2),'-','linewidth',1.5,'color',roi_map(counter,:));

	set(0,'CurrentFigure',save_fig);
	plot(h(:,1),h(:,2),'-','linewidth',1.5,'color',roi_map(counter,:));
	

	set(0,'CurrentFigure',overview_fig);

	% what's inside of the ROI?  this could also be used to normalize fluorescene per 
	% Svoboda et al. (take an annulus surrounding the ROI)

	roi=inpolygon(xi,yi,h(:,1),h(:,2));
	[idx]=find(roi);

	% xi=columns yi=rows

	% collect the roi

	EXTRACTED_ROI{end+1}=[yi(idx) xi(idx)];
	counter=counter+1; % increment the colormap


end

fb_multi_fig_save(save_fig,save_dir,'rois','tiff','res',100);
close([save_fig]);
