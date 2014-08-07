function [MAXDFF]=fb_maxdff(DATA,varargin)
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
lims=2; % contrast prctile limits
roi_map='lines';
save_dir='auto_roi';
per=5; % baseline percentile (0 for min)
ave_scale=80; % for adaptive threshold, scale for averaging filter
resize_correct=1; % correction of parameters for resized movies

background_subtract=0;
% parameters used for morphological opening

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
		case 'ave_scale'
			ave_scale=varargin{i+1};
		case 'per'
			per=varargin{i+1};
		case 'lims'
			lims=varargin{i+1};
	end
end


im_resize=1; % if im_resize does not exist as a variable, the data has not been resized!

%TODO: allow for multiple file selection, concatenate movies

if resize_correct & im_resize~=1

	disp('Correcting parameters since file has been downsampled...');
	ave_scale=round(ave_scale.*im_resize);
	filt_rad=round(filt_rad.*im_resize);
	filt_alpha=filt_alpha.*im_resize;

end

disp('Filtering images, this may take a minute...');

[rows,columns,frames]=size(DATA);
h=fspecial('gaussian',filt_rad,filt_alpha);
%h=fspecial('disk',filt_rad);

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for j=1:frames
	fprintf(1,formatstring,round((j/frames)*100));	
	DATA(:,:,j)=imfilter(DATA(:,:,j),h,'circular');
end

fprintf(1,'\n');


baseline=repmat(prctile(DATA,per,3),[1 1 frames]);
dff=((DATA-baseline)./baseline).*100;

if background_subtract

	dff_mu=zeros(size(dff));

	h=fspecial('average',ave_scale);

	[nblanks formatstring]=fb_progressbar(100);
	fprintf(1,['Progress:  ' blanks(nblanks)]);

	for j=1:frames
		fprintf(1,formatstring,round((j/frames)*100));	
		dff_mu(:,:,j)=imfilter(dff(:,:,j),h,'circular');
	end

	fprintf(1,'\n');

	dff=dff-dff_mu;

end

MAXDFF=max(dff,[],3);
%MAXDFF=mean(dff,3);
clims=prctile(MAXDFF(:),[lims 100-lims]);

% basic morphological operators before conncomp 
%
%
MAXDFF=min(MAXDFF,clims(2));
MAXDFF=max(MAXDFF,clims(1))-clims(1);

MAXDFF=round(MAXDFF./(clims(2)-clims(1)).*255);
