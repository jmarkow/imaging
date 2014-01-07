function EXTRACTED_ROI=fb_auto_roi(DIR,varargin)
%fb_select_roi selects an arbitrary number of roi's for plotting
%
%
%
%


% select file to load

nparams=length(varargin);
baseline=2; % 0 for mean, 1 for median, 2 for trimmed mean
filt_rad=30; % disk filter radius
filt_alpha=20;
lims=1;
trim_per=20;
save_dir='proc';
roi_map=colormap('lines');
save_dir='auto_roi';

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

mov_filt=ones(1,3)*1/3;

%for i=1:rows
%	i
%	for j=1:columns
%		mov_data(i,j)=conv(mov_data(i,j),mov_filt,'same');
%	end
%end

mov_data=imfilter(mov_data,h);
baseline=repmat(prctile(mov_data,10,3),[1 1 frames]);

mov_data=((mov_data-baseline)./baseline).*100;
max_proj=max(mov_data,[],3);
max_proj=medfilt2(max_proj,[10 10]); % get rid of salt and pepper noise
max_proj=max(max_proj./max(max_proj(:)),0); % convert to [0,1]

% need to smooth in time as well

EXTRACTED_ROI={};

[rows,columns]=size(max_proj);

[xi,yi]=meshgrid(1:columns,1:rows); % collect all coordinates into xi and yi
exit_flag=0;
counter=1;

slmin=min(max_proj(:));
slmax=max(max_proj(:));

slider_fig=figure();
imagesc(max_proj>0);

hsl = uicontrol(slider_fig,'Style','slider','Min',slmin,'Max',slmax,...
          'SliderStep',[1e-3 1e-3],'Value',slmin,...
             'Position',[20 10 200 20]);
set(hsl,'Callback',{@slider_callback,slider_fig,max_proj})

pause();
threshold=get(hsl,'value');

% convert threshold to [0,1)

threshold=threshold;
new_image=im2bw(max_proj,threshold);
conn_comp=bwconncomp(new_image);

for i=1:length(conn_comp.PixelIdxList);
	[xi,yi]=ind2sub(size(max_proj),conn_comp.PixelIdxList{i});
	EXTRACTED_ROI{i}=[xi(:) yi(:)];
end



end

function slider_callback(hObject,eventdata,fig,max_proj)

val=get(hObject,'Value');
set(0,'CurrentFigure',fig)
imagesc(max_proj>val)
title(val);
setappdata(fig,'threshold',val);

end
