function [DFF_MAT,DFF_COMBINE]=fb_batch_maxdff(varargin)
%fb_select_roi selects an arbitrary number of roi's for plotting
%
%
%
%


% select file to load

nparams=length(varargin);
filt_rad=12; % gauss filter radius
filt_alpha=4; % gauss filter alpha
lims=2; % contrast prctile limits
save_dir='dff';
per=2; % baseline percentile (0 for min)
resize_correct=1; % correction of parameters for resized movies
activity_colormap='gray'; % colormap for activity
resize=1;
combine='max';

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'activity_colormap'
			activity_colormap=varargin{i+1};
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
		case 'label_color'
			label_color=varargin{i+1};
		case 'combine'
			combine=varargin{i+1};
	end
end


resize_skip=0;
im_resize=1; % if im_resize does not exist as a variable, the data has not been resized!

[filename,pathname]=uigetfile({'*.mat';'*.tif'},'Pick mat files to extract the image data from',pwd,'Multiselect','on');
%[path,file,ext]=fileparts(filename);

if ~iscell(filename)
	tmp=filename;
	clear filename;
	filename{1}=tmp;
end

nfiles=length(filename);
mkdir(save_dir);

dff_store={};

for ii=1:length(filename)

	[path,file,ext]=fileparts(filename{ii});

	if strcmp(ext,'.mat')
		load(fullfile(pathname,filename{ii}),'mov_data','im_resize');
	else
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

	mov_data=single(mov_data); % convert to single to save memory here
	baseline=repmat(prctile(mov_data,per,3),[1 1 frames]);
	
	dff=((mov_data-baseline)./baseline).*100;
	clear mov_data;
	
	dff=max(dff,[],3);
	dff_clims=prctile(dff(:),[lims 100-lims]);

	% convert to uint8, faster loading for animation

	dff=min(dff,dff_clims(2)); % clip to max
	dff=max(dff-dff_clims(1),0); % clip min
	dff=dff./(dff_clims(2)-dff_clims(1)); % normalize to [0,1]
	dff=uint8(dff.*255);

	imwrite(dff,gray(256),fullfile(save_dir,[ file '_maxdff.tiff' ]),'tiff');

	% write dff to image, save to mat, and we're done

	dff_store{ii}=dff;

end

DFF_MAT=cat(3,dff_store{:});
clear dff_store;

DFF_COMBINE.ave=mean(DFF_MAT,3);
DFF_COMBINE.max=max(DFF_MAT,[],3);

imwrite(DFF_COMBINE.ave,gray(256),fullfile(save_dir,[ 'ave_maxdff.tiff' ]),'tiff');
imwrite(DFF_COMBINE.max,gray(256),fullfile(save_dir,[ 'max_maxdff.tiff' ]),'tiff');

save(fullfile(save_dir,[ 'dff_data.mat' ]),'DFF_MAT','DFF_COMBINE','filename','-v7.3');
