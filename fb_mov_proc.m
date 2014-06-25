function fb_mov_proc(DIR,varargin)
% fb_mov_proc processes movie data and creates a series of movies and images for
% further analysis (NOTE:  this requires the VideoWriter class, not available 
% in older versions of MATLAB)
%

nparams=length(varargin);
fs=24.414e3; % frame rate for TDT data
ave_playback_fs=60; % frame rate of average movie
playback_fs=20; % frame fs for single trial movies
debug=0; % show debug image?
baseline=3; % 0 for mean, 1 for median, 2 for trimmed mean
filt_rad=15; % Gauss filter size
filt_alpha=5; % Gauss filter variance
lims=1; % prctile limits for contrast enhancement (trims <lims and >lims prctile)
trim_per=20; % prctile for trimmed mean (trim_per/2 from each end is trimmed)
save_dir='proc'; % save directory
activity_map='gray'; % colormap for activity
cb_height=.03; % colorbar height (in normalized units)
motion_correction=1; % 0 for no correction, 1 for correction
per=0; % percentile to use for baseline compution (robust minimum)
outlier_detect=1; % 0 for no outlier detection, 1 for detection
outlier_mads=6; % number of median absolute deviations from the median a pixel must be
		% to be designated an outlier
outlier_frac=.1; % fraction of pixels that must be outliers in a single frame to reject a
		 % movie
junk_dir='junk';
resize_correct=1; % correct parameters if movie has been downsampled
im_resize=1; % change to resize the motion-corrected movies (similar to im_resize in fb_template_match.m
motion_auto=1;

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})	
		case 'fs'
			fs=varargin{i+1};
		case 'playback_fs'
			playback_fs=varargin{i+1};
		case 'debug'
			debug=varargin{i+1};
		case 'baseline'
			baseline=varargin{i+1};
		case 'filt_rad'
			filt_rad=varargin{i+1};
		case 'lims'
			lims=varargin{i+1};
		case 'trim_per'
			trim_per=varargin{i+1};
		case 'ave_playback_fs'
			ave_playback_fs=varargin{i+1};
		case 'save_dir'
			save_dir=varargin{i+1};
		case 'filt_alpha'
			filt_alpha=varargin{i+1};
		case 'motion_correction'
			motion_correction=varargin{i+1};
		case 'outlier_detect'
			outlier_detect=varargin{i+1};
		case 'outlier_mads'
			outlier_mads=varargin{i+1};
		case 'outlier_frac'
			outlier_frac=varargin{i+1};
		case 'junk_dir'
			junk_dir=varargin{i+1};
		case 'resize_correct'
			resize_correct=varargin{i+1};
		case 'im_resize'
			im_resize=varargin{i+1};
		case 'motion_auto'
			motion_auto=varargin{i+1};
	end
end

outlier_flag=0;

if nargin<1 | isempty(DIR), DIR=pwd; end

if ~exist(fullfile(DIR,save_dir),'dir') 
	mkdir(fullfile(DIR,save_dir)); 
else
	rmdir(fullfile(DIR,save_dir),'s');
	mkdir(fullfile(DIR,save_dir));
end

if ~exist(fullfile(DIR,junk_dir),'dir')
	mkdir(fullfile(DIR,junk_dir));
else
	rmdir(fullfile(DIR,junk_dir),'s');
	mkdir(fullfile(DIR,junk_dir));
end

mov_listing=dir(fullfile(DIR,'*.mat'));
mov_listing={mov_listing(:).name};

load(fullfile(DIR,mov_listing{1}),'mov_data');

[height,width,frames]=size(mov_data);

min_time=inf;
max_time=-inf;

% get the first and last frame of the average movie, loop through movie files

disp('Determining first and last frame of the average movie...');

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for i=1:length(mov_listing)

	fprintf(1,formatstring,round((i/length(mov_listing))*100));

	load(fullfile(DIR,mov_listing{i}),'frame_idx')

	min_tmp=min(frame_idx);
	max_tmp=max(frame_idx);

	if min_tmp<min_time, min_time=min_tmp; end
	if max_tmp>max_time, max_time=max_tmp; end

end

fprintf(1,'\n');


if debug
	fig=figure('visible','on','renderer','zbuffer','PaperPositionMode','auto');
else
	fig=figure('visible','off','renderer','zbuffer','PaperpositionMode','auto');
end

old_save_dir=save_dir;
motion_flag=0;

ave_time=min_tmp:1/ave_playback_fs*fs:max_tmp;
ave_frames=length(ave_time);

filt_rad1=filt_rad;
filt_alpha1=filt_alpha;

dff={};
im_resize_store=im_resize;

for i=1:length(mov_listing)

	save_dir=old_save_dir;
	outlier_flag=0;

	disp(['Processing file ' num2str(i) ' of ' num2str(length(mov_listing))]);
	disp([mov_listing{i}]);

	[path,file,ext]=fileparts(mov_listing{i});

	% compute sonogram of mic_data

	load(fullfile(DIR,mov_listing{i}),'mov_data','mov_idx','frame_idx','mic_data','fs','movie_fs','im_resize');

	[rows,columns,frames]=size(mov_data);
	
	if im_resize_store~=1

		[rows,columns]=size(imresize(mov_data(:,:,1),im_resize_store));
		new_mov_data=zeros(rows,columns,frames);

		for j=1:frames
			new_mov_data(:,:,j)=imresize(mov_data(:,:,j),im_resize_store);
			im_resize=im_resize*im_resize_store;
		end

		mov_data=new_mov_data;
		clear new_mov_data;

	end

	if resize_correct & im_resize~=1
		disp('Correcting parameters since file has been downsampled...');
		filt_rad=max(round(filt_rad1.*im_resize),1);
		filt_alpha=max(filt_alpha1.*im_resize,1);
	end
		
	h = fspecial('gaussian',filt_rad, filt_alpha); % could use a gaussian as well 

	[template_image,f,t]=fb_pretty_sonogram(double(mic_data),fs,'low',1.5,'zeropad',1024,'N',2048,'overlap',2000);

	% first do outlier detection
	
	if outlier_detect

		disp('Entering outlier detection...');

		% get the median

		mov_med=median(mov_data,3);

		% get the mad

		mov_mad=mad(mov_data,1,3);

		% do any frames exceed the limit

		mask1=mov_med-mov_mad*outlier_mads;
		mask2=mov_med+mov_mad*outlier_mads;

		len=rows*columns;

		for j=1:frames

			flag1=mov_data(:,:,j)<mask1;
			flag2=mov_data(:,:,j)>mask2;

			flag1=sum(flag1(:)==1)>len*outlier_frac;
			flag2=sum(flag2(:)==1)>len*outlier_frac;	

			if flag1 || flag2
				warning('Outlier frame %g detected',j);
				outlier_flag=1;
				save_dir=junk_dir;
				break;
			end


		end

	end

	if motion_correction && ~outlier_flag

		disp('Performing motion correction...');

		% have user select file for motion correction template
		
		if ~motion_flag

			% correction template, select coordinates
			% mean works for now, could get more complicated in the future

			corr_roi=mean(mov_data,3);
		
			% for now use entire image, could consider cropping

			y_segment=1:rows;
			x_segment=1:columns;
	
			correction_fft=fft2(corr_roi);

			motion_flag=1;

		end

		[nblanks formatstring]=fb_progressbar(100);
		fprintf(1,['Progress:  ' blanks(nblanks)]);

		for j=1:frames

			fprintf(1,formatstring,round((j/frames)*100));	
			
			% uncomment to test for motion correction, introduces random x,y shifts

			%tmp=circshift(mov_data(y_segment,x_segment,j),[randi(100) randi(100)]);
			
			tmp=mov_data(y_segment,x_segment,j);
	
			% last argument is upsample factor

			[output Greg]=dftregistration(correction_fft,fft2(tmp),100);
			mov_data(:,:,j)=abs(ifft2(Greg)); % recover corrected image

		end

		% save the motion corrected data

		save(fullfile(save_dir,[file '_motioncorrected.mat']),'mov_data','mov_idx',...
			'frame_idx','mic_data','fs','movie_fs','im_resize','motion_correction','-v7.3');
		
		fprintf(1,'\n');

		clear mic_data;

	end

	mov_filt=zeros(size(mov_data));

	disp('Filtering image data...');

	[nblanks formatstring]=fb_progressbar(100);
	fprintf(1,['Progress:  ' blanks(nblanks)]);

	for j=1:frames
		fprintf(1,formatstring,round((j/frames)*100));	
		mov_filt(:,:,j)=imfilter(mov_data(:,:,j),h);
	end

	fprintf(1,'\n');

	clear mov_data;
	mov_filt=single(mov_filt);

	if baseline==0
		norm_fact=mean(mov_filt,3);
	elseif baseline==1
		norm_fact=median(mov_filt,3);
	elseif baseline==2
		norm_fact=trimmean(mov_filt,trim_per,'round',3);
	elseif baseline==3
		norm_fact=prctile(mov_filt,per,3);
	else
		norm_fact=min(mov_filt,[],3);
	end

	norm_fact=repmat(norm_fact,[1 1 frames]);

	disp('Computing df/f...');
	
	mov_norm=((mov_filt-norm_fact)./norm_fact).*100;

	clearvars hline;

	set(0,'currentfigure',fig);
    
	subplot(3,1,1);cla;
	imagesc(t,f,template_image);axis xy;hold on;
	colormap(hot);freezeColors();

	% normalize each pixel by its mean over time (median is probably better given
	% the few frames we're working with)

	% average the prctiles across frames, otherwise too costly an operation (would require flattening the whole
	% whole volume, sucking up way too much RAM)

	disp('Writing maximum projection');

	% also save max projection

	max_projection=max(mov_norm,[],3); % maximum df/f across time

	max_projection_clims=prctile(max_projection(:),[lims 100-lims]);
	max_projection=min(max_projection,max_projection_clims(2)); % clip to max
	max_projection=max(max_projection-max_projection_clims(1),0); % clip min
	max_projection=max_projection./(max_projection_clims(2)-max_projection_clims(1)); % normalize to [0,1]
	max_projection=uint8(max_projection.*255);

	imwrite(max_projection,gray(256),fullfile(save_dir,[file '_maxproj.tiff']),'tiff');

	norm_lim=prctile(mov_norm(:),[lims 100-lims]);
	raw_lim=prctile(mov_filt(:),[lims 100-lims]);

	% dff movie

	movie_t=zeros(1,frames);
	movie_idx=zeros(1,frames);

	for j=1:frames
		[~,idx]=min(abs(frame_idx(j)-ave_time));
		movie_t(j)=ave_time(idx)./fs;
		movie_idx(j)=idx;
	end

	mov_filt=min(mov_filt,raw_lim(2)); % clip to max
	mov_filt=max(mov_filt-raw_lim(1),0); % clip min
	mov_filt=mov_filt./(raw_lim(2)-raw_lim(1)); % normalize to [0,1]
	mov_filt=uint8(mov_filt.*255);

	mov_norm=min(mov_norm,norm_lim(2)); % clip to max
	mov_norm=max(mov_norm-norm_lim(1),0); % clip min
	mov_norm=mov_norm./(norm_lim(2)-norm_lim(1)); % normalize to [0,1]
	mov_norm=uint8(mov_norm.*255);

	if ~outlier_flag
		dff{end+1}=max_projection;
	end

	writer_obj=VideoWriter(fullfile(DIR,save_dir,[file '_dff.avi']));
	writer_obj.FrameRate=playback_fs;

	open(writer_obj);	

	disp('Writing df/f movie...');

	[nblanks formatstring]=fb_progressbar(100);
	fprintf(1,['Progress:  ' blanks(nblanks)]);

	for j=1:frames

		fprintf(1,formatstring,round((j/frames)*100));

		set(0,'currentfigure',fig);

		if exist('hline','var'), delete(hline); end

		frame_t=movie_t(j);

		subplot(3,1,1);
		hline=plot([frame_t;frame_t],[f(end);f(1)],'r--','linewidth',2);
		title(['Time ' num2str(frame_t)]);

		ax=subplot(3,1,2:3);
		image(mov_norm(:,:,j));
		colormap([ activity_map '(256)' ]);freezeColors();
		axis off;

		pos=get(ax,'pos');
		set(ax,'pos',[pos(1) pos(2)*1.05 pos(3) pos(4)*0.95]);
		pos=get(ax,'pos');

		hc=colorbar('location','southoutside','position',[pos(1) pos(2)-cb_height*2 pos(3) cb_height]);
		set(hc,'xaxisloc','bottom');
		set(hc,'xtick',[1 256],'xticklabel',[ norm_lim ],'TickLength',[0 0]);

		%frame=getframe(fig);
		frame=im2frame(hardcopy(fig,'-Dzbuffer','-r0'));
		writeVideo(writer_obj,frame);

		if debug==1 pause(.1); end


	end

	fprintf(1,'\n');

	close(writer_obj);

	writer_obj=VideoWriter(fullfile(DIR,save_dir,[file '_raw.avi']));
	writer_obj.FrameRate=playback_fs;

	open(writer_obj);

	% raw movie

	disp('Writing raw movie...');

	[nblanks formatstring]=fb_progressbar(100);
	fprintf(1,['Progress:  ' blanks(nblanks)]);

	for j=1:frames

		fprintf(1,formatstring,round((j/frames)*100));


		set(0,'currentfigure',fig);

		if exist('hline','var'), delete(hline); end

		frame_t=movie_t(j);

		subplot(3,1,1);
		hline=plot([frame_t;frame_t],[f(end);f(1)],'r--','linewidth',2);
		title(['Time ' num2str(frame_t)]);

		ax=subplot(3,1,2:3);
		image(mov_filt(:,:,j));
		colormap([ activity_map '(256)' ]);freezeColors();
		axis off;

		pos=get(ax,'pos');
		set(ax,'pos',[pos(1) pos(2)*1.05 pos(3) pos(4)*0.95]);
		pos=get(ax,'pos');

		hc=colorbar('location','southoutside','position',[pos(1) pos(2)-cb_height*2 pos(3) cb_height]);
		set(hc,'xaxisloc','bottom');
		set(hc,'xtick',[1 256],'xticklabel',[ raw_lim ],'TickLength',[0 0]);

		frame=im2frame(hardcopy(fig,'-Dzbuffer','-r0'));
		writeVideo(writer_obj,frame);

		if debug==1 pause(.1); end

	end

	close(writer_obj);
	fprintf(1,'\n');

	% maybe save the average as a 4-d array so you can remove outliers
	% save the processed data nad 


end

% for average, pre-flight check to make sure camera is on

% establish 20-30 hz time-frame from min to max timestamp, each frame from each file counts once
% (counts toward the closets time point)

DFF_MAT=cat(3,dff{:});
clear dff;

DFF_COMBINE.ave=mean(DFF_MAT,3);
DFF_COMBINE.max=max(DFF_MAT,[],3);

imwrite(DFF_COMBINE.ave,gray(256),fullfile(save_dir,[ 'ave_maxdff.tiff' ]),'tiff');
imwrite(DFF_COMBINE.max,gray(256),fullfile(save_dir,[ 'max_maxdff.tiff' ]),'tiff');

save(fullfile(save_dir,[ 'dff_data.mat' ]),'DFF_MAT','DFF_COMBINE','-v7.3');

end
