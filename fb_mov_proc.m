function mov_ave=fb_mov_proc(DIR,varargin)
% fb_mov_proc processes movie data and creates a series of movies and images for
% further analysis (NOTE:  this requires the VideoWriter class, not available 
% in older versions of MATLAB)
%

nparams=length(varargin);
fs=24.414e3;
ave_playback_fs=60; % frame rate of average movie
playback_fs=20;
debug=0; % show debug image?
baseline=3; % 0 for mean, 1 for median, 2 for trimmed mean
filt_rad=15; % disk filter radius
filt_alpha=5;
lims=1;
trim_per=20;
save_dir='proc';
high_pass=0;
activity_map='gray';
cb_height=.03;
motion_correction=1; % 0 for no correction, 1 for correction
motion_crop=20; % allowable crop for motion correction (deleted for later movies)
per=8;
outlier_detect=1;

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
		case 'high_pass'
			high_pass=varargin{i+1};
		case 'save_dir'
			save_dir=varargin{i+1};
		case 'filt_alpha'
			filt_alpha=varargin{i+1};
		case 'motion_correction'
			motion_correction=varargin{i+1};
		case 'outlier_detect'
			outlier_detect=varargin{i+1};
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

ave_time=min_tmp:1/ave_playback_fs*fs:max_tmp;
ave_frames=length(ave_time);
ave_counter=zeros(1,length(ave_time));
mov_ave=zeros(height,width,ave_frames);

h = fspecial('gaussian',filt_rad, filt_alpha); % could use a gaussian as well 

disp('Processing image data');


if debug
	fig=figure('visible','on','renderer','zbuffer','PaperPositionMode','auto');
else
	fig=figure('visible','off','renderer','zbuffer','PaperpositionMode','auto');
end

old_save_dir=save_dir;
motion_flag=0;

for i=1:length(mov_listing)

	save_dir=old_save_dir;
	outlier_flag=0;

	disp(['Processing file ' num2str(i) ' of ' num2str(length(mov_listing))]);
	disp([mov_listing{i}]);

	[path,file,ext]=fileparts(mov_listing{i});

	% compute sonogram of mic_data

	load(fullfile(DIR,mov_listing{i}),'mov_data','mov_idx','frame_idx','mic_data','fs','movie_fs');

	[template_image,f,t]=fb_pretty_sonogram(double(mic_data),fs,'low',1.5,'zeropad',1024,'N',2048,'overlap',2040);
	[rows,columns,frames]=size(mov_data);
	
	mov_filt=imfilter(mov_data,h);

	% first do outlier detection

	
	if outlier_detect

		if ~exist(save_dir,'dir'); mkdir(save_dir); end

		% get the median

		mov_med=median(mov_filt,3);

		% get the mad

		mov_mad=mad(mov_filt,1,3);

		% do any frames exceed the limit

		mask1=mov_med-mov_mad*6;
		mask2=mov_med+mov_mad*6;

		len=rows*columns;

		for j=1:frames

			flag1=mov_filt(:,:,j)<mask1;
			flag2=mov_filt(:,:,j)>mask2;

			flag1=sum(flag1(:)==1)>len*.1;
			flag2=sum(flag2(:)==1)>len*.1;	

			if flag1 || flag2
				warning('Outlier frame %g detected',j);
				outlier_flag=1;
				save_dir='junk';
			end


		end

	end

	if motion_correction && ~outlier_flag

		new_mov_filt=zeros(rows-2*motion_crop,columns-2*motion_crop,frames);

		disp('Performing motion correction...');

		%[nblanks formatstring]=fb_progressbar(100);
		%fprintf(1,['Progress:  ' blanks(nblanks)]);

		% have user select file for motion correction template
		
		if ~motion_flag

			%[filename,pathname]=uigetfile({'*.mat'},'Pick a movie file to extract the template from',pwd);
			%template=load(fullfile(pathname,filename),'mov_data');

			% correction template, select coordinates

			%corr_tmp=imfilter(template.mov_data(:,:,:),h);
			%corr_tmp=mean(corr_tmp,3); % average

			%corr_tmp=corr_tmp(:,:,1); % take the first frame

			disp('Motion correction region selection...');
			disp('Drag the rectangle to enclose the template region and double-click to finish...');

			corr_tmp=mov_filt(:,:,1);

			%corr_tmp(corr_tmp>clim(2))=clim(2);
			%corr_tmp=max(corr_tmp-clim(1),0);
			%corr_tmp=corr_tmp./max(corr_tmp(:));

			overview_fig=figure('Toolbar','none','Menubar','none');
			overview_img=imshow(corr_tmp./max(corr_tmp(:)));
			hold on;	

			overview_scroll=imscrollpanel(overview_fig,overview_img);
			imoverview(overview_img);
			contrast_handle=imcontrast(overview_img);

			uiwait(contrast_handle);
			clims=caxis();

			rect_handle=imrect(get(overview_fig,'CurrentAxes'));

			rect_pos=wait(rect_handle);

			y_segment=round([rect_pos(2):rect_pos(2)+rect_pos(4)]);
			y_segment(y_segment>rows)=[];

			x_segment=round([rect_pos(1):rect_pos(1)+rect_pos(3)]);
			x_segment(x_segment>columns)=[];

			corr_roi=corr_tmp(y_segment,x_segment);

			clim(1)=prctile(corr_roi(:),30);
			clim(2)=prctile(corr_roi(:),70);

			corr_roi(corr_roi>clim(2))=clim(2);
			corr_roi=corr_roi-clim(1);
			corr_roi=max(corr_roi,0);
			corr_roi=corr_roi./clim(2);
			
			corr_roi=corr_roi-mean(corr_roi(:));
			corr_roi=corr_roi./std(corr_roi(:));

			% normalize for motion correction
			
			correction_fft=fft2(corr_roi);
			close(overview_fig);

			motion_flag=1;

		end


		for j=1:frames

			%fprintf(1,formatstring,round((j/frames)*100));	
			tmp=mov_filt(y_segment,x_segment,j);
		
			clim(1)=prctile(tmp(:),30);
			clim(2)=prctile(tmp(:),70);

			tmp(tmp>clim(2))=clim(2);
			tmp=tmp-clim(1);
			tmp=max(tmp,0);
			tmp=tmp./clim(2);

			tmp=tmp-mean(tmp(:));
			tmp=tmp./std(tmp(:));

			[output Greg]=dftregistration(correction_fft,fft2(tmp),100);	
			output

			shift=round(output(3:4));
			shift(shift>motion_crop)=motion_crop;

			mov_filt(:,:,j)=circshift(mov_filt(:,:,j),[shift]);

		end

		new_mov_filt=mov_filt(motion_crop:rows-motion_crop,motion_crop:columns-motion_crop,:);
		mov_filt=new_mov_filt;

		% save the motion corrected data

		save(fullfile(save_dir,[file '_motioncorrected.mat']),'mov_data','mov_idx','frame_idx','mic_data','fs','movie_fs','motion_crop');


	end

	[b,a]=ellip(4,.2,40,[1]/(20/2),'high');

	% warning: high_pass will throw off df/f, would advise not to use for now

	if high_pass==1
		for j=1:rows
			for k=1:columns
				mov_filt(j,k,:)=filtfilt(b,a,squeeze(mov_filt(j,k,:)));
			end
		end

	end




	writer_obj=VideoWriter(fullfile(DIR,save_dir,[file '_dff.avi']));
	writer_obj.FrameRate=playback_fs;

	open(writer_obj);	

	mov_norm=mov_filt;

	if baseline==0
		norm_fact=mean(mov_filt,3);
	elseif baseline==1
		norm_fact=median(mov_filt,3);
	elseif baseline==2
		norm_fact=trimmean(mov_filt,trim_per,'round',3);
	else
		norm_fact=prctile(mov_filt,per,3);
	end

	std_fact=std(mov_filt,[],3);
	clearvars hline;

	set(0,'currentfigure',fig);
    
	subplot(3,1,1);cla;
	imagesc(t,f,template_image);axis xy;hold on;
	colormap(hot);freezeColors();

	if high_pass
		for j=1:frames
			mov_norm(:,:,j)=(mov_norm(:,:,j)./norm_fact)*100;
		end
	else

		% normalize each pixel by its mean over time (median is probably better given
		% the few frames we're working with)

		for j=1:frames
			mov_norm(:,:,j)=((mov_norm(:,:,j)-norm_fact)./norm_fact).*100; % df/f (in percent)
		end
	end

	norm_lim(1)=prctile(mov_norm(:),lims);
	norm_lim(2)=prctile(mov_norm(:),100-lims);

	raw_lim(1)=prctile(mov_filt(:),lims);
	raw_lim(2)=prctile(mov_filt(:),100-lims);

	% dff movie

	movie_t=zeros(1,frames);
	movie_idx=zeros(1,frames);

	for j=1:frames
		[~,idx]=min(abs(frame_idx(j)-ave_time));
		movie_t(j)=ave_time(idx)./fs;
		movie_idx(j)=idx;
	end


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
		imagesc(mov_norm(:,:,j));caxis(norm_lim);
		colormap(activity_map);freezeColors();
		axis off;

		pos=get(ax,'pos');
		set(ax,'pos',[pos(1) pos(2)*1.05 pos(3) pos(4)*0.95]);
		pos=get(ax,'pos');

		hc=colorbar('location','southoutside','position',[pos(1) pos(2)-cb_height*2 pos(3) cb_height]);
		set(hc,'xaxisloc','bottom');

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
		imagesc(mov_filt(:,:,j));caxis(raw_lim);
		colormap(activity_map);freezeColors();
		axis off;

		pos=get(ax,'pos');
		set(ax,'pos',[pos(1) pos(2)*1.05 pos(3) pos(4)*0.95]);
		pos=get(ax,'pos');

		hc=colorbar('location','southoutside','position',[pos(1) pos(2)-cb_height*2 pos(3) cb_height]);
		set(hc,'xaxisloc','bottom');

		frame=im2frame(hardcopy(fig,'-Dzbuffer','-r0'));
		writeVideo(writer_obj,frame);

		if debug==1 pause(.1); end

	end

	fprintf(1,'\n');

	disp('Writing maximum projection');

	close(writer_obj);

	% also save max projection

	max_projection=max(mov_norm,[],3); % maximum df/f across time

	max_lim(1)=prctile(max_projection(:),lims);
	max_lim(2)=prctile(max_projection(:),100-lims);

	newfig=figure('visible','off');
	imagesc(max_projection);caxis([max_lim]);
	colorbar();colormap(gray);

	fb_multi_fig_save(newfig,fullfile(DIR,save_dir),[file '_maxproj'],'png');
	close([newfig]);

	% maybe save the average as a 4-d array so you can remove outliers
	% save the processed data nad 

end

% for average, pre-flight check to make sure camera is on

% establish 20-30 hz time-frame from min to max timestamp, each frame from each file counts once
% (counts toward the closets time point)

end
