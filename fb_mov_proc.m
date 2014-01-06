function mov_ave=fb_mov_proc(DIR,varargin)
% fb_mov_proc processes movie data and creates a series of movies and images for
% further analysis (NOTE:  this requires the VideoWriter class, not available 
% in older versions of MATLAB)
%

nparams=length(varargin);
fs=24.414e3;
ave_movie_fs=60; % frame rate of average movie
movie_fs=5;
debug=0; % show debug image?
baseline=0; % 0 for mean, 1 for median, 2 for trimmed mean
filt_rad=30; % disk filter radius
filt_alpha=25;
lims=1;
trim_per=20;
save_dir='proc';
high_pass=0;
activity_map='gray';
cb_height=.03;

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})	
		case 'fs'
			fs=varargin{i+1};
		case 'movie_fs'
			movie_fs=varargin{i+1};
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
		case 'ave_movie_fs'
			ave_movie_fs=varargin{i+1};
		case 'high_pass'
			high_pass=varargin{i+1};
		case 'save_dir'
			save_dir=varargin{i+1};
		case 'filt_alpha'
			filt_alpha=varargin{i+1};
	end
end

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

ave_time=min_tmp:1/ave_movie_fs*fs:max_tmp;
ave_frames=length(ave_time);
ave_counter=zeros(1,length(ave_time));
mov_ave=zeros(height,width,ave_frames);

h = fspecial('gaussian',filt_rad, filt_alpha); % could use a gaussian as well 

disp('Processing image data');

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);


if debug
	fig=figure('visible','on');
else
	fig=figure('visible','off');
end


for i=1:length(mov_listing)

	[path,file,ext]=fileparts(mov_listing{i});

	writer_obj=VideoWriter(fullfile(DIR,save_dir,[file '_dff.avi']));
	writer_obj.FrameRate=movie_fs;

	open(writer_obj);

	fprintf(1,formatstring,round((i/length(mov_listing))*100));

	% compute sonogram of mic_data

	load(fullfile(DIR,mov_listing{i}),'mov_data','mov_idx','frame_idx','mic_data');

	[template_image,f,t]=fb_pretty_sonogram(double(mic_data),fs,'low',1.5,'zeropad',1024,'N',2048,'overlap',2040);	
	[rows,columns,frames]=size(mov_data);

	mov_filt=imfilter(mov_data,h);
	%mov_filt=mov_data;

	[b,a]=ellip(4,.2,40,[1]/(20/2),'high');

	% warning: high_pass will throw off df/f

	if high_pass==1
		for j=1:rows
			for k=1:columns
				mov_filt(j,k,:)=filtfilt(b,a,squeeze(mov_filt(j,k,:)));
			end
		end

	end

	mov_norm=mov_filt;

	if baseline==0
		norm_fact=mean(mov_filt,3);
	elseif baseline==1
		norm_fact=median(mov_filt,3);
	else
		norm_fact=trimmean(mov_filt,trim_per,'round',3);
	end

	std_fact=std(mov_filt,[],3);
	clearvars hline;
	
	set(0,'currentfigure',fig);
	
	figure(fig);
	subplot(3,1,1);cla;
	imagesc(t,f,template_image);axis xy;hold on;
	colormap(hot);freezeColors();

	%mov_filt=hist_eq(mov_filt,255);

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

	for j=1:frames

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

		frame=getframe(fig);
		writeVideo(writer_obj,frame);

		if debug==1 pause(.1); end


	end

	close(writer_obj);

	writer_obj=VideoWriter(fullfile(DIR,save_dir,[file '_raw.avi']));
	writer_obj.FrameRate=movie_fs;

	open(writer_obj);

	% raw movie

	for j=1:frames

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

		frame=getframe(fig);
		writeVideo(writer_obj,frame);

		if debug==1 pause(.1); end

	end

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

	for j=1:frames
		idx=movie_idx(j);
		mov_ave(:,:,idx)=mov_ave(:,:,idx)+mov_filt(:,:,j);
		ave_counter(idx)=ave_counter(idx)+1;
	end
end

fprintf(1,'\n');

for i=1:size(mov_ave,3)
	mov_ave(:,:,i)=mov_ave(:,:,i)./ave_counter(i);
end


% for average, pre-flight check to make sure camera is on

% establish 20-30 hz time-frame from min to max timestamp, each frame from each file counts once
% (counts toward the closets time point)

end

function new_movie = hist_eq(matrix,max_val)
	%
	%
%

matrix=uint8(matrix);
vals=matrix(:);
bins=0:max_val;

density=hist(vals,bins);
density=density./sum(density);

cdf=cumsum(density);

new_movie=cdf(matrix+1);

max(new_movie)
min(new_movie)

new_movie=uint8(new_movie*max_val);

end
