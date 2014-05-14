function roi_ave=fb_plot_roi(ROIS,varargin)
%fb_select_roi selects an arbitrary number of roi's for plotting
%
%
%
%
%

colors=colormap('lines');
sono_colormap='hot';
baseline=3;
ave_fs=100;
save_dir='roi';
template=[];
fs=24.414e3;
per=2;
max_row=5;
filt_rad=30;
filt_alpha=10;
min_f=0;
max_f=9e3;
lims=1;
dff_scale=20;
t_scale=2;

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})	
		case 'colors'
			colors=varargin{i+1};
		case 'sono_colormap'
			sono_colormap=varargin{i+1};
		case 'baseline'
			baseline=varargin{i+1};
		case 'ave_fs'
			ave_fs=varargin{i+1};
		case 'save_dir'
			save_dir=varargin{i+1};
		case 'template'
			template=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'per'
			per=vargin{i+1};
		case 'max_row'
			max_row=varargin{i+1};
		case 'filt_alpha'
			filt_alpha=varargin{i+1};
		case 'filt_rad'
			filt_rad=varargin{i+1};
		case 'dff_scale'
			dff_scale=varargin{i+1};
		case 't_scale'
			t_scale=varargin{i+1};
	end
end


mkdir(save_dir);

% ROIS is a cell array of image indices returned by fb_select_roi
%

% first convert ROIS to row and column indices, then average ROI and plot
% the time course

% TODO save image with ROIs (using color scheme that's used for time plots)

mov_listing=dir(fullfile(pwd,'*.mat'));
mov_listing={mov_listing(:).name};
roi_n=length(ROIS);

load(fullfile(pwd,mov_listing{1}),'mov_data','mic_data','fs');
[rows,columns,frames]=size(mov_data);

ave_time=0:1/ave_fs:length(mic_data)/fs;

% need to interpolate the average onto a new time bases

roi_ave.raw={};
roi_ave.interp_dff=zeros(roi_n,length(ave_time),length(mov_listing));
roi_ave.interp_raw=zeros(roi_n,length(ave_time),length(mov_listing));
disp('Generating single trial figures...');

for i=1:length(mov_listing)

	disp(['Processing file ' num2str(i) ' of ' num2str(length(mov_listing))]);
	load(fullfile(pwd,mov_listing{i}),'mov_data','mov_idx','frame_idx','mic_data','fs');

	[path,file,ext]=fileparts(mov_listing{i});
	save_file=[ file '_roi' ];

	[b,a]=ellip(5,.2,80,[500]/(fs/2),'high');

	[song_image,f,t]=fb_pretty_sonogram(filtfilt(b,a,double(mic_data)),fs,'low',1.5,'zeropad',1024,'N',2048,'overlap',2000);	

	% roi_traces

	[rows,columns,frames]=size(mov_data);
	roi_t=zeros(roi_n,frames);

	timevec=frame_idx./fs;

	disp('Computing ROI averages...');

	[nblanks formatstring]=fb_progressbar(100);
	fprintf(1,['Progress:  ' blanks(nblanks)]);

	% unfortunately we need to for loop by frames, otherwise
	% we'll eat up too much RAM for large movies

	for j=1:roi_n
		fprintf(1,formatstring,round((j/roi_n)*100));

		for k=1:frames
			tmp=mov_data(ROIS{j}(:,1),ROIS{j}(:,2),k);
			roi_t(j,k)=mean(tmp(:));
		end
	end

	fprintf(1,'\n');

	dff=zeros(size(roi_t));

	for j=1:roi_n

		tmp=roi_t(j,:);

		if baseline==0
			norm_fact=mean(tmp,3);
		elseif baseline==1
			norm_fact=median(tmp,3);
		elseif baseline==2
			norm_fact=trimmean(tmp,trim_per,'round',3);
		else
			norm_fact=prctile(tmp,per);
		end

		dff(j,:)=((tmp-norm_fact)./norm_fact).*100;

		yy=interp1(timevec,dff(j,:),ave_time,'spline');
		yy2=interp1(timevec,tmp,ave_time,'spline');

		roi_ave.interp_dff(j,:,i)=yy;
		roi_ave.interp_raw(j,:,i)=yy2;

	end

	detrended=fb_roi_detrend(dff,timevec,'normalize',0);
	save_fig=figure('visible','off');box off;
	set(save_fig,'paperpositionmode','auto','position',[100 100 450 900])

	ax(1)=subplot(7,1,1:2);
	imagesc(t,f./1e3,song_image);axis xy;
	xlabel('Time (in s)');
	ylabel('Freq. (kHz)');
	ylim([min_f/1e3 max_f/1e3]);
	set(gca,'TickDir','out');
	colormap(sono_colormap);freezeColors();

	ax(2)=subplot(7,1,3:7);
	fb_plot_roi_stackplot(detrended,timevec,'spacing',8,'colors',colors);axis off;
	line([-.25 -.25],[0 dff_scale],'clipping','off','linewidth',2.5);
	line([0 t_scale],[-5 -5],'clipping','off','linewidth',2.5);	
	linkaxes(ax,'x');

	% loop through each ROI, average data

	% break into a grid of columns and rows to plot >5 rois

	% plot each roi in a subplot, perhaps use subaxis to keep things close together
	
	fb_multi_fig_save(save_fig,save_dir,save_file,'eps,png,fig','res',100);
	save(fullfile(save_dir,[save_file '.mat']),'roi_t','frame_idx','fs','timevec');

	roi_ave.raw{i}=roi_t; % store for average
	roi_ave.filename{i}=mov_listing{i};
	% 1d interpolate all rois to common frame

end

save(fullfile(save_dir,['ave_roi.mat']),'ave_time','roi_ave');
disp('Generating average ROI figure...');

% plot the averages with confidence intervals

%timevec=ave_time;

% if template is passed use the template mic trace, otherwise use the last song

%roi_mu=mean(roi_ave.interp_dff,3);
%roi_sem=std(roi_ave.interp_dff,[],3)./sqrt(size(roi_ave.interp_dff,3));

%if ~isempty(template)
%	[song_image,f,t]=fb_pretty_sonogram(double(template),fs,'low',1.5,'zeropad',1024,'N',2048,'overlap',2040);	
%end

%fb_multi_fig_save(save_fig,save_dir,'ave_roi','eps,png,fig','res',100);

%%%%%%%%%%%%% CELL MASK MATCHED TO ROI
% plot cell masks color-matched to their ROIs 
