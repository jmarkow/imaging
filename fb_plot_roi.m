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
ave_fs=80;
save_dir='roi';
template=[];
fs=24.414e3;
per=2;

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
roi_ave.interp=zeros(roi_n,length(ave_time),length(mov_listing));

disp('Generating single trial figures...');

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for i=1:length(mov_listing)

	fprintf(1,formatstring,round((i/length(mov_listing))*100));

	load(fullfile(pwd,mov_listing{i}),'mov_data','mov_idx','frame_idx','mic_data','fs');

	[path,file,ext]=fileparts(mov_listing{i});
	save_file=[ file '_roi' ];

	[song_image,f,t]=fb_pretty_sonogram(double(mic_data),fs,'low',1.5,'zeropad',1024,'N',2048,'overlap',2040);	

	mov_norm=mov_data;

	%for j=1:frames
	%	mov_norm(:,:,j)=((mov_norm(:,:,j)-norm_fact)./norm_fact).*100; % df/f (in percent)
	%end

	% roi_traces

	[rows,columns,frames]=size(mov_data);
	roi_t=zeros(roi_n,frames);

	for j=1:frames
		for k=1:roi_n
			tmp=mov_norm(ROIS{k}(:,1),ROIS{k}(:,2),j);
			roi_t(k,j)=mean(tmp(:));
		end
	end

	% plot each roi in a subplot, perhaps use subaxis to keep things close together

	nplots=roi_n+1;
	timevec=frame_idx./fs;

	ax=[];
	save_fig=figure('visible','off');
	set(save_fig,'position',[100 100 500 200*nplots],'PaperPositionMode','auto')
	clf;
	
	ax(1)=subaxis(nplots,1,1,1,1,1,'spacingvert',.012,'marginbottom',.15);

	imagesc(t,f./1e3,song_image);axis xy;box off;
	colormap(sono_colormap);
	ylabel('Fs (kHz)');

	set(gca,'TickDir','out','linewidth',1,'FontSize',12);
	set(gca,'xcolor',get(gcf,'color'),'xtick',[]);

	for j=1:nplots-1

		size(timevec)
		size(roi_t(j,:))

		ax(j+1)=subaxis(nplots,1,1,j+1,1,1,'spacingvert',.012,'marginbottom',.15);
		set(gca,'TickDir','out','linewidth',1,'FontSize',12);
		plot(timevec,roi_t(j,:),'color',colors(j,:));
		box off;

		if j<nplots-1
			set(gca,'xcolor',get(gcf,'color'),'xtick',[]);
		end
		
		axis tight;

	end

	linkaxes(ax,'x');box off;
	xlim([timevec(1) timevec(end)]);
	xlabel('Time (in s)');

	fb_multi_fig_save(save_fig,save_dir,save_file,'eps,png,fig','res',100);
	close([save_fig]);
	save(fullfile(save_dir,[save_file '.mat']),'roi_t','frame_idx');

	roi_ave.raw{i}=roi_t; % store for average

	% 1d interpolate all rois to common frame

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

		dff=((tmp-norm_fact)./norm_fact).*100;

		yy=interp1(frame_idx./fs,dff,ave_time,'spline');
		roi_ave.interp(j,:,i)=yy;
	end

	% loop through each ROI, average data

end


fprintf(1,'\n');
disp('Generating average ROI figure...');

% plot the averages with confidence intervals

nplots=roi_n+1;
timevec=ave_time;

% if template is passed use the template mic trace, otherwise use the last song

roi_mu=mean(roi_ave.interp,3);
roi_sem=std(roi_ave.interp,[],3)./sqrt(size(roi_ave.interp,3));

if ~isempty(template)
	[song_image,f,t]=fb_pretty_sonogram(double(mic_data),fs,'low',1.5,'zeropad',1024,'N',2048,'overlap',2040);	
end


save_fig=figure('visible','off');
set(save_fig,'position',[100 100 500 200*nplots],'PaperPositionMode','auto','visible','off')
clf;

ax=[];
ax(1)=subaxis(nplots,1,1,1,1,1,'spacingvert',.012,'marginbottom',.15);

imagesc(t,f./1e3,song_image);axis xy;box off;
colormap(sono_colormap);

ylabel('Fs (Hz)');

set(gca,'TickDir','out','linewidth',1,'FontSize',12);
set(gca,'xcolor',get(gcf,'color'),'xtick',[]);

for i=1:nplots-1

	ax(i+1)=subaxis(nplots,1,1,i+1,1,1,'spacingvert',.012,'marginbottom',.15);
	set(gca,'TickDir','out','linewidth',1,'FontSize',12);

	% plot confidence bars
	%

	plot(timevec,roi_mu(i,:),'color',colors(i,:));
	hold on;
	plot(timevec,roi_mu(i,:)+roi_sem(i,:),'k--','color',colors(i,:));
	plot(timevec,roi_mu(i,:)-roi_sem(i,:),'k--','color',colors(i,:));

	box off;

	if i<nplots-1
		set(gca,'xcolor',get(gcf,'color'),'xtick',[]);
	end	

	axis tight;

end

box off;
linkaxes(ax,'x');
xlabel('Time (in s)');
xlim([timevec(1) timevec(end)]);
fb_multi_fig_save(save_fig,save_dir,'ave_roi','eps,png,fig','res',100);
save(fullfile(save_dir,['ave_roi.mat']),'ave_time','roi_ave','roi_mu','roi_sem');

close([save_fig]);

