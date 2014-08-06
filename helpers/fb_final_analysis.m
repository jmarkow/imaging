%%%%%% final analysis of organized data

load custom_colormaps;
mkdir('analysis');

%%%% assume mic_data is loaded in, along with fs, movie_fs and rois

um_per_px=.6250; % um per pixel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% aggregating data

if ~exist(fullfile('analysis','agg_data.mat'),'file')
	[agg_detrended,agg_ispeak]=fb_aggregate_trials(pwd,'thresh',-inf);


	% get the ca data for each ROI only on trials with peaks


	agg_pkonly=[];
	include=[];
	for i=1:length(agg_ispeak)
		if ~isempty(agg_ispeak{i})
			agg_pkonly=[agg_pkonly mean(agg_ispeak{i},2)];
			include=[include i];
		end
	end

	%mergelist=fb_merge_peaks(agg_peak_locs,agg_peak_vals,'win',2,'thresh',-inf);
	%mergepeaks={};
	%
	%for i=1:length(mergelist)
	%	mergepeaks{i}=cellfun(@mean,mergelist{i})./movie_fs;
	%end
	%

	newrois=rois;
	newrois.coordinates=newrois.coordinates(include);
	newrois.stats=newrois.stats(include);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% new peak selection
	[mergepeaks,mergevals]=fb_compute_peak_simple(agg_pkonly,...
		'thresh_t',.2,'debug',1,'onset_only',0,'thresh_hi',1,'thresh_int',8,'thresh_dist',.2); % thresh_int previously 5

	save(fullfile('analysis','agg_data.mat'),'agg_detrended','agg_ispeak','agg_pkonly','mergepeaks','newrois');
else
	load(fullfile('analysis','agg_data.mat'));
end

% remove peaks in the padded regions

[nsamples,nrois,ntrials]=size(agg_detrended);
startpt=padding(1);
stoppt=padding(1)+length(TEMPLATE.data)/fs;

for i=1:length(mergepeaks)
	mergepeaks{i}=mergepeaks{i}/movie_fs;
	mergepeaks{i}(mergepeaks{i}<startpt|mergepeaks{i}>stoppt)=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% normalizing and peak sorting

pkinclude=find(cellfun(@length,mergepeaks)>0);
plotrois=newrois;
plotrois.coordinates=plotrois.coordinates(pkinclude);
plotrois.stats=plotrois.stats(pkinclude);

agg_pkonly=agg_pkonly(:,pkinclude);

movie_idx=round([startpt stoppt].*movie_fs);

[maxca maxloc]=max(agg_pkonly(movie_idx(1):movie_idx(2),:));

% normalize ca traces to [0,1]

[minca]=min(agg_pkonly);

norm_pkonly=agg_pkonly;

lbound=repmat(minca,[nsamples 1]);
ubound=repmat(maxca,[nsamples 1]);

norm_pkonly=(norm_pkonly-lbound)./(ubound-lbound); 

[nframes,nrois]=size(norm_pkonly);
%[~,sortidx]=sort(maxloc);
%[~,sortidx]=sort(pktime);

interp_fs=20;
interp_x=[0:1/interp_fs:nframes]/movie_fs;
movie_x=[0:nframes-1]/movie_fs;
interp_pkonly=interp1(1:nframes,agg_pkonly,[1:1/interp_fs:nframes]);

[interp_nsamples,~]=size(interp_pkonly);

[minca]=min(interp_pkonly);
[maxca maxloc]=max(interp_pkonly);

norm_interp=interp_pkonly;

lbound=repmat(minca,[interp_nsamples 1]);
ubound=repmat(maxca,[interp_nsamples 1]);

norm_interp=(norm_interp-lbound)./(ubound-lbound); 

[~,sortidx]=sort(maxloc);

% bin by prctile
%

interp_pkonly_bins=zeros(size(interp_pkonly));

for i=1:nrois

	% get prctile cutoffs
	%
	
	bins=prctile(interp_pkonly(:,i),[0:25:100]);
	[n,binned]=histc(interp_pkonly(:,i),bins);

	binned
	length(binned)
	interp_pkonly_bins(:,i)=binned;


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% spatial correlations

[dist_t,dist_s,dist_y,dist_bins]=fb_plot_spatialcorrelation_times(newrois,mergepeaks,'bins',0:200:2e3); % good results no pad, 125 um spacing
save(fullfile('analysis','spatial_corr.mat'),'dist_t','dist_s','dist_y','dist_bins');

bin_c=[];
for i=1:length(dist_bins)-1
	bin_c(i)=mean(dist_bins(i:i+1));
end

dist_y=dist_y(1:end-1);
len=cellfun(@length,dist_y);
%cutoff=min(find(len(2:end)<5)+1);

dist_y(len<5)=[];
bin_c(len<5)=[];

if length(dist_y)>3
	[lambda]=fb_plot_symbols([bin_c],dist_y);

	xlabel('Distance (um)');
	ylabel('Peak separation (s)');

	set(gca,'TickDir','out');
	set(gcf,'Position',[457 883 263 223]);
	set(gcf,'PaperPositionMode','auto');

	multi_fig_save(gcf,'analysis','spatial_correlation','eps,png,fig');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sonogram

[b,a]=ellip(5,.2,40,[700]/(fs/2),'high');
align_mic_data(1:startpt*fs)=0;
align_mic_data(stoppt*fs:end)=0;
[s,f,t]=fb_pretty_sonogram(filtfilt(b,a,align_mic_data),fs,'n',4096,'overlap',4090,'zeropad',0,'low',1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sonogram with peak sorted traces

ax=[];
fig1=figure();
ax(1)=subplot(2,1,1);
imagesc(t,f,s);
axis xy;
colormap(fee_map);
box off;
set(gca,'TickDir','out');
ylim([0 9e3]);
xlabel('Time (s)');
ylabel('Fs (Hz)');
freezeColors();

ax(2)=subplot(2,1,2);
imagesc(interp_x,[],norm_interp(:,sortidx)');
set(gca,'xtick',[]);
ylabel('Cell');
colormap(hot);
h=colorbar();
caxis([0 1]);
set(h,'position',[.92 .13 .025 .3])
% add colorbar

set(gcf,'position',[64 708 311 359]);
set(gcf,'PaperPositionMode','auto');
linkaxes(ax,'x');
%xlim([startpt stoppt]);

multi_fig_save(gcf,'analysis','sorted_peaks','eps,png,fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sonogram with roi marked by peak time

for i=1:length(mergepeaks)
	if ~isempty(mergepeaks{i})
		minpeaks{i}=min(mergepeaks{i});
	else
		minpeaks{i}=[];
	end
end

fig2=figure();
subplot(6,1,1:2);
imagesc(t,f,s);
axis xy;
colormap(fee_map);
box off;
set(gca,'TickDir','out');
ylim([0 9e3]);
xlim([startpt stoppt]);

%xlabel('Time (s)');
%ylabel('Fs (Hz)');

freezeColors();

subplot(6,1,3:5);
fb_visualize_rois(newrois,'weights',minpeaks,'filled',1,'weights_map','jet',...
	'ref_image',[],'fig_num',fig2,'weights_scale',[startpt stoppt])
axis on;set(gca,'xtick',[],'ytick',[]);

subplot(6,1,6);
pts=linspace(startpt,stoppt,10e3);
imagesc(pts,[],pts);
colormap(jet)
set(gca,'YTick',[]);

set(gcf,'position',[190 607 309 408]);
set(gcf,'PaperPositionMode','auto');
multi_fig_save(gcf,'analysis','roi_weighted','eps,png,fig');


