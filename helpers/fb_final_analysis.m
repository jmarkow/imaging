%%%%%% final analysis of organized data

load custom_colormaps;
mkdir('analysis');

%%%% assume mic_data is loaded in, along with fs, movie_fs and rois

um_per_px=.6250; % um per pixel
timevec=[1:length(align_mic_data)]/fs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% aggregating data

[agg_detrended,agg_peak_locs,agg_peak_vals,agg_ispeak]=fb_aggregate_trials(pwd,'thresh',-inf);

% get the ca data for each ROI only on trials with peaks

[nsamples,nrois,ntrials]=size(agg_detrended);

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
[mergepeaks,mergevals]=fb_compute_peak(agg_pkonly,'thresh_t',.2,'debug',1,'onset_only',0,'thresh_hi',1,'thresh_int',5);

save(fullfile('analysis','agg_data.mat'),'agg_detrended','agg_peak_locs','agg_peak_vals','agg_ispeak','agg_pkonly','mergepeaks');

% remove peaks in the padded regions

startpt=padding(1)-.25;
stoppt=padding(1)+length(TEMPLATE.data)/fs+.25;

for i=1:length(mergepeaks)
	mergepeaks{i}=mergepeaks{i}/movie_fs;
	mergepeaks{i}(mergepeaks{i}<startpt|mergepeaks{i}>stoppt)=[];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% normalizing and peak sorting


[maxca maxloc]=max(agg_pkonly);

% normalize ca traces to [0,1]

[minca]=min(agg_pkonly);

norm_pkonly=agg_pkonly;

lbound=repmat(minca,[nsamples 1]);
ubound=repmat(maxca,[nsamples 1]);

norm_pkonly=(norm_pkonly-lbound)./(ubound-lbound); 

[~,sortidx]=sort(maxloc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% spatial correlations

[dist_t,dist_s,dist_y,dist_bins]=fb_plot_spatialcorrelation_times(newrois,mergepeaks,'bin_spacing',50);
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
[lambda]=fb_plot_symbols([bin_c].*um_per_px,dist_y);

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
[s,f,t]=fb_pretty_sonogram(filtfilt(b,a,align_mic_data),fs,'n',4096,'overlap',4090,'zeropad',0,'low',1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sonogram with peak sorted traces

fig1=figure();
subplot(2,1,1);
imagesc(t,f,s);
axis xy;
colormap(fee_map);
box off;
set(gca,'TickDir','out');
ylim([0 10e3]);
xlabel('Time (s)');
ylabel('Fs (Hz)');
freezeColors();

subplot(2,1,2);
imagesc(norm_pkonly(:,sortidx)');
set(gca,'xtick',[]);
ylabel('Cell');
colormap(hot);
h=colorbar();
caxis([0 1]);
set(h,'position',[.92 .13 .025 .3])
% add colorbar

set(gcf,'position',[64 708 311 359]);
set(gcf,'PaperPositionMode','auto');
multi_fig_save(gcf,'analysis','sorted_peaks','eps,png,fig');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sonogram with roi marked by peak time

left_edge=0;
right_edge=ceil(timevec(end)*10)/10;

fig2=figure();
subplot(6,1,1:2);
imagesc(t,f,s);
axis xy;
colormap(fee_map);
box off;
set(gca,'TickDir','out');
ylim([0 10e3]);
%xlabel('Time (s)');
%ylabel('Fs (Hz)');
freezeColors();

subplot(6,1,3:5);
fb_visualize_rois(newrois,'weights',mergepeaks,'filled',1,'weights_map','jet',...
	'ref_image',0,'fig_num',fig2,'weights_scale',[startpt stoppt])
axis on;set(gca,'xtick',[],'ytick',[]);

subplot(6,1,6);
pts=linspace(startpt,stoppt,10e3);
imagesc(pts,[],pts);
colormap(jet)
set(gca,'YTick',[]);

set(gcf,'position',[190 607 309 408]);
set(gcf,'PaperPositionMode','auto');
multi_fig_save(gcf,'analysis','roi_weighted','eps,png,fig');


