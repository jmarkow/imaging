%%%%%% final analysis of organized data

mkdir('analysis');

%%%% assume mic_data is loaded in, along with fs, movie_fs and rois

um_per_px=.7; % um per pixel
timevec=[1:length(align_mic_data)]/fs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% aggregating data

[agg_detrended,agg_peak_locs,agg_peak_vals,agg_ispeak]=fb_aggregate_trials(pwd);

% get the ca data for each ROI only on trials with peaks

[nsamples,nrois,ntrials]=size(agg_detrended);

agg_pkonly=[];
for i=1:length(agg_ispeak)
	if ~isempty(agg_ispeak{i})
		agg_pkonly=[agg_pkonly mean(agg_ispeak{i},2)];
	end
end

mergelist=fb_merge_peaks(agg_peak_locs,agg_peak_vals,'win',2,'thresh',-inf);
mergepeaks={};

for i=1:length(mergelist)
	mergepeaks{i}=cellfun(@mean,mergelist{i})./movie_fs;
end

save(fullfile('analysis','agg_data.mat'),'agg_detrended','agg_peak_locs','agg_peak_vals','agg_ispeak','mergepeaks');

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

[dist_t,dist_s,dist_y,dist_bins]=fb_plot_spatialcorrelation_times(rois,mergepeaks,'bin_spacing',100);
bin_c=[];
for i=1:length(dist_bins)-1
	bin_c(i)=mean(dist_bins(i:i+1));
end

len=cellfun(@length,dist_y);
cutoff=min(find(len(2:end)<5))+1;
[lambda]=fb_plot_symbols([bin_c(1:cutoff)].*um_per_px,dist_y(1:cutoff));

xlabel('Distance (um)');
ylabel('Peak separation (s)');

save(fullfile('analysis','spatial_corr.mat'),'dist_t','dist_s','dist_y','dist_bins');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sonogram

[b,a]=ellip(5,.2,40,[300]/(fs/2),'high');
[s,f,t]=fb_pretty_sonogram(filtfilt(b,a,align_mic_data),fs,'n',4096,'overlap',4090,'zeropad',0);

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
caxis([0 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sonogram with roi marked by peak time

left_edge=0;
right_edge=ceil(timevec(end)*10)/10;

fig2=figure();
subplot(8,1,1:2);
imagesc(t,f,s);
axis xy;
colormap(fee_map);
box off;
set(gca,'TickDir','out');
ylim([0 10e3]);
%xlabel('Time (s)');
%ylabel('Fs (Hz)');
freezeColors();

subplot(8,1,3:7);
fb_visualize_rois(rois,'weights',mergepeaks,'filled',1,'weights_map','jet',...
	'ref_image',1,'fig_num',fig2,'weights_scale',[0 right_edge])
axis off

subplot(8,1,8);
pts=linspace(left_edge,right_edge,1e3);
imagesc(pts,[],pts);
colormap(jet)
set(gca,'YTick',[]);


