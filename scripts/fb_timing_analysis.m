%%%% script summarizing actions to analyze long calcium extractions

% clean the rois

% useful data:
%	1) LW76, 3/09/14 files 7 and 63
%

%%%%%%%%%%%%%%%%% user parameters %%%%%%%%%%%%%%%%%

load custom_colormaps;
fs=24.414e3; % tdt sampling rate
tmp=ROI_LW76; % rois to use
cut=25; % frames to cut out of front

% mic_filter

[b,a]=ellip(6,.2,80,[500]/(fs/2),'high');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

extractions=fb_quick_template_match(mic_data,'fs',fs);

% assume template matches are saved as the variable extractions (in time)
% get the movie data

[mov_data,frame_idx]=fb_retrieve_mov();
correction=frame_idx(cut-1)./fs;
movie_fs=round((1)./((frame_idx(2)-frame_idx(1))/fs));

% make sure we've removed all duplicates

[curr_roi,include]=fb_clean_roi(tmp);

% assume we've loaded in mic_data from the appropriate file

[roi_traces]=fb_extract_roi_traces(curr_roi,mov_data,frame_idx);

% detrend and make dff traces

dff_detrended=fb_roi_detrend(roi_traces.raw(cut:end,:),'fs',movie_fs,'dff',1,'win',.6);

% compute peaks using fit routine

sample_vec=frame_idx(cut:end);
dff_peaks=fb_compute_peak(dff_detrended,'method','p','debug',1,'onset_only',0,'thresh_hi',1);

corrected_peaks=dff_peaks;
for i=1:length(corrected_peaks)
	if ~isempty(corrected_peaks{i})
		corrected_peaks{i}=sample_vec(corrected_peaks{i})/fs;
	end
end

% TODO: merge files effectively
% grab peaks in each extraction, combine for further analysis

[s,f,t]=pretty_sonogram(filtfilt(b,a,double(mic_data)),fs,'low',1);

figure();imagesc(t,f,s);
set(gca,'TickDir','out')
xlabel('Time (s)');
ylabel('Fs (Hz)');
axis xy;
colormap(fee_map);

for i=1:size(extractions,1)
	fb_visualize_rois(curr_roi,'weights',corrected_peaks,'weights_map','jet','weights_range',...
		[extractions(i,1)/fs-.05 extractions(i,2)/fs-.05],'filled',1,'weights_correction',0);
end







