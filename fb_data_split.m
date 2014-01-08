function fb_data_split(MIC_SIGNAL,SYNC_SIGNAL,TTL_SIGNAL,varargin)
%fb_data_split organizes a directory of movie files using the movie sync 
%signal.  
%
%	fb_data_split(MIC_SIGNAL,SYNC_SIGNAL,TTL_SIGNAL,varargin)
%
%	MIC_SIGNAL
%	microphone signal
%
%	SYNC_SIGNAL
%	movie sync signal
%
%	TTL_SIGNAL
%	custom file creation TTL signal (not currently used)
%	
%	the following parameters may be specified using parameter/value pairs:
%
%		high
%		value for positive-going threshold crossing (float, default: 200)
%
%		low 
%		value for negative-going threshold crossing (float, default: 200)
%
%		gap_thresh
%		gap in sync signal that indicates a recording stop (int, default: 10 samples)
%
%		fs
%		sampling frequency of TDT (float, default: 24.414e3)
%		
%		mat_dir
%		where to store mat files with mic and image data (string, default: 'mat_dir')
%
%		gif_dir
%		where to store sonograms (string, default:'git_dir')
%
%
%%%%%%

% logic:
% 1) find "long" gaps in the sync signal, this specifies points where
%    song has not been detected, hence the camera is off
% 2) the threshold for what constitutes a long gap has been determined empirically.
%    When the camera is on, gaps are on the order of 250 µs (6-8 samples at 24.4414 kHz)
% 3) take the rising edge after a long gap as the onset of the file
% 4) take the falling edge before a long gap as the offset of the file
% 5) make sure the number of file onsets matches the number of .tif files in the directory
% 6) match the onset number to the .tif file number
% 7) split the data along the onset points

% parameter collection

nparams=length(varargin);

high=200;
low=200;
gap_thresh=10; % how long of a gap in the sync signal for file splitting?
fs=24.414e3;
mat_dir='mat';
gif_dir='gif';

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'high'
			high=varargin{i+1};
		case 'low'
			low=varargin{i+1};
		case 'gap_thresh'
			gap_thresh=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
		case 'mat_dir'
			mat_dir=varargin{i+1};
		case 'gif_dir'
			gif_dir=varargin{i+1};
	end
end

mkdir(mat_dir);
mkdir(gif_dir);

% TTL signal isn't used as the moment, just put in a dummy signal in case
% we add it later

if nargin<3 | isempty(TTL_SIGNAL), TTL_SIGNAL=zeros(size(SYNC_SIGNAL)); end 

idx=1:length(SYNC_SIGNAL)-1;

rising_edges=find(SYNC_SIGNAL(idx)<low&SYNC_SIGNAL(idx+1)>high);
rising_edges=rising_edges(:);

falling_edges=find(SYNC_SIGNAL(idx)>high&SYNC_SIGNAL(idx+1)<low);
falling_edges=unique([1;falling_edges(:);length(SYNC_SIGNAL)]);

size(rising_edges)
size(falling_edges)

falling_edges=falling_edges(1:length(rising_edges)+1);

% first rising edge is the first data split

disp('Detecting gaps in sync signal...');

gap_lengths=rising_edges-falling_edges(1:end-1);

file_onsets=find(gap_lengths>gap_thresh);

file_offsets=falling_edges(file_onsets);
file_onsets=rising_edges(file_onsets); % where are the file onsets?

file_offsets=[file_offsets;length(SYNC_SIGNAL)];

tif_list=dir(fullfile(pwd,'*.tif')); % list of files

% if number of tif>file_onsets then movie software was left on, if vice versa
% then TDT was left on

if length(tif_list)~=length(file_onsets)
	warning('Number of files (%g) should equal the number of detected intervals (%g)',...
		length(tif_list),length(file_onsets));
end

% parse trial number from filenames

trial_num=zeros(length(tif_list),1);

% should deal with all filename variants Bill's thrown at me so far...

for i=1:length(tif_list)

	tmp=tif_list(i).name;

	tokens=regexp(tmp,'\-','split');

	if length(tokens)==1
		tokens=regexp(tmp,'\_','split');
		trial_num(i)=str2num(tokens{3});
	else

		tokens2=regexp(tokens{2},'\_','split');

		if isempty(str2num(tokens2{1}))
			tokens2=regexp(tokens{1},'\_','split');
			trial_num(i)=str2num(tokens2{3});
		else
			trial_num(i)=str2num(tokens2{1});
		end

	end

end

% TODO check to see if we can align data by shifting trial numbers by 
% 1-5 in either direction

[val,idx]=sort(trial_num,'ascend');

to_del=find(val>length(file_onsets));

val(to_del)=[]; % if we have more trial numbers than file onsets, trim to prevent errors
idx(to_del)=[];

filenames={tif_list(:).name};
filenames=filenames(idx);

disp('Checking data integrity...');

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for i=1:length(filenames)

	fprintf(1,formatstring,round((i/length(filenames))*100));

	[path,file,ext]=fileparts(filenames{i});
	idxs_segment=file_onsets(val(i)):file_offsets(val(i)+1);

	% number of frames in sync signal from onset to offset

	sync_data=SYNC_SIGNAL(idxs_segment);
	mic_data=MIC_SIGNAL(idxs_segment);
	ttl_data=TTL_SIGNAL(idxs_segment);

	% how many rising edges in the sync signal?

	idx=1:length(sync_data)-1;

	rising=find(sync_data(idx)<low&sync_data(idx+1)>high);
	falling=find(sync_data(idx)>high&sync_data(idx+1)<low);

	% for each rising edge is there a nearby falling edge?

	if isempty(falling)
		continue;
	end

	% get the number of frames to check data integrity

	image_info=imfinfo(filenames{i});
	frame_num=length(image_info);

	% save the rising edges and falling edges (save as point process, 1 indicates
	% rising edge or falling edge)

	rising_data=zeros(size(sync_data));
	rising_data(rising)=1:length(rising);
	
	falling_data=zeros(size(sync_data));
	falling_data(falling)=1:length(falling);


	if frame_num~=length(rising)
	
		warning('Frame number (%g) not equal to number of rising edges (%g) in file %s',...
			frame_num,length(rising),filenames{i});

		% save data to junk directory or skip if the frame number does not match

	else

		% if the frame numbers match the number of rising edges then
		% filter the mic data, generate a sonogram, and store everything

		[b,a]=ellip(5,.2,80,[500]/(fs/2),'high');
		plot_data=mic_data./abs(max(mic_data));

		[s,f,t]=fb_pretty_sonogram(filtfilt(b,a,mic_data./abs(max(mic_data))),fs,'low',1.5,'zeropad',0);

		minpt=1;
		maxpt=min(find(f>=10e3));

		imwrite(flipdim(s(minpt:maxpt,:),1),hot,fullfile(gif_dir,[file '.gif']),'gif');
		save(fullfile(mat_dir,[file '.mat']),'sync_data','mic_data','ttl_data','rising_data','falling_data','fs');

	end

end

fprintf('\n');



