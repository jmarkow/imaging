function [DETRENDED,BASELINE]=fb_roi_detrend(ROI,T,varargin)
%
%
%
%
%
%

% parameter collection

nparams=length(varargin);
per=0;
normalize=1;
cutoff=.2;

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'per'
			per=varargin{i+1};
		case 'normalize'
			normalize=varargin{i+1};
	end
end

ave_fs=1./(diff(T(1:2)));

% filter out slow drift

BASELINE=prctile(ROI,per,2);
[b,a]=butter(5,[cutoff]/(ave_fs/2),'high');
%[b,a]=ellip(6,.2,50,[cutoff]/(ave_fs/2),'high');

% assume each row is a new roi, each column is a timepoint

DETRENDED=filtfilt(b,a,ROI')';
if normalize
	DETRENDED=DETRENDED./repmat(BASELINE,[1 size(DETRENDED,2)]);
end


