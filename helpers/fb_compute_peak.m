function PEAKS=fb_compute_peak(CA_DATA,varargin)
% Computes the peaks for a calcium trace or series of calcium traces
%
%
%
% algorithm:  schmitt trigger, double exponential fit (for now)
%
%

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})	
		case 'roi_map'
			roi_map=varargin{i+1};

	end
end

% ensure formatting is correct

if isvector(CA_DATA), CA_DATA=CA_DATA(:); end
[samples,nrois]=size(CA_DATA);

PEAKS={};

for i=1:nrois
	
	curr_trial=CA_DATA(:,i);
	z_trial=zscore(curr_trial);
	PEAKS{i}=find(z_trial>2);

	% schmitt trigger, integrate, find t0 then fit double exponential (subtract???)

end
