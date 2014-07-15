function PEAKS=fb_compute_peak(CA_DATA,varargin)
% Computes the peaks for a calcium trace or series of calcium traces
%
%
%
% algorithm:  schmitt trigger, double exponential fit (for now)
%
%

thresh_hi=2;
thresh_lo=1;
thresh_t=100;
fs=200;

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})	
		case 'roi_map'
			roi_map=varargin{i+1};
		case 'thresh_hi'
			thresh_hi=varargin{i+1};
		case 'thresh_lo'
			thresh_lo=varargin{i+1};
		case 'thresh_t'
			thresh_t=varargin{i+1};
		case 'fs'
			fs=varargin{i+1};
	end
end

% ensure formatting is correct

if isvector(CA_DATA), CA_DATA=CA_DATA(:); end
[samples,nrois]=size(CA_DATA);

PEAKS={};

idx=1:samples-1;
a=optimset('MaxFunEvals',1e3);

% options for optimization algorithms
% not all options are used for all algorithms
maxIter=1e3;
options.Display = 'off';
options.MaxIter = maxIter;
options.MaxIter = Inf;
options.UseParallel = 'always';
options.ObjectiveLimit = 0;
% options.TimeLimit = 10; % in s / default is Inf

% experimental
options.MeshAccelerator = 'on'; % off by default
options.TolFun = 1e-9; % default is 1e-6
options.TolMesh = 1e-9; % default is 1e-6
options.TolX = 1e-9; % default is 1e-6

for i=1:nrois

	% find first threshold crossing

	% breakpoints in .4 steps

	curr_roi=CA_DATA(:,i);
	curr_roi=curr_roi-smooth(curr_roi,round(1.5*fs));
	pos_crossing=find(curr_roi(idx)<thresh_hi&curr_roi(idx+1)>thresh_hi);

	% take the first crossing

	if isempty(pos_crossing)
		continue;
	end

	init_guess=pos_crossing(1);

	schmitt_flag=0;
	if init_guess+thresh_t<samples
		schmitt_flag=all(curr_roi(init_guess:init_guess+thresh_t)>thresh_lo);
	end


	% attempt to fit the double exponential model, fit A, onset time, and tau

	figure(1);plot(curr_roi);

	schmitt_flag
	init_guess

	init_guess=init_guess./fs;
	
	%[x]=fminsearch(@(x) obj_function(x,curr_roi,fs),[ 1 init_guess .03 .2 ],a);
	x=simulannealbnd(@(x) obj_function(x,curr_roi,fs),... 
		[ 1 init_guess .03 .2 ] ,[ 0 init_guess-.04 .002 .05 ],[ 5 init_guess+.05 .03 .3 ],options);

	time_vec=[1:length(curr_roi)]./fs;

	y=calcium_model(x(1),x(2),x(3),x(4),time_vec);

	hold on;	
	plot(y,'r-');

	pause();
	cla;


end

end

function res = obj_function(x,dff,fs)
%
%
%
%

A=x(1);
t_o=x(2);
t_on=x(3);
t_1=x(4);

time_vec=[1:length(dff)]./fs;

y=calcium_model(A,t_o,t_on,t_1,time_vec);
res=sum((y(:)-dff(:)).^2)

end




function y = calcium_model(A,t_o,t_on,t_1,t)
%
%
%
%

y = A.*(1-exp(-(t-t_o)./t_on)).*exp(-(t-t_o)./t_1);
y(t<t_o)=0;

end
