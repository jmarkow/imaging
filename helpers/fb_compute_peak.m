function PEAKS=fb_compute_peak(CA_DATA,varargin)
% Computes the peaks for a calcium trace or series of calcium traces
%
%
%
% algorithm:  schmitt trigger, double exponential fit (for now)
%
%

thresh_hi=1.5;
thresh_lo=0;
thresh_t=40;
fs=200;
method='f'; % f-min, simulated annealing, pattern search, etc.
max_iter=1e3; % maximum iterations for optimization
t_1=.07
spk_delta=.04;

onset_init_guess= [ 1 .03 ];
onset_lbound= [ 0 .002  ];
onset_hbound= [ 10 .04 ];

full_init_guess= [ 1 1 .1 .1 ];
full_lbound= [ 0 0 .03 .03 ];
full_hbound= [ 10 10 .5 .5 ];

onset_only=1;

debug=1;

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
		case 'method'
			method=varargin{i+1};
		case 't_1'
			t_1=varargin{i+1};
		case 'max_iter'
			max_iter=varargin{i+1};
		case 'onset_init_guess'
			onset_init_guess=varargin{i+1};
		case 'onset_lbound'
			onset_lbound=varargin{i+1};
		case 'onset_only'
			onset_only=varargin{i+1};
		case 'debug'
			debug=varargin{i+1};
	end
end

% ensure formatting is correct

if isvector(CA_DATA), CA_DATA=CA_DATA(:); end
[samples,nrois]=size(CA_DATA);

PEAKS={};

idx=1:samples-1;
a=optimset('MaxFunEvals',1e3);

maxIter=1e3;
options.Display = 'off';
options.MaxIter = max_iter;
options.UseParallel = 'always';
options.ObjectiveLimit = 0;

[b,a]=butter(2,.5/(fs/2),'low');

for i=1:nrois

	PEAKS{i}=[];

	% find first threshold crossing

	% breakpoints in .4 steps

	curr_roi=CA_DATA(:,i);

	filt_size=round(1.5*fs);

	if exist('smooth')==2
		smooth_roi=smooth(curr_roi,filt_size);
	else	
		smooth_roi=conv(curr_roi,ones(1,filt_size)./filt_size,'same');
	end

	curr_roi=curr_roi-smooth_roi;

	% get the positive threshold crossings

	pos_crossing=find(curr_roi(idx)<thresh_hi&curr_roi(idx+1)>thresh_hi);

	% take the first crossing

	if isempty(pos_crossing)
		continue;
	end

	if curr_roi(1)>thresh_hi
		pos_crossing=[ 1;pos_crossing ];
	end

	schmitt_flag=zeros(1,length(pos_crossing));
	for j=1:length(pos_crossing)

		init_guess=pos_crossing(j);

		% do we stay above the low threshold for a sufficient amount of time?

		if init_guess+thresh_t<samples
			schmitt_flag(j)=all(curr_roi(init_guess:init_guess+thresh_t)>thresh_lo);
		end

		% attempt to fit the double exponential model, fit A, onset time, and tau

	end

	pos_crossing=pos_crossing(schmitt_flag==1);

	i
	for j=1:length(pos_crossing)

		spk_t=pos_crossing(j)./fs;	
		time_vec=[1:length(curr_roi)]./fs;

		switch lower(method(1))

			case 'f'
				[x]=fminsearch(@(x) obj_function_onset(x,curr_roi,fs,t_1),...
					[ onset_init_guess spk_t ],options);
			case 's'
				x=simulannealbnd(@(x) obj_function_onset(x,curr_roi,fs,t_1),... 
					[ onset_init_guess spk_t ] ,...
					[ onset_lbound spk_t-spk_delta ],...
					[ onset_hbound spk_t+spk_delta ],...
					options);
			otherwise
				error('Did not understand optimization method');
		
		end

		A=x(1);
		t_on=x(2);
		t_0=x(3);

		PEAKS{i}(end+1)=t_0;
		
		if onset_only
			continue;
		end

		switch lower(method(1))

			case 'f'
				[x]=fminsearch(@(x) obj_function_full(x,curr_roi,fs,t_0,t_on),...
					[ full_init_guess spk_t ],options);
			case 's'
				x=simulannealbnd(@(x) obj_function_full(x,curr_roi,fs,t_0,t_on),... 
					[ full_init_guess ] ,...
					[ full_lbound ],...
					[ full_hbound ],...
					options);
			otherwise
				error('Did not understand optimization method');
		
		end
		
		A_1=x(1);
		A_2=x(2);
		t_1=x(3);
		t_2=x(4);
	
		y1=calcium_model_onset(A,t_on,t_0,t_1,time_vec);
		y2=calcium_model_full(t_0,t_on,A_1,A_2,t_1,t_2,time_vec);


		if debug

			time_vec=1:length(curr_roi);
			figure(1);plot(time_vec,curr_roi);

			hold on;	
			plot(time_vec,y1,'r-');
			plot(time_vec,y2,'g-');

			title(['ROI:  ' num2str(i)]);

			pause();
			cla;
		end

	end

end

end

function res = obj_function_onset(x,dff,fs,t_1)
%
%
%
%

A=x(1);
t_on=x(2);
t_0=x(3);

time_vec=[1:length(dff)]./fs;

y=calcium_model_onset(A,t_on,t_0,t_1,time_vec);
res=sum((y(:)-dff(:)).^2);

end

function res = obj_function_full(x,dff,fs,t_0,t_on)

A_1=x(1);
A_2=x(2);
t_1=x(3);
t_2=x(4);

time_vec=[1:length(dff)]./fs;

y=calcium_model_full(t_0,t_on,A_1,A_2,t_1,t_2,time_vec);
res=sum((y(:)-dff(:)).^2);

end

function y = calcium_model_onset(A,t_on,t_0,t_1,t)
%
%
%
%

y = A.*(1-exp(-(t-t_0)./t_on)).*exp(-(t-t_0)./t_1);
y(t<t_0)=0;

end

function y = calcium_model_full(t_0,t_on,A_1,A_2,t_1,t_2,t)

y=(1-exp(-(t-t_0)./t_on)).*(A_1.*exp(-(t-t_0)./t_1)+A_2.*exp(-(t-t_0)./t_2));
y(t<t_0)=0;

end
