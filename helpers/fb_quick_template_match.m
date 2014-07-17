function [MATCHES,SCORE]=fb_quick_template_match(MIC_DATA,varargin)
%
%
%
%
%

fs=24.414e3;
n=1024;
overlap=1e3;
down_factor=5;

nparams=length(varargin);

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})	
		case 'fs'
			fs=varargin{i+1};
	end
end


[b,a]=ellip(5,.2,80,[500]/(fs/2),'high');
MIC_DATA=filtfilt(b,a,double(MIC_DATA));

temp_mat=[];
TEMPLATE.data=fb_spectro_navigate(MIC_DATA);
TEMPLATE.features=fb_smscore(TEMPLATE.data,fs);
file_features=fb_smscore(MIC_DATA,fs);
file_length=size(file_features{1},2)
template_length=size(TEMPLATE.features{1},2)-1;

for j=1:length(file_features)
	
	score_temp{j}=[];

	for k=1:file_length-template_length
		score_temp{j}=[score_temp{j} sum(sum(abs( file_features{j}(:,k:k+template_length)-TEMPLATE.features{j} )))];
	end

	% keep the raw scores for further analysis

	raw_temp{j}=score_temp{j};

	score_temp{j}=score_temp{j}-mean(score_temp{j});
	score_temp{j}=score_temp{j}/std(score_temp{j});
	score_temp{j}(score_temp{j}>0)=0;
	score_temp{j}=abs(score_temp{j});

end

attributes=length(score_temp);
product_score=score_temp{1};

for j=2:attributes, product_score=product_score.*score_temp{j}; end

[pks,locs]=findpeaks(product_score,'MINPEAKHEIGHT',.5,'MINPEAKDISTANCE',template_length);

SCORE=score_temp;
MATCHES(:,1)=(locs*(n-overlap)*down_factor)-n;
MATCHES(:,2)=MATCHES(:,1)+length(TEMPLATE.data);
