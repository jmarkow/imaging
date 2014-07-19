%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% aggregate spatial data

[~,result]=system('find . -type f -name spatial_corr.mat');

files=regexp(result,'\n','split');
files(end)=[];

load(files{1},'dist_y','dist_bins');

nbins=length(dist_y);

for i=1:nbins
	agg_spatial{i}=[];
end

for i=1:length(files)

	load(files{i},'dist_y');

	for j=1:length(dist_y)
		agg_spatial{j}=[agg_spatial{j} dist_y{j}];
	end


end
