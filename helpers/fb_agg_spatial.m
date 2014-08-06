%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% aggregate spatial data

[~,result]=system('find . -type f -name spatial_corr.mat');

um_per_px=.6250;
bird_id={};

files=regexp(result,'\n','split');
files(end)=[];

load(files{1},'dist_y','dist_bins');

nbins=length(dist_y);

range=1:5;

agg_spatial={};
for i=1:nbins
	agg_spatial{i}=[];
end

for i=1:length(files)

	load(files{i},'dist_y');
	for j=1:length(dist_y)
		agg_spatial{j}=[agg_spatial{j} dist_y{j}];
	end


end

%%%% put together labels, bin center is used for measuring lambda
%

labels={};
bin_c=[];

dist_bins=dist_bins.*um_per_px;

for i=1:length(dist_bins)-1
	labels{i}=[ sprintf('%i-%i',dist_bins(i),dist_bins(i+1)) ];
	%labels{i}=[ num2str(dist_bins(i));num2str(dist_bins(i+1)) ];
	bin_c(i)=mean(dist_bins(i:i+1));
end

[lambda]=fb_plot_symbols([bin_c(range)],agg_spatial(range),'xlabels',labels(range),...
	'markerfacecolor',[1 1 1],'markeredgecolor',[0 0 0],'error_color',[0 0 0],'linespec','go','w',15,'linewidth',1);
%axis tight;
set(gca,'TickDir','out','TickLength',[.025 .025],'FontSize',9,'FontName','Helvetica');
set(gcf,'Position',[680 826 240 176],'paperpositionmode','auto');
xticklabel_rotate([],45);
title([num2str(lambda)]);
