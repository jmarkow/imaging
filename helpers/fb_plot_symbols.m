function fig=fb_plot_symbols(XDATA,YDATA,varargin)
%
%
%
%



to_del=find(cellfun(@length,YDATA)<2);

XDATA(to_del)=[];
YDATA(to_del)=[];

mu=cellfun(@mean,YDATA)
sem=cellfun(@std,YDATA)./sqrt(cellfun(@length,YDATA))

to_del=find(isnan(mu));

mu(to_del)=[];
sem(to_del)=[];
XDATA(to_del)=[];

spacefun=@(M,lambda,C,x) M*exp(-x/(lambda))+C;
g=fittype(@(M,lambda,C,x) M*exp(-x/(lambda))+C);

f1=fit(XDATA(:),mu(:),g);
coeffs=coeffvalues(f1)
%ci=predint(f1,uniq);
%funx=1:10:1e3;
funx=linspace(XDATA(1),XDATA(end),100);
funeval=spacefun(coeffs(1),coeffs(2),coeffs(3),funx);

w=10;
color=[1 0 0];
w=.25;
figure();
for i=1:length(mu)
	line([XDATA(i) XDATA(i)],[mu(i)-sem(i) mu(i)+sem(i)],'color',color);
	hold on;
	line([XDATA(i)-w XDATA(i)+w],[mu(i)-sem(i) mu(i)-sem(i)],'color',color);
	line([XDATA(i)-w XDATA(i)+w],[mu(i)+sem(i) mu(i)+sem(i)],'color',color);
end

plot(XDATA,mu,'r^','color',[1 0 0],'markersize',9,'markerfacecolor',[1 1 1],'markeredgecolor',[1 0 0]);
plot(funx,funeval,'r-','linewidth',1.5);
