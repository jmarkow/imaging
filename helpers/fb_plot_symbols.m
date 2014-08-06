function [LAMBDA]=fb_plot_symbols(XDATA,YDATA,varargin)
%
%
%
%


% parameter collection

w=5;
symbol_color=[1 0 0];
linespec='r^';
error_color=[1 0 0];
xlabels='';
markersize=10;
linewidth=1.25;
nparams=length(varargin);

markerfacecolor=[1 1 1];
markeredgecolor=[1 0 0];

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'w'
			w=varargin{i+1};
		case 'symbol_color'
			symbol_color=varargin{i+1};
		case 'linespec'
			linespec=varargin{i+1};
		case 'error_color'
			error_color=varargin{i+1};
		case 'xlabels'
			xlabels=varargin{i+1};
		case 'markersize'
			markersize=varargin{i+1};
		case 'markerfacecolor'
			markerfacecolor=varargin{i+1};
		case 'markeredgecolor'
			markeredgecolor=varargin{i+1};
		case 'linewidth'
			linewidth=varargin{i+1};
	end
end


options.MaxIter=inf;
options.TolMesh=1e-9;
options.TolX=1e-9;
options.TolFun=1e-9;

mu=cellfun(@mean,YDATA);
st_dev=cellfun(@std,YDATA);
len=cellfun(@length,YDATA);

sem=st_dev./(sqrt(len));

%[x] = patternsearch(@(x) obj_function(x,mu,XDATA),[ 0 1 mu(1) ],[],[],[],[],[],[],[],options);

spacefun=@(M,lambda,C,x) M*(1-exp(-x/(lambda)))+C;
g=fittype(@(M,lambda,C,x) M*(1-exp(-x/(lambda)))+C);
f1=fit(XDATA(:),mu(:),g,'startpoint',[ -20 100 -5 ]);
x=coeffvalues(f1)

LAMBDA=x(2);

%funx=1:10:1e3;
funx=linspace(XDATA(1),XDATA(end),100);
funeval=spacefun(x(1),x(2),x(3),funx);

figure();
for i=1:length(mu)
	line([XDATA(i) XDATA(i)],[mu(i)-sem(i) mu(i)+sem(i)],'color',error_color,'linewidth',linewidth);
	hold on;
	line([XDATA(i)-w XDATA(i)+w],[mu(i)-sem(i) mu(i)-sem(i)],'color',error_color,'linewidth',linewidth);
	line([XDATA(i)-w XDATA(i)+w],[mu(i)+sem(i) mu(i)+sem(i)],'color',error_color,'linewidth',linewidth);
end

plot(XDATA,mu,linespec,'color',symbol_color,'markersize',markersize,...
	'markerfacecolor',markerfacecolor,'markeredgecolor',markeredgecolor,'linewidth',linewidth);
plot(funx,funeval,'r-','linewidth',1.5);

if ~isempty(xlabels)
	set(gca,'XTick',XDATA,'XTickLabel',xlabels);
end

end

function res=obj_function(x,ydata,t)

M=x(1);
lambda=x(2);
C=x(3);

funeval=M*(1-exp(-t/(lambda)))+C;
res=sum((ydata(:)-funeval(:)).^2);

res

end

%
%
%
