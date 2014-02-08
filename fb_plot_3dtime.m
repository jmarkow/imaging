function fb_plot_3dtime(DATA,T,varargin)
%
%
%
%
%
%
%

nparams=length(varargin);

linewidth=1.2;
colors='winter';
spacing=10;
fs=24.414e3;
dims=1:3;

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'colors'
			colors=varargin{i+1};
		case 'linewidth'
			linewidth=varargin{i+1};
		case 'dims'
			dims=varargin{i+1};
	end
end

plotdata=DATA(:,dims);

x=plotdata(:,1);
y=plotdata(:,2);
z=plotdata(:,3);
cvec=T;

x=[x';x'];
y=[y';y'];
z=[z';z'];
col=[cvec;cvec];

surface(x,y,z,col,'facecol','no','edgecol','interp','linew',1.5);

%surface(x,y,z,col,'facecol','no','edgecol','interp','linew',2);
%h=waterfall(x,y,col);
%set(h,'linewidth',2,'facecolor','interp');
%CD = get (h, 'CData');
%CD(1,:) = nan;
%CD(end-2:end,:) = nan; % remove the "falling boundary"
%set (h, 'CData', CD)


grid off;

%axis([-.12 .12 -.12 .12]);
%colormap(flipud(hot(len)));



