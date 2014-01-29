function [MOV_DATA,FRAME_IDX]=fb_retrieve_mov(FILE,varargin)
%
%
%
%
%
%

nparams=length(varargin);

im_resize=1;

% parameter collection

if mod(nparams,2)>0
	error('Parameters must be specified as parameter/value pairs');
end

for i=1:2:nparams
	switch lower(varargin{i})
		case 'im_resize'
			im_resize=varargin{i+1};
	end
end

if nargin<1 | isempty(FILE)
	[filename,pathname]=uigetfile({'*.mat'},'Pick a .mat file to retrieve movie data for',fullfile(pwd));
end

load(fullfile(pathname,filename),'rising_data','fs');

[path,file,ext]=fileparts(filename);

% movie filename one directory down

mov_filename=fullfile('..',[file '.tif']);
image_info=imfinfo(mov_filename);

FRAME_IDX=find(rising_data>0);
frame_val=rising_data(FRAME_IDX);

first_frame=frame_val(1);
last_frame=frame_val(end);

mov_idx=first_frame:last_frame;

width=image_info.Width;
height=image_info.Height;

if ~isempty(im_resize)
	width=width*im_resize;
	height=height*im_resize;
end

MOV_DATA=zeros(height,width,length(mov_idx));

disp('Retrieving movie data...');

counter=1;
for i=first_frame:last_frame
	imdata=imread(mov_filename,i);

	if ~isempty(im_resize)
		imdata=imresize(imdata,im_resize);
	end

	MOV_DATA(:,:,counter)=imdata;
	counter=counter+1;
end

% FRAME_IDX is used to determine the closest point in time for the onset
% of a given frame
