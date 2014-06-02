function ROIS=fb_convert_imagej_roi(FILENAME,varargin)
%
%
%
%
%

% wrapper for ReadImageJROI

% pull in data, if cell, then map to multiple ROIs

rois=ReadImageJROI(FILENAME);

if iscell(rois)
	for i=1:length(rois)
		% return in x,y pair format (all coordinates inside ROI)

		tmp=rois{i}.vnRectBounds; % top, left, bottom, right
		[y x]=meshgrid(tmp(1):tmp(3),tmp(2):tmp(4));
		ROIS{i}=[y(:) x(:)];
		
	end
else

	tmp=rois.vnRectBounds; % top, left, bottom, right
	[y x]=meshgrid(tmp(1):tmp(3),tmp(2):tmp(4));
	ROIS{1}=[y(:) x(:)];

end


