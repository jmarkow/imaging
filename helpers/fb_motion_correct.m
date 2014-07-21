function CORRECT=fb_motion_correct(MOV_DATA,TEMPLATE)
%motion corrects the data in mov_data to match a template image using xcorr
%
%
%


if ~motion_flag

	% correction template, select coordinates
	% mean works for now, could get more complicated in the future

	corr_roi=mean(mov_data,3);

	% for now use entire image, could consider cropping

	y_segment=1:rows;
	x_segment=1:columns;

	correction_fft=fft2(corr_roi);

	motion_flag=1;

end

[nblanks formatstring]=fb_progressbar(100);
fprintf(1,['Progress:  ' blanks(nblanks)]);

for j=1:frames

	fprintf(1,formatstring,round((j/frames)*100));	

	% uncomment to test for motion correction, introduces random x,y shifts

	%tmp=circshift(mov_data(y_segment,x_segment,j),[randi(100) randi(100)]);

	tmp=mov_data(y_segment,x_segment,j);

	% last argument is upsample factor

	[output Greg]=dftregistration(correction_fft,fft2(tmp),100);
	mov_data(:,:,j)=abs(ifft2(Greg)); % recover corrected image

end


