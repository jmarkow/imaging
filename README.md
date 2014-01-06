# Welcome

This package is designed to align movie files to corresponding experimental data using the sync signal.

===

# Quick Start


If you want to get running right away.  First, ensure that the movie files are all in one directory alone with the file that contains your aligned sync and experimental data.  Now, load the sync and experimental data. 

First, you need to *resplit* the movie files with the aligned microphone data.    

```
>> cd ~/newdata; % cd into the dir with movie files
>> mic_data=wavDat{1}.dat; # load in the microphone data
>> sync_data=wavDat{2}.dat; # load in the sync data
>> wav_fs=24.414e3; % sampling frequency of collected data
>> fb_data_split(mic_data,sync_data,'fs',wav_fs)
```
The data should be split, you may see some warnings (ignore them for now).  If you get a lot of "frame number not equal to number of rising edges" errors, then you may need to change the *gap_thresh* parameter.   Type `help fb_data_split` for more info on changing this. 

Next, we need to perform our template alignments.  This involves computing spectral features for a user-selected template, then finding matches to the template using a running Euclidean distance score.  First, compute the spectral features of *split* MATLAB files.

```
cd mat/; % cd into the newly created mat dir
matlabpool open 4; % open 4 parfor workers
fb_compute_features; % compute the features for all files
matlbapool close all; % close the parfor workers
```

Now the template alignment. The following command will start the template matching process and extract an additional 150 msecs before and after the match. 

```
fb_template_match([],'padding',[.15 .15])
```

If you haven't performed an extraction before, this will call a GUI to first select a file to get the template from (check the gif directory to find a good example).  Then draw a box around the template, double click, and confirm in command window (type d for done).  You'll see a lot of status updates.  Then, you will need to perform a manual cluster cut of the matches.  Draw a polygon around your cluster, type enter to confirm, then be sure to save the appropriate cluster (by default this is cluster 2).  Next, to generate pretty movies,

```
cd extraction/
fb_mov_proc();
```

Make sure your screensaver is off before running this command, since getframe (UGH) managed to draw anything on the screen in the vicinity if the figure window.  This command will generate a *raw* movie, *df/f* movie and a maximum projection of *df/f* in the folder *proc*.  To proceed with roi analysis,

```
cd mov/
rois=fb_select_roi();
fb_plot_roi(rois);
```

To select ROIs, you select a file in extraction/mov to get the ROI from. Then, you need to adjust the contrast of the image to see the ROI appropriately.  Close the contrast adjustment tool when finished, then drag an ellipsis around each ROI and double click to confirm.  When you're finished, close the main figure window.  Your ROIs are saved to the environment variables rois.  Feed this to fb_plot_roi, and you'll see some plots of the ROI time-course aligned to song in the folder *roi*.


===
