%% Test script for MS2 tracking. 


clear; clc;
close all;

sigma = 2;          % blurring for FINE nuclei identification and tracking [px], 

minSize = 1000;   % minimum size for nuclei identification [px^2]
maxSize=18000; %This is also useful because sometimes two cells are attached (SZ)

nclusters = 3;  % Clusters of the image used for the rought segmentation. The threshold is the mean of the centroids of the intensities from the lower to the ncluster-1 max of the value
treshFactor_FINE = 1;  % threshold factor for nuclei identification
TolArea=0.25; 

parCellTrack = [sigma, minSize, nclusters, treshFactor_FINE, maxSize,TolArea]; % Parameters for cell tracking



% Set parameters for MS-2 spot identification and tracking

hold off
close all
loBP = 1;   % lower limit for bandpass filter [px]
hiBP = 5;   % higher limit for band pass filtering [px]
Ns = 1.25; % Ntimes value for the spot quantification
windowSz = 7; % Size of window for centroid quantification
ExcludeThresh = 15; % Exclude threshold distance for tracking of MS2 spot. 6 worked well but might be too stringent.
sizemaxdot=25; 
parTrackSpot = [loBP,hiBP, Ns, windowSz,ExcludeThresh,sizemaxdot];  % Parameters for MS-2 spot tracking.


pathInput=strcat(pwd,'\');
pathOutput=strcat(pwd,'\');
pathCalculations=strcat(pwd,'\');

fileName='Short_TestMS2.tif'

Results=functionMS2SpotsTracking_v8( parCellTrack, parTrackSpot, fileName, pathInput, pathOutput,pathCalculations,true);


 
load('QUANTv8_Short_TestMS2.mat') 


fileName=Results.nameFileInput;

pathName= Results.nameDirInput;

%%
if fileName ~= 0
    [Stack, nFrames] = TIFread([pathName, fileName]);
end

OUTfile=Results.OUT;

Intensities=Results.SPOT_INTENSITY;

[nFrames,nCells]=size(Intensities);

figure(2)
plot_MS2Results_v8_PEAK(Results, Stack, nCells, nFrames);