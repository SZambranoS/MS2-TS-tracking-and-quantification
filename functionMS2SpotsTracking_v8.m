function Results = functionMS2SpotsTracking_v8_testnoise(parCellTrack, parTrackSpot, fileName, pathInput, pathOutput,pathCalculations,I_wanna_plot)%% 

%May 2017: uses clustering 

parCellSeg=parCellTrack(1:5);
TolArea=parCellTrack(6);



% Parameters for MS-2 spot identification and tracking

loBP=parTrackSpot(1);
hiBP=parTrackSpot(2);
Ns=parTrackSpot(3); 
windowSz=parTrackSpot(4);
ExcludeThresh=parTrackSpot(5); 
sizemaxdot=parTrackSpot(6);



currentdir=pwd; 

if fileName ~= 0
    [Stack, nFrames] = TIFread([pathInput, fileName]);
end

%% Find nuclei in each frame
nFrames = length(Stack);


cd(pathCalculations);



for i = 1:nFrames;
    cellROIs_refined(i) = segment_fromImage_refined_cluster_NOWATERSHED_NOBORDER(Stack(i).data, parCellSeg);
    disp(i)
end
clear i



%% Track Cells based on Nearest Neighbor
[LabelsMap, OUT, nCells,TotalintensityTrack,TotalintensityQuant,ObjectsPerFrame] = track_and_quantify_2channels(cellROIs_refined, Stack, Stack, nFrames, TolArea);



disp('Cells tracked'); 

%% Find MS-2 spots.


% bandpass filtering:

for i = 1:nFrames

filterStack(i).data = ...
        bpass(Stack(i).data, loBP, hiBP);    

end

% find peaks and compute peak centroid.

% Calculate Particles centroids
%SpotData = findParticlesv3(filterStack, Ns, hiBP,windowSz,sizemaxdot);

%SpotData = findParticlesv4(filterStack, LabelsMap, Ns, hiBP,windowSz,sizemaxdot);


SpotData = findParticlesv5(filterStack, Stack, LabelsMap, Ns, hiBP,windowSz,sizemaxdot);




for i = 1:nFrames

filterStack(i).data = ...
        bpass(Stack(i).data, loBP, hiBP);    
end






% finf
  

SpotDataFinal=SpotData; 


%% Assign the identified MS-2 spots to the different cells.
for i = 1:length(SpotData(:,1));
    SpotDataFinal(i,7) = LabelsMap(SpotData(i,6)).data(round(SpotData(i,2)), round(SpotData(i,1)));
end

idx = find(SpotDataFinal(:,7) ==0);
SpotDataFinal(idx,:) = [];



clear idx
clear i
 

%% Compute trace of the MS-2 Spot.


OUT = trackSpot_MS2_v8(OUT,SpotDataFinal,Stack,nCells, ExcludeThresh,LabelsMap);



 for iCell = 1:nCells  
   SPOT_INTENSITY(:,iCell) = [OUT{iCell}.Trajectory(:,5); zeros(nFrames-OUT{iCell}.maxFrame,1)];%Intensity of the spot
   SPOT_BGINTENSITY(:,iCell) = [OUT{iCell}.Trajectory(:,6); zeros(nFrames-OUT{iCell}.maxFrame,1)];%Mean of a crown around the spot (inside the cell)
   SPOT_STDBGINTENSITY(:,iCell) = [OUT{iCell}.Trajectory(:,7); zeros(nFrames-OUT{iCell}.maxFrame,1)];%Std of a crown around the spot
   SPOT_MAXINTENSITY(:,iCell) = [OUT{iCell}.Trajectory(:,8); zeros(nFrames-OUT{iCell}.maxFrame,1)];%Intensity of the spot
   
  
   
   
   NUCLEARINTENSITY(:,iCell)=[OUT{iCell}.TotalIntensityTrack;zeros(nFrames-OUT{iCell}.maxFrame,1)];
   NUCLEARAREAS(:,iCell)=[OUT{iCell}.Area;zeros(nFrames-OUT{iCell}.maxFrame,1)];

 end
clear iCell



 
%% Plot output cell
    

outputname=strcat('QUANTv8_',fileName,'.mat');
outputname=strrep(outputname,'.tif','');

Results.OUT=OUT; 
Results.SPOT_INTENSITY=SPOT_INTENSITY;
Results.SPOT_BGINTENSITY=SPOT_BGINTENSITY;
Results.SPOT_STDBGINTENSITY=SPOT_STDBGINTENSITY;
Results.SPOT_MAXINTENSITY=SPOT_MAXINTENSITY;


Results.NUCLEARINTENSITY=NUCLEARINTENSITY;
Results.NUCLEARAREAS=NUCLEARAREAS;


Results.parCellSeg=parCellSeg;
Results.parTrackSpot=parTrackSpot;

Results.nameDirInput=pathInput;

Results.nameFileInput=fileName;


cd(pathOutput);
save(outputname,'Results');
cd(currentdir);


% if I_wanna_plot==true
% plot_MS2Results_v4(OUT, Stack, nCells, nFrames);
% end; 


%  
% if I_wanna_plot==true
% plot_MS2Results_v8_PEAK(Results, Stack, nCells, nFrames);
% end; 





    
    
    
    







