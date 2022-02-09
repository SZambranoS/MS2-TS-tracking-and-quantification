function plot_MS2Results_v8_PEAK(Results, Stack, nCells, nFrames)

figure(1)

OUT=Results.OUT;

SPOTminusBG=Results.SPOT_MAXINTENSITY-Results.SPOT_BGINTENSITY; 

stdBG=3*Results.SPOT_STDBGINTENSITY;

maxaxis=max(max(SPOTminusBG));


hStack = axes('Position', [.05, .3, .45, .6]);

A=double(Stack(1).data);

minimum=mean(A(:))-std(A(:));

maximum=mean(A(:))+5*std(A(:)); 

limitstouse=[minimum, maximum]; 

imagesc(Stack(1).data,limitstouse);
colormap('jet')

hPlot = axes('Position', [.55, .3, .40, .6]);

plot(hPlot, [1:nFrames], SPOTminusBG(:,1),...
    [1,1],[0,maxaxis],[1:nFrames], stdBG(:,1),'r');


setappdata(hStack, 'iFrame', 1);
setappdata(hStack, 'iCell', 1);




hFrameSlider = uicontrol( ...
    'Style', 'slider', ...
    'Value', 1,...
    'Min', 1, 'Position',[30,64,500,23],...
    'Max',  nFrames, ...
    'SliderStep', [1/nFrames, 1/nFrames],...
    'Callback', @(src, evt) changeFrame(hStack, hPlot, get(src, 'Value'), ...
     OUT, Stack, nFrames,limitstouse,SPOTminusBG,maxaxis,stdBG));


hCellSlider = uicontrol( ...
    'Style', 'slider', ...
    'Value', 1,...
    'Min', 1,'Position',[30,24,500,23],...
    'Max',  nCells, ...
    'SliderStep', [1/nCells, 1/nCells],...
    'Callback', @(src, evt) changeCell(hStack, hPlot, get(src, 'Value'),...
     OUT, Stack, nFrames,limitstouse,SPOTminusBG,maxaxis,stdBG));









function changeFrame(hStack, hPlot, iFrame, ...
    OUT, Stack, nFrames ,limitstouse,SPOTminusBG,maxaxis,stdBG)



iFrame = round(iFrame);
iCell = getappdata(hStack, 'iCell');


imagesc(Stack(iFrame).data, 'Parent', hStack,limitstouse);


axis([max(OUT{iCell}.Trajectory(iFrame,1)-100,1), min(OUT{iCell}.Trajectory(iFrame,1)+100,1024), max(OUT{iCell}.Trajectory(iFrame,2)-100,1),min(OUT{iCell}.Trajectory(iFrame,2)+100,1024)]);



if iFrame <= OUT{iCell}.maxFrame
hold(hStack, 'on');
plot(hStack,OUT{iCell}.Trajectory(iFrame,1), OUT{iCell}.Trajectory(iFrame,2),'or',...
    'MarkerSize', 9);

title(num2str(iCell))
hold(hStack, 'off');
end



plot(hPlot, [1:nFrames], SPOTminusBG(:,iCell),...
    [iFrame,iFrame],[0,maxaxis],[1:nFrames], stdBG(:,iCell),'r');
title(num2str(iCell))
 
 


setappdata(hStack, 'iFrame', iFrame);




function changeCell(hStack, hPlot,  iCell, ...
    OUT, Stack, nFrames,limitstouse,SPOTminusBG,maxaxis,stdBG)

iFrame = getappdata(hStack, 'iFrame');
iCell = round(iCell);



imagesc(Stack(iFrame).data, 'Parent', hStack,limitstouse);


axis([max(OUT{iCell}.Trajectory(iFrame,1)-100,1), min(OUT{iCell}.Trajectory(iFrame,1)+100,1024), max(OUT{iCell}.Trajectory(iFrame,2)-100,1),min(OUT{iCell}.Trajectory(iFrame,2)+100,1024)]);

if iFrame <= OUT{iCell}.maxFrame
hold(hStack, 'on');
plot(hStack,OUT{iCell}.Trajectory(iFrame,1), OUT{iCell}.Trajectory(iFrame,2),'or',...
    'MarkerSize', 9);
hold(hStack, 'off');
end

plot(hPlot, [1:nFrames], SPOTminusBG(:,iCell),...
    [iFrame,iFrame],[0,maxaxis],[1:nFrames], stdBG(:,iCell),'r');
title(num2str(iCell))
setappdata(hStack, 'iCell', iCell);




