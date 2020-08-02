function [OUT] = trackSpot_MS2_v8(IN,SpotData,Stack,nCells,ExcludeThresh,LabelsMap)


OUT = IN;

% Identify how many spots are in each cell


[M,N]=size(Stack(1).data); 

% figure(3000)
% imagesc(Stack(1).data)
% colormap('gray');
% hold on;
% plot(SpotData(:,1),SpotData(:,2),'go');
% hold on;
% title('First detected spots')

for iCell = 1:nCells;
    iSpot = find(SpotData(:,7) == iCell);
    
    OUT{iCell}.nDetectedSpots = length(iSpot);
    OUT{iCell}.DetectedSpots = SpotData(iSpot,:);
    
    
    if OUT{iCell}.nDetectedSpots == 0 % If no spot is found in iCell
        
        OUT{iCell}.Trajectory = OUT{iCell}.Baricenter;
        OUT{iCell}.Trajectory(:,3) = [1:1:OUT{iCell}.maxFrame];
        OUT{iCell}.Trajectory(:,4) = 0;
        OUT{iCell}.frameswithspot=[];
        % Assign the baricenter of the cell as the trajectory of the MS-2
        % spot.
        
    else % if there are spots in the cells
        
        % Correct track for Cell Baricenter movement.
        BarCorrected = zeros(OUT{iCell}.nDetectedSpots,3);
        
        
        Intensities=[];
        
        for ispot = 1:OUT{iCell}.nDetectedSpots % Cycle on the spots found the iCell
            iFrame = OUT{iCell}.DetectedSpots(ispot,6);
            BarCorrected(ispot,1) = OUT{iCell}.DetectedSpots(ispot,1)- ...
                OUT{iCell}.Baricenter(iFrame,1);
            BarCorrected(ispot,2) = OUT{iCell}.DetectedSpots(ispot,2)- ...
                OUT{iCell}.Baricenter(iFrame,2);
            BarCorrected(ispot,3) = iFrame;
            Intensities(ispot)=OUT{iCell}.DetectedSpots(ispot,5); 
                        
        end
        
        
        % New in this version: once we find the barcorrected, we just keep
        % the barcorrected spot that are in the BRIGHTEST group -with higher maximum-  that falls
        % within ExcludeThresh: Change in v5_ALT
        
        
        BarCorrected=FindGroupSpotsWithBrighterSpot_v2(BarCorrected,ExcludeThresh,Intensities); %Change in v5_ALT
        
%figure(3000)
%[Mb,Nb]=size(BarCorrected);

% for nb=1:Mb
% plot(BarCorrected(nb,1)+OUT{iCell}.Baricenter(BarCorrected(nb,3),1), BarCorrected(nb,2)+OUT{iCell}.Baricenter(BarCorrected(nb,3),2),'r+')
% end;
%pause(0.25)
        
        if OUT{iCell}.nDetectedSpots > 1
            CoM = [median(BarCorrected(:,1)),median(BarCorrected(:,2))];
        else
            CoM = BarCorrected(:,1:2);
        end
        

        BarCorrected(:,4)=1;
                
        if ~isempty(BarCorrected)
          [order,indexesordered]=sort(BarCorrected(:,3));
        end; 

        disp('This is the baricenter corrected \n');

        BarCorrected=BarCorrected(indexesordered,:);
        
        
        
        
        
        disp('This is by getting rid of too jumpy pixels')
        
        BarCorrected  = DiscardJumpsBarcorrected( BarCorrected,10);
        
   
        
        [Nelements,Ncolumns]=size(BarCorrected);
        
        matCoM=[CoM(1)*ones(Nelements,1),CoM(2)*ones(Nelements,1)];
        
        Distance = BarCorrected(:,1:2) - matCoM;
        Distance = sqrt(Distance(:,1).^2 + Distance(:,2).^2);
        
        
        
        
        
        
        
        
        % Find non-unique spots in each frames and keep the closest to the
        % center of mass of detected spots
        OUT{iCell}.Trajectory = [];
        
        
        %It might happen that all spots were eliminated, so we need this
        %if.
        if ~isempty(unique(BarCorrected(:,3)'))
           
            
            for iFrame = unique(BarCorrected(:,3)');
                
                idx = find(BarCorrected(:,3) == iFrame);
                if length(idx) == 1;
                    OUT{iCell}.Trajectory = [OUT{iCell}.Trajectory; BarCorrected(idx,:)];
                elseif length(idx) > 1;
                    [~, iMin] = min(Distance(idx));
                    OUT{iCell}.Trajectory = [OUT{iCell}.Trajectory; BarCorrected(idx(iMin),:)];
                    tempExcluded = BarCorrected(idx,:);
                    tempExcluded(iMin, :) = [];
                    tempExcluded(:,4) = -2;

                end
                
            end
            
            
            OUT{iCell}.frameswithspot=unique(BarCorrected(:,3));
            
            
            % if the first spot in cell is  found after  frame 1.
            iFirst = OUT{iCell}.Trajectory(1,3);
            if iFirst > 1;
                tempTrajectory = repmat(OUT{iCell}.Trajectory(1,1:2),iFirst - 1, 1);
                tempTrajectory(:,3) = 1: iFirst - 1;
                tempTrajectory(:,4) = 2;
                OUT{iCell}.Trajectory = [tempTrajectory; OUT{iCell}.Trajectory];
            end;
            
            % if the last spot in cell is  found before the end of the tracked cell.
            iLast = OUT{iCell}.Trajectory(end,3);
            if iLast < OUT{iCell}.maxFrame
                tempTrajectory = repmat(OUT{iCell}.Trajectory(end,1:2),...
                    OUT{iCell}.maxFrame - iLast, 1);
                
                tempTrajectory(:,3) = iLast + 1:OUT{iCell}.maxFrame;
                tempTrajectory(:,4) = 2;
                
                OUT{iCell}.Trajectory = [OUT{iCell}.Trajectory; tempTrajectory];
            end
            
            
            % Find gaps in trajectory and fill them by interpolation.
            
            
            track = OUT{iCell}.Trajectory;
            
            %find and fill gaps in tracks
            idx_gaps = find(track(2:end,3)- track(1:end-1,3) > 1); % find gaps
            if ~isempty(idx_gaps)
                
                for j = 1:length(idx_gaps);
                    k = idx_gaps(j);
                    N_steps = track(k+1,3)-track(k,3) + 1;
                    trackTemp = [];
                    trackTemp(:,1) = linspace(track(k,1),track(k+1,1), N_steps);
                    trackTemp(:,2) = linspace(track(k,2),track(k+1,2), N_steps);
                    trackTemp(:,3) = [track(k,3):track(k+1,3)]';
                    trackTemp(:,4) = 3;
                    trackTemp = trackTemp(2:end-1,:);
                    track = [track;trackTemp];
                end
                track = sortrows(track,3);
                OUT{iCell}.Trajectory = track;
            end
            % Add baricenter position back to trajectory
            OUT{iCell}.Trajectory(:,1) =  OUT{iCell}.Trajectory(:,1) + ...
                OUT{iCell}.Baricenter(:,1);
            
            OUT{iCell}.Trajectory(:,2) =  OUT{iCell}.Trajectory(:,2) + ...
                OUT{iCell}.Baricenter(:,2);
            
            %In the unlikely case in which the Trajectory goes out of the
            %limits (this might happen)
            
            iwrong1=[];  
            iwrong1=find(OUT{iCell}.Trajectory(:,1)<1); 
            OUT{iCell}.Trajectory(iwrong1,1)=1;   
           
            
            iwrong2=[];
            iwrong2=find(OUT{iCell}.Trajectory(:,1)>M); 
            OUT{iCell}.Trajectory(iwrong2,1)=M; 
            iwrong3=[];
            iwrong3=find(OUT{iCell}.Trajectory(:,2)<1); 
            OUT{iCell}.Trajectory(iwrong3,2)=1;            
            iwrong4=[];
            iwrong4=find(OUT{iCell}.Trajectory(:,2)>N); 
            OUT{iCell}.Trajectory(iwrong4,2)=N; 
            
            
            
            for m=1:OUT{iCell}.maxFrame
                
                maskcell=(LabelsMap(m).data==iCell);
                
                    
                if (maskcell(round(OUT{iCell}.Trajectory(m,2)),round(OUT{iCell}.Trajectory(m,1)))==0)
                 
                    %disp('Problem outside');
                    %iCell
                    %figure(800)
                    %imagesc(maskcell)
                   
                    %hold on;
                    %plot(round(OUT{iCell}.Trajectory(m,1)),round(OUT{iCell}.Trajectory(m,2)),'wo')
                   
                   
                    [I,J]=find(maskcell);
                    [indexmin,~] = knnsearch([I,J],[round(OUT{iCell}.Trajectory(m,2)),round(OUT{iCell}.Trajectory(m,1))]);
                    %plot(J(indexmin), I(indexmin), 'ko')
                    %hold on;
                    %axis([min(J)-10, max(J)+10, min(I)-10, max(I)+10])     
                    %figure(801)
                    %imagesc(Stack(m).data)
                    %hold on;
                    %plot(round(OUT{iCell}.Trajectory(m,1)),round(OUT{iCell}.Trajectory(m,2)),'ro')
                    %plot(J(indexmin), I(indexmin), 'go')
                    %hold on;
                    %axis([min(J)-10, max(J)+10, min(I)-10, max(I)+10])
                    %pause(0.5)
                    

                    OUT{iCell}.Trajectory(m,2)=I(indexmin(1)); 
                    OUT{iCell}.Trajectory(m,1)=J(indexmin(1)); 
                    
                    
                   % hold off;
                    
                    
                    
                    
                end;
                
            end; 
            
            
            
            
            
            
            
            
        else
            
            OUT{iCell}.Trajectory = OUT{iCell}.Baricenter;
            OUT{iCell}.Trajectory(:,3) = [1:1:OUT{iCell}.maxFrame];
            OUT{iCell}.Trajectory(:,4) = 0;
            
            
            
            OUT{iCell}.frameswithspot=[];
            
        end;
        
    end
    
    
 
% 
% figure(3000)  
% plot(OUT{iCell}.Trajectory(:,1), OUT{iCell}.Trajectory(:,2),'r+')
% title(num2str(iCell))
% pause(0.25);
% hold on;
%     
    
    
end


vectorInitialIntensities=[];

for iCell=1:nCells
    
    vectorInitialIntensities=[vectorInitialIntensities,OUT{iCell}.  TotalIntensityTrack(1)/OUT{iCell}.Area(1)];
    
end;

meanInitialIntensity=mean(vectorInitialIntensities); 


% Compute intensities along the trajectory.
for iCell = 1:nCells;

   
    
    for iFrame = 1:OUT{iCell}.maxFrame
        
        
        
        maskcell=(LabelsMap(iFrame).data==iCell);
        
        Mask = zeros(size(Stack(1).data));
        x_center = round(OUT{iCell}.Trajectory(iFrame,1));
        y_center = round(OUT{iCell}.Trajectory(iFrame,2));
        ymin=max(1, y_center-2);
        ymax=min(N, y_center+2);
        xmin=max(1, x_center-2);
        xmax=min(M, x_center+2);
  
        Mask(ymin: ymax, xmin: xmax) = 1;
        areadot=length(find(Mask>=1));
        if sum(size(Mask) > size(Stack(1).data)) > 0
            Mask = zeros(size(Stack(1).data));       
        end
        
        
        OUT{iCell}.Trajectory(iFrame,5) = sum(sum(double(Stack(iFrame).data).*Mask));
             
      
        ymin=max(1, y_center-4);
        ymax=min(N, y_center+4);
        xmin=max(1, x_center-4);
        xmax=min(M, x_center+4);     
        
        Mask(ymin: ymax, xmin: xmax) = 1;
        
        ymin=max(1, y_center-2);
        ymax=min(N, y_center+2);
        xmin=max(1, x_center-2);
        xmax=min(M, x_center+2);
  
        Mask(ymin: ymax, xmin: xmax) = 0;
        
        Mask=Mask.*maskcell;
        
%         figure(800)
%         imagesc(maskcell+Mask)
%         xlabel(num2str(iFrame));
%         ylabel(num2str(iCell));
%         pause(0.5)
        
   
        if sum(size(Mask) > size(Stack(1).data)) > 0
            Mask = zeros(size(Stack(1).data));
        end        
        %OUT{iCell}.Trajectory(iFrame,6) =
        %sum(sum(double(Stack(iFrame).data).*Mask));
        MatrixMaskBG=double(Stack(iFrame).data).*Mask;
        elementsBG=MatrixMaskBG(MatrixMaskBG>0);
        
        %This is the BG AROUND THE SPOT
        
        %OUT{iCell}.Trajectory(iFrame,6) = median(elementsBG);
        
        %Having selected the points inside the cell, we calculate the mean

        
        OUT{iCell}.Trajectory(iFrame,6) = mean(elementsBG);
        
         
        OUT{iCell}.Trajectory(iFrame,7) = std(elementsBG);
        
        
        OUT{iCell}.Trajectory(iFrame,8) =  Stack(iFrame).data(y_center,x_center);
        

     
    end
    
    
end

