function out = findParticlesv5(StackFiltered, Stack, LabelsMap, Ntimes, hiBP,windowSz,szmax)

% findparticles uses the functions pkfind and cntr to identify the
% centroids of peaks in a stack of bandpass filtered images. The output is 
% Nx6 array, where N is the number of particles identified,and:
%
% out(:,1) = x-coordinate
% out(:,2) = y-coordinate
% out(:,3) = brightness of particle
% out(:,4) = giration radius
% out(:,5) = peak intensity of particle
% out(:,6) = frame at which particle has been identified

nFrames = size(StackFiltered,2);

out = [];



%SZ, it uses a threshold that is Ntimes the 99% pctile of the filtered image. 
% Three times its value (Ntimes=3) seems to work. 

%v5: One uses the filtered stack for the detection of peaks, but the true
%stack for the quantification of the peak intensity. We considered the
%possibility of using a mask but this detects "nucleoli"



%[Mim,Nim]=size(StackFiltered(imageIx).data);

for imageIx = 1:nFrames  

    
    matrixfiltered=StackFiltered(imageIx).data;
    %disp('Threshold for spots n') 
    
    %figure(800)
    %pause(0.5)
    %imagesc(matrixfiltered.*(LabelsMap(imageIx).data>0))
    
    peaks = pkfnd_mod_LabelMaps(matrixfiltered, LabelsMap(imageIx).data, Ntimes, szmax, hiBP);
    if ~isempty(peaks)
    centroids = cntrd(StackFiltered(imageIx).data,peaks,windowSz);
    else
    disp('Warning: no spots in frame...')
    disp(num2str(imageIx))
    %centroids = [0 0 0 0];
    centroids = [];
    isempty(centroids)
    end;
    
    if ~isempty(centroids)
        for i = 1: length(centroids(:,1))
        
            x = round(centroids(i,1));
            y = round(centroids(i,2));
            
            
%         ymin=max(1, y_center-2);
%         ymax=min(Nim, y_center+2);
%         xmin=max(1, x_center-2);
%         xmax=min(Mim, x_center+2);
%   
%         Mask=zeros(Mim,Nim); 
%         Mask(ymin: ymax, xmin: xmax) = 1;
%       
             
            
            
            
            
            if x ~= 0 && y ~= 0;
                %centroids(i,5) =  sum(sum(double(Stack(imageIx).data).*Mask));
                 centroids(i,5) = Stack(imageIx).data(y,x);
            else
                centroids(i,5) = 0;            
            end
        end
    centroids(:,6) = imageIx;
    end
    out = [out; centroids];
   
     
            
end

out(1,:);
    
