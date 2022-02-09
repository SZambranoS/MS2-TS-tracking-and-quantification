function [ BarOutput, BarExcluded] = FindGroupSpotsWithBrighterSpot_v2(BarCorrected,ExcludeThresh,Intensities)

[idsclose,D]= rangesearch(BarCorrected(:,1:2),BarCorrected(:,1:2),ExcludeThresh);






numberelements=[];

for n=1:length(idsclose)    
    numberelements(n)=length(idsclose{n});    
   intensities(n)=max(Intensities(idsclose{n}));
    
    
    BarCluster(n,:)=[mean(BarCorrected(idsclose{n},1)),mean(BarCorrected(idsclose{n},2))];
    
    
    
end;


disp('Total number of elements in cluster \n')
sum(numberelements)
disp('Total number points \n')
length(BarCorrected(:,2))


% 
% if iCell==23
%     
%     ExcludeThresh
%     numberelements
%     intensities
% end; 

indexesaccepted=find(numberelements>=1); 

if length(indexesaccepted)>0
[valuemax, indexmaxvalue]=max(intensities(indexesaccepted));

indexmax=indexesaccepted(indexmaxvalue);
else

[valuemax, indexmax]=max(intensities);
    
end; 

%disp('Number of clusters with maximum equal to valuemax \n'); 
indexcandidates=(find(intensities==valuemax))
ncandidates=numberelements(indexcandidates)

[maxn, indexn]=max(ncandidates);


indexmax=indexcandidates(indexn)

%disp('Number of elements in the cluster \n'); 
%numberelements(intensities==valuemax)    





spotssaved=idsclose{indexmax};

BarOutput=BarCorrected(spotssaved,:);






end

