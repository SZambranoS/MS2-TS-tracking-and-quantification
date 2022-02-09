function outfinal=pkfnd_mod_LabelMaps(Filteredimage, LabelsMap, Ntimes, szmax,sz)


%This one applies a threshold which is mean + Ntimes SD of each mask that
%one can obtain from LabelsMap;

outfinal=[]; 

labels=unique(LabelsMap); 

labels=labels(find(labels));


for m=1:length(labels)

    Mask=(LabelsMap==labels(m));
    
    Masked_Filteredimage=Mask.*Filteredimage; 
    
    intensityvalues=Masked_Filteredimage(:);
    
    intensityvalues=intensityvalues(find(intensityvalues));
    
    th=Ntimes*prctile(intensityvalues,99); 
    
    %figure(m+10)
    %hist(intensityvalues,100,100);
    %prctile(intensityvalues,99)
    %Ntimes
    %title(num2str(th));
    %pause(0.5);
    close
    


    out=pkfnd_mod(Masked_Filteredimage,th,szmax,sz);
    
%     display('Cell:')
%     labels(m)
%     display('Number of spots')
%     length(out)
%     
%     
%     xlabel(num2str(length(out))); 
%     
    outfinal=[outfinal;out]; 


end; 










end

