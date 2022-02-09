function [ BarCorrected ] = DiscardJumpsBarcorrected(BarCorrected, pixelsperiter)

[M,N]=size(BarCorrected);

indextodiscard=[]; 

if M>1

posinitial=BarCorrected(1,1:2); 
tinitial=BarCorrected(1,3);

posfinal=BarCorrected(2,1:2); 
tfinal=BarCorrected(2,3);

if norm(posinitial-posfinal)>(tfinal-tinitial)*pixelsperiter
    indextodiscard=[indextodiscard, 1]; 
end; 

if M>=3
    
    for n=2:M-1
        
        pos=BarCorrected(n,1:2);
        t0=BarCorrected(n,3);
        
        pospre=BarCorrected(n-1,1:2); 
        tpre=BarCorrected(n-1,3);
        
        

        pospost=BarCorrected(n+1,1:2); 
        tpost=BarCorrected(n+1,3);
        
        if  (norm(pos-pospre)>(t0-tpre)*pixelsperiter)||(norm(pos-pospost)>(tpost-t0)*pixelsperiter)
         indextodiscard=[indextodiscard, n]; 
        end;
        
        
    end; 
    
    
posinitial=BarCorrected(M-1,1:2); 
tinitial=BarCorrected(M-1,3);

posfinal=BarCorrected(M,1:2); 
tfinal=BarCorrected(M,3); 
    
    
if norm(posinitial-posfinal)>(tfinal-tinitial)*pixelsperiter
    indextodiscard=[indextodiscard, M]; 
end; 





end; 





BarCorrected(indextodiscard,:)=[]; 


end

