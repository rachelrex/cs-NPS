
function [events]=reanalysis_find_events(pulsesog,threshptsc)

szp=length(pulsesog(:,1));
szt=length(threshptsc(:,1));
events=zeros(szp,6);
found=0;
k=1;
for i=1:szp
    
       found=0;
 
        %find the threshpt (that is going down=0) taht is nearest to the
        %pulsesog(i,col) 
        while k<=szt && threshptsc(k,1)~=0 && found ==0
            if pulsesog(i,1)>threshptsc(k,1)
                k=k+1;
            else
                %is it closer to the threshpt at k, or k-1? 
                if threshptsc(k,2)==1
                    k=k-1;
                end
                if k<=2
                    events(i,1)=threshptsc(k,1);
                    events(i,2)=threshptsc(k+1,1);
                    k=k+1;
                    found=1;
                elseif abs(threshptsc(k,1)-pulsesog(i,1))<abs(threshptsc(k-2,1)-pulsesog(i,1))
                    events(i,1)=threshptsc(k,1);
                    events(i,2)=threshptsc(k+1,1);
                    k=k+1;
                    found=1;
                else
                    events(i,1)=threshptsc(k-2,1);
                    events(i,2)=threshptsc(k-1,1);
                    k=k-1;
                    found=1;
                end
            end
        end
        %MAKE IT MOVE ONTO THE NEXT EVENTS LINE 

  
end
