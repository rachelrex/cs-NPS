function editedthresh=findnearestthresh(new_threshptsc,deletethresh,whichaxes)
%find the threshpts that are nearest to the cursor selection, specifically
%find the threshpts taht is "going up" (1) that is nearest, and the
%threshpt taht is "going down" (0) after that 

%this wont work if you click on the first tiled axes 
zone_threshptsc=new_threshptsc(:,(whichaxes-1)*2-1:(whichaxes-1)*2);
thresh_diff=abs(zone_threshptsc(:,1)-deletethresh);
[~,imin]=min(thresh_diff);
if zone_threshptsc(imin,2)==0
    %it will either be the one before or after this one
    if thresh_diff(imin-1,1)<=thresh_diff(imin+1,1)
        deletethreshi=imin-1;
    else
        deletethreshi=imin+1;
    end
else
    deletethreshi=imin;
end

zone_threshptsc(deletethreshi:deletethreshi+1,:)=[];
zone_threshptsc=[zone_threshptsc;zeros(2,2)];
new_threshptsc(:,(whichaxes-1)*2-1:(whichaxes-1)*2)=zone_threshptsc;
editedthresh=new_threshptsc;
        





    