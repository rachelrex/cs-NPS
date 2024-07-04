function [events]=aretheytogether(threshptsc,eventlength,mc)
%eventlength=the amount of data points necessary to encompass an entire
%cell event (all zones), dependent on samplerate, geometry, pressure 
%mc= the number of rows in the threshptsc matrix (the number of corrected
%threshold crossings in whichever zone has the most) 

%events is a matrix where col1-2, col3-4, col5-6 are for indices from each
%zone respectively, each row is a separate cell event
threshz1=threshptsc(:,1:2);
threshz2=threshptsc(:,3:4);
threshz3=threshptsc(:,5:6);

%this assumes that the first threshold crossing of zone 1 is actually an
%event
events=zeros(1,6);
events(1,1)=threshz1(1,1);
events(1,2)=threshz1(2,1);


stopit=0;
i=3;%identifies which row of thresh we are on 
j=2;%identifies which row of events we are on

%this identifies event starts in ZONE 1
while i<=mc && stopit==0 
   if threshz1(i,1)==0  %at the end of the thresh there will just be zeros, so we know we've come to the end
       stopit=1;
   else
       if threshz1(i,1)-events(j-1,1)<eventlength %if there is a pulse within the "eventlength" of an identified pulse then we ignore the second pulse 
           i=i+2;
       else % if the next pulse is after the last events entire "eventlength" then it is the start of a new event
           events(j,1)=threshz1(i,1);
           events(j,2)=threshz1(i+1,1);
           i=i+2;
           j=j+1;
       end
   end
end

%this identifies events in zone 2 that match with those found in zone 1 
[me,~]=size(events);
i=1;%identifies which row of thresh we are on 
j=1;%identifies which row of events we are on
stopit=0;
while i<=mc && stopit==0 && j<=me 
   if threshz2(i,1)==0 %when we've reached the end of threshptsc for this zone there will only be zeros --> stop the while loop 
       stopit=1;
   else
       if threshz2(i,1)<=events(j,2) && threshz2(i+1,1)<= events(j,2) % does the zone 2 pulse start and end before the zone 1 pulse ends for this event? (pulse means crossing and crossing back) 
           i=i+2; %if yes, skip this pulse 
       elseif j+1<=me %is this the last event?  
           % we have determined that the zone 2 pulse starts or ends after the zone 1 pulse ends 
           if threshz2(i+1,1)>=events(j+1,1) %does the zone 2 pulse end after the next event's zone 1 pulse starts? 
               j=j+1; %if yes skip to next event 
           else %we know the zone 2 pulse starts or ends after the zone 1 pulse ends, and ends before the next event's zone 1 pulse starts 
             events(j,3)=threshz2(i,1); %therefore this zone 2 pulse is part of this event 
             events(j,4)=threshz2(i+1,1);
             i=i+2;%go to next pulse
             j=j+1;%go to next event 
           end
       else % this is the last event, so this zone 2 pulse MUST go in this event 
           events(j,3)=threshz2(i,1);
           events(j,4)=threshz2(i+1,1);
           i=i+2;
           j=j+1;
       end
   end
end

%this identifies events in zone 3 that match with those found and matched
%in zone 1 and 2
i=1;%identifies which row of thresh we are on 
j=1;%identifies which row of events we are on
stopit=0;
[mt,~]=size(threshz3);
while i<mc && stopit==0 && j<=me %used to be i<=mc but that lead to an error **
   if threshz3(i,1)==0 %when we've reached the end of threshptsc for this zone there will only be zeros --> stop the while loop 
       stopit=1;
   else
       if threshz3(i,1)<=events(j,4) && threshz3(i+1,1)<=events(j,4)  % does the zone 3 pulse start and end before the zone 2 pulse ends for this event? 
           i=i+2; %skip to next pulse 
       elseif threshz3(i,1)<=events(j,2) && threshz3(i+1,1)<=events(j,2) % does the zone 3 pulse start and end before the zone 2 pulse ends for this event?
           %this is the case when thresholds in zone2 = 0 
           i=i+2;
       elseif j+1<=me %threshold are nonzero in zone 2 of this event, and the zone 3 pulse starts or ends after the zone 2 pulse ends 
           if threshz3(i+1,1)>=events(j+1,1) %does the zone 3 pulse end after the next events zone 1 pulse ends? 
               j=j+1; %go to next event
           else %the zone 3 pulse starts or ends after the zone 2 pulse ends but before the next event's zone 1 pulse ends - this zone 3 pulse is part of this event! 
             events(j,5)=threshz3(i,1); 
             events(j,6)=threshz3(i+1,1);
             i=i+2; %go to next pulse
             j=j+1;  %go to next event 
           end
       else %this is the last event 
           events(j,5)=threshz3(i,1);
           if i+1<=mt
            events(j,6)=threshz3(i+1,1);
           end
           i=i+2;
           j=j+1;
       end
   end
end

    
    