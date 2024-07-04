f1= figure('Name','Window of Interest');
a4=subplot(2,1,1);

a6=subplot(2,1,2);
wp=6000;
dinterest=10;
Deff=12.78153763; %effective diameter 
L=12650;%overall length of device 
dinterest_line=(((dinterest^3)/((Deff^2)*L))*(1/(1-(0.8*((dinterest/Deff)^3)))))*yasls;
dinterest_line=dinterest_line+yasls;

lengthpost=length(ydetrend(:,1));
xaxispost = linspace(1,lengthpost,lengthpost)';
ynz=ydetrend+yasls;
yasdetnz=yas2det+yasls;
sz=3;
newsize=zeros(length(pulses(:,1)),2);
squeezesizing=zeros(1,4);
if length(pulses(1,:))<13
    pulses(:,13)=zeros(length(pulses(:,1)),1);
end
for i=1:length(pulses(:,1))
%i=1;
for j=1:1
    zonenum=1;
    switch zonenum
        case 1
             a=pulses(i,1)-200;
             if pulses(i,13)==0
                 b=a+wp+1000;
             else
                b=pulses(i,13)+100;
             end
             pulses2plot=[pulses(i,1:4) pulses(i,12:13)];
        case 2
            
             a=pulses(i,5)-100;
             if pulses(i,15)==0
                 b=a+wp;
             else
                b=pulses(i,15)+100;
             end
             pulses2plot=[pulses(i,5:6) pulses(i,14:15)];
        case 3
             a=pulses(i,7)-100;
             if pulses(i,17)==0
                 b=a+wp;
             else
                b=pulses(i,17)+100;
             end
             pulses2plot=[pulses(i,7:8) pulses(i,16:17)];
    end
   
      cla(a4)
    reset(a4)
   
    cla(a6)
    reset(a6)
   
    plot(a4,xaxispost(a:b),ydetrend(a:b,j),'r');
       xlim(a4,[a,b]) 
        plot(a6,xaxispost(a-wp:b+wp),ynz(a-wp:b+wp,j),'r');
        yline(a4,0,'k')
        
        hold(a4,'on')
       
        hold(a6,'on')
        
        plot(a6,xaxispost(a-wp:b+wp),yasls(a-wp:b+wp,j),'k')
        %plot(a6,xaxispost(a-wp:b+wp),yaslsnew(a-wp:b+wp,j),'m')
        %plot(a6,xaxispost(a-wp:b+wp),yaslsog(a-wp:b+wp,j),'g')
        
        plot(a6,xaxispost(a-wp:b+wp),yasdetnz(a-wp:b+wp,j),'b','LineWidth',2)
        plot(a6,xaxispost(a-wp:b+wp),dinterest_line(a-wp:b+wp,j),'Color',"#77AC30",'LineWidth',2)      
      
        title(a4,strcat(' zone number=',num2str(zonenum)));
        title(a6,strcat(' diameter=',num2str(pulses(i,11))));
        scatter(a4,xaxispost(a:b),ydetrend(a:b,j),sz,'r');

        
        plot(a4,xaxispost(a:b),yas2det(a:b,j),'b','LineWidth',1.25)
        
        for k=1:length(pulses2plot)
            if pulses2plot(1,k)~=0
                xline(a4,pulses2plot(1,k),'k','LineWidth',1.25)
                %xline(a6,pulses2plot(1,k),'k','LineWidth',1.25)
                if mod(k,2)==1
                    xline(a4,pulses2plot(1,k)-50,'r','LineWidth',1.25)
                else
                    xline(a4,pulses2plot(1,k)+50,'r','LineWidth',1.25)
                end
            end
        end
        if j==1
            deltaR=mean(ydetrend(pulses(i,1):pulses(i,2),1))
            deltRog=pulses(i,9)
            Rbase=mean(yasls(pulses(i,1):pulses(i,2),1))
            Rbaseog=pulses(i,10)
            deltaRR=deltaR/Rbase
            size=calcsize(Deff,L,deltaRR)
            newsize(i,1)=pulses(i,11);
            newsize(i,2)=size;
            
            
        end
end %ADD DEBUG HERE!!
end