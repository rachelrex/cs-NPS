f1= figure('Name','Window of Interest');

a4=subplot(2,1,1);

a6=subplot(2,1,2);
f2=figure('Name','Big Picture');
a7=axes(f2);
f3=figure('Name','All Zones');
a8=axes(f3);
wp=2000;
wp2=10000;
dinterest=10;
Deff=14.73; %effective diameter 
L=9050;%overall length of device 
dinterest_line=(((dinterest^3)/((Deff^2)*L))*(1/(1-(0.8*((dinterest/Deff)^3)))))*yasls;
dinterest_line=dinterest_line+yasls;

lengthpost=length(ydetrend(:,1));
xaxispost = linspace(1,lengthpost,lengthpost)';
ynz=ydetrend+yasls;
yasdetnz=yas2det+yasls;
sz=3;
newsize=zeros(length(pulses(:,1)),4);
squeezesizing=zeros(1,4);
if length(pulses(1,:))<13
    pulses(:,13:17)=zeros(length(pulses(:,1)),5);
end
for i=1:length(pulses(:,1))
     cla(a8)
    reset(a8)
%i=1;
for j=1:3
    zonenum=j;
    switch zonenum
        case 1
             a=pulses(i,1)-200;
             if pulses(i,13)==0
                 b=a+wp+1000;
             else
                b=pulses(i,13)+300;
             end
             pulses2plot=[pulses(i,1:4) pulses(i,12:13)];
             if en==1
                 if pulses(i,21)~=0
                     pulses2plot=[pulses2plot pulses(i,21:22)];
                     b=pulses(i,22)+100;
                 end
             end
        case 2
            if pulses(i,5)==0
                'no second zone recorded'
                continue
            end 
             a=pulses(i,5)-100;
             if pulses(i,15)==0
                 b=a+wp;
             else
                b=pulses(i,15)+300;
             end
             pulses2plot=[pulses(i,5:6) pulses(i,14:15)];
             if en==1
                 if pulses(i,23)~=0
                     pulses2plot=[pulses2plot pulses(i,23:24)];
                     b=pulses(i,24)+100;
                 end
             end
        case 3
            if pulses(i,7)==0
                'no third zone recorded'
                continue
            end 
             a=pulses(i,7)-100;
             if pulses(i,17)==0
                 b=a+wp;
             else
                b=pulses(i,17)+300;
             end
             
             pulses2plot=[pulses(i,7:8) pulses(i,16:17)];
             if en==1
                 if pulses(i,25)~=0
                     pulses2plot=[pulses2plot pulses(i,25:26)];
                     b=pulses(i,26)+100;
                 end
             end
    end
   
      cla(a4)
    reset(a4)
   
    cla(a6)
    reset(a6)
    cla(a7)
    reset(a7)
    plot(a4,xaxispost(a:b),ydetrend(a:b,j),'r');
       xlim(a4,[a,b]) 
        plot(a6,xaxispost(a-wp:b+wp),ynz(a-wp:b+wp,j),'r');
        yline(a4,0,'k')
        
        hold(a4,'on')
       
        hold(a6,'on')
        hold(a7,'on')
        plot(a6,xaxispost(a-wp:b+wp),yasls(a-wp:b+wp,j),'k')
        %plot(a6,xaxispost(a-wp:b+wp),yaslsnew(a-wp:b+wp,j),'m')
        %plot(a6,xaxispost(a-wp:b+wp),yaslsog(a-wp:b+wp,j),'g')
        
        plot(a6,xaxispost(a-wp:b+wp),yasdetnz(a-wp:b+wp,j),'b','LineWidth',1)
        %plot(a4,xaxispost(a-wp:b+wp),ydetrendnew(a-wp:b+wp,j),'k','LineWidth',2)
        plot(a7,xaxispost(a-wp:b+wp),ynz(a-wp:b+wp,j).*(-1),'r');
        plot(a6,xaxispost(a-wp:b+wp),dinterest_line(a-wp:b+wp,j),'Color',"#77AC30",'LineWidth',2)      
      
        
        plot(a7,xaxispost(a-10000:b+10000),yasls(a-10000:b+10000,j).*(-1),'k')
        plot(a7,xaxispost(a-10000:b+10000),ynz(a-10000:b+10000,j).*(-1),'r');
        
        title(a4,strcat(' zone number=',num2str(zonenum)));
        title(a6,strcat(' diameter=',num2str(pulses(i,11))));
        scatter(a4,xaxispost(a:b),ydetrend(a:b,j),sz,'r');

        
        plot(a4,xaxispost(a:b),yas2det(a:b,j),'b','LineWidth',1.25)
        
        for k=1:length(pulses2plot)
            if pulses2plot(1,k)~=0
                xline(a4,pulses2plot(1,k),'k','LineWidth',1.25)
                %xline(a6,pulses2plot(1,k),'k','LineWidth',1.25)
%                 if mod(k,2)==1
%                     xline(a4,pulses2plot(1,k)-50,'r','LineWidth',1.25)
%                 else
%                     xline(a4,pulses2plot(1,k)+50,'r','LineWidth',1.25)
%                 end
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
            newsize(i,3)=deltaR;
            newsize(i,4)=Rbase;
            squeezesizing(i,1)=pulses(i,11);
            if pulses(i,3)>1
                squeezesizing(i,2)=mean(ydetrend(pulses(i,3):pulses(i,4),1))/mean(yasls(pulses(i,3):pulses(i,4),1));
            end
            if pulses(i,5)>1
                squeezesizing(i,3)=mean(ydetrend(pulses(i,5):pulses(i,6),2))/mean(yasls(pulses(i,5):pulses(i,6),2));
            end
            if pulses(i,7)>1
                squeezesizing(i,4)=mean(ydetrend(pulses(i,7):pulses(i,8),3))/mean(yasls(pulses(i,7):pulses(i,8),3));
            end 
        end
end %ADD DEBUG HERE!!
end