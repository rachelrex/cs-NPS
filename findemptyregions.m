emptyregions=zeros(2,2);

%%

wp=8200;
dinterest=10;
Deff=16.964379152595436416572891977727; %effective diameter 
L=9050;%overall length of device 
dinterest_line=(((dinterest^3)/((Deff^2)*L))*(1/(1-(0.8*((dinterest/Deff)^3)))))*yasls;
linecolor=["#D95319";"#7E2F8E";"#77AC30"];
linecolor2=["#f09973";"#d574e8";"#b0f05d"];

lengthpost=length(ydetrend(:,1));
xaxispost = linspace(1,lengthpost,lengthpost)';

numevents=length(events(:,1));
f1=figure; 
ax=axes;
i=1;
k=1;
while i<=numevents
    cla(ax)
    reset(ax)
    endreg=events(i,1)-100;
    if i==1
        startreg=201;
    else
       if events(i-1,6)==0
            startreg=events(i-1,1)+wp;
        else
            startreg=events(i-1,6)+100;
        end
    end
    if endreg-startreg>=10000
         plot(xaxispost(startreg-200:endreg+200,1),dinterest_line(startreg-200:endreg+200,1),'Color',"#77AC30",'LineWidth',2)  
         hold(ax,'on')
        for j=1:3
            plot(xaxispost(startreg-200:endreg+200,1),ydetrend(startreg-200:endreg+200,j),'Color',linecolor(j,1))
            plot(xaxispost(startreg-200:endreg+200,1),yas2det(startreg-200:endreg+200,j),'Color',linecolor2(j,1),'LineWidth',2)

        end
        xline(ax,startreg,'k')
        xline(ax,endreg,'k')
        useranswer=input('region good?, leave empty for yes, type something for no');
        if isempty(useranswer)
            emptyregions(k,1)=startreg;
            emptyregions(k,2)=endreg;
            k=k+1;
        else
            useranswer2=input('do you want to repick start empty for no, 1 for start');
            if isempty(useranswer2)==0
                d=datacursormode(f1);
                 d.Enable='on';
             d.DisplayStyle='window';
                 startstdi=input('click where you want region to start, then press enter');
                if isempty(startstdi)==1
                 vals=getCursorInfo(d);
                 startreg=vals.Position(1,1);
                end
             d.Enable='off';
            end
            useranswer3=input('do you want to repick end empty for no, 1 for end');
            if isempty(useranswer3)==0
                d=datacursormode(f1);
                 d.Enable='on';
             d.DisplayStyle='window';
                 startstdi=input('click where you want region to start, then press enter');
                if isempty(startstdi)==1
                 vals=getCursorInfo(d);
                 endreg=vals.Position(1,1);
                end
             d.Enable='off';
            end
            if isempty(useranswer2)==0 || isempty(useranswer3)==0
                emptyregions(k,1)=startreg;
                emptyregions(k,2)=endreg;
                k=k+1;
            end
        end
    else
        'skipped'
        i
    end
    i=i+1;
end


%%
for i=1:length(emptyregions(:,1))
    for j=1:3
        emptyregions(i,j+2)=std(ydetrend(emptyregions(i,1):emptyregions(i,2),j));
    end
end
figure
scatter(emptyregions(:,1),emptyregions(:,3))

%% find "moving average std" 
movingstd=zeros(2,4);
k=1;
fwait=waitbar(0,'Calculating moving std');
for i=1:length(emptyregions(:,1))
    waitbar(i/length(emptyregions(:,1)),fwait,'Calculating moving std');
    a=emptyregions(i,1);
    b=emptyregions(i,2);
    window=10000;
     if b-window-1>=a
        for j=a:b-window-1

            movingstd(k,2)=std(ydetrend(j:j+window-1,1));
            movingstd(k,3)=std(ydetrend(j:j+window-1,2));
            movingstd(k,4)=std(ydetrend(j:j+window-1,3));
            movingstd(k,1)=j;
            k=k+1;
        end
     end
end
figure
scatter(movingstd(:,1),movingstd(:,2),5)
hold on
scatter(movingstd(:,1),movingstd(:,3),5)
scatter(movingstd(:,1),movingstd(:,4),5)


%getting std of region where window is x, starting pt is a and ending pt is
%b, first window is a:a+x-1, last window is b-x-1:b
%imagining 1-10 with x=4,  1-4;2-5;3-6;4-7;5-8;6-9;7-10

%% find moving average noise (ish) 
%moving average over 200,000

movavgstd=movmean(movingstd(:,2),200000);
figure
scatter(movingstd(:,1),movingstd(:,2),5)
hold on
plot(movingstd(:,1),movavgstd,'LineWidth',3)

%% find min and maxes of movingstd and movavgstd

minmovstd=min(movingstd(:,2))
maxmovstd=max(movingstd(:,2))
minavgstd=min(movavgstd)
maxavgstd=max(movavgstd)




