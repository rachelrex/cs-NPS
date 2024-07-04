pulsesnew=pulses(pulses(:,3)~=0,:);
lp=length(pulsesnew(:,1));
endzero=lp;
for i=1:lp
    if pulsesnew(i,1)==0
        endzero=i-1;
        break
    end
end
lp=endzero;

%make sure the data you want to look at is in your workspace as "pulses" 
%currdata: col1=sizing start, col2=time in sizing channel, col3=diameter 
currdata=zeros(lp,3);
pulsesnew=pulsesnew(1:lp,:);
currdata(:,1)=pulsesnew(:,1);
currdata(:,2)=(pulsesnew(:,2)-pulsesnew(:,1))./10;
currdata(:,3)=pulsesnew(:,11);

%% load reference speed data 
filename="/Users/rachelrex/Documents/drive-download-20220825T191243Z-001/analyzed_data/swssnps/hl60_nb4_9_13_23.mat";
load(filename,"blue","red","green","purple","orange");
load(filename,"cell_hl60","cell_nb4");



%% low speed graph
diameterx=linspace(7,16,19)';
%non-empty are: 
%10l- 80,90,100
%10r-70,80,90,100,110
%11r-70,90,100




figure 
hold on 
scatter(cell_hl60{2,2}.Diameter,cell_hl60{2,2}.timeSizemsec,'filled','MarkerFaceColor',	orange(1,1))
scatter(cell_nb4{3,2}.Diameter,cell_nb4{3,2}.timeSizemsec,'filled','MarkerFaceColor',green(2,1))
scatter(currdata(:,3),currdata(:,2),'filled','MarkerFaceColor',blue(2,1))

%scatter(cell_nb4{3,2}.Diameter,cell_nb4{3,2}.timeSizemsec,'filled','MarkerFaceColor',green(3,1))
ylabel("Time in Sizing Channel (msec)")
xlabel("Diameter (microns)") 
legend(strcat(string(cell_hl60{1,2}),string(cell_hl60{2,1})),strcat(string(cell_nb4{1,2}),string(cell_nb4{3,1})),'9R 80 kas')


%% high speed graph
diameterx=linspace(7,16,19)';
%non-empty are: 
%10l- 80,90,100
%10r-70,80,90,100,110
%11r-70,90,100




figure 
hold on 
scatter(cell_hl60{4,2}.Diameter,cell_hl60{4,2}.timeSizemsec,'filled','MarkerFaceColor',	orange(1,1))
scatter(cell_nb4{6,2}.Diameter,cell_nb4{6,2}.timeSizemsec,'filled','MarkerFaceColor',green(2,1))
%scatter(currdatakasumi110_10523(:,3),currdatakasumi110_10523(:,2),'filled','MarkerFaceColor',blue(2,1))
scatter(currdata(:,3),currdata(:,2),'filled','MarkerFaceColor',purple(2,1))


%scatter(cell_nb4{3,2}.Diameter,cell_nb4{3,2}.timeSizemsec,'filled','MarkerFaceColor',green(3,1))
ylabel("Time in Sizing Channel (msec)")
xlabel("Diameter (microns)") 
legend(strcat(string(cell_hl60{1,2}),string(cell_hl60{4,1})),strcat(string(cell_nb4{1,2}),string(cell_nb4{6,1})),'9R 110 kas')