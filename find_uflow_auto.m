function [slopesmatrix,slopebtwnpt,slopesmatrixcurve,kx,kynew,kfit2]=find_uflow_auto(cell_array,celltypes,pressure,devicename,mindiameter,diameter_range,ch_width,ch_height)
%this code will be good for debugging when there are issues with find_uflow_auto_curve_2   

%cell_array should contain a "cell array" with all Cell types you want
    %to include in your uflow calculation, (Cell= biological cell) (cell= 
    %an element of a cell matrix) one Cell type per cell of the array in
    %separate columns
    % the format of the each  Cell type will be a cell array with the first
    % row as the device name (where replicate devices end with a
    % 1,2,3..etc) and the first column as a pressue
    
    %celltypes= strings of the cell type names, in the order they appear in
    %the "cell_array" 
    
    %pressure= the pressure at which you want to calculate the uflow
    
    %devicename= the device name of which you want to calculate the uflow
    %(replicate devices will be included in the calculation as long as they
    %have identical names ending with different numbers
    
    %fixedslope=3.5
    
    %mindiameter=9.5 microns (the minimum cell diameter at which you will
    %calculate wCDI, should be the same across all devices used in a study! 
    %this code has only been tested with this set to 9.5 microns so change at
    % your own peril, must be a whole number or in increments of 0.5) 
    
    %diameter_range= [smallestdiamter biggestdiameter] the range of cell diameters you measured, the minimum
    %can be greater than your "mindiameter"
    
    
    %(ex), pressure=110, devicename='9enR' - may include '9enR1' and
    %'9enR2'
    blue=["#0000FF" "#4DBEEE" "#0072BD"];
    red=["#A2142F" "#FF1B00"  "#E74C3C"];
    green=["#196F3D" "#2ECC71" "#719B02"];
    purple=["#4A235A" "#9B59B6" "#D2B4DE"];
    orange=["#F39C12"];
    colors=[blue;red;green;purple];
    
    num_cell_type=length(cell_array(1,:));
    diameterx=linspace(diameter_range(1,1),diameter_range(1,2),(diameter_range(1,2)-diameter_range(1,1))*2+1)';
    diameterxmin=linspace(mindiameter,diameter_range(1,2),(diameter_range(1,2)-mindiameter)*2+1)';
    of_interest=cell(num_cell_type,4);
    num_datasets=num_cell_type;
    i=1; %the row number of of_interest that we're on 
    k=1; %the col number of celltypes that we're on 
    %create a cell array containing all the tables of interest 
    while i<=num_datasets
        current=cell_array{1,k}; %the current cell type we're looking at
        
        %find the device(s) we want to look at (will be in the first row,
        %this cannot handle>9 replicate)
        correct_col=0; %will contain the columns of the device(s) we want to look at 
        for j=2:length(current(1,:)) %starts at 2 becasue the first col should be empty
            dev=char(current{1,j});
            %check if there are replicate devices
            devrep=dev(1:end-1);
            
            
            devicename=char(devicename);
            if devicename(end)~='R' && devicename(end)~='L'
                startname=length(devicename);
                if startname>2 
                    if devicename(end-1)=='R' || devicename(end-1)=='L'
                        devrep=dev;
                    else
                       devrep=dev(1:startname);
                    end
                else
                    
                    devrep=dev(1:startname);
                end
            end

            
            if string(dev)==string(devicename)
                %not a replicate, but we are in the correct column
                correct_col=j;
            elseif string(devrep)==string(devicename)
                %it is a replicate and we are in one of the correct columns
                if correct_col==0
                    correct_col=j;
                else
                    correct_col=[correct_col;j];
                end
            end
        end
        if correct_col==0
            'device not found'
        end
        
        %find the pressures we want to look at 
        correct_row=0;
        current_pressures=cell2mat(current(:,1));
        for j=1:length(current_pressures(:,1))
            if current_pressures(j,1)==pressure
                correct_row=j+1;
                break
            end
        end
        if correct_row==0
            'pressure not found'
        end
        
        if correct_row(1,1)~=0 && correct_col(1,1)~=0
            %check to make sure there is actually data in this cell type at the
            %given device and pressure
            numcorrectcol=length(correct_col);
            j=1;
            while j<=numcorrectcol
                if isempty(current{correct_row,correct_col(j)})
                    'data not found in one of the cells'
                    correct_col(j)=[];
                    numcorrectcol=numcorrectcol-1;
                else
                    j=j+1;
                end
            end
            if numcorrectcol==0
                %of_interest(num_datasets,:)=[];
                num_datasets=num_datasets-1;
            else
            %now check how many device replicates there are 
                num_rep=length(correct_col);
                num_datasets=num_datasets+num_rep-1;
                for j=1:num_rep 
                    of_interest{i,1}=current{correct_row,correct_col(j)};
                    of_interest{i,2}=celltypes(1,k);
                    i=i+1;
                end
            end
        else
            of_interest(num_datasets,:)=[];
            num_datasets=num_datasets-1;
            
        end
        
        k=k+1; %move onto the next cell type 
        
        
    end
    numthingsofint=length(of_interest(:,1));
for i=1:numthingsofint
    if isempty(of_interest{i,1})
        of_interest(i,:)=[];
        numthingsofint=numthingsofint-1;
    end
end

totalsamp=0;    
all=0;
k=1;
j=1;
figure
hold on
for i=1:length(of_interest(:,1))
    dev=unique(of_interest{i,1}.deviceName);
    dateofsamp=unique(of_interest{i,1}.date);
    of_interest{i,3}=strcat(string(pressure),string(of_interest{i,2}),string(dev));
    of_interest{i,1}=of_interest{i,1}(of_interest{i,1}.Diameter<=diameter_range(1,2),:);
    of_interest{i,4}=strcat(string(of_interest{i,2}),string(dev)," ",string(dateofsamp));
    newset=[of_interest{i,1}.Diameter of_interest{i,1}.timeSizemsec];
    if all(1,1)==0
        all=newset;
        plotcolor=colors(1,1);
    else
        all=[all;newset];
        if string(of_interest{i,2})~=string(of_interest{i-1,2})
            k=k+1;
            j=1;
        else
            j=j+1;
        end
        plotcolor=colors(k,j);
    end
    
    scatter(of_interest{i,1}.Diameter,of_interest{i,1}.timeSizemsec,'filled','MarkerFaceColor',	plotcolor) 
    totalsamp=totalsamp+height(of_interest{i,1});
end
legend(of_interest(:,3))

ylabel("Time in Sizing Channel (msec)")
xlabel("Diameter (microns)") 


% start of calculating uflow 

if totalsamp>= 150
    al=2;
elseif totalsamp<=25
    al=10
else
    al=5;
end

shp=alphaShape(all,al,'RegionThreshold',1);
[bf,P]=boundaryFacets(shp);
line_handle=plot(P(:,1),P(:,2));

if isempty(selfintersect(P(:,1),P(:,2)))~=1
    while isempty(selfintersect(P(:,1),P(:,2)))~=1
        delete(line_handle)
        'alpha adjusted to prevent self-intersection of boundary'
        al=al+.01
        shp=alphaShape(all,al,'RegionThreshold',1);
        [bf,P]=boundaryFacets(shp);
        line_handle=plot(P(:,1),P(:,2));
    end
    
end

tf=inShape(shp,all(:,1),all(:,2));
[outx,outy]=find(~tf);

scatter(all(outx,1),all(outx,2),100,"x",'k')

kx=P(:,1);
ky=P(:,2);
%find lower boundary 
[~,minindex]=min(ky);
%start at the minimum time in pore (should be smallest on-axis) 
if minindex==length(kx)
    kx=[kx(end,1);kx(1:end-1)];
    ky=[ky(end,1);ky(1:end-1)];
end
[~,minindex]=min(ky);
if minindex~=1
    kx=[kx(minindex:end,1);kx(1:minindex-1)];
    ky=[ky(minindex:end,1);ky(1:minindex-1)];
end

tooearly=1;
endofchosen=0;
while tooearly==1
    slopebtwnpt=zeros(length(kx)-1,1);
    for i=1:(length(kx)-1)
        slopebtwnpt(i,1)=(ky(i+1,1)-ky(i,1))/(kx(i+1,1)-kx(i,1));
        
        if slopebtwnpt(i,1)<0 && slopebtwnpt(i,1)> -50
            
            endofchosen=1;
        elseif i>1
            if kx(i+1)<kx(i) && ky(i+1)<ky(i)% was i<i-1
            endofchosen=1;
            end
        end
        if endofchosen==1
            if i<=2 && length(kx)>3
                'error probs wrong'
                kx(1)=[];
                ky(1)=[];
                endofchosen=0;
            else
                tooearly=0;
            end
            break
        end
    end
end

kx=kx(1:i);
ky=ky(1:i);



scatter(kx,ky,100,"x",'r')



[minDindex,~]=find(diameterx==mindiameter);


slopesmatrix{1,1}=devicename;
 slopesmatrix{1,2}=pressure;
 slopesmatrix{1,3}=length(all(:,1));
 slopesmatrix{1,4}=al;

kfitunfix=polyfit(kx,ky,1);
klineunfix=polyval(kfitunfix,diameterx);
plot(diameterx,klineunfix,'--','Color',orange(1,1),'LineWidth',1.5)

slopesmatrix{1,5}=kfitunfix(1,1);
slopesmatrix{1,6}=klineunfix(minDindex,1);

%b_fit=mean(ky-fixedslope*kx);
%kfit=[fixedslope b_fit];
%kline=polyval(kfit,diameterx);
%plot(diameterx,kline,'Color',orange(1,1),'LineWidth',1.5)
 
%slopesmatrix{1,7}=kfit(1,1);
%slopesmatrix{1,8}=kline(minDindex,1);
%slopesmatrix{1,9}=round(800/(round(kline(minDindex,1),1)),1);

scatter(P(:,1),P(:,2),'filled','k')

figure
hold on 
k=1;
j=1;
for i=1:length(of_interest(:,1))
    if i==1
        plotcolor=colors(1,1);
        sampleinfo=of_interest{i,4};
    else
      sampleinfo=strcat(sampleinfo,"; ",of_interest{i,4});
        if string(of_interest{i,2})~=string(of_interest{i-1,2})
            k=k+1;
            j=1;
        else
            j=j+1;
        end
        plotcolor=colors(k,j);
    end
    
    scatter(of_interest{i,1}.Diameter,800./of_interest{i,1}.timeSizemsec,'filled','MarkerFaceColor',	plotcolor) 
    
    
end
legend(of_interest(:,3))

ylabel("Particle Velocity (mm/sec)")
xlabel("Diameter (microns)") 

plot(P(:,1),800./P(:,2))
scatter(all(outx,1),800./all(outx,2),100,"x",'k')
kynew=800./ky;
scatter(kx,kynew,100,"x",'r')

syms vf
% w0=14;
% h=19.45;


D=sqrt(ch_width*ch_height*4/pi());
%D=18.4;
solution=zeros(length(kx),1);
for i=1:length(kx)
    d=kx(i);
    vp=kynew(i);
eqnright=(-2/(3*(D^2)))*vf*(d^2)+vf;
solution(i)=vpasolve(eqnright==vp,vf);
end

p1=(-2/(3*(D^2)))*mean(solution);
p2=0;
p3=mean(solution);


slopesmatrixcurve{1,1}=devicename;
 slopesmatrixcurve{1,2}=pressure;
 slopesmatrixcurve{1,3}=length(all(:,1));
 slopesmatrixcurve{1,4}=al;
 
kfit2=[p1 p2 p3];
kline2=polyval(kfit2,diameterx);
plot(diameterx,kline2,'Color',orange(1,1),'LineWidth',1.5); 
slopesmatrixcurve{1,5}=round(p3,1);
slopesmatrixcurve{1,6}=round(kline2(minDindex,1),2);
slopesmatrixcurve{2,6}=mindiameter;
xlim(diameter_range)
% kfit2unfix=polyfit(kx,kynew,2);
% kline2unfix=polyval(kfit2unfix,diameterx);
% plot(diameterx,kline2unfix,'--','Color',orange(1,1),'LineWidth',1.5)
% at0=polyval(kfit2unfix,linspace(0,14,29)');
% slopesmatrixcurve{1,7}=at0(1,1);

% for i=1:(length(kx)-1)
%     slopebtwnpt(i,1)=(ky(i+1,1)-ky(i,1))/(kx(i+1,1)-kx(i,1));
% end
allvel=800./all(:,2);
smallestvel=min(allvel);

ylim([smallestvel-1 p3])
%klinevel=800./kline;
klineunfixvel=800./klineunfix


%give info on samples in slopesmatrixcurve

slopesmatrixcurve{1,7}=sampleinfo;




% start of calculating new stuff

for j=1:length(of_interest(:,1))
    all=[of_interest{j,1}.Diameter of_interest{j,1}.timeSizemsec];
    totalsamp=height(of_interest{j,1});
    if totalsamp>= 150
        al=2;
    elseif totalsamp<=25
        al=10;
    else
        al=5;
    end
    figure 
    hold on 
    scatter(of_interest{j,1}.Diameter,800./of_interest{j,1}.timeSizemsec,'filled')
    shp=alphaShape(all,al,'RegionThreshold',1);
    [bf,P]=boundaryFacets(shp);
   % line_handle=plot(P(:,1),P(:,2));

    if isempty(selfintersect(P(:,1),P(:,2)))~=1
        while isempty(selfintersect(P(:,1),P(:,2)))~=1
            delete(line_handle)
            'alpha adjusted to prevent self-intersection of boundary'
            al=al+.01
            shp=alphaShape(all,al,'RegionThreshold',1);
            [bf,P]=boundaryFacets(shp);
            %line_handle=plot(P(:,1),P(:,2));
        end

    end

    tf=inShape(shp,all(:,1),all(:,2));
    [outx,outy]=find(~tf);

    scatter(all(outx,1),all(outx,2),100,"x",'k')

    kx=P(:,1);
    ky=P(:,2);
    %find lower boundary 
    [~,minindex]=min(ky);
    %start at the minimum time in pore (should be smallest on-axis) 
    if minindex==length(kx)
        kx=[kx(end,1);kx(1:end-1)];
        ky=[ky(end,1);ky(1:end-1)];
    end
    [~,minindex]=min(ky);
    if minindex~=1
        kx=[kx(minindex:end,1);kx(1:minindex-1)];
        ky=[ky(minindex:end,1);ky(1:minindex-1)];
    end

    tooearly=1;
    endofchosen=0;
    while tooearly==1
        slopebtwnpt=zeros(length(kx)-1,1);
        for i=1:(length(kx)-1)
            slopebtwnpt(i,1)=(ky(i+1,1)-ky(i,1))/(kx(i+1,1)-kx(i,1));

            if slopebtwnpt(i,1)<0 && slopebtwnpt(i,1)> -50

                endofchosen=1;
            elseif i>1
                if kx(i+1)<kx(i) && ky(i+1)<ky(i)% was i<i-1
                endofchosen=1;
                end
            end
            if endofchosen==1
                if i<=2 && length(kx)>3
                    'error probs wrong'
                    kx(1)=[];
                    ky(1)=[];
                    endofchosen=0;
                else
                    tooearly=0;
                end
                break
            end
        end
    end

    kx=kx(1:i);
    ky=ky(1:i);



%     scatter(kx,ky,100,"x",'r')



    [minDindex,~]=find(diameterx==mindiameter);

% 
%     slopesmatrix{1,1}=devicename;
%      slopesmatrix{1,2}=pressure;
%      slopesmatrix{1,3}=length(all(:,1));
%      slopesmatrix{1,4}=al;

    kfitunfix=polyfit(kx,ky,1);
    klineunfix=polyval(kfitunfix,diameterx);
%     plot(diameterx,klineunfix,'--','Color',orange(1,1),'LineWidth',1.5)
% 
%     slopesmatrix{1,5}=kfitunfix(1,1);
%     slopesmatrix{1,6}=klineunfix(minDindex,1);

%     b_fit=mean(ky-fixedslope*kx);
%     kfit=[fixedslope b_fit];
%     kline=polyval(kfit,diameterx);
%     %plot(diameterx,kline,'Color',orange(1,1),'LineWidth',1.5)

%     slopesmatrix{1,7}=kfit(1,1);
%     slopesmatrix{1,8}=kline(minDindex,1);
%     slopesmatrix{1,9}=round(800/(round(kline(minDindex,1),1)),1);
sampleinfo=of_interest{j,4};
%     scatter(P(:,1),P(:,2),'filled','k')
% 
%     figure
%     hold on 
%     k=1;
%     j=1;
%     for i=1:length(of_interest(:,1))
%         if i==1
%             plotcolor=colors(1,1);
%             sampleinfo=of_interest{i,4};
%         else
%           sampleinfo=strcat(sampleinfo,"; ",of_interest{i,4});
%             if string(of_interest{i,2})~=string(of_interest{i-1,2})
%                 k=k+1;
%                 j=1;
%             else
%                 j=j+1;
%             end
%             plotcolor=colors(k,j);
%         end
% 
%         scatter(of_interest{i,1}.Diameter,800./of_interest{i,1}.timeSizemsec,'filled','MarkerFaceColor',	plotcolor) 
% 
% 
%     end
    legend(of_interest(j,3))
% 
     ylabel("Particle Velocity (mm/sec)")
    xlabel("Diameter (microns)") 
% 
     plot(P(:,1),800./P(:,2))
     scatter(all(outx,1),800./all(outx,2),100,"x",'k')
     kynew=800./ky;
     scatter(kx,kynew,100,"x",'r')

    syms vf
    % w0=14;
    % h=19.45;


    D=sqrt(ch_width*ch_height*4/pi());
    %D=18.4;
    solution=zeros(length(kx),1);
    for i=1:length(kx)
        d=kx(i);
        vp=kynew(i);
    eqnright=(-2/(3*(D^2)))*vf*(d^2)+vf;
    solution(i)=vpasolve(eqnright==vp,vf);
    end

    p1=(-2/(3*(D^2)))*mean(solution);
    p2=0;
    p3=mean(solution);


    slopesmatrixcurve{j+1,1}=devicename;
     slopesmatrixcurve{j+1,2}=pressure;
     slopesmatrixcurve{j+1,3}=length(all(:,1));
     slopesmatrixcurve{j+1,4}=al;

    kfit2=[p1 p2 p3];
    kline2=polyval(kfit2,diameterx);
    plot(diameterx,kline2,'Color',orange(1,1),'LineWidth',1.5); 
    slopesmatrixcurve{j+1,5}=round(p3,1);
    slopesmatrixcurve{j+1,6}=round(kline2(minDindex,1),2);
   % slopesmatrixcurve{2,6}=mindiameter;
    xlim(diameter_range)
 
    allvel=800./all(:,2);
    smallestvel=min(allvel);

    ylim([smallestvel-1 p3])
    %klinevel=800./kline;
    klineunfixvel=800./klineunfix


    %give info on samples in slopesmatrixcurve

    slopesmatrixcurve{j+1,7}=sampleinfo;
    
end
% plot(diameterx,klinevel,'Color',"#DC7633",'LineWidth',1.5)
% 
% plot(diameterx,klineunfixvel,'--','Color',"#DC7633",'LineWidth',1.5)
%  scatter(P(:,1),800./P(:,2),'filled','k')
%  xline(9.5)
