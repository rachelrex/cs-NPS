%this code takes data saved after cs-NPS experimental runs 

%% 1. calculate Resistance

numZones=6; %num zones recorded
numZonesi=3; %num zones of interest 
I = transpose(data(1:numZones+1:end) * ampsPerVolt);

% 2. Retrieve Raw Voltage Signals -----------------------------------------

% intialize raw voltage matrix
V = zeros(length(data)/(numZones+1), numZones);

% populate matrix with raw voltage signal for pore and each zone
for i = 1:numZones
    V(:,i) = transpose(data(i+1:numZones+1:end));
end

% create total voltage vector
Vtotal = sum(V(:,1:numZones),2);

% 3. Compute Time Vector At Original Sample Rate --------------------------

t = (0:size(I,1)-1).' ./ sampleRate; 
 
%4. Calculate Resistance Signals -----------------------------------------




figure 
hold on
for i = 1:numZones
    plot(V(:,i))
end

R = V ./ I;



%% 2. Filter Resistance Signals --------------------------------------------
fpass = 50; % passband frequency of the lowpass filter (in Hz)
window = 31; % averaging window width (in units of samples) 
window2= 71;
polynomial=3;
polynomial2=1;
polynomial3=2;

numsgo=10;
Rfalt1= sgolayfilt(lowpass(R, fpass, sampleRate),polynomial,window);  
Rfalt2= sgolayfilt(lowpass(R, fpass, sampleRate),polynomial,window2);  
Rfalt3= sgolayfilt(lowpass(R, fpass, sampleRate),polynomial2,window2);  
Rfalt4= sgolayfilt(lowpass(R, fpass, sampleRate),polynomial3,101);  
If= sgolayfilt(lowpass(I, fpass, sampleRate),polynomial3,101);  
for i=numsgo
    Rfalt1= sgolayfilt(Rfalt1,polynomial,window); 
    Rfalt2= sgolayfilt(Rfalt2,polynomial,window2); 
    Rfalt3= sgolayfilt(Rfalt3,polynomial2,window2); 
    Rfalt4= sgolayfilt(Rfalt4,polynomial3,101);
    If= sgolayfilt(If,polynomial3,101);
end
Rlow=lowpass(R,50,sampleRate);
Rlow=Rlow(window2:end-window2,:);

Rfalt1=Rfalt1(window2:end-window2,:);
Rfalt2=Rfalt2(window2:end-window2,:);
Rfalt3=Rfalt3(window2:end-window2,:);
Rfalt4=Rfalt4(window2:end-window2,:);
If=If(window2:end-window2,:);



f1=figure;
tiledlayout(3,1)
ax1=nexttile; 
plot(Rlow(:,1),'LineWidth',2.0)
hold on
plot(Rfalt1(:,1),'LineWidth',2.0)
plot(Rfalt2(:,1),'LineWidth',2.0)
plot(Rfalt3(:,1),'LineWidth',2.0)
plot(Rfalt4(:,1),'LineWidth',2.0)
ax2=nexttile; plot(Rfalt1(:,2),'k')
ax3=nexttile; plot(Rfalt1(:,3),'r')
%ax3=nexttile; plot(If,'r')
linkaxes([ax1 ax2 ax3],'x')

Rf=Rfalt1;




%% 3. plotting



f1=figure;
tiledlayout(3,1)
ax1=nexttile; 
plot(Rlow(:,1),'LineWidth',2.0)
hold on
plot(Rfalt1(:,1),'LineWidth',2.0)
plot(Rfalt2(:,1),'LineWidth',2.0)
plot(Rfalt3(:,1),'LineWidth',2.0)
plot(Rfalt4(:,1),'LineWidth',2.0)
ax2=nexttile; plot(Rfalt1(:,2),'k')
ax3=nexttile; plot(Rfalt1(:,3),'r')
%ax3=nexttile; plot(If,'r')
linkaxes([ax1 ax2 ax3],'x')

Rf=Rfalt1;



%% 4. Get stderror

 %pick region for stdev calculation --> uncomment when you are doing the
 stderror=zeros(1,3);
   
        d=datacursormode(f1);
          d.Enable='on';
         d.DisplayStyle='window';
         startstdi=input('click where you want region to start, then press enter');
         if isempty(startstdi)==1
             vals=getCursorInfo(d);
             startstd=vals.Position(1,1);
         end
         d.Enable='off';
         d.Enable='on';
         endstdi=input('click where you want region to end, then press enter');
         if isempty(endstdi)==1
             vals=getCursorInfo(d);
             endstd=vals.Position(1,1);
         end
         for i=1:3
            stderror(1,i)=std(Rf(startstd:endstd,i));
         end

%% 5. and 9. fit baseline
%for asls if not downsampled 
lambda1 = 1e12; %larger means smoother background 
% lambda2=1e5;
%for asls if downsampled
% lambda2=5e3; %50
% lambda3=1e5; %1000
p=0; % less than 0.5 means negative is strongly supressed %0.02
maxiter=20; % maximum iterations 
noisez1=stderror(1,1);
noisez2=stderror(1,2);
noisez3=stderror(1,3);

aslsparamz1=struct('lambda', lambda1, 'p', p, 'max_iter', maxiter, 'noise_margin', noisez1);
aslsparamz2=struct('lambda', lambda1, 'p', p, 'max_iter', maxiter, 'noise_margin', noisez2);
aslsparamz3=struct('lambda', lambda1, 'p', p, 'max_iter', maxiter, 'noise_margin', noisez3);
aslsparam=[aslsparamz1;aslsparamz2; aslsparamz3];
% asls2param=aslsparam;
% asls3param=aslsparam;
% for i=1:numZonesi
%     asls2param(i).lambda=lambda2;
%     asls3param(i).lambda=lambda3;
% end



yasls=zeros(length(Rf(:,1)),numZonesi);
yasls2=yasls;
yasls3=yasls;
ydetrend=yasls;
yas2det=yasls;
yas3det=yasls;
for i=1:numZonesi
       yasls(:,i)=ASLS2(Rf(:,i),aslsparam(i,1));
       %yasls2(:,i)=ASLS2(Rf(:,i),asls2param(i,1));
       %yasls3(:,i)=ASLS2(Rf(:,i),asls3param(i,1));
end

load chirp 
y=y(1:2458,:); 
sound(y,Fs);
%% plot yasls results - OPTIONAL - just for visualization

figure
tiledlayout(3,1)
ax1=nexttile; plot(yasls(:,1))
ax2=nexttile; plot(yasls(:,2))
ax3=nexttile; plot(yasls(:,3))
axhand=[ax1;ax2;ax3];


for i=1:numZonesi
    hold(axhand(i,1),'on')
    axes(axhand(i,1));
    plot(Rf(:,i))
    %plot(yasls2(:,i))
    %plot(yasls3(:,i))
    plot(yasls(:,i))
    plot(Rfalt3(:,i))
end
linkaxes([ax1 ax2 ax3],'x')
% Rfsmooth=smoothdata(Rf,'movmean',91);
% plot(Rfsmooth(:,1))
 %yas3det will follow the sizing pulse, but MAY NOT go into the sq pulse depending on transit time, rely on yas2det for that  


    
 %% 6. and 10. subtract baseline_replot
f1=figure;
hold on
for i=1:numZonesi 
    ydetrend(:,i)=Rf(:,i)-yasls(:,i);
    yas2det(:,i)=Rfalt2(:,i)-yasls(:,i);
    yas3det(:,i)=Rfalt3(:,i)-yasls(:,i);
    plot(ydetrend(:,i))
    plot(yas2det(:,i))
    plot(yas3det(:,i),'LineWidth',2)
end

%stop here and save your pre-processed data 
%% 7a. calculate noise again- use same region as before
%you can use either 7a or 7b to recalculate the "zeroed" noise- this one is
%less work 

         for i=1:3
            stderror(2,i)=std(ydetrend(startstd:endstd,i));
   
         end

%% 7b. calculate noise again
%pick region for stdev calculation 

       d=datacursormode(f1);
          d.Enable='on';
         d.DisplayStyle='window';
         startstdi=input('click where you want region to start, then press enter');
         if isempty(startstdi)==1
             vals=getCursorInfo(d);
             startstd=vals.Position(1,1);
         end
         d.Enable='off';
        

         d.Enable='on';
         endstdi=input('click where you want region to end, then press enter');
         if isempty(endstdi)==1
             vals=getCursorInfo(d);
             endstd=vals.Position(1,1);
         end
          std4errorplot=zeros(endstd-startstd+1,3);
         for i=1:3
            stderror(3,i)=std(ydetrend(startstd:endstd,i));
            std4errorplot(:,i)=ydetrend(startstd:endstd,i);
         end



%% 8. reset stderror

stderror(1,:)=stderror(2,:)

% go back to (9)"fit baseline" section after this and proceed until
% (10)"subtract baseline_replot", then save your workspace 