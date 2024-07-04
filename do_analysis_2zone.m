%this code runs your analysis script 

if exist('dmin','var') ~=1
    dmin=input('dmin is not saved in your workspace, type what you want dmin to be right now');
end
if exist('plot_thresh','var') ~=1
    plot_thresh=0;
end
if exist('en','var')~=1
    en=input('is this an extra node device? if so type 1');
    if isempty(en)==1
        en=0;
    elseif en~=1
        en=0;
    end
end
%if en =1 then there is an extra node in the recovery section 
% run if no figures are up (automatically starts from the beginning) 
zones_of_interest=2;
[pulses,threshpts,fig_handles]=mznps_pickpulse_semiauto_nodownsamp_zone(ydetrend,yasls,yas2det,dmin,stderror,yas3det,plot_thresh,0,0,0,en,zones_of_interest);

%% run if your figures are NOT already up- but you want to start from where you left off 
if exist('dmin','var') ~=1
    dmin=input('dmin is not saved in your workspace, type what you want dmin to be right now');
end
if exist('plot_thresh','var') ~=1
    plot_thresh=0;
end
if exist('en','var')~=1
    en=input('is this an extra node device? if so type 1');
    if isempty(en)==1
        en=0;
    elseif en~=1
        en=0;
    end
end
%if you copied and pasted the pulses as a new variable there needs to be at
%least one row of zeros at the end - so lets add one..
numpulses=length(pulses(:,1));
pulses(numpulses+1,:)=zeros(1,length(pulses(1,:)));
numpulses=numpulses+1;
i=2;
while i<=numpulses
    if pulses(i,1)==0
        i=i-1;
        skipto=pulses(i,1);
        break
    else
        i=i+1;
    end
end
[pulses,threshpts,fig_handles]=mznps_pickpulse_semiauto_nodownsamp_zone(ydetrend,yasls,yas2det,dmin,stderror,yas3det,plot_thresh,skipto,0,0,en,zones_of_interest);


%% run if your figures are already up and you want to start from the beginning 
use_old_fig=1;
skipto=0;
[pulses,threshpts,fig_handles]=mznps_pickpulse_semiauto_nodownsamp_zone(ydetrend,yasls,yas2det,dmin,stderror,yas3det,plot_thresh,skipto,use_old_fig,fig_handles,en,zones_of_interest);

%% run if your figures are already up and you want to start from where you left off
use_old_fig=1;
i=2;
while i<=length(pulses(:,1))
    if pulses(i,1)==0
        i=i-1;
        skipto=pulses(i,1);
        break
    else
        i=i+1;
    end
end

[pulses,threshpts,fig_handles]=mznps_pickpulse_semiauto_nodownsamp_zone(ydetrend,yasls,yas2det,dmin,stderror,yas3det,plot_thresh,skipto,use_old_fig,fig_handles,en,zones_of_interest);

%% run for REANALYSIS of old data, figures already up
use_old_fig=1;
if exist('dmin','var') ~=1
    dmin=input('dmin is not saved in your workspace, type what you want dmin to be right now');
end
if exist('plot_thresh','var') ~=1
    plot_thresh=0;
end
if exist('en','var')~=1
    en=input('is this an extra node device? if so type 1');
    if isempty(en)==1
        en=0;
    elseif en~=1
        en=0;
    end
end

skipto=0;
[pulses,threshpts,fig_handles]=mznps_pickpulse_semiauto_nodownsamp_zone(ydetrend,yasls,yas2det,dmin,stderror,yas3det,plot_thresh,skipto,use_old_fig,fig_handles,en,zones_of_interest,pulsesog);


%% run for REANALYSIS of old data, figures already up - START from where you left off 
use_old_fig=1;
if exist('dmin','var') ~=1
    dmin=input('dmin is not saved in your workspace, type what you want dmin to be right now');
end
if exist('plot_thresh','var') ~=1
    plot_thresh=0;
end
if exist('en','var')~=1
    en=input('is this an extra node device? if so type 1');
    if isempty(en)==1
        en=0;
    elseif en~=1
        en=0;
    end
end
i=2;
while i<=length(pulses(:,1))
    if pulses(i,1)==0
        i=i-1;
        skipto=pulses(i,1);
        break
    else
        i=i+1;
    end
end

[pulses,threshpts,fig_handles]=mznps_pickpulse_semiauto_nodownsamp_zone(ydetrend,yasls,yas2det,dmin,stderror,yas3det,plot_thresh,skipto,use_old_fig,fig_handles,en,zones_of_interest,pulsesog);




%% find fig handles if figures are up, but handles aren't saved
%this might happen if it was your first run of the session and it ended in an error 

fa=findobj('type','figure','Name','Data Set');
fb=findobj('type','figure','Name','Entire Event');
f4=findobj('type','figure','Name','Window of Interest');
fig_handles=[fa;fb;f4];
