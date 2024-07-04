
                    
                    
function   [line_handles, new_threshptsc]=redo_threshpts(newdmin,line_handles,axesbounds,Deff,L,yasls,xaxis,yas2det_for_thresh,p,axar,mininterval,linecolor,linew)                 
                    %inputs = newdmin, line_handles,axesbounds, Deff, L, yasls, p,
                    %axar, xaxis,yas2det_for_thresh,mininterval,
                    %linecolor,linew
                    
                    new_minexpected=(((newdmin^3)/((Deff^2)*L))*(1/(1-(0.8*((newdmin/Deff)^3)))))*yasls;
                    new_minexp_for_thresh=-1*new_minexpected;
                    %recalculate threshpts with "new_minexpected" as the
                    %threshold
                    new_threshpts=zeros(2,6);
                    new_threshptsc=zeros(2,6);
                    for j1=1:3
                         start_thresh=axesbounds(p,j1*2+1)-1;
                         end_thresh=axesbounds(p,j1*2+2);
                         if start_thresh<=1 || end_thresh<=1 
                             continue
                         end
                         new_threshptszone=zeros(2,2);
                         line_handles(1,j1)=plot(axar(j1,1),xaxis(start_thresh:end_thresh,1),new_minexpected(start_thresh:end_thresh,j1),'k');
                         new_threshptszone=classicthresholding(yas2det_for_thresh(start_thresh:end_thresh,j1),new_minexp_for_thresh(start_thresh:end_thresh,j1),new_threshptszone,2,1);
                         [m,~]=size(new_threshptszone);
                         new_threshptszonec=fixthreshptserrormz(new_threshptszone,mininterval,1,m);
                         [mc,~]=size(new_threshptszonec);
    %                      new_threshpts(1:m,j1*2-1)=new_threshptszone(1:m,1)+start_thresh-1;
    %                      new_threshpts(1:m,j1*2)=new_threshptszone(1:m,2);
                         if isempty(new_threshptszonec)==1 || new_threshptszonec(1,1)==0
                             'no threshpts found'
                         else
                             new_threshptsc(1:mc,j1*2-1)=new_threshptszonec(1:mc,1)+start_thresh-1;
                             new_threshptsc(1:mc,j1*2)=new_threshptszonec(1:mc,2);
                         end 
                         %now plot new thresh lines
                         for x=1:mc
                             line_handles(x+1,j1)=xline(axar(j1,1),new_threshptsc(x,j1*2-1),'Color',linecolor(j1,1),'LineWidth',linew);
                         end
                    end
                    
                    
                    
                