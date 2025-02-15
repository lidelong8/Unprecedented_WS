%%%%%%%%%%%%%%
clear;clc

[~,R]=geotiffread('D:\RecoverTime\dem_05dgr.tif');%%      
[row col]=latlon2pix(R,[435 35.83],[-0.847 -0.9]);% first lat then lon


col=round(abs ((-180 -(-0.847)) / .5) + 1) %lon 
row=round(abs ((90 - (35.83)) / .5) + 1) %lat
col=round(abs ((-180 -(3.0)) / .5) + 1) %lon 
row=round(abs ((90 - (36.5)) / .5) + 1)-1 %lat

%%%%%%%%%
ghms={'h08','cwatm'};
gcms={'gfdl-esm4','ipsl-cm6a-lr','mpi-esm1-2-hr','mri-esm2-0','ukesm1-0-ll'}; % no atotuse

names=cell(length(ghms)*length(gcms),1);
c=0;
for i=1:length(ghms)
for j=1:length(gcms)
%[k,i,j,]
c=c+1
names{c,1}=sprintf('%s_%s',ghms{i},gcms{j}); 
end
end


%%%%%%%%%%%%%%%%%%  
%
%i=1,j=1,k=1
hist_end = (2015-1850+1)*12;
future_start=hist_end+1;
future_end = (2100-1850+1)*12;
future_time = future_start:future_end; % 未来时期的年份


ave_window = 20;
freq_window = 20;
min_duration=6; 
min_month=1;
threshold=0.4;
max_gap=3;
threshold_f=1/freq_window;
%%%%%%%%%%%%%%%%%%  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%the tested WSI series
x=linspace(1850,2100,3012);
rand('seed',10);
y1=rand(3012,1);
y1=y1*0.5;
x0=linspace(0,1,3012)';
y2=-(x0-0.3).*(x0-0.35).*(x0-0.45).*(x0-0.5).*(x0-0.55)*10;
y3=0.30;
y=y1+y2+y3;
y=flipud(y);

y_ori=y;
y=movmean(y,[20-1 0],1,"omitnan"); 
y=y';

%%%%%%%frequency
wsi_ori=y_ori';
wss=wsi_ori;%water shortage state, first cal ws;culate to save time!
wss(wss <= threshold) = 0; % set non-water scarcity to 0
wss(wss >  threshold) = 1; % set water scarcity to 1
wss_ori=wss;
%%%%cal ES frequency,here do not consider the short periods of non-WS months between lengthy periods of WS months
temp = wss; % whole period
transitions = diff([temp,temp(:,end)],1,2); %diff(X,n,dim)
transitions(transitions~=1)=0;
freq=filter(ones(1,freq_window)/freq_window,1,transitions,[],2);
freq=movmean(freq,[freq_window-1 0],2,"omitnan"); %window:2+0+1
freq(1:ave_window-1)=nan;
%%%%%%%frequency
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%the tested WSI series


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
figure('unit','centimeters','Position',[2 2 30 25]);

%%%%%%%%%%plot1
subplot_tight(4,1,1,[0.04,0.04]);

%xi = x(1):1:x(end);
%yi = interp1(x,y,xi,'linear'); %only show 1850-2100
xi = x;
yi =y;

hold on
h2=area(xi,bsxfun(@max, yi, max(y(1:hist_end))),max(y(1:hist_end)),'LineStyle','none','ShowBaseLine','off','FaceColor','r','EdgeColor','none','FaceAlpha',0.4);%
h4=plot(x(1:(2014-1850)*12+12),y(1:(2014-1850)*12+12),'Linestyle','-','Marker','none','color','k','Linewidth',1);
hold on
h5=plot(x((2015-1850)*12+1:(2100-1850)*12+12),y((2015-1850)*12+1:end),'Linestyle','-','Marker','none','color','b','Linewidth',1);
%%%%%%%%%%
hold on;
xline(2014,'k--','2014','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
hold on;
yline(max(y(1:hist_end)),'r-',['Maximum intensity: ',num2str(max(y(1:hist_end)),'%.2f')],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'center','Fontsize',12,'fontweight','bold');
hold on;
yline(0.4,'r--','water scarcity threshold: 0.40','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'center','Fontsize',12,'fontweight','bold');
fill([x(1)+0.1,x(1)+0.1,x(end)-0.1,x(end)-0.1],[0,max(y(1:hist_end)),max(y(1:hist_end)),0],[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.4);
%
set(gca,'Xlim',[1850 2100],'xtick',[1850:50:2100],'xticklabel',[])
set(gca,'ylim',[0,1.0],'ytick',0:0.2:1,'yticklabel',{'0.0','0.2','0.4','0.6','0.8','1.0'})%
set(gca,'XMinorTick','on')%
set(gca,'YMinorTick','on')%

ax = gca;
ax.XAxis.TickDirection='out';
ax.YAxis.TickDirection='in';
ax.XAxis.TickLength = [0.005 0.015]; %
%ax.YAxis.TickLength = [0.005 0.015]; % 

ylabel ('Intensity','Fontsize',12,'fontweight','bold'); 
%xlabel ('Year','Fontsize',12,'fontweight','bold');
set(gca,'FontSize',12,'FontName','Times New Roman','fontweight','bold')
%%%%%processing...
[max_i_h,max_i_f,max_d_h,max_d_f,max_f_h,max_f_f,TFE_i,TFE_d,TFE_f,TFE_w] = cal_TFE_mon_demo(y_ori',ave_window,freq_window,threshold,threshold_f,min_duration,max_gap,min_month,hist_end)

hold on;
%xline(TFE_i,'k-',['TFE_i: ',num2str(TFE_i+1),'   (>=5 yr)'],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment', 'left','Fontsize',12,'fontweight','bold'); 
xline(TFE_i,'b-','','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment', 'left','Fontsize',12,'fontweight','bold'); 
%text('Units','normalized','Position',[0.70, 0.90],'String',['TFE_i: ',num2str(TFE_i),' (WSI > ',num2str(max_i_h,'%.2f'),')'],'rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');
text('Units','normalized','Position',[0.54, 0.20],'String',{['TFE_i: ',num2str(TFE_i),' (WSI > ',num2str(max_i_h,'%.2f'),') and Kolmogorov-Smirnov test'];' between X1 and X2 shows significant differences'},'rotation',0,'Color','b','FontName','Arial','FontSize',12,'fontweight','bold');


%
temp = y(1,1:hist_end); % only account for historical period
%%historical maximum 
[thres_i,loc_i_h]=nanmax(temp,[],2);
base_series=temp(1,loc_i_h-ave_window+1:loc_i_h);
base_x=x((loc_i_h-ave_window+1):loc_i_h);
%%     
temp = y(future_time); % 
thres_i=max(y(1:hist_end));
%
temp(temp <= thres_i) = 0; % 
temp(temp >  thres_i) = 1; % 
%
run_lengths = diff([0, temp, 0]); % 
start_index = find(run_lengths == 1); % 
end_index = find(run_lengths == -1) - 1; % 
durations = end_index - start_index + 1; % 
start_x=x(hist_end+start_index);
end_x=x(hist_end+end_index);      

      
  
%%%%%processing...
%% Add labels
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
%x1=0.02+(ToE-xlim(1))/diff(xlim);
%x2=0.02+(ToD-xlim(1))/diff(xlim);
%
%annotation('doublearrow',[0.02 (2014-1850)/(2100-1850)],[0.65 0.65]); 
%text(1900,1.4,'The historical period','rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');
%text(TFE_w,1.6,'(>=5 yr)','rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');
%
ct=1;
for j=1:length(durations)
if durations(j)>=6&ct==1;
line([start_x(j) end_x(j)+0],[thres_i-0.10 thres_i-0.10],'linestyle','-','color','r','Linewidth',1.0)
line([start_x(j) start_x(j)],[thres_i thres_i-0.10],'linestyle','-','color','r','Linewidth',1.0)
line([end_x(j)+0 end_x(j)+0],[thres_i thres_i-0.10],'linestyle','-','color','r','Linewidth',1.0)
text((start_x(j)+end_x(j))/2,thres_i-0.20,num2str(durations(j)),'rotation',0,'Color','r','FontName','Arial','FontSize',12,'fontweight','bold','HorizontalAlignment', 'center');
ct=ct+1;
future_series=y(1,start_x(j)-ave_window+1:start_x(j));
future_x=x((hist_end+start_index(j)-ave_window+1):(hist_end+start_index(j)));
end
end

hold on
fill([base_x(1),base_x(1),base_x(end),base_x(end)],[min(ylim),max(base_series),max(base_series),min(ylim)],[0.2,0.2,0.2],'EdgeColor','none','FaceAlpha',0.3);
%fill([base_x(1),base_x(1),base_x(end),base_x(end)],[min(ylim),max(ylim),max(ylim),min(ylim)],[0.2,0.2,0.2],'EdgeColor','none','FaceAlpha',0.3);
hold on
%fill([future_x(1),future_x(1),future_x(end),future_x(end)],[min(future_series),max(future_series),max(future_series),min(future_series)],[0.2,0.2,0.2],'EdgeColor','none','FaceAlpha',0.3);
fill([future_x(1),future_x(1),future_x(end),future_x(end)],[min(ylim),max(future_series),max(future_series),min(ylim)],[0.2,0.2,0.2],'EdgeColor','none','FaceAlpha',0.3);

hold on
plot(base_x,base_series,'Linestyle','-','Marker','none','color','c','Linewidth',2);
hold on
plot(future_x,future_series,'Linestyle','-','Marker','none','color','y','Linewidth',2);

text((base_x(1)+base_x(end))/2,max(base_series)+0.1,'x1','FontSize',12,'fontweight','bold','color','b','HorizontalAlignment', 'center');
text((future_x(1)+future_x(end))/2,max(future_series)+0.1,'x2','FontSize',12,'fontweight','bold','color','b','HorizontalAlignment', 'center');



lgd=legend([h2,h4,h5,],{'Exceed the threshold','WSI (1850-2014)','WSI (2015-2100)'},'Location','northwest','FontSize',14,'FontWeight','bold','NumColumnsMode','manual','NumColumns',3);
lgd.Box='off';
%%
box on
%%%
box off
ax2 = axes('Position',get(gca,'Position'),...
                 'XAxisLocation','top',...
                 'YAxisLocation','right',...
                 'Color','none',...
                 'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on
%%%
text('Units','normalized','Position',[0.005, 1.09],'String','(a)','rotation',0,'Color','k','FontName','Arial','FontSize',16,'fontweight','bold');
%
%title(['Representative Grid (36.5°N, 3.0°E) in {\color{blue}watergap2-2e' '\color[rgb]{0 0 0}\_' '\color[rgb]{1 0 0}mpi-esm1-2-hr}'],'FontWeight','bold','Fontsize',14,'FontName','Arial')
title('Unprecedented Intensity of Water Scarcity','FontName','Arial','FontSize',14,'fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold','FontName','Times new Roman');  
%%%%%%%%%%plot1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%plot2
subplot_tight(4,1,2,[0.04,0.04]);
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% duration
wsi=y_ori';
wss=wsi;%water shortage state, first cal ws;culate to save time!
wss(wss <= threshold) = 0; % set non-water scarcity to 0
wss(wss >  threshold) = 1; % set water scarcity to 1
wss_ori=wss;

temp = wss_ori;% 
%%%%
%%%%%%%%%% cal ES frequency,here do not consider the short periods of non-WS months between lengthy periods of WS months
  %%%%
%%1: the maximum gap of Non-WS duration between two WS evetns
run_lengths1 = diff([1, temp, 1]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
start_index1 = find(run_lengths1 == -1); % 
end_index1 = find(run_lengths1 == 1) - 1; 
durations1 = end_index1 - start_index1 + 1; 
%find origional data
id_min_dur1=find(durations1<max_gap);
widths1=durations1(id_min_dur1);
indices1=start_index1(id_min_dur1);
%%
id_min1 = [];
for k = 1:length(indices1)
   id_min1 = [id_min1, indices1(k):indices1(k)+widths1(k)-1];
end
%%Repalce origional data
temp(id_min1)= ~temp(id_min1);
t=wsi(id_min1);
wsi(id_min1)=0.4+(0.4-t);
%% The short periods of non-WS months between lengthy periods of WS months were considered WS due to the pooling effect.
%%%%%%
%%%%%%
%%2: the minimum duration of WS between two Non-WS evetns
run_lengths = diff([0, temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
start_index = find(run_lengths == 1); % 
end_index = find(run_lengths == -1) - 1; 
durations = end_index - start_index + 1; 
%find origional data
id_min_dur=find(durations<min_duration);
widths=durations(id_min_dur);
indices=start_index(id_min_dur);
%%
id_min = [];
for k = 1:length(indices)
   id_min = [id_min, indices(k):indices(k)+widths(k)-1];
end
%%Repalce origional data
temp(id_min)= ~temp(id_min);
t=wsi(id_min);
wsi(id_min)=0.4+(0.4-t);
%% The duration of a negligible WS period shorter than x-months
%%%%%%
%%%%%%   
%%3: detect the WS evetns ignoring the min non-WS events or min WS events
%% In whole peried, find if the year 2015 is within the WS duaration and dates
run_lengths = diff([0, temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
start_index = find(run_lengths == 1); % 
end_index = find(run_lengths == -1) - 1; 
durations = end_index - start_index + 1; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           temp_h=temp(1,1:hist_end);% only account for historical period
           run_lengths = diff([0, temp_h, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
           start_index = find(run_lengths == 1); % 
           end_index = find(run_lengths == -1) - 1; 
           durations = end_index - start_index + 1; 
           start_x1= x(start_index);
           end_x1= x(end_index);
           durations1=durations;
           [max_d_h,id_max_dur_h] = nanmax(durations,[],2); % find historical maximum duration    

           temp_f=temp(1,future_time);% only account for future period
           run_lengths = diff([0, temp_f, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
           start_index = find(run_lengths == 1); % 
           end_index = find(run_lengths == -1) - 1; 
           durations = end_index - start_index + 1;  
           durations2=durations; 
           start_x2=x(hist_end+start_index);
           end_x2=x(hist_end+end_index);               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xi = 1850:0.1:2100;
yi = interp1(x,wsi,xi,'linear'); %only show 1850-2100
h1=area(x,bsxfun(@max, wsi, 0.4),0.4,'LineStyle','none','ShowBaseLine','off','FaceColor','r','EdgeColor','none','FaceAlpha',0.4);%
hold on
h2=area(x,bsxfun(@min, wsi, 0.4),0.4,'LineStyle','none','ShowBaseLine','off','FaceColor','b','EdgeColor','none','FaceAlpha',0.4);%
hold on
xline(2014,'k--','2014','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
hold on
yline(0.4)


%%%%%%%%%%
set(gca,'Xlim',[1850 2100],'xtick',[1850:50:2100],'xticklabel',[])
set(gca,'ylim',[0,1.0],'ytick',0:0.2:1,'yticklabel',{'0.0','0.2','0.4','0.6','0.8','1.0'})%
set(gca,'XMinorTick','on')%
set(gca,'YMinorTick','on')%
% 
ax = gca;
ax.XAxis.TickDirection='out';
% 
ax.YAxis.TickDirection='in';
ax.XAxis.TickLength = [0.005 0.015]; %  
% 
%ax.YAxis.TickLength = [0.005 0.015]; %  

ylabel ('Intensity','Fontsize',12,'fontweight','bold'); 
%xlabel ('Year','Fontsize',12,'fontweight','bold');
set(gca,'FontSize',12,'FontName','Times New Roman','fontweight','bold')
%%%%%processing...


j=id_max_dur_h(1);
fill([start_x1(j)-0.1,start_x1(j)-0.1,end_x1(j)+0.1,end_x1(j)+0.1],[0,ylim(2),ylim(2),0],[0.4 0.4 0.4],'EdgeColor','none','FaceAlpha',0.4);
j=id_max_dur_h(1);
%xline((start_x1(j)+end_x1(j))/2,'w--',['Past maximum duration: ',num2str(max_d_h+0,'%.0f'),' months'],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
%annotation('textarrow',[(start_x1(j)-1850)/(2100-1850+1)+0.2,(end_x1(j)-1850)/(2100-1850+1)+0.2],[0.60 0.60],'String',['Past maximum duration: ',num2str(max_d_h+0,'%.0f'),' months'],'FontSize',12,'fontweight','bold','color','r');
text((start_x1(j)+end_x1(j))/2,0.80,['Past maximum duration: ',num2str(max_d_h+0,'%.0f'),' months'],'FontSize',12,'fontweight','bold','color','r','HorizontalAlignment', 'center');

%annotation('doublearrow',[0.02 (2014-1850)/(2100-1850)],[0.4 0.4]); %
%text(1900,1.4,'The historical period','rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');
%text(TFE_w,1.6,'(>=5 yr)','rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');
%
thres=0.4;
ct=1;
for j=1:length(durations1)
if mod(j,2)==1
line([start_x1(j) end_x1(j)+0],[thres-0.10 thres-0.10],'linestyle','-','color','r','Linewidth',1.0)
line([start_x1(j) start_x1(j)],[thres thres-0.10],'linestyle','-','color','r','Linewidth',1.0)
line([end_x1(j)+0 end_x1(j)+0],[thres thres-0.10],'linestyle','-','color','r','Linewidth',1.0)
text((start_x1(j)+end_x1(j))/2,thres-0.20,num2str(durations1(j)),'rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold','HorizontalAlignment', 'center');
else
line([start_x1(j) end_x1(j)+0],[thres+0.10 thres+0.10],'linestyle','-','color','r','Linewidth',1.05)
line([start_x1(j) start_x1(j)],[thres thres+0.10],'linestyle','-','color','r','Linewidth',1.0)
line([end_x1(j)+0 end_x1(j)+0],[thres thres+0.10],'linestyle','-','color','r','Linewidth',1.0)
text((start_x1(j)+end_x1(j))/2,thres+0.20,num2str(durations1(j)),'rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold','HorizontalAlignment', 'center');
end
end

for j=1:length(durations2)
  if mod(j,2)==1
    line([start_x2(j) end_x2(j)+0],[thres-0.10 thres-0.10],'linestyle','-','color','r','Linewidth',1.0)
    line([start_x2(j) start_x2(j)],[thres thres-0.10],'linestyle','-','color','r','Linewidth',1.0)
    line([end_x2(j)+0 end_x2(j)+0],[thres thres-0.10],'linestyle','-','color','r','Linewidth',1.0)
    if durations2(j)>max_d_h
        text((start_x2(j)+end_x2(j))/2,thres-0.20,num2str(durations2(j)),'rotation',0,'Color','r','FontName','Arial','FontSize',12,'fontweight','bold','HorizontalAlignment', 'center');
    else
        text((start_x2(j)+end_x2(j))/2,thres-0.20,num2str(durations2(j)),'rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold','HorizontalAlignment', 'center');
    end
  else
    line([start_x2(j) end_x2(j)+0],[thres+0.10 thres+0.10],'linestyle','-','color','r','Linewidth',1.0)
    line([start_x2(j) start_x2(j)],[thres thres+0.10],'linestyle','-','color','r','Linewidth',1.0)
    line([end_x2(j)+0 end_x2(j)+0],[thres thres+0.10],'linestyle','-','color','r','Linewidth',1.0)
    if durations2(j)>max_d_h
        text((start_x2(j)+end_x2(j))/2,thres+0.20,num2str(durations2(j)),'rotation',0,'Color','r','FontName','Arial','FontSize',12,'fontweight','bold','HorizontalAlignment', 'center');
    else
        text((start_x2(j)+end_x2(j))/2,thres+0.20,num2str(durations2(j)),'rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold','HorizontalAlignment', 'center');
    end
  end %if mod
end %end for
box on

lgd=legend([h1,h2,],{'WS','Non-WS'},'Location','northwest','FontSize',14,'FontWeight','bold','NumColumnsMode','manual','NumColumns',2);
lgd.Box='off';

%%%
box off
ax2 = axes('Position',get(gca,'Position'),...
                 'XAxisLocation','top',...
                 'YAxisLocation','right',...
                 'Color','none',...
                 'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on
%%%

%
hold on;
if TFE_d<2100
xline(TFE_d,'b-',['TFE_d: ',num2str(TFE_d),' (Duration > ',num2str(max_d_h),' yr)'],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','left','Fontsize',12,'fontweight','bold'); 
else
text('Units','normalized','Position',[0.70, 0.90],'String',['No  TFE_d: (Max Duration_f < ',num2str(max_d_h),' months)'],'rotation',0,'Color','b','FontName','Arial','FontSize',12,'fontweight','bold');
end

%% Add labels
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
%x1=0.02+(ToE-xlim(1))/diff(xlim);
%x2=0.02+(ToD-xlim(1))/diff(xlim);

text('Units','normalized','Position',[0.005, 1.09],'String','(b)','rotation',0,'Color','k','FontName','Arial','FontSize',16,'fontweight','bold');
title('Unprecedented Duration of Water Scarcity','FontName','Arial','FontSize',14,'fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold','FontName','Times new Roman');  
%%%%%%%%%%plot2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%plot3
subplot_tight(4,1,3,[0.04,0.04]);
%%%%
xi = x(1):0.01:x(end);
yi = interp1(x,freq,xi,'linear'); %only show 1850-2100
hold on
h2=area(xi,bsxfun(@max, yi, max(freq(1:hist_end))),max(freq(1:hist_end)),'LineStyle','none','ShowBaseLine','off','FaceColor','r','EdgeColor','none','FaceAlpha',0.4);%
h4=plot(x(1:(2014-1850)*12+12),freq(1:(2014-1850)*12+12),'Linestyle','-','Marker','none','color','k','Linewidth',1);
hold on
h5=plot(x((2015-1850)*12+1:(2100-1850)*12+12),freq((2015-1850)*12+1:end),'Linestyle','-','Marker','none','color','b','Linewidth',1);
%%%%%%%%%%
hold on;
xline(2014,'k--','2014','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
hold on;
yline(max(freq(1:hist_end)),'r-',['Maximum intensity: ',num2str(max(freq(1:hist_end)),'%.2f')],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'center','Fontsize',12,'fontweight','bold');
hold on;
%yline(0.4,'r--','water scarcity threshold: 0.40','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'center','Fontsize',12,'fontweight','bold');
fill([x(1)+0.1,x(1)+0.1,x(end)-0.1,x(end)-0.1],[0,max(freq(1:hist_end)),max(freq(1:hist_end)),0],[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.4);
%
set(gca,'Xlim',[1850 2100],'xtick',[1850:50:2100],'xticklabel',[])
set(gca,'ylim',[0,0.5],'ytick',0:0.1:0.5,'yticklabel',{'0.0','0.1','0.2','0.3','0.4','0.5'})%
%set(gca,'ylim',[0,0.5],'ytick',[0,0.5])%
set(gca,'XMinorTick','on')%
set(gca,'YMinorTick','on')%
% 
ax = gca;
ax.XAxis.TickDirection='out';
% 
ax.YAxis.TickDirection='in';
ax.XAxis.TickLength = [0.005 0.015]; %  
% 
%ax.YAxis.TickLength = [0.005 0.015]; %  

ylabel ('Frequency','Fontsize',12,'fontweight','bold'); 
%xlabel ('Year','Fontsize',12,'fontweight','bold');
set(gca,'FontSize',12,'FontName','Times New Roman','fontweight','bold')
%%%%%processing...
[max_i_h,max_i_f,max_d_h,max_d_f,max_f_h,max_f_f,TFE_i,TFE_d,TFE_f,TFE_w] = cal_TFE_mon_demo(y_ori',ave_window,freq_window,threshold,threshold_f,min_duration,max_gap,min_month,hist_end)

hold on;
%xline(TFE_f+1,'k-',['TFE_f: ',num2str(TFE_f+1),'   (>=5 yr)'],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment', 'left','Fontsize',12,'fontweight','bold'); 
if TFE_f<2100
%xline(TFE_f,'k-',['TFE_f: ',num2str(TFE_f+1),'   (>=5 yr)'],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment', 'left','Fontsize',12,'fontweight','bold'); 
text('Units','normalized','Position',[0.70, 0.90],'String',['TFE_f: ',num2str(TFE_f),' (Freq > ',num2str(max_f_h,'%.2f'),')'],'rotation',0,'Color','b','FontName','Arial','FontSize',12,'fontweight','bold');
else
text('Units','normalized','Position',[0.70, 0.90],'String',['No TFE_f: (Max Freq_f < ',num2str(max_f_h,'%.2f'),')'],'rotation',0,'Color','b','FontName','Arial','FontSize',12,'fontweight','bold');
end
text('Units','normalized','Position',[0.70, 0.78],'String','No X2, no Kolmogorov-Smirnov test','rotation',0,'Color','b','FontName','Arial','FontSize',12,'fontweight','bold');
%
temp = freq(1,1:hist_end); % only account for historical period
%%historical maximum 
[thres_f,loc_f_h]=nanmax(temp,[],2);
base_series=temp(1,loc_f_h-ave_window+1:loc_f_h);
base_x=x((loc_f_h-ave_window+1):loc_f_h);
        
temp = freq(future_time); % 
thres_i=max(freq(1:hist_end));
%
temp(temp <= thres_i) = 0; % 
temp(temp >  thres_i) = 1; % 
%
run_lengths = diff([0, temp, 0]); % 
start_index = find(run_lengths == 1); % 
end_index = find(run_lengths == -1) - 1; % 
durations = end_index - start_index + 1; % 
start_x=x(hist_end+start_index);
end_x=x(hist_end+end_index);        
%%%%%processing...
%% Add labels
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
%x1=0.02+(ToE-xlim(1))/diff(xlim);
%x2=0.02+(ToD-xlim(1))/diff(xlim);
%
%annotation('doublearrow',[0.02 (2014-1850)/(2100-1850)],[0.65 0.65]); %
%text(1900,1.4,'The historical period','rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');
%text(TFE_w,1.6,'(>=5 yr)','rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');
%
ct=1;
future_x=[];
future_series=[];
for j=1:length(durations)
if durations(j)>=6&ct==1;
line([start_x(j) end_x(j)+1],[thres_i-0.10 thres_i-0.10],'linestyle','-','color','r','Linewidth',2.0)
line([start_x(j) start_x(j)],[thres_i thres_i-0.10],'linestyle','-','color','r','Linewidth',2.0)
line([end_x(j)+1 end_x(j)+1],[thres_i thres_i-0.10],'linestyle','-','color','r','Linewidth',2.0)
text((start_x(j)+end_x(j))/2,thres_i-0.20,num2str(durations(j)),'rotation',0,'Color','r','FontName','Arial','FontSize',12,'fontweight','bold');
ct=ct+1;
future_series=y(1,start_x(j)-ave_window+1:start_x(j));
future_x=x((hist_end+start_index(j)-ave_window+1):(hist_end+start_index(j)));
end
end

hold on
fill([base_x(1),base_x(1),base_x(end),base_x(end)],[min(ylim),max(base_series),max(base_series),min(ylim)],[0.2,0.2,0.2],'EdgeColor','none','FaceAlpha',0.3);
hold on
if future_x
fill([future_x(1),future_x(1),future_x(end),future_x(end)],[min(ylim),max(future_series),max(future_series),min(ylim)],[0.2,0.2,0.2],'EdgeColor','none','FaceAlpha',0.3);
end

hold on
plot(base_x,base_series,'Linestyle','-','Marker','none','color','c','Linewidth',2);
hold on
if future_x
plot(future_x,future_series,'Linestyle','-','Marker','none','color','y','Linewidth',2);
end

text((base_x(1)+base_x(end))/2,max(base_series)+0.05,'x1','FontSize',12,'fontweight','bold','color','b','HorizontalAlignment', 'center');
if future_x
text((future_x(1)+future_x(end))/2,max(future_series)+0.02,'x2','FontSize',12,'fontweight','bold','color','b','HorizontalAlignment', 'center');
end


lgd=legend([h2,h4,h5,],{'Exceed the threshold','Freq (1850-2014)','Freq (2015-2100)'},'Location','northwest','FontSize',14,'FontWeight','bold','NumColumnsMode','manual','NumColumns',3);
lgd.Box='off';
%%
box on
%%%
box off
ax2 = axes('Position',get(gca,'Position'),...
                 'XAxisLocation','top',...
                 'YAxisLocation','right',...
                 'Color','none',...
                 'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on
%%%
text('Units','normalized','Position',[0.005, 1.09],'String','(c)','rotation',0,'Color','k','FontName','Arial','FontSize',16,'fontweight','bold');

title('Unprecedented Frequency of Water Scarcity','FontName','Arial','FontSize',14,'fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold','FontName','Times new Roman');  
%%%%%%%%%%plot3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%plot4  
subplot_tight(4,1,4,[0.06,0.04]);        
%%
xi = x;
yi = y; %only show 1850-2100
hold on
h2=area(xi,bsxfun(@max, yi, 0.4),0.4,'LineStyle','none','ShowBaseLine','off','FaceColor','r','EdgeColor','none','FaceAlpha',0.4);%
h4=plot(x(1:hist_end),y(1:hist_end),'Linestyle','-','Marker','none','color','k','Linewidth',1);
hold on
h5=plot(x(hist_end+1:end),y(hist_end+1:end),'Linestyle','-','Marker','none','color','b','Linewidth',1);
%%%%%%%%%%
hold on;
xline(2014,'k--','2014','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
hold on;
yline(0.4,'r--','water scarcity threshold: 0.40','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
%fill([x(2),x(2),x(end-1),x(end-1)],[0,0.4,0.4,0],[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.4);
%
set(gca,'Xlim',[1850 2100],'xtick',[1850:50:2100],'xticklabel',[1850:50:2100])
set(gca,'ylim',[0,1.0],'ytick',0:0.2:1,'yticklabel',{'0.0','0.2','0.4','0.6','0.8','1.0'})%
set(gca,'XMinorTick','on')%
set(gca,'YMinorTick','on')%
%set(gca,'TickDir','out')%
%set(gca,'TickLength',[0.01 0.01])
% 
ax = gca;
ax.XAxis.TickDirection='out';
% 
ax.YAxis.TickDirection='in';
ax.XAxis.TickLength = [0.005 0.015]; %  
% 
%ax.YAxis.TickLength = [0.005 0.015]; %  



ylabel ('Intensity','Fontsize',12,'fontweight','bold'); 
xlabel('Year','Fontsize',12,'fontweight','bold');
pos=get(get(gca,'xlabel'),'Position');
pos(2)=-0.2;
set(get(gca,'xlabel'),'Position',pos);
set(gca,'FontSize',12,'FontName','Times New Roman','fontweight','bold')
%%%%%%%%%%%%
%%%%%processing...
[n_i,n_d,n_f,nt_i,nt_d,nt_f,max_i,max_d,max_f,max_i_h,max_i_f,max_d_h,max_d_f,max_f_h,max_f_f,TFE_i,TFE_d,TFE_f,TFE_w,i_f,d_f,f_f] = calculate_TFE_yr(wsi,freq,threshold,threshold_f,min_duration,hist_end);

hold on;
%xline(TFE_w,'k-',num2str(TFE_w),'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'center','Fontsize',12,'fontweight','bold'); 
xline(TFE_w,'b-','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'center','Fontsize',12,'fontweight','bold'); 
text('Units','normalized','Position',[0.70, 0.90],'String',['TFE_{w}: ',num2str(TFE_w+0),' (WSI > ',num2str(0.4,'%.2f'),')'],'rotation',0,'Color','b','FontName','Arial','FontSize',12,'fontweight','bold');
%% Add labels
xlim = get(gca,'XLim');
ylim = get(gca,'YLim');
%x1=0.02+(ToE-xlim(1))/diff(xlim);
%x2=0.02+(ToD-xlim(1))/diff(xlim);
%
%text(TFE_w,1.6,'(>=5 yr)','rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');
annotation('doublearrow',[0.04 (2014-1850)/(2100-1850+1)-0.01],[0.18 0.18]); %
%text(1900,0.97,'The historical period','rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');
text('Units','normalized','Position',[0.31, 0.65],'String','The historical period','rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');
%
thres=0.4;
if 0
for j=1:length(durations2)
line([start_x2(j) end_x2(j)+1],[thres-0.10 thres-0.10],'linestyle','-','color','r','Linewidth',2.0)
line([start_x2(j) start_x2(j)],[thres thres-0.10],'linestyle','-','color','r','Linewidth',2.0)
line([end_x2(j)+1 end_x2(j)+1],[thres thres-0.10],'linestyle','-','color','r','Linewidth',2.0)
text((start_x2(j)+end_x2(j))/2,thres-0.20,num2str(durations2(j)),'rotation',0,'Color','r','FontName','Arial','FontSize',12,'fontweight','bold');
end
end
%%%
box off
ax2 = axes('Position',get(gca,'Position'),...
                 'XAxisLocation','top',...
                 'YAxisLocation','right',...
                 'Color','none',...
                 'XColor','k','YColor','k');
set(ax2,'YTick', []);
set(ax2,'XTick', []);
box on
%%%
%%%%%%
text('Units','normalized','Position',[0.005, 1.09],'String','(d)','rotation',0,'Color','k','FontName','Arial','FontSize',16,'fontweight','bold');
%title('Representative Grid (36.5°N, 3.0°E) in atergap2-2e\_mpi-esm1-2-hr','FontWeight','bold','Fontsize',14,'FontName','Arial');
%title(['\fontsize{14}Representative Grid (36.5°N, 3.0°E) in {\color{blue}watergap2-2e ' '\_|' '\color[rgb]{0 .5 .5 }mpi-esm1-2-hr}'])
%
%text('Units','normalized','Position',[0.05, 0.50],'String',['Events: ',num2str(length(durations2))],'rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');
%text('Units','normalized','Position',[0.15, 0.50],'String',['Maximum Duration: ',num2str(max(durations2)),' yr'],'rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');
%text('Units','normalized','Position',[0.32, 0.50],'String',['Cumulative occurring time: ',num2str(sum(durations2)),' yr'],'rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');
lgd=legend([h2,h4,h5,],{'Exceed the threshold','WSI (1850-2014)','WSI (2015-2100)'},'Location','northwest','FontSize',14,'FontWeight','bold','NumColumnsMode','manual','NumColumns',3);
lgd.Box='off';

title('Conventional Water Scarcity','FontName','Arial','FontSize',14,'fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold','FontName','Times new Roman');  
%%%%%%%%%%plot4



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save pic
%exportgraphics(gcf,'D:\ws_fluctuation\figures\Fig.s1.tif','Resolution',350);% no white 
exportgraphics(gcf,'D:\ws_fluctuation\figures\Fig.s1.jpg','Resolution',350);% no white 
%%%%%%%%%%%%%%%%%%%%%