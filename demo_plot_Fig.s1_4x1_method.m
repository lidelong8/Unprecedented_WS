%%%%%%%%%%%%%%
%% the reference time is fixed in the year of 2014.
% wsi: monthly water shortage indices (1850-2100), with dimensions 360x720x3012
% pop: yearly population in each grid from 2015-2100
% areas:Areas of each grid
% ave_window: smooth wiondows
% freq_window: wiondows to calculate frequency
% threshold: 0.4
% threshold_d: if historical duration is less than 12 months, in the future, WS with at least 12 months will be identified as UDWS.
% threshold_f: at least 1 times in 60 months.
% min_month: UWS last at least for 5 months.
% min_duration: minimum duration of consecutive exceedance in months (i.e., 5 months)
% max_gap: the non-WS events are seen as WS shorter than three months
% pre_industria: pre-industria periods 1850-1900
% hist_end: index of last month of historical period (must be less than or equal to 2014)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
i=1

ave_window=60;
freq_window=60;
threshold=0.4;
threshold_d=12;
threshold_f=1/60;
min_duration=5;
max_gap=3;
min_month=5;
pre_industrial=(1900-1850)*12;
hist_end=(2014-1850+1)*12;
future_start=hist_end+1;
future_end = (2100-1850+1)*12;
future_time = future_start:future_end; % 
%%%%%%%%%%%%%%%%%%




%%%%%%%%%
rng(123); % seed 123
wsi=rand(1,3012)+0.2;
wsi(1,1:pre_industrial)=wsi(1,1:pre_industrial)*0.8;
wsi(1,2000:end)=wsi(1,2000:end)*1.2;
wsi=wsi-0.32;

wsi_ori=wsi;
wss=wsi;%water shortage state, first cal ws;culate to save time!
wss(wss <= threshold) = 0; % set non-water scarcity to 0
wss(wss >  threshold) = 1; % set water scarcity to 1
wss_ori=wss;
%%smooth
wsi_s=movmean(wsi,[ave_window-1 0],2,"omitnan"); %window:2+0+1
wsi_s(wsi_s<=0)=0;
wss=wsi_s;
wss(wss <= threshold) = 0; % set non-water scarcity to 0
wss(wss >  threshold) = 1; % set water scarcity to 1
wss_s=wss;
%%%%%%%%%smooth
y_ori=wsi_ori;
y=wsi_s;
y(1:ave_window)=nan;%;%set the first ave_window-1 value to nan
x=1850:2100;
x=linspace(1850,2100,3012);
%%%%%%%%%%%%%%

% Calculate regional average water scarcity indices
[~,future_end]=size(wsi);
future_start = hist_end+1;
future_time = future_start:future_end; %
all_time=1:future_end;
yr=length(future_time)/12;
A=1;
P=ones(1,86);%2025-2100


A_i = zeros(1,yr); %
A_d = zeros(1,yr); %
A_f = zeros(1,yr); %
A_w = zeros(1,yr); %

P_i = zeros(1,yr); %
P_d = zeros(1,yr); %
P_f = zeros(1,yr); %
P_w = zeros(1,yr); %

mon_i = zeros(1,yr); %
mon_d = zeros(1,yr); %
mon_f = zeros(1,yr); %
mon_w = zeros(1,yr); %

max_i_h= zeros(1,1); %
max_d_h = zeros(1,1); %
max_f_h = zeros(1,1); %

max_i_f= zeros(1,1); %
max_d_f = zeros(1,1); %
max_f_f = zeros(1,1); %

TFE_i= 2101*ones(1,1);
TFE_d= 2101*ones(1,1);
TFE_f= 2101*ones(1,1);
TFE_w= 2101*ones(1,1);
%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 1;
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TFE_i=2101;
TFE_d=2101;
TFE_f=2101;
TFE_w=2101;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%cal Traditional WS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%cal Traditional WS
temp = wss_s(i,future_time); %  future time

ws=squeeze(nansum(reshape(temp,12,[]),1));%save
%%
tmp=ws;
tmp(tmp<1)=0;
tmp(tmp>=1)=1;
%%
mon_w(i,:)=ws;
A_w(i,:)=A*tmp;
P_w(i,:)=P(i,:).*tmp;
%%
idw=find(ws>=1);
if idw
    TFE_w=2015+idw(1)-1;%Time of first emergence
end
%%%%%%%%%%%%%
%%Traditional WS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%cal Unprecedented intensity
%% intensity
%%%%%%%%%%%%%%%his
temp = wsi_s(i,1:hist_end); % only account for historical period
x_h=x(i,1:hist_end);
Noise_i=nanstd(temp(1:pre_industrial),0,2);
%%historical maximum
[thres_i,loc_i_h]=nanmax(temp,[],2);
thres_i=thres_i+1*Noise_i;% important
max_i_h=thres_i;

if thres_i<threshold
    thres_i=single(threshold);%if there is water scarcity (<0.4), set 0
end
%%historical maximum base period
if loc_i_h>ave_window
    base_series=temp(1,loc_i_h-ave_window+1:loc_i_h);
    base_x=x_h(1,loc_i_h-ave_window+1:loc_i_h);
else
    base_x=x_h(1,1:loc_i_h);
    base_series=temp(1,1:loc_i_h);
end
%% intensity

%%TFE of intensity
future_time = future_start:future_end; % future time
%%Unprecedented intensity of WS
temp = wsi_s(i,future_time); %  future time
x_f= x(i,future_time);
wsi_s_f = wsi_s(i,future_time); %  future time

%%
temp(temp <= thres_i) = 0; % if there is unprecedented water scarcity, set to 0
temp(temp >  thres_i) = 1; % set unprecedented water scarcity to 1
%%
%%calculate the duration of UIWS
run_lengths = diff([0, temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
start_index = find(run_lengths == 1); %
end_index = find(run_lengths == -1) - 1;
durations = end_index - start_index + 1;
%%mask the event less than min_month
id_i_dur=find(durations>=min_month);%last at leat for 5 months
widths=durations(id_i_dur);
%
if length(widths)>0
    indices=start_index(id_i_dur)+future_start-1;%convert to all series
    %find origional data
    %%%the index of UIWS at leat 5 months
    id_mask = [];
    for k = 1:length(indices)
        id_mask = [id_mask, indices(k):indices(k)+widths(k)-1];
    end
    %%Extract future series
    temp0=zeros(1,length(future_time));

    toe=0;
    for k = 1:length(id_mask)  %这可能错误，一个事件只有一个标记
        loc_i_f=id_mask(k);
        future_series=wsi_s(i,loc_i_f-ave_window+1:loc_i_f);
        future_x=x(i,loc_i_f-ave_window+1:loc_i_f);
        %%significant test:  two-sample Kolmogorov-Smirnov test
        [h,p]=my_kstest2(base_series',future_series');%null hypothesis vectors x1 and x2 are from the same continuous, 0:same; 1,significant difference.
        if h==1 %significant difference
            %%Repalce origional data
            temp0(1,loc_i_f-future_start+1)=1;
            toe=toe+1;
            if toe==1
                max_i_f=wsi_s_f(1,loc_i_f-future_start+1);
                1850+floor((loc_i_f-1)/12)
                break
            end
        end
    end  %end for k

    uiws=squeeze(nansum(reshape(temp0,12,[]),1));


    %%
    tmp=uiws;
    tmp(tmp>=1)=1;
    %%
    mon_i(i,:)=uiws;
    A_i(i,:)=A*tmp;
    P_i(i,:)=P(i,:).*tmp;

    %%UIWS
    %
    idi=find(uiws>=1);
    if idi
        TFE_i=2015+idi(1)-1;%Time of first emergence
    end

end %end if length(widths)>0
%%Unprecedented intensity
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
figure('unit','centimeters','Position',[2 2 32 25]);

%%%%%%%%%%plot1
subplot_tight(4,1,1,[0.04,0.04]);

%xi = x(1):1:x(end);
%yi = interp1(x,y,xi,'linear'); %only show 1850-2100
xi = x;
yi =y;

hold on
h2=area(xi,bsxfun(@max, yi,thres_i),thres_i,'LineStyle','none','ShowBaseLine','off','FaceColor','r','EdgeColor','none','FaceAlpha',0.9);%
h4=plot(x(1:(2014-1850)*12+12),y(1:(2014-1850)*12+12),'Linestyle','-','Marker','none','color','k','Linewidth',1);
hold on
h5=plot(x((2015-1850)*12+1:(2100-1850)*12+12),y((2015-1850)*12+1:end),'Linestyle','-','Marker','none','color','b','Linewidth',1);
%%%%%%%%%%
hold on;
xline(2014,'k--','2014','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
hold on;
xline(1900,'k--','1900','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
hold on;
%yline(max(y(1:hist_end)),'k--',['Maximum intensity: ',num2str(max(y(1:hist_end)),'%.2f')],'Linewidth',0.01,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'left','Fontsize',12,'fontweight','bold');
line([x(1) x(hist_end)],[max(y(1:hist_end)),max(y(1:hist_end))],'linestyle','--','color','r','Linewidth',1.0)
hold on;
yline(thres_i,'r-',['Dynamic threshold: ',num2str(thres_i,'%.2f')],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment', 'center','Fontsize',12,'fontweight','bold');
hold on;
%yline(0.4,'k--','WS threshold: 0.40','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'left','Fontsize',12,'fontweight','bold');
fill([x(1)+0.1,x(1)+0.1,x(end)-0.1,x(end)-0.1],[0,thres_i,thres_i,0],[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.4);
%
line([x(1) x(pre_industrial)],[max(y(1:pre_industrial)) max(y(1:pre_industrial))],'linestyle','--','color','r','Linewidth',1.0)
line([x(1) x(pre_industrial)],[min(y(1:pre_industrial)) min(y(1:pre_industrial))],'linestyle','--','color','r','Linewidth',1.0)
fill([x(1),x(1),x(pre_industrial), x(pre_industrial)],[min(y(1:pre_industrial)),max(y(1:pre_industrial)),max(y(1:pre_industrial)),min(y(1:pre_industrial))],[0.2 0.2 0.2],'EdgeColor','none','FaceAlpha',0.4);
text(1901,(max(y(1:pre_industrial))+min(y(1:pre_industrial)))/2,['Noise (1SD): ',num2str(Noise_i,'%.2f')],'rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');%,'HorizontalAlignment', 'right'
text(1855,0.92*max(y(1:hist_end)),['Maximum intensity: ',num2str(max(y(1:hist_end)),'%.2f')],'rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');
%
set(gca,'Xlim',[1850 2100],'xtick',[1850:50:2100],'xticklabel',[])
set(gca,'ylim',[0,0.8],'ytick',0:0.2:0.8,'yticklabel',{'0.0','0.2','0.4','0.6','0.8'})%
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


hold on;
%xline(TFE_i,'k-',['TFE_i: ',num2str(TFE_i+1),'   (>=5 yr)'],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment', 'left','Fontsize',12,'fontweight','bold');
xline(TFE_i,'b-','','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment', 'left','Fontsize',12,'fontweight','bold');
%text('Units','normalized','Position',[0.70, 0.90],'String',['TFE_i: ',num2str(TFE_i),' (WSI > ',num2str(max_i_h,'%.2f'),')'],'rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');
text('Units','normalized','Position',[0.70, 0.20],'String',['TFE_i: ',num2str(TFE_i)],'rotation',0,'Color','b','FontName','times new roman','FontSize',12,'fontweight','bold');
text('Units','normalized','Position',[0.68, 0.10],'String',['(WSI>',num2str(max_i_h,'%.2f'),', KS(X1,X2)=1, at least 5 months)'],'rotation',0,'Color','b','FontName','times new roman','FontSize',12,'fontweight','bold');

%
temp = y(1,1:hist_end); % only account for historical period
%%historical maximum

base_series=temp(1,loc_i_h-ave_window+1:loc_i_h);
base_x=x((loc_i_h-ave_window+1):loc_i_h);
%%
temp = y(future_time); %
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

ct=1;
for j=1:length(durations)
    if durations(j)>=5&ct==1;
        line([start_x(j) end_x(j)+0],[thres_i-0.10 thres_i-0.10],'linestyle','-','color','r','Linewidth',1.0)
        line([start_x(j) start_x(j)],[thres_i thres_i-0.10],'linestyle','-','color','r','Linewidth',1.0)
        line([end_x(j)+0 end_x(j)+0],[thres_i thres_i-0.10],'linestyle','-','color','r','Linewidth',1.0)
        text((start_x(j)+end_x(j))/2,thres_i-0.20,num2str(durations(j)),'rotation',0,'Color','r','FontName','times new roman','FontSize',12,'fontweight','bold','HorizontalAlignment', 'center');
        ct=ct+1;
    end
end


%%%%%processing...
%% Add labels
xlims = get(gca,'XLim');
ylims = get(gca,'YLim');
%x1=0.02+(ToE-xlims(1))/diff(xlims);
%x2=0.02+(ToD-xlims(1))/diff(xlims);
%
%annotation('doublearrow',[0.02 (2014-1850)/(2100-1850)],[0.65 0.65]);
%text(1900,1.4,'The historical period','rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');
%text(TFE_w,1.6,'(>=5 yr)','rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');
%

hold on
fill([base_x(1),base_x(1),base_x(end),base_x(end)],[min(ylims),max(base_series),max(base_series),min(ylims)],[0.2,0.2,0.2],'EdgeColor','none','FaceAlpha',0.3);
%fill([base_x(1),base_x(1),base_x(end),base_x(end)],[min(ylims),max(ylims),max(ylims),min(ylims)],[0.2,0.2,0.2],'EdgeColor','none','FaceAlpha',0.3);
hold on
%fill([future_x(1),future_x(1),future_x(end),future_x(end)],[min(future_series),max(future_series),max(future_series),min(future_series)],[0.2,0.2,0.2],'EdgeColor','none','FaceAlpha',0.3);
fill([future_x(1),future_x(1),future_x(end),future_x(end)],[min(ylims),max(future_series),max(future_series),min(ylims)],[0.2,0.2,0.2],'EdgeColor','none','FaceAlpha',0.3);

hold on
plot(base_x,base_series,'Linestyle','-','Marker','none','color','c','Linewidth',2);
hold on
plot(future_x,future_series,'Linestyle','-','Marker','none','color','y','Linewidth',2);

text((base_x(1)+base_x(end))/2,max(base_series)+0.01,'x1','FontSize',12,'fontweight','bold','color','b','HorizontalAlignment', 'center');
text((future_x(1)+future_x(end))/2,max(future_series)+0.01,'x2','FontSize',12,'fontweight','bold','color','b','HorizontalAlignment', 'center');



lgd=legend([h2,h4,h5,],{'UIWS','WSI (1850-2014)','WSI (2015-2100)'},'Location','northwest','FontSize',14,'FontWeight','bold','NumColumnsMode','manual','NumColumns',3);
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
text('Units','normalized','Position',[0.005, 1.09],'String','(a)','rotation',0,'Color','k','FontName','times new roman','FontSize',16,'fontweight','bold');
%
%title(['Representative Grid (36.5°N, 3.0°E) in {\color{blue}watergap2-2e' '\color[rgb]{0 0 0}\_' '\color[rgb]{1 0 0}mpi-esm1-2-hr}'],'FontWeight','bold','Fontsize',14,'FontName','times new roman')
title('Unprecedented Intensity of Water Scarcity','FontName','times new roman','FontSize',14,'fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold','FontName','Times new Roman');
%%%%%%%%%%plot1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%Unprecedented duration

temp = wss_ori(i,:); % whole period
%%1:detect the Non-WS duration between two WS evetns
run_lengths1 = diff([1, temp, 1]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
start_index1 = find(run_lengths1 == -1); %
end_index1 = find(run_lengths1 == 1) - 1;
durations1 = end_index1 - start_index1 + 1;
%find origional data
id_min_dur1=find(durations1<=max_gap);
widths1=durations1(id_min_dur1);
indices1=start_index1(id_min_dur1);
%%
id_min1 = [];
for k = 1:length(indices1)
    id_min1 = [id_min1, indices1(k):indices1(k)+widths1(k)-1];
end
%%Repalce origional data
temp(id_min1)= ~temp(id_min1);
%% The short periods of non-WS months between lengthy periods of WS months were considered WS due to the pooling effect.
%%%%%%
%%2:discard the WS duration between two Non-WS evetns
run_lengths2 = diff([0, temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
start_index2 = find(run_lengths2 == 1); %
end_index2 = find(run_lengths2 == -1) - 1;
durations2 = end_index2 - start_index2 + 1;
%find origional data
id_min_dur2=find(durations2<threshold_d);% ignore the events less than two years. eg. 12 mongths
widths2=durations2(id_min_dur2);
indices2=start_index2(id_min_dur2);
%%
id_min2 = [];
for k = 1:length(indices2)
    id_min2 = [id_min2, indices2(k):indices2(k)+widths2(k)-1];
end
%%Repalce origional data
temp(id_min2)= ~temp(id_min2);
%% The duration of a negligible WS period shorter than x-months
%%%%%%
%%%%%%
%%%%%%
%%3: detect the WS evetns
run_lengths = diff([0, temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
start_index = find(run_lengths == 1); %
end_index = find(run_lengths == -1) - 1;
durations = end_index - start_index + 1;
%%%%%%
%%%%%%
%%4: %%Assign the length of water scarcity events to each period
idx_ones = find(temp == 1);% find all locations with WS
% the boundary with consecutive WS
diff_idx = [1,diff(idx_ones)]; % Calculate the difference between adjacent WS
boundaries = find(diff_idx > 1); % Find the boundaries of different continuous WS
% Calculate the length of each continuous WS
starts = [1, boundaries]; % Index at the beginning of consecutive WS
ends = [boundaries - 1,length(idx_ones)]; % Index at the end of consecutive WS
segment_lengths = ends - starts + 1; % The length of each continuous WS
%
group_ids = cumsum([1,diff_idx > 1]); % the group ID of each WS
group_ids = group_ids(1,2:end);%Because one value was added, the initial redundant values were deleted
length_d = zeros(size(temp)); % Initialize output
if length(idx_ones) >0
    length_d(idx_ones) = segment_lengths(group_ids); % Allocate length to each original position  %下标赋值维度不匹配 当转为mex函数时,且idx_ones为空报错，在普通函数正常
end
%
%%eg. temp = [0, 1, 1, 0, 1, 1, 1, 0, 1, 0];% 1 denote WS
%% length_d= [0, 2, 2, 0, 3, 3, 3, 0, 1, 0];
%%4: %%Assign the length of water scarcity events to each period
%%%%%%%%%%%%%%%
%%%%%%
if length(durations)>0 %if no WS duration, skip...
    %% find if the year 2015 is within the WS duaration and dates
    %% In whole peried, Detect if the 2015 is within the duration of WS events


    if temp(hist_end)==1 %if year 2015 is with the duration of WS events, it is very complex
        %%Detect the start and end date of WS events corss 2015
        id_st=hist_end; %Not fully defined in some execution paths
        id_ed=hist_end; %Not fully defined in some execution paths
        loc=hist_end;
        while loc<future_end&temp(loc+1)
            id_ed=all_time(loc);
            loc=loc+1;
        end
        %%
        loc=hist_end;
        while loc>=2&temp(loc-1)
            loc=loc-1;
            id_st=all_time(loc);
        end
        %%Detect the start and end date of WS events corss 2015
        max_d2=id_ed-id_st+1;
        max_d21=hist_end-id_st+1;
        max_d22=id_ed-hist_end;
        %%%%%%%%%%%%
        %%1：id_st-1,id_st:id_ed,id_ed+1:end
        temp_h1=temp(1:id_st-1);% only account for historical first period
        run_lengths = diff([0, temp_h1, 0]); % % calculate the transition form non-water scarcity (0) to water scarcity (1)
        start_index = find(run_lengths == 1); %
        end_index = find(run_lengths == -1) - 1;
        durations = end_index - start_index + 1;

        %%%%%%%% find historical maximum duration
        if length(durations)>0 %
            [max_d1,id_max_dur_h] = nanmax(durations,[],2); % find historical maximum duration of all previous WS event
            thres_d = max([max_d1 max_d21]);
        else
            thres_d=max_d21;
        end
        if thres_d<threshold_d
            thres_d=threshold_d;%%very important, or use the original value
        end
        max_d_h=thres_d;
        %%%%%%%%%% find historical maximum duration

        %%%%%%%%%% judge if there will have WS in the future
        temp_f2=temp(1,id_ed+1:end);% only account for future second period
        run_lengths = diff([0, temp_f2, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
        start_index = find(run_lengths == 1); %
        end_index = find(run_lengths == -1) - 1;
        durations = end_index - start_index + 1;
        if length(durations)>0 %Or no WS if no value
            [max_d3,~] = nanmax(durations,[],2); % find future maximum duration
        else
            max_d3=0;% no WS
        end
        %%%%%%%%%%

        %%muti scenarios
        [~,id] = max([thres_d max_d2 max_d3]);
        if id==1 %No TFE_d in the future
            TFE_d=2101;%after 2100
            max_d_f=max([max_d2,max_d3]);
        else% Have TFE_d in the future
            temp_f=temp(id_st+1:end);% account for future period plus the second historical period
            length_f=length_d(id_st+1:end);
            run_lengths = diff([0,temp_f,0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
            start_index = find(run_lengths == 1); %
            end_index = find(run_lengths == -1) - 1;
            durations = end_index - start_index + 1;
            id_min_dur=find(durations>thres_d);%%duration less than the thres_d is not regarded as udws
            widths=durations(id_min_dur);
            indices=start_index(id_min_dur);
            %% Repalce origional data
            temp0=zeros(1,length(temp_f));
            toe=0;
            for k = 1:length(indices)
                temp0(1,indices(k)+thres_d:indices(k)+widths(k)-1)=1;%ignore the length less than thres_d, must plus +thres_d !!!
                toe=toe+1;
                if toe==1
                    max_d_f=length_f(1,indices(k)+thres_d);
                    2015+floor((indices(k)+thres_d-max_d21)/12)
                end
            end
            %%temp0(id_ed+1:end)=0; %%the current is the maximum duration,after this event,no UDWS
            temp_f1 =temp0(1,hist_end-id_st+1:end);%% need to delete the second historical period
            length_f1 =length_f0(1,hist_end-id_st+1:end);
            udws=squeeze(nansum(reshape(temp_f1,12,[]),1));
            length_d1=squeeze(nanmean(reshape(length_f1,12,[]),1));
            %%
            tmp=udws;
            tmp(tmp>=1)=1;
            %%

            mon_d(i,:)=udws;
            A_d(i,:)=A*tmp;
            P_d(i,:)=P(i,:).*tmp;
            %%
            idd=find(udws>=1);%%
            if idd
                TFE_d=2015+idd(1)-1;%Time of first emergence
            end
            %%Unprecedented durations
        end %%%end muti scnerios
        %%muti scenarios

        %%%%%%%%%%%%
        %%%%%%%%%%%%
        %%%%%%
    else %if 2015 not in ,just seperately excute
        temp_h=temp(1,1:hist_end);% only account for historical period
        run_lengths = diff([0, temp_h, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
        start_index = find(run_lengths == 1); %
        end_index = find(run_lengths == -1) - 1;
        durations = end_index - start_index + 1;
        if length(durations)>0 %
            [thres_d,id_max_dur_h] = nanmax(durations,[],2); % find historical maximum duration
        else
            thres_d=threshold_d;
        end
        if thres_d<threshold_d
            thres_d=threshold_d;
        end
        max_d_h=thres_d;

        %%find future unprecedented duration of WS
        temp_f=temp(1,future_time);% only account for future period
        length_f=length_d(1,future_time);% only account for future period
        run_lengths = diff([0,temp_f,0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
        start_index = find(run_lengths == 1); %
        end_index = find(run_lengths == -1) - 1;
        durations = end_index - start_index + 1;
        id_min_dur=find(durations>thres_d);
        widths=durations(id_min_dur);
        indices=start_index(id_min_dur);
        %%Repalce origional data
        temp0=zeros(1,length(future_time));

        toe=0;
        for k = 1:length(indices)
            temp0(1,indices(k)+thres_d:indices(k)+widths(k)-1)=1;
            toe=toe+1;
            if toe==1
                max_d_f=length_f(1,indices(k)+thres_d);
                2015+floor((indices(k)+thres_d)/12)
            end
        end
        udws=squeeze(nansum(reshape(temp0,12,[]),1));
        %%
        tmp=udws;
        tmp(tmp>=1)=1;
        %%
        mon_d(i,:)=udws;
        A_d(i,:)=A*tmp;
        P_d(i,:)=P(i,:).*tmp;
        %%
        %%Unprecedented durations
        idd=find(udws>=1);%
        if idd
            TFE_d=2015+idd(1)-1;%Time of first emergence
        end
    end %end   if temp(hist_end)==1
end %end if length(durations)>0
%%%%%%
%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%plot2
subplot_tight(4,1,2,[0.04,0.04]);
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% duration
wsi=y_ori;
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
id_min_dur1=find(durations1<=max_gap);
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
id_min_dur=find(durations<threshold_d);
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
yi = interp1(x,wsi_ori,xi,'linear'); %only show 1850-2100
h1=area(x,bsxfun(@max, wsi_ori, 0.4),0.4,'LineStyle','none','ShowBaseLine','off','FaceColor','r','EdgeColor','none','FaceAlpha',0.4);%
hold on
h2=area(x,bsxfun(@min, wsi_ori, 0.4),0.4,'LineStyle','none','ShowBaseLine','off','FaceColor','b','EdgeColor','none','FaceAlpha',0.4);%
hold on
xline(2014,'k--','2014','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
hold on
yline(0.4)

%
hold on;
if TFE_d<2100
    hold on
    xline(TFE_d,'b-',['TFE_d: ',num2str(TFE_d),' (Duration > ',num2str(max_d_h),' months)'],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
else
    text('Units','normalized','Position',[0.90, 0.90],'String',['No  TFE_d: (Max Duration_f < ',num2str(max_d_h),' months)'],'rotation',0,'Color','b','FontName','times new roman','FontSize',12,'fontweight','bold');
end
%%%%%%%%%%
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


j=id_max_dur_h(1);
fill([start_x1(j)-0.1,start_x1(j)-0.1,end_x1(j)+0.1,end_x1(j)+0.1],[0,ylims(2),ylims(2),0],[0.4 0.4 0.4],'EdgeColor','none','FaceAlpha',0.4);
j=id_max_dur_h(1);
%xline((start_x1(j)+end_x1(j))/2,'w--',['Dynamic threshold: ',num2str(max_d_h+0,'%.0f'),' months'],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
%annotation('textarrow',[(start_x1(j)-1850)/(2100-1850+1)+0.2,(end_x1(j)-1850)/(2100-1850+1)+0.2],[0.60 0.60],'String',['Past maximum duration: ',num2str(max_d_h+0,'%.0f'),' months'],'FontSize',12,'fontweight','bold','color','r');
text((start_x1(j)+end_x1(j))/2,0.80,['Dynamic threshold: ',num2str(thres_d+0,'%.0f'),' months'],'FontSize',12,'fontweight','bold','color','r','HorizontalAlignment', 'center');

%annotation('doublearrow',[0.02 (2014-1850)/(2100-1850)],[0.4 0.4]);
%text(1900,1.4,'The historical period','rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');
%
thres=0.4;
ct=1;
up=1;
dw=1;
for j=1:length(durations1)
    line([start_x1(j) end_x1(j)+0],[thres+0.10 thres+0.10],'linestyle','-','color','k','Linewidth',1.05)
    line([start_x1(j) start_x1(j)],[thres thres+0.10],'linestyle','-','color','k','Linewidth',1.0)
    line([end_x1(j)+0 end_x1(j)+0],[thres thres+0.10],'linestyle','-','color','k','Linewidth',1.0)
    if durations1(j)>50
        text((start_x1(j)+end_x1(j))/2,thres+0.20,num2str(durations1(j)),'rotation',0,'Color','k','FontName','times new roman','FontSize',10,'fontweight','bold','HorizontalAlignment', 'center');
    end
    up=up+1;
end

for j=1:length(durations2)
    if durations2(j)>thres_d
        line([start_x2(j) end_x2(j)+0],[thres+0.10+0.01 thres+0.10+0.01],'linestyle','-','color','r','Linewidth',1.0)
        line([start_x2(j) start_x2(j)],[thres thres+0.10+0.01],'linestyle','-','color','r','Linewidth',1.0)
        line([end_x2(j)+0 end_x2(j)+0],[thres thres+0.10+0.01],'linestyle','-','color','r','Linewidth',1.0)
        text((start_x2(j)+end_x2(j))/2,thres+0.20,num2str(durations2(j)),'rotation',0,'Color','r','FontName','times new roman','FontSize',10,'fontweight','bold','HorizontalAlignment', 'center');
        hold on
        h8=fill([start_x2(j)+thres_d/12,start_x2(j)+thres_d/12,end_x2(j),end_x2(j)],[thres,thres+0.10,thres+0.10,thres],[255,1,1]/255,'EdgeColor','none','FaceAlpha',1);%[255,153,153]/255
        %
        line([start_x2(j) start_x2(j)+thres_d/12],[thres-0.10 thres-0.10],'linestyle','-','color','k','Linewidth',1.0)
        line([start_x2(j) start_x2(j)],[thres thres-0.10],'linestyle','-','color','k','Linewidth',1.0)
        line([start_x2(j)+thres_d/12 start_x2(j)+thres_d/12],[thres thres-0.10],'linestyle','-','color','k','Linewidth',1.0)
        text((start_x2(j)+start_x2(j)+thres_d/12)/2,thres-0.20,num2str(thres_d),'rotation',0,'Color','k','FontName','times new roman','FontSize',10,'fontweight','bold','HorizontalAlignment', 'center');
    else
        line([start_x2(j) start_x2(j)],[thres+0.10 thres+0.10],'linestyle','-','color','k','Linewidth',1.0)
        line([start_x2(j) start_x2(j)],[thres thres+0.10],'linestyle','-','color','k','Linewidth',1.0)
        line([end_x2(j)+0 end_x2(j)+0],[thres thres+0.10],'linestyle','-','color','k','Linewidth',1.0)
        %text((start_x2(j)+end_x2(j))/2,thres+0.20,num2str(durations2(j)),'rotation',0,'Color','k','FontName','times new roman','FontSize',10,'fontweight','bold','HorizontalAlignment', 'center');
    end
end %end for
box on

lgd=legend([h8,h1,h2,],{'UDWS','WS','Non-WS'},'Location','northwest','FontSize',14,'FontWeight','bold','NumColumnsMode','manual','NumColumns',3);
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



%% Add labels
xlims = get(gca,'XLim');
ylims = get(gca,'YLim');
%x1=0.02+(ToE-xlims(1))/diff(xlims);
%x2=0.02+(ToD-xlims(1))/diff(xlims);

text('Units','normalized','Position',[0.005, 1.09],'String','(b)','rotation',0,'Color','k','FontName','times new roman','FontSize',16,'fontweight','bold');
title('Unprecedented Duration of Water Scarcity','FontName','times new roman','FontSize',14,'fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold','FontName','Times new Roman');
%%%%%%%%%%plot2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%
%%Unprecedented frenquency
% Identify the maximum value during the historical baseline period

temp = wss_ori(i,:); % whole period
%%1:detect the Non-WS duration between two WS evetns
run_lengths1 = diff([1, temp, 1]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
start_index1 = find(run_lengths1 == -1); %
end_index1 = find(run_lengths1 == 1) - 1;
durations1 = end_index1 - start_index1 + 1;
%find origional data
id_min_dur1=find(durations1<=max_gap);
widths1=durations1(id_min_dur1);
indices1=start_index1(id_min_dur1);
%%
id_min1 = [];
for k = 1:length(indices1)
    id_min1 = [id_min1, indices1(k):indices1(k)+widths1(k)-1];
end
%%Repalce origional data
temp(id_min1)= ~temp(id_min1);%% The short periods of non-WS months between lengthy periods of WS months were considered WS due to the pooling effect.
%%%%%%
%%%%%%
%%2:discard the WS duration between two Non-WS evetns
run_lengths2 = diff([0, temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
start_index2 = find(run_lengths2 == 1); %
end_index2 = find(run_lengths2 == -1) - 1;
durations2 = end_index2 - start_index2 + 1;
%find origional data
id_min_dur2=find(durations2<min_month);
widths2=durations2(id_min_dur2);
indices2=start_index2(id_min_dur2);
%%
id_min2 = [];
for k = 1:length(indices2)
    id_min2 = [id_min2, indices2(k):indices2(k)+widths2(k)-1];
end
%%Repalce origional data
temp(id_min2)= ~temp(id_min2);%% The duration of a negligible WS period shorter than x-months
%%%%%%
%%%%%%
%%3: cal ES frequency,here consider the short periods of non-WS months between lengthy periods of WS months
transitions = diff([temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
transitions(transitions~=1)=0;
freq=filter(ones(1,freq_window)/freq_window,1,transitions,[],2);%the frequency of local and previous ave_window-1 values
freq(1,1:freq_window-1)=nan;%the first ave_window-1 value is nan
freq=movmean(freq,[ave_window-1 0],2,"omitnan"); %ave_window:2+0+1 the mean of local and previous ave_window-1 values
%%%%cal WS frequency

%%%%%%%%%%%%%%%%%%%% Identify the maximum value during the historical baseline period
temp = freq(1,1:hist_end); % only account for historical period
Noise_f=nanstd(temp(1:pre_industrial),0,2);
%Noise_f=0;
%%historical maximum
[thres_f,loc_f_h]=nanmax(temp,[],2);
thres_f=thres_f+1*Noise_f;% important
max_f_h=thres_f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%plot3
subplot_tight(4,1,3,[0.04,0.04]);
%%%%
xi = x(1):0.01:x(end);
yi = interp1(x,freq,xi,'linear'); %only show 1850-2100
hold on
h2=area(xi,bsxfun(@max, yi, thres_f),thres_f,'LineStyle','none','ShowBaseLine','off','FaceColor','r','EdgeColor','none','FaceAlpha',0.9);%
h4=plot(x(1:(2014-1850)*12+12),freq(1:(2014-1850)*12+12),'Linestyle','-','Marker','none','color','k','Linewidth',1);
hold on
h5=plot(x((2015-1850)*12+1:(2100-1850)*12+12),freq((2015-1850)*12+1:end),'Linestyle','-','Marker','none','color','b','Linewidth',1);
%%%%%%%%%%
hold on;
xline(1900,'k--','1900','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
hold on;
xline(2014,'k--','2014','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
hold on;
%yline(max(freq(1:hist_end)),'k--',['Maximum frequency: ',num2str(max(freq(1:hist_end)),'%.2f')],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'left','Fontsize',12,'fontweight','bold');
hold on
line([x(1) x(hist_end)],[max(freq(1:hist_end)),max(freq(1:hist_end))],'linestyle','--','color','r','Linewidth',1.0)
hold on
yline(thres_f,'r-',['Dynamic threshold: ',num2str(thres_f,'%.2f')],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment', 'center','Fontsize',12,'fontweight','bold');
hold on;
%yline(0.4,'r--','water scarcity threshold: 0.40','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'center','Fontsize',12,'fontweight','bold');
fill([x(1)+0.1,x(1)+0.1,x(end)-0.1,x(end)-0.1],[0,thres_f,thres_f,0],[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.4);
%
line([x(1) x(pre_industrial)],[max(freq(1:pre_industrial)) max(freq(1:pre_industrial))],'linestyle','--','color','r','Linewidth',1.0)
line([x(1) x(pre_industrial)],[min(freq(1:pre_industrial)) min(freq(1:pre_industrial))],'linestyle','--','color','r','Linewidth',1.0)
fill([x(1),x(1),x(pre_industrial), x(pre_industrial)],[min(freq(1:pre_industrial)),max(freq(1:pre_industrial)),max(freq(1:pre_industrial)),min(freq(1:pre_industrial))],[0.2 0.2 0.2],'EdgeColor','none','FaceAlpha',0.4);
text(1901,(max(freq(1:pre_industrial))+min(freq(1:pre_industrial)))/2,['Noise (1SD): ',num2str(Noise_f,'%.2f')],'rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');%,'HorizontalAlignment', 'right'
text(1960,0.92*max(freq(1:hist_end)),['Maximum frequency: ',num2str(max(freq(1:hist_end)),'%.2f')],'rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');

%
%
set(gca,'Xlim',[1850 2100],'xtick',[1850:50:2100],'xticklabel',[])
set(gca,'ylim',[0,0.1],'ytick',0:0.02:0.1,'yticklabel',{'0.0','.02','.04','.06','.08',0.10''})%
%set(gca,'ylim',[0,0.2],'ytick',[0,0.5])%
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

hold on;
%xline(TFE_f+1,'k-',['TFE_f: ',num2str(TFE_f+1),'   (>=5 yr)'],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment', 'left','Fontsize',12,'fontweight','bold');
if TFE_f<2100
    %xline(TFE_f,'k-',['TFE_f: ',num2str(TFE_f+1),'   (>=5 yr)'],'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment', 'left','Fontsize',12,'fontweight','bold');
    text('Units','normalized','Position',[0.70, 0.90],'String',['TFE_f: ',num2str(TFE_f),' (Freq > ',num2str(max_f_h,'%.2f'),')'],'rotation',0,'Color','b','FontName','times new roman','FontSize',12,'fontweight','bold');
else
    text('Units','normalized','Position',[0.70, 0.90],'String',['No TFE_f: (Max Freq_f < ',num2str(max_f_h,'%.2f'),')'],'rotation',0,'Color','b','FontName','times new roman','FontSize',12,'fontweight','bold');
end
%text('Units','normalized','Position',[0.70, 0.78],'String','No X2, no Kolmogorov-Smirnov test','rotation',0,'Color','b','FontName','times new roman','FontSize',12,'fontweight','bold');
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
run_lengths = diff([0, temp, 0]); % 0  [0 1 1 1 0 1 0] 0  dif:[0 1 0 0 -1 1 -1 0]
start_index = find(run_lengths == 1); %  %[2,6]
end_index = find(run_lengths == -1) - 1; % %[5,7]
durations = end_index - start_index + 1; % %[3,1]
start_x=x(hist_end+start_index);
end_x=x(hist_end+end_index);
%%%%%processing...
%% Add labels
xlims = get(gca,'XLim');
ylims = get(gca,'YLim');
%x1=0.02+(ToE-xlims(1))/diff(xlims);
%x2=0.02+(ToD-xlims(1))/diff(xlims);
%
%annotation('doublearrow',[0.02 (2014-1850)/(2100-1850)],[0.65 0.65]); %X 和 Y 值必须介于 0 与 1 之间
%text(1900,1.4,'The historical period','rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');
%text(TFE_w,1.6,'(>=5 yr)','rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');
%
ct=1;
future_x=[];
future_series=[];
for j=1:length(durations)
    if durations(j)>=6&ct==1;
        line([start_x(j) end_x(j)+1],[thres_i-0.10 thres_i-0.10],'linestyle','-','color','r','Linewidth',2.0)
        line([start_x(j) start_x(j)],[thres_i thres_i-0.10],'linestyle','-','color','r','Linewidth',2.0)
        line([end_x(j)+1 end_x(j)+1],[thres_i thres_i-0.10],'linestyle','-','color','r','Linewidth',2.0)
        text((start_x(j)+end_x(j))/2,thres_i-0.20,num2str(durations(j)),'rotation',0,'Color','r','FontName','times new roman','FontSize',12,'fontweight','bold');
        ct=ct+1;
    end
end


lgd=legend([h2,h4,h5,],{'UFWS','Freq (1850-2014)','Freq (2015-2100)'},'Location','northwest','FontSize',14,'FontWeight','bold','NumColumnsMode','manual','NumColumns',3);
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
text('Units','normalized','Position',[0.005, 1.09],'String','(c)','rotation',0,'Color','k','FontName','times new roman','FontSize',16,'fontweight','bold');

title('Unprecedented Frequency of Water Scarcity','FontName','times new roman','FontSize',14,'fontweight','bold')
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

%h1=plot(x,y_ori,'Linestyle','-','Marker','none','color','c','Linewidth',0.01);
hold on
h2=area(xi,bsxfun(@max, yi, 0.4),0.4,'LineStyle','none','ShowBaseLine','off','FaceColor','r','EdgeColor','none','FaceAlpha',0.4);%
h4=plot(x(1:hist_end),y(1:hist_end),'Linestyle','-','Marker','none','color','k','Linewidth',1);
hold on
h5=plot(x(hist_end+1:end),y(hist_end+1:end),'Linestyle','-','Marker','none','color','b','Linewidth',1);
%%%%%%%%%%
hold on;
xline(1900,'k--','1900','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
hold on;
xline(2014,'k--','2014','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','bottom','LabelHorizontalAlignment','center','Fontsize',12,'fontweight','bold');
hold on;
yline(0.4,'r-','WS threshold: 0.40','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','top','LabelHorizontalAlignment','left','Fontsize',12,'fontweight','bold');
fill([x(2),x(2),x(end-1),x(end-1)],[0,0.4,0.4,0],[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.4);
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

hold on;
%xline(TFE_w,'k-',num2str(TFE_w),'Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'center','Fontsize',12,'fontweight','bold');
xline(TFE_w,'b-','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'center','Fontsize',12,'fontweight','bold');
text('Units','normalized','Position',[0.70, 0.90],'String',['TFE_{w}: ',num2str(TFE_w+0),' (WSI > ',num2str(0.4,'%.2f'),')'],'rotation',0,'Color','b','FontName','times new roman','FontSize',12,'fontweight','bold');
%% Add labels
xlims = get(gca,'XLim');
ylims = get(gca,'YLim');
%x1=0.02+(ToE-xlims(1))/diff(xlims);
%x2=0.02+(ToD-xlims(1))/diff(xlims);
%
%text(TFE_w,1.6,'(>=5 yr)','rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');
annotation('doublearrow',[0.04 (2014-1850)/(2100-1850+1)-0.01],[0.18 0.18]); %X 和 Y 值必须介于 0 与 1 之间
%text(1900,0.97,'The reference period','rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');
text('Units','normalized','Position',[0.31, 0.65],'String','Reference period','rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');
%
annotation('doublearrow',[0.04 (1900-1850)/(2100-1850+1)+0.025],[0.085 0.085]);
text('Units','normalized','Position',[0.025, 0.10],'String','Pre-industrial period','rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');
%
thres=0.4;
if 0
    for j=1:length(durations2)
        line([start_x2(j) end_x2(j)+1],[thres-0.10 thres-0.10],'linestyle','-','color','r','Linewidth',2.0)
        line([start_x2(j) start_x2(j)],[thres thres-0.10],'linestyle','-','color','r','Linewidth',2.0)
        line([end_x2(j)+1 end_x2(j)+1],[thres thres-0.10],'linestyle','-','color','r','Linewidth',2.0)
        text((start_x2(j)+end_x2(j))/2,thres-0.20,num2str(durations2(j)),'rotation',0,'Color','r','FontName','times new roman','FontSize',12,'fontweight','bold');
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
text('Units','normalized','Position',[0.005, 1.09],'String','(d)','rotation',0,'Color','k','FontName','times new roman','FontSize',16,'fontweight','bold');
%title('Representative Grid (36.5°N, 3.0°E) in atergap2-2e\_mpi-esm1-2-hr','FontWeight','bold','Fontsize',14,'FontName','times new roman');
%title(['\fontsize{14}Representative Grid (36.5°N, 3.0°E) in {\color{blue}watergap2-2e ' '\_|' '\color[rgb]{0 .5 .5 }mpi-esm1-2-hr}'])
%
%text('Units','normalized','Position',[0.05, 0.50],'String',['Events: ',num2str(length(durations2))],'rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');
%text('Units','normalized','Position',[0.15, 0.50],'String',['Maximum Duration: ',num2str(max(durations2)),' yr'],'rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');
%text('Units','normalized','Position',[0.32, 0.50],'String',['Cumulative occurring time: ',num2str(sum(durations2)),' yr'],'rotation',0,'Color','k','FontName','times new roman','FontSize',12,'fontweight','bold');
%lgd=legend([h2,h4,h5,h1],{'WS','WSI (1850-2014)','WSI (2015-2100)','WSI (Original)'},'Location','northwest','FontSize',14,'FontWeight','bold','NumColumnsMode','manual','NumColumns',4);
lgd=legend([h2,h4,h5],{'WS','WSI (1850-2014)','WSI (2015-2100)'},'Location','northwest','FontSize',14,'FontWeight','bold','NumColumnsMode','manual','NumColumns',4);
lgd.Box='off';


title('Conventional Water Scarcity','FontName','times new roman','FontSize',14,'fontweight','bold')
set(gca,'fontsize',12,'fontweight','bold','FontName','Times new Roman');
%%%%%%%%%%plot4



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save pic
%exportgraphics(gcf,'D:\ws_fluctuation\figures\Fig.s1.tif','Resolution',350);% no white
exportgraphics(gcf,'D:\ws_fluctuation\figures\Fig.s1.jpg','Resolution',350);% no white
%%%%%%%%%%%%%%%%%%%%%
