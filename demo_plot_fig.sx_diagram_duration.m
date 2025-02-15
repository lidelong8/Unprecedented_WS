
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%the tested WSI series
x=linspace(1850,2100,3012);
rand('seed',10);
y1=rand(3012,1);
y1=y1*0.5;
x0=linspace(0,1,3012)';
y2=-(x0-0.3).*(x0-0.35).*(x0-0.45).*(x0-0.5).*(x0-0.55)*10;
y3=0.30;
y=y1+y2+y3;
y=flipud(y);
y(y<0)=0;
%%%%%%%%%%%%
y_ori=y;
wsi=y_ori';
y2=y+0.3;
y3=y-0.4;y3(y3<0)=0;
%%%%%%%%%%%%the tested WSI series


%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%  
%
%i=1,j=1,k=1
hist_end = (2015-1850+1)*12;
future_start=hist_end+1;
future_end = (2100-1850+1)*12;
future_time = future_start:future_end; % 


ave_window = 20;
freq_window = 20;
min_duration=6; 
min_month=1;
threshold=0.4;
max_gap=3;
threshold_f=1/freq_window;
%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
close all
figure('unit','centimeters','Position',[2 2 25 15]);

%%%%%%%%%%plot1
subplot_tight(3,1,1,[0.06,0.06])
%yyaxis left
h1=area(x,bsxfun(@max, y_ori, 0.4),0.4,'LineStyle','none','ShowBaseLine','off','FaceColor','r','EdgeColor','none','FaceAlpha',0.4);%
hold on
h2=area(x,bsxfun(@min, y_ori, 0.4),0.4,'LineStyle','none','ShowBaseLine','off','FaceColor','b','EdgeColor','none','FaceAlpha',0.4);%
hold on
xline(2014,'k--','2014','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'center','Fontsize',12,'fontweight','bold');
hold on
yline(0.4)

yt=get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',arrayfun(@(x) num2str(x,'%1.1f'),yt,'UniformOutput',false), 'Fontname', 'Times new Roman', 'Fontsize',12, 'Fontname', 'Times new Roman', 'Fontsize',12)
set(gca,'YMinorTick','on')
set(gca,'XTickLabel', []); 
set(gca,'fontsize',12,'fontweight','bold','FontName','Times new Roman');  
ylabel ('WSI','fontsize',12,'fontweight','bold','FontName','Times new Roman'); 
text('Units','normalized','Position',[0.00, 1.05],'String','(a)','Color','k','fontsize',14,'fontweight','bold','rotation',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
strs=['Events: ',num2str(length(durations1)),'   (1 month: ',num2str(length(durations1(durations1==1))),'; 2 months: ',num2str(length(durations1(durations1==2))),'; 3 months: ',num2str(length(durations1(durations1==3))),'; > 3 months: ',num2str(length(durations1(durations1>3))),')']
text('Units','normalized','Position',[0.05, 0.90],'String',strs,'rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');

lg=legend([h1,h2],{'WS','Non-WS'},'Location', 'northeast','NumColumnsMode','manual','NumColumns',2,'fontsize',12,'fontweight','bold','Fontname','Times new Roman'); 
legend('boxoff')

title('Original Time Series of WSI','FontName','Arial','FontSize',14,'fontweight','bold')

%%%%%%%%%%plot1







%%%%%%%%%%plot2
subplot_tight(3,1,2,[0.06,0.06])


h1=area(x,bsxfun(@max, wsi, 0.4),0.4,'LineStyle','none','ShowBaseLine','off','FaceColor','r','EdgeColor','none','FaceAlpha',0.4);%
hold on
h2=area(x,bsxfun(@min, wsi, 0.4),0.4,'LineStyle','none','ShowBaseLine','off','FaceColor','b','EdgeColor','none','FaceAlpha',0.4);%
hold on
xline(2014,'k--','2014','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'center','Fontsize',12,'fontweight','bold');
hold on
yline(0.4)

yt=get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',arrayfun(@(x) num2str(x,'%1.1f'),yt,'UniformOutput',false), 'Fontname', 'Times new Roman', 'Fontsize',12, 'Fontname', 'Times new Roman', 'Fontsize',12)
set(gca,'YMinorTick','on')
set(gca,'XTickLabel', []);  
set(gca,'fontsize',12,'fontweight','bold','FontName','Times new Roman');  
ylabel ('WSI','fontsize',12,'fontweight','bold','FontName','Times new Roman'); 


text('Units','normalized','Position',[0.00, 1.05],'String','(b)','Color','k','fontsize',14,'fontweight','bold','rotation',0);

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

strs=['Events: ',num2str(length(durations)),'   (1 month: ',num2str(length(durations(durations==1))),'; 2 months: ',num2str(length(durations(durations==2))),'; 3 months: ',num2str(length(durations(durations==3))),'; > 3 months: ',num2str(length(durations(durations>3))),')']
text('Units','normalized','Position',[0.05, 0.90],'String',strs,'rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');

lg=legend([h1,h2],{'WS','Non-WS'},'Location', 'northeast','NumColumnsMode','manual','NumColumns',2,'fontsize',12,'fontweight','bold','Fontname','Times new Roman'); 
legend('boxoff')

title('Treat intervals shorter than 4 months between two WS as a single event','FontName','Arial','FontSize',14,'fontweight','bold')
%%%%%%%%%%plot2




%%%%%%%%%%plot3
subplot_tight(3,1,3,[0.06,0.06])

%%%%%%   
%%3: detect the WS evetns ignoring the min non-WS events or min WS events
%% In whole peried, find if the year 2015 is within the WS duaration and dates
run_lengths = diff([0, temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
start_index = find(run_lengths == 1); % 
end_index = find(run_lengths == -1) - 1; 
durations = end_index - start_index + 1; 
%%%%%%

h1=area(x,bsxfun(@max, wsi, 0.4),0.4,'LineStyle','none','ShowBaseLine','off','FaceColor','r','EdgeColor','none','FaceAlpha',0.4);%
hold on
h2=area(x,bsxfun(@min, wsi, 0.4),0.4,'LineStyle','none','ShowBaseLine','off','FaceColor','b','EdgeColor','none','FaceAlpha',0.4);%
hold on
xline(2014,'k--','2014','Linewidth',1,'LabelOrientation','horizontal','LabelVerticalAlignment','middle','LabelHorizontalAlignment', 'center','Fontsize',12,'fontweight','bold');
hold on
yline(0.4)

yt=get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',arrayfun(@(x) num2str(x,'%1.1f'),yt,'UniformOutput',false), 'Fontname', 'Times new Roman', 'Fontsize',12, 'Fontname', 'Times new Roman', 'Fontsize',12)
set(gca,'YMinorTick','on')
set(gca,'fontsize',12,'fontweight','bold','FontName','Times new Roman');  
ylabel ('WSI','fontsize',12,'fontweight','bold','FontName','Times new Roman'); 


text('Units','normalized','Position',[0.00, 1.05],'String','(c)','Color','k','fontsize',14,'fontweight','bold','rotation',0);
%set(gca,'XTickLabel', []);  

strs=['Events: ',num2str(length(durations)),'   (< 6 month: ',num2str(length(durations(durations<6))),'; >= 6 months: ',num2str(length(durations(durations>3))),')']
text('Units','normalized','Position',[0.05, 0.90],'String',strs,'rotation',0,'Color','k','FontName','Arial','FontSize',12,'fontweight','bold');


lg=legend([h1,h2],{'WS','Non-WS'},'Location', 'northeast','NumColumnsMode','manual','NumColumns',2,'fontsize',12,'fontweight','bold','Fontname','Times new Roman'); 
legend('boxoff')

title('Ignore water scarcity events lasting less than 6 months','FontName','Arial','FontSize',14,'fontweight','bold')
%%%%%%%%%%plot3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exportgraphics(gcf,'D:\ws_fluctuation\figures\Fig.sx.diagram_duration.jpg','Resolution',350);% no white 
%%%%%%%%%%%%%%%%%%%%%
