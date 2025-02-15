

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
future_time = future_start:future_end; % 未来时期的年份


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
figure('unit','centimeters','Position',[2 2 25 25]);

%%%%%%%%%%plot1
subplot_tight(3,1,1,[0.04,0.04])

yyaxis left
h1=plot(x,y_ori,'c','Linewidth',1);
yline(0.4)
yt=get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',arrayfun(@(x) num2str(x,'%1.1f'),yt,'UniformOutput',false), 'Fontname', 'Times new Roman', 'Fontsize',12, 'Fontname', 'Times new Roman', 'Fontsize',12)
set(gca,'YMinorTick','on')
ylabel ('WSI','fontsize',12,'fontweight','bold','FontName','Times new Roman'); 
%%%%%%%%%%%%%%%%frequency
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
freq1=movmean(freq,[freq_window-1 0],2,"omitnan"); %window:2+0+1
freq1(1:ave_window-1)=nan;

yyaxis right
h2=plot(x(ave_window:end),freq(ave_window:end),'-r','Linewidth',1);
hold on
h3=plot(x(ave_window:end),freq1(ave_window:end),'-b','Linewidth',2);
yt=get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',arrayfun(@(x) num2str(x,'%1.2f'),yt,'UniformOutput',false), 'Fontname', 'Times new Roman', 'Fontsize',12, 'Fontname', 'Times new Roman', 'Fontsize',12)
set(gca,'YMinorTick','on')
set(gca,'XTickLabel', []);  
set(gca,'fontsize',12,'fontweight','bold','FontName','Times new Roman'); 
ylabel ('Frequency','fontsize',12,'fontweight','bold','FontName','Times new Roman'); 
text('Units','normalized','Position',[0.00, 1.05],'String','(a)','Color','k','fontsize',14,'fontweight','bold','rotation',0);
lg=legend([h1,h2,h3],{'Original Time Series of WSI','Original frequency of WSI','Smoothed frequency of WSI'},'Location', 'north','NumColumnsMode','manual','NumColumns',3,'fontsize',12,'fontweight','bold','Fontname','Times new Roman'); 
legend('boxoff')
title('Fluctuating frequently around the water scarcity threshold','FontName','Arial','FontSize',14,'fontweight','bold')
%%%%%%%%%%plot1


%%%%%%%%%%plot2
subplot_tight(3,1,2,[0.04,0.04])

yyaxis left
h1=plot(x,y3,'c','Linewidth',1);
yline(0.4)
yt=get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',arrayfun(@(x) num2str(x,'%1.1f'),yt,'UniformOutput',false), 'Fontname', 'Times new Roman', 'Fontsize',12, 'Fontname', 'Times new Roman', 'Fontsize',12)
set(gca,'YMinorTick','on')
ylabel ('WSI','fontsize',12,'fontweight','bold','FontName','Times new Roman'); 
%%%%%%%%%%%%%%%%frequency
wsi_ori=y3';
wss=wsi_ori;%water shortage state, first cal ws;culate to save time!
wss(wss <= threshold) = 0; % set non-water scarcity to 0
wss(wss >  threshold) = 1; % set water scarcity to 1
wss_ori=wss;
%%%%cal ES frequency,here do not consider the short periods of non-WS months between lengthy periods of WS months
temp = wss; % whole period
transitions = diff([temp,temp(:,end)],1,2); %diff(X,n,dim)
transitions(transitions~=1)=0;
freq=filter(ones(1,freq_window)/freq_window,1,transitions,[],2);
freq1=movmean(freq,[freq_window-1 0],2,"omitnan"); %window:2+0+1
freq1(1:ave_window-1)=nan;

yyaxis right
h2=plot(x(ave_window:end),freq(ave_window:end),'r','Linewidth',1);
hold on
h3=plot(x(ave_window:end),freq1(ave_window:end),'-b','Linewidth',2);
yt=get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',arrayfun(@(x) num2str(x,'%1.2f'),yt,'UniformOutput',false), 'Fontname', 'Times new Roman', 'Fontsize',12, 'Fontname', 'Times new Roman', 'Fontsize',12)
set(gca,'YMinorTick','on')
set(gca,'XTickLabel', []);  
set(gca,'fontsize',12,'fontweight','bold','FontName','Times new Roman'); 
ylabel ('Frequency','fontsize',12,'fontweight','bold','FontName','Times new Roman'); 

text('Units','normalized','Position',[0.00, 1.05],'String','(b)','Color','k','fontsize',14,'fontweight','bold','rotation',0);
lg=legend([h1,h2,h3],{'Original Time Series of WSI','Original frequency of WSI','Smoothed frequency of WSI'},'Location', 'north','NumColumnsMode','manual','NumColumns',3,'fontsize',12,'fontweight','bold','Fontname','Times new Roman');
legend('boxoff')
title('Almost always in a state of water safety','FontName','Arial','FontSize',14,'fontweight','bold')
%%%%%%%%%%plot2




%%%%%%%%%%plot3
subplot_tight(3,1,3,[0.04,0.04])

yyaxis left
h1=plot(x,y2,'c','Linewidth',1);
yline(0.4)
yt=get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',arrayfun(@(x) num2str(x,'%1.1f'),yt,'UniformOutput',false), 'Fontname', 'Times new Roman', 'Fontsize',12, 'Fontname', 'Times new Roman', 'Fontsize',12)
set(gca,'YMinorTick','on')

ylabel ('WSI','fontsize',12,'fontweight','bold','FontName','Times new Roman'); 
%%%%%%%%%%%%%%%%frequency
wsi_ori=y2';
wss=wsi_ori;%water shortage state, first cal ws;culate to save time!
wss(wss <= threshold) = 0; % set non-water scarcity to 0
wss(wss >  threshold) = 1; % set water scarcity to 1
wss_ori=wss;
%%%%cal ES frequency,here do not consider the short periods of non-WS months between lengthy periods of WS months
temp = wss; % whole period
transitions = diff([temp,temp(:,end)],1,2); %diff(X,n,dim)
transitions(transitions~=1)=0;
freq=filter(ones(1,freq_window)/freq_window,1,transitions,[],2);
freq1=movmean(freq,[freq_window-1 0],2,"omitnan"); %window:2+0+1
freq1(1:ave_window-1)=nan;

yyaxis right
h2=plot(x(ave_window:end),freq(ave_window:end),'r','Linewidth',1);
hold on
h3=plot(x(ave_window:end),freq1(ave_window:end),'-b','Linewidth',2);
yt=get(gca,'YTick');
set(gca,'YTick',yt,'YTickLabel',arrayfun(@(x) num2str(x,'%1.2f'),yt,'UniformOutput',false), 'Fontname', 'Times new Roman', 'Fontsize',12, 'Fontname', 'Times new Roman', 'Fontsize',12)
set(gca,'YMinorTick','on')
set(gca,'fontsize',12,'fontweight','bold','FontName','Times new Roman'); 
ylabel ('Frequency','fontsize',12,'fontweight','bold','FontName','Times new Roman'); 

text('Units','normalized','Position',[0.00, 1.05],'String','(c)','Color','k','fontsize',14,'fontweight','bold','rotation',0);

lg=legend([h1,h2,h3],{'Original Time Series of WSI','Original frequency of WSI','Smoothed frequency of WSI'},'Location', 'north','NumColumnsMode','manual','NumColumns',3,'fontsize',12,'fontweight','bold','Fontname','Times new Roman');
legend('boxoff')
title('Almost always in a state of water scarcity','FontName','Arial','FontSize',14,'fontweight','bold')
%%%%%%%%%%plot3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
exportgraphics(gcf,'D:\ws_fluctuation\figures\Fig.sx.diagram_frequency.jpg','Resolution',350);% no white 
%%%%%%%%%%%%%%%%%%%%%
