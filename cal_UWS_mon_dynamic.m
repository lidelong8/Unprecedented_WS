function [A_i,A_d,A_f,A_w,P_i,P_d,P_f,P_w,mon_i,mon_d,mon_f,mon_w,max_i_h,max_i_f,max_d_h,max_d_f,max_f_h,max_f_f,TFE_i,TFE_d,TFE_f,TFE_w] = cal_UWS_mon_dynamic(atots,qtots,pop,areas,ave_window,freq_window,threshold,threshold_d,threshold_f,min_duration,max_gap,min_month,pre_industrial,hist_end)
%% the reference time is dynamic.
% wsi: monthly water shortage indices (1850-2100), with dimensions 360x720x3012
% pop: yearly population in each grid from 2015-2100
% areas:Areas of each grid 
% ave_window: smooth wiondows 
% freq_window: windows to extract time series of X1 and X2; and windows to calculate frequency
% threshold: ws threshold, eg.0.4
% threshold_d: if historical duration is less than one year, in the future, WS with at least 24 months will be identified as UDWS. 
% threshold_f: at least 5 times in 60 months.
% min_month: the minimum duration of consecutive exceedance in months (i.e., 5 months) when identifying UIWS and UFWS.
% min_duration: the minimum duration of an indenpendent WS (i.e., 3,4,5,6 months)
% max_gap: the non-WS events are seen as WS shorter than three months (i.e., 3,4,5,6 months)
% pre_industria: pre-industria periods 1850-1900
% hist_end: index of last month of historical period (must be less than or equal to 2014)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculate regional average water scarcity indices
[row,future_end0]=size(qtots);
hist_end0=hist_end;
future_start0 = hist_end0+1;
future_time0 = future_start0:future_end0; % 
all_time=1:future_end0;
yr=length(future_time0)/12;
A=areas;
P=pop(:,hist_end0/12+1:end);

A_i = zeros(row,yr); % 
A_d = zeros(row,yr); % 
A_f = zeros(row,yr); % 
A_w = zeros(row,yr); %

P_i = zeros(row,yr); % 
P_d = zeros(row,yr); % 
P_f = zeros(row,yr); %
P_w = zeros(row,yr); %

mon_i = zeros(row,yr); % 
mon_d = zeros(row,yr); % 
mon_f = zeros(row,yr); %
mon_w = zeros(row,yr); %

max_i_h= zeros(row,1); % 
max_d_h = zeros(row,1); % 
max_f_h = zeros(row,1); % 

max_i_f= zeros(row,1); % 
max_d_f = zeros(row,1); % 
max_f_f = zeros(row,1); % 

TFE_i= 2101*ones(row,1); 
TFE_d= 2101*ones(row,1); 
TFE_f= 2101*ones(row,1); 
TFE_w= 2101*ones(row,1); 
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%smooth
qtots(qtots<=1)=1;
atots(atots<=1)=1;
wsi=atots./qtots;
wsi_ori=wsi;
wss=wsi;%water shortage state, first cal ws;culate to save time!
wss(wss <= threshold) = 0; % set non-water scarcity to 0
wss(wss >  threshold) = 1; % set water scarcity to 1
wss_ori=wss;

%%smooth
wsi_s=movmean(wsi,[ave_window-1 0],2,"omitnan"); %ave_window:2+0+1
wsi_s(wsi_s<=0)=0;
wss=wsi_s;
wss(wss <= threshold) = 0; % set non-water scarcity to 0
wss(wss >  threshold) = 1; % set water scarcity to 1
wss_s=wss;
%%%%%%%%%smooth


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%loop for pixel...
for i = 1:row %
    %%%%%%%%%%%%%%%%%%%%%
    TFE_i(i)=2101;
    TFE_w(i)=2101;
    TFE_d(i)=2101;
    TFE_f(i)=2101;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%cal Traditional WS
    temp = wss_s(i,future_time0); %  future time
    if length(find(~isnan(temp)))/length(temp)<0.9 | all(temp==0);% if less than 90% have value or never WS, then skip...
      continue
    end
    ws=squeeze(nansum(reshape(temp,12,[]),1));%save
    %%
    tmp=ws;
    tmp(tmp<1)=0;
    tmp(tmp>=1)=1;
    %%
    mon_w(i,:)=ws;
    A_w(i,:)=A(i)*tmp;
    P_w(i,:)=P(i,:).*tmp;        
    %%
    idw=find(ws>=1);
    if idw
       TFE_w(i)=2015+idw(1)-1;%Time of first emergence  
    end
    %%%%%%%%%%%%%
    %%Traditional WS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%Unprecedented intensity
    %%%%%%%%%%%%%%%%%%%
    temp0=zeros(1,length(future_time0));
    %
    %% historical intensity
    temp = wsi_s(i,1:hist_end0); % only account for historical period
    Noise_i=nanstd(temp(1:pre_industrial),0,2);
    %Noise_i=0;
    N_fut=length(future_time0);  
    
    %%%%%To overcome the unnecessary loop, we skip the pixel with no UIWS
    temp_his = wsi_s(i,1:hist_end0); % only account for historical period
    temp_fut = wsi_s(i,future_time0); %  future time
    [thres_i,~]=nanmax(temp_his,[],2);
    thres_i=thres_i+1*Noise_i;% important
    if thres_i<threshold
        thres_i=single(threshold);%if there is water scarcity (<0.4), set 0
    end 
    temp_fut(temp_fut <= thres_i) = 0; % if there is unprecedented water scarcity, set to 0
    temp_fut(temp_fut >  thres_i) = 1; % set unprecedented water scarcity to 1  
    run_lengths = diff([0, temp_fut, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
    start_index = find(run_lengths == 1); % 
    end_index = find(run_lengths == -1) - 1; 
    durations = end_index - start_index + 1;  
    widths=durations(durations>=min_month);
    %%%%%To overcome the unnecessary loop, we skip the pixel with no UIWS
    
    if ~isempty(widths)
    %%%%%
    toe=0; 
    for j=1:N_fut-min_month
        %%Extract future series  
        hist_end=future_time0(j)-1;%dynamic historical time bin
        loc_i_f=future_time0(j);
        %%
        temp_his = wsi_s(i,1:hist_end); % only account for historical period
        
        %%historical maximum 
        [thres_i,loc_i_h]=nanmax(temp_his,[],2);
        thres_i=thres_i+1*Noise_i;% important
                
        if thres_i<threshold
           thres_i=single(threshold);%if there is water scarcity (<0.4), set 0
        end
        %%historical maximum base period
        if loc_i_h>freq_window
            base_series=temp(1,loc_i_h-freq_window+1:loc_i_h);
        else 
            base_series=temp(1,1:loc_i_h);
        end
        %% historical intensity
        
        %%TFE of intensity
        %%Unprecedented intensity of WS
        temp = wsi_s(i,all_time); %  all time     
        %%
        temp(temp <= thres_i) = 0; % if there is unprecedented water scarcity, set to 0
        temp(temp >  thres_i) = 1; % set unprecedented water scarcity to 1
        %%
       future_series=wsi_s(i,loc_i_f-freq_window+1:loc_i_f);
       %%significant test:  two-sample Kolmogorov-Smirnov test
       [h,~]=my_kstest2(base_series',future_series');%my_kstest2 null hypothesis vectors x1 and x2 are from the same continuous, 0:same; 1,significant difference. 
       if j<=N_fut-min_month
          loc_i_f_out=loc_i_f+min_month-1;
        else
          loc_i_f_out=N_fut;
       end
        
       if h==1 && all(temp(loc_i_f:loc_i_f_out)==1) %%  %significant difference and last at least 5 months 
           %%Repalce origional data
           temp0(1,j)=1;
           toe=toe+1; 
           if toe==1 %
              max_i_h(i)=thres_i;
              max_i_f(i)=wsi_s(1,loc_i_f);%only mark the wsi with the first UIWS
           end
       end
   end  %end for j

        uiws=squeeze(nansum(reshape(temp0,12,[]),1));
        %%
        tmp=uiws;
        tmp(tmp>=1)=1;
        %%
        mon_i(i,:)=uiws;
        A_i(i,:)=A(i)*tmp;
        P_i(i,:)=P(i,:).*tmp;
        
        %%UIWS            
        %
        idi=find(uiws>=1);%               
        if idi
          TFE_i(i)=2015+idi(1)-1;%Time of first emergence  
        end  
        
   end %if length(widths)>0 
   %%%%%%%%%%%%%%%%      
   %%Unprecedented intensity
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%Unprecedented duration 
   %%%%%%%%%%%%%%%%  
   %%%%%%%%%1: fill the gap and identify the independent evets
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
   temp(id_min1)= ~temp(id_min1);%% Here the short periods of non-WS months between lengthy periods of WS months were considered WS due to the pooling effect.
   %%%%%%
   %%%%%%
   %%%%%%%%%1:  fill the gap and identify the independent evets
   %%%%%%
   %%%%%%
   %%2:discard the WS duration between two Non-WS evetns
   run_lengths2 = diff([0, temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
   start_index2 = find(run_lengths2 == 1); %
   end_index2 = find(run_lengths2 == -1) - 1;
   durations2 = end_index2 - start_index2 + 1;
   %find origional data
   id_min_dur2=find(durations2<min_duration);% ignore the WS events less than min_duration. eg. 5 months
   widths2=durations2(id_min_dur2);
   indices2=start_index2(id_min_dur2);
   %%
   id_min2 = [];
   for k = 1:length(indices2)
       id_min2 = [id_min2, indices2(k):indices2(k)+widths2(k)-1];
   end
   %%Repalce origional data
   temp(id_min2)= ~temp(id_min2);
   %%2: The duration of a negligible WS period shorter than 12-months
   %%%%%%
   
   %%%%%%
   %%%%%%3: detect the WS evetns ignoring the min non-WS events or min WS events
   run_lengths = diff([0, temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
   start_index = find(run_lengths == 1); %
   end_index = find(run_lengths == -1) - 1;
   durations = end_index - start_index + 1;
   %%%%%%%%%%%%%%%
   %%%%%%
   %%%%%%
   %%%%4: Assign the length of water scarcity events to each period
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
   if ~isempty(idx_ones)
       length_d(idx_ones) = segment_lengths(group_ids); % Allocate length to each original position  %下标赋值维度不匹配 当转为mex函数时,且idx_ones为空报错，在普通函数正常
   end
   %
   %%eg. temp = [0, 1, 1, 0, 1, 1, 1, 0, 1, 0];% 1 denote WS
   %% length_d= [0, 2, 2, 0, 3, 3, 3, 0, 1, 0];
   %%%% 4:ssign the length of water scarcity events to each period
   %%%%%%%%%%%%%%%
   %%%%%%

 
   
   %%%%%%%%% loop to identify the UDWS evets
   if ~isempty(durations) %if have WS duration, then excute...
   temp0=zeros(1,length(future_time0));
   length_d0=length_d(1,future_time0);
   toe=0; 
   for j=1:length(future_time0)
        %%Extract future series  
        hist_end=future_time0(j)-1;%%dynamic historical time bin
        loc_d_f=future_time0(j);       
        
        %%%%%%
       %%%%%%%%%%%%%%%%%%% 
        %% In whole peried, find if the target month is within the WS duaration
        if temp(loc_d_f)==1 %if the target month is within the duration of WS events, just skip....
          id_st=loc_d_f; 
          c=loc_d_f;
          while c>1 && temp(c-1)
            c=c-1;
            id_st=all_time(c);
          end  
          %%Detect the start and end date of WS events corss loc_d_f  
          max_d21=loc_d_f-id_st+1; 
          %%%%%%%%%%%%
          %%1：id_st-1,id_st:id_ed,id_ed+1:end 
           temp_h1=temp(1:id_st-1);% only account for historical first period
           run_lengths = diff([0, temp_h1, 0]); % % calculate the transition form non-water scarcity (0) to water scarcity (1)
           start_index = find(run_lengths == 1); % 
           end_index = find(run_lengths == -1) - 1; 
           durations = end_index - start_index + 1;     
           %%%%%%
           %%%%%%%%%%%%
           %%  
               
           if ~isempty(durations) %
               [max_d1,~] = nanmax(durations,[],2); % find historical maximum duration
               if max_d1 < max_d21 && max_d21 >= threshold_d %exceed the longest duration of all previous WS event and at least last for 12 months if historical maximum duration less than 24 months.
                 temp0(1,j)=1;
                 toe=toe+1;
                 if toe==1 %save the durtion of TFE_d
                     max_d_h(i)=max([max_d1;max_d21]);
                     max_d_f(i)=length_d0(1,j); %only mark the wsi with the first UDWS
                 end %end if toe==1
               end
           else
               if max_d21 >= threshold_d %No historical WS events, at least last for 12 months.
                   temp0(1,j)=1;
                   toe=toe+1;
                   if toe==1 %save the durtion of TFE_d
                       max_d_h(i)=max_d21;
                       max_d_f(i)=length_d0(1,j);%only mark the duration with the first UDWS
                   end %end if toe==1
               end           
           end %length(durations)>0
           %%%%%%%%%
         end %end if within WS,   temp(loc_d_f)==1
   end %end for j
           %%%%%%%%%%%%%
        %%%%%%
           
           %%%%%%%%%%%%%TFE_d 
           udws=squeeze(nansum(reshape(temp0,12,[]),1));
           tmp=udws;
           tmp(tmp>=1)=1;
                %%
           mon_d(i,:)=udws;
           A_d(i,:)=A(i)*tmp;
           P_d(i,:)=P(i,:).*tmp;
           %%%%Unprecedented durations        
           idd=find(udws>=1); 
           if idd
                TFE_d(i)=2015+idd(1)-1;%Time of first emergence  
           end 
           %%%%%%
    end %end if have WS
    %%Unprecedented duration  
    %%%%%%       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%
    %%Unprecedented frenquency
    
    %%%%%%%%%%%%%%%%%%% fill the gap and identify the independent evets
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
    %%%%%%
    %%2:ignore the WS duration < 5 months between two Non-WS evetns
    run_lengths2 = diff([0, temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
    start_index2 = find(run_lengths2 == 1); %
    end_index2 = find(run_lengths2 == -1) - 1;
    durations2 = end_index2 - start_index2 + 1;
    %find origional data
    id_min_dur2=find(durations2<min_duration);
    widths2=durations2(id_min_dur2);
    indices2=start_index2(id_min_dur2);
    %%
    id_min2 = [];
    for k = 1:length(indices2)
        id_min2 = [id_min2, indices2(k):indices2(k)+widths2(k)-1];
    end
    %%Repalce origional data
    temp(id_min2)= ~temp(id_min2);% %% The duration of a negligible WS period shorter than x-months
    %%%%%%%%%%%%%%%%%%% fill the gap and identify the independent evets
    %%%%%%
    %%%%%%
    %%3: cal ES frequency,here consider the short periods of non-WS months between lengthy periods of WS months
    transitions = diff([temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
    transitions(transitions~=1)=0;
    %a=1:10;filter(ones(1,3)/3,1,a,[],2),movmean(a,[3-1 0],2,"omitnan")
    freq=filter(ones(1,freq_window)/freq_window,1,transitions,[],2); % the frequency of local and previous ave_window-1 values
    freq(1,1:freq_window-1)=nan; %%the first ave_window-1 value is nan
    freq=movmean(freq,[ave_window-1 0],2,"omitnan"); %ave_window:2+0+1  the mean of local and previous ave_window-1 values
    %%%%cal WS frequency
    %%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%% Identify the maximum value during the historical baseline period
    temp0=zeros(1,length(future_time0));
    %
    %% historical frequency
    temp = freq(1,1:hist_end); % only account for historical period
    Noise_f=nanstd(temp(1:pre_industrial));
    %Noise_f=0;
    N_fut=length(future_time0);
    
     %%%%%To overcome the unnecessary loop, we skip the pixel with no UFWS
    temp_his =  freq(1,1:hist_end0); % only account for historical period
    temp_fut =  freq(1,future_time0); %  future time
    [thres_f,~]=nanmax(temp_his,[],2);
    thres_f=thres_f+1*Noise_f;% important
    if thres_f<threshold_f
        thres_f=single(threshold_f);%
    end 
    temp_fut(temp_fut <= thres_f) = 0; % if there is unprecedented water scarcity, set to 0
    temp_fut(temp_fut >  thres_f) = 1; % set unprecedented water scarcity to 1  
    run_lengths = diff([0, temp_fut, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
    start_index = find(run_lengths == 1); % 
    end_index = find(run_lengths == -1) - 1; 
    durations = end_index - start_index + 1;  
    id_min_dur=durations>=min_month;
    widths=durations(id_min_dur);
    %%%%%To overcome the unnecessary loop, we skip the pixel with no UFWS
    if ~isempty(widths)
    toe=0; 
    for j=1:N_fut-min_month
        %%Extract future series  
        future_time=future_time0(j:end);
        hist_end=future_time0(j)-1;%dynamic historical time bin
        loc_f_f=future_time0(j);
    
        %%historical maximum 
        temp_his = freq(1,1:hist_end); % only account for historical period
        [thres_f,~]=nanmax(temp_his,[],2);
        thres_f=thres_f+1*Noise_f;% important
                
        %%historical maximum
        if thres_f<threshold_f
          thres_f=single(threshold_f);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%TFE of  frequency   
        temp=freq(1,all_time);
        freq_f=freq(1,future_time);

        %%
        temp(temp <= thres_f) = 0; % if there is unprecedented water scarcity, set to 0
        temp(temp > thres_f) = 1; % set unprecedented water scarcity to 1
        %%   
        if j<=N_fut-min_month
          loc_f_f_out=loc_f_f+min_month-1;
        else
          loc_f_f_out=N_fut;
        end
        
  
        if all(temp(loc_f_f:loc_f_f_out)==1) %%  %significant difference and last at least 5 months 
              temp0(1,j)=1;
              toe=toe+1; 
              if toe==1
                  max_f_h(i)=thres_f;
                  max_f_f(i)=freq(1,loc_f_f);%only mark the frequency with the first UFWS
              end
        end

    end  %end for j
         %%Unprecedented frequency of WS
        %
        ufws=squeeze(nansum(reshape(temp0,12,[]),1));
        %%
        tmp=ufws;
        tmp(tmp>=1)=1;
        %%
        mon_f(i,:)=ufws;
        A_f(i,:)=A(i)*tmp;
        P_f(i,:)=P(i,:).*tmp;
        %%UFWS            
        %
        idf=find(ufws>=1);                    
        if idf
           TFE_f(i)=2015+idf(1)-1;%Time of first emergence  
        end    
        
    end % if length(widths)>0    
   %%Unprecedented frequency  
%%%%%%%%%%%%%%%  
end %end loops rows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
