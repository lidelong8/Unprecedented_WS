function [A_i,A_d,A_f,A_w,P_i,P_d,P_f,P_w,mon_i,mon_d,mon_f,mon_w,max_i_h,max_i_f,max_d_h,max_d_f,max_f_h,max_f_f,TFE_i,TFE_d,TFE_f,TFE_w] = cal_UWS_mon_fixed(atots,qtots,pop,areas,ave_window,freq_window,threshold,threshold_d,threshold_f,min_duration,max_gap,min_month,pre_industrial,hist_end)
%% the reference ime is fixed in year of 2014.
% wsi: monthly water shortage indices (1850-2100), with dimensions 360x720x3012
% pop: yearly population in each grid from 2015-2100
% areas:Areas of each grid 360x720
% ave_window: windows to extract time series or smooth
% freq_window: windows to calculate frequency of WS
% threshold: the water scarcity thresholds (i.e., 0.4)
% threshold_d: if historical duration is less than one year, in the future, WS with at least 12 months will be identified as UDWS.
% threshold_f: at least 1 times in 60 months.
% min_month: minimum duration of consecutive exceedance in months (i.e., 5 months)
% min_duration: minimum duration of a single WS events (i.e., 5 months)
% max_gap: the non-WS events are seen as WS shorter than three months (i.e., 3 months)
% pre_industria: pre-industria periods 1850-1900
% hist_end: index of last month of historical period (must be less than or equal to 2014)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate regional average water scarcity indices
[row,future_end]=size(qtots);
future_start = hist_end+1;
future_time = future_start:future_end; %
all_time=1:future_end;
yr=length(future_time)/12;
A=areas;
P=pop(:,hist_end/12+1:end);

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
%wsi_ori=wsi;
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:row %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%cal Traditional WS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%cal Traditional WS
    temp = wss_s(i,future_time); %  future time
    if length(find(~isnan(temp)))/length(temp)<0.9 || all(temp==0)% if less than 90% have value or never WS, then skip...
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
        TFE_w(i)=1850+hist_end/12+idw(1)-1;%Time of first emergence
    end
    %%%%%%%%%%%%%
    %%Traditional WS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%cal Unprecedented intensity
    %% intensity
    %%%%%%%%%%%%%%%his
    temp = wsi_s(i,1:hist_end); % only account for historical period
    Noise_i=nanstd(temp(1:pre_industrial),0,2);
    %%historical maximum
    [thres_i,loc_i_h]=nanmax(temp,[],2);
    thres_i=thres_i+1*Noise_i;% important
    max_i_h(i)=thres_i;
    
    if thres_i<threshold
        thres_i=single(threshold);%if there is water scarcity (<0.4), set 0
    end
    %%historical maximum base period
    if loc_i_h>freq_window
        base_series=temp(1,loc_i_h-freq_window+1:loc_i_h);
    else
        base_series=temp(1,1:loc_i_h);
    end
    %% intensity
    
    %%TFE of intensity
    %%Unprecedented intensity of WS
    temp = wsi_s(i,future_time); %  future time
    wsi_fut = wsi_s(i,future_time); %  future time
    
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
    if ~isempty(widths)
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
            future_series=wsi_s(i,loc_i_f-freq_window+1:loc_i_f);
            %%significant test:  two-sample Kolmogorov-Smirnov test
            [h,~]=my_kstest2(base_series',future_series');%null hypothesis vectors x1 and x2 are from the same continuous, 0:same; 1,significant difference.
            if h==1 %significant difference
                %%Repalce origional data
                temp0(1,loc_i_f-future_start+1)=1;
                toe=toe+1;
                if toe==1
                    max_i_f(i)=wsi_fut(1,loc_i_f-future_start+1);
                end
            end
        end  %end for k
        
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
        idi=find(uiws>=1);
        if idi
            TFE_i(i)=1850+hist_end/12+idi(1)-1;%Time of first emergence
        end
        
    end %end if length(widths)>0
    %%Unprecedented intensity
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
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
    %%%%%%
    %%2:discard the WS duration between two Non-WS evetns
    run_lengths2 = diff([0, temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
    start_index2 = find(run_lengths2 == 1); %
    end_index2 = find(run_lengths2 == -1) - 1;
    durations2 = end_index2 - start_index2 + 1;
    %find origional data
    id_min_dur2=find(durations2<min_duration);% ignore the events less than one years. eg. 12 mongths
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
    %%3: detect the WS evetns ignoring the min non-WS events or min WS events
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
    if ~isempty(idx_ones)
        length_d(idx_ones) = segment_lengths(group_ids); % Allocate length to each original position  %下标赋值维度不匹配 当转为mex函数时,且idx_ones为空报错，在普通函数正常
    end
    %
    %%eg. temp = [0, 1, 1, 0, 1, 1, 1, 0, 1, 0];% 1 denote WS
    %% length_d= [0, 2, 2, 0, 3, 3, 3, 0, 1, 0];
    %%4: %%Assign the length of water scarcity events to each period
    %%%%%%%%%%%%%%%
    %%%%%%
    if ~isempty(durations) %if no WS duration, skip...
        %% find if the year 2015 is within the WS duaration and dates
        %% In whole peried, Detect if the 2015 is within the duration of WS events
        
        
        if temp(hist_end)==1 %if year 2015 is with the duration of WS events, it is very complex
            %%Detect the start and end date of WS events corss 2015
            id_st=hist_end; %Not fully defined in some execution paths
            id_ed=hist_end; %Not fully defined in some execution paths
            loc=hist_end;
            while loc < future_end && temp(loc+1)  
                loc=loc+1;
                id_ed=all_time(loc);
            end
            %%
            loc=hist_end;
            while loc > 1 && temp(loc-1)
                loc=loc-1;
                id_st=all_time(loc);
            end
            %%Detect the start and end date of WS events corss 2015
            max_d2=id_ed-id_st+1;
            max_d21=hist_end-id_st+1;
            % max_d22=id_ed-hist_end;
            %%%%%%%%%%%%
            %%1：id_st-1,id_st:id_ed,id_ed+1:end
            temp_h1=temp(1:id_st-1);% only account for historical first period
            run_lengths = diff([0, temp_h1, 0]); % % calculate the transition form non-water scarcity (0) to water scarcity (1)
            start_index = find(run_lengths == 1); %
            end_index = find(run_lengths == -1) - 1;
            durations = end_index - start_index + 1;
            
            %%%%%%%% find historical maximum duration
            if ~isempty(durations) %
                [max_d1,~] = nanmax(durations,[],2); % find historical maximum duration of all previous WS event
                thres_d = max([max_d1 max_d21]);
            else
                thres_d=max_d21;
            end
            if thres_d<threshold_d
                thres_d=threshold_d;%%very important, or use the original value
            end
            max_d_h(i)=thres_d;
            %%%%%%%%%% find historical maximum duration
            
            %%%%%%%%%% judge if there will have WS in the future
            temp_f2=temp(1,id_ed+1:end);% only account for future second period
            run_lengths = diff([0, temp_f2, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
            start_index = find(run_lengths == 1); %
            end_index = find(run_lengths == -1) - 1;
            durations = end_index - start_index + 1;
            if ~isempty(durations) %Or no WS if no value
                [max_d3,~] = nanmax(durations,[],2); % find future maximum duration
            else
                max_d3=0;% no WS
            end
            %%%%%%%%%%
            
            %%muti scenarios
            [~,id] = max([thres_d max_d2 max_d3]);
            if id==1 %No TFE_d in the future
                TFE_d(i)=2101;%after 2100
                max_d_f(i)=max([max_d2,max_d3]);
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
                        max_d_f(i)=length_f(1,indices(k)+thres_d);
                    end
                end
       
                temp_f1 =temp0(1,hist_end-id_st+1:end);%% need to delete the second historical period
                udws=squeeze(nansum(reshape(temp_f1,12,[]),1));
                
                %%
                tmp=udws;
                tmp(tmp>=1)=1;
                %%
                
                mon_d(i,:)=udws;
                A_d(i,:)=A(i)*tmp;
                P_d(i,:)=P(i,:).*tmp;
                %%
                idd=find(udws>=1);%%
                if idd
                    TFE_d(i)=1850+hist_end/12+idd(1)-1;%Time of first emergence
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
            if ~isempty(durations) %
                [thres_d,~] = nanmax(durations,[],2); % find historical maximum duration
            else
                thres_d=threshold_d;
            end
            if thres_d<threshold_d
                thres_d=threshold_d;
            end
            max_d_h(i)=thres_d;
            
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
                    max_d_f(i)=length_f(1,indices(k)+thres_d);
                end
            end
            udws=squeeze(nansum(reshape(temp0,12,[]),1));
            %%
            tmp=udws;
            tmp(tmp>=1)=1;
            %%
            mon_d(i,:)=udws;
            A_d(i,:)=A(i)*tmp;
            P_d(i,:)=P(i,:).*tmp;
            %%
            %%Unprecedented durations
            idd=find(udws>=1);%
            if idd
                TFE_d(i)=1850+hist_end/12+idd(1)-1;%Time of first emergence
            end
            
        end %end   if temp(hist_end)==1
    end %end if length(durations)>0
    %%%%%%
    %%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
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
    id_min_dur2=find(durations2<min_duration);
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
    %%historical maximum
    [thres_f,~]=nanmax(temp,[],2);
    thres_f=thres_f+1*Noise_f;% important
    max_f_h(i)=thres_f;
    %%historical maximum base period
    %%historical maximum
    if thres_f<threshold_f
        thres_f=single(threshold_f);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%TFE of  frequency
    temp=freq(1,future_time);
    freq_f=freq(1,future_time);
    
    temp(temp <= thres_f) = 0; % if there is unprecedented water scarcity, set to 0
    temp(temp > thres_f) = 1; % set unprecedented water scarcity to 1
    %%
    %%
    %%calculate the duration of UFWS
    run_lengths = diff([0, temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
    start_index = find(run_lengths == 1); %
    end_index = find(run_lengths == -1) - 1;
    durations = end_index - start_index + 1;
    id_f_dur=find(durations>=min_month);%at leat x months
    widths=durations(id_f_dur);
    %widths=durations;  %not consider the duaration for frequency
    %
    if ~isempty(widths)
        indices=start_index(id_f_dur)+future_start-1;%convert to full series
        %find origional data
        %%%the index of UIWS at leat x months
        id_mask = [];
        for k = 1:length(indices)
            id_mask = [id_mask, indices(k):indices(k)+widths(k)-1];
        end
        %%Extract future series
        temp0=zeros(1,length(future_time));
        toe=0;
        for k = 1:length(id_mask)
            loc_f_f=id_mask(k);
            temp0(1,loc_f_f-future_start+1)=1;
            toe=toe+1;
            if toe==1
                max_f_f(i)=freq_f(1,loc_f_f-future_start+1);
            end
        end  %end for k
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
            TFE_f(i)=1850+hist_end/12+idf(1)-1;%Time of first emergence
        end
        
    end %end if length(widths)>0
    %%Unprecedented frequency
    %%%%%%%%%%%%%%%
end %end loops rows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
