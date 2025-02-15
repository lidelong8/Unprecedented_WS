function [max_i_h,max_i_f,max_d_h,max_d_f,max_f_h,max_f_f,TFE_i,TFE_d,TFE_f,TFE_w] = cal_TFE_mon_demo(wsi,window,freq_window,threshold,threshold_f,min_duration,max_gap,min_month,hist_end)
% wsi: monthly water shortage indices (1850-2100), with dimensions 360x720x3012
% pop: yearly population in each grid from 2015-2100
% areas:Areas of each grid 360x720
% window: smooth windows for WSI (i.e., 20 months)
% freq_window: windows to calculate frequency (i.e., 20 months)
% threshold: the water scarcity thresholds (i.e., 0.4 months)
% threshold_f: at least 1 times in 20 months if no WS in historical period.
% min_month: at least 1 months in each year.
% min_duration: minimum duration of consecutive exceedance in months (i.e., 6 months)
% max_gap: the non-WS events are seen as WS shorter than three months (i.e., 3 months)
% hist_end: index of last month of historical period (must be less than or equal to (i.e., 2014-1850+1)*12)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate regional average water scarcity indices
[row,future_end]=size(wsi);
future_start = hist_end+1;
future_time = future_start:future_end; % 
all_time=1:future_end;
yr=length(future_time)/12;


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
wsi_ori=wsi;
wss=wsi;%water shortage state, first cal ws;culate to save time!
wss(wss <= threshold) = 0; % set non-water scarcity to 0
wss(wss >  threshold) = 1; % set water scarcity to 1
wss_ori=wss;
%%%%cal ES frequency,here do not consider the short periods of non-WS months between lengthy periods of WS months
temp = wss; % whole period
transitions = diff([temp,temp(:,end)],1,2); %diff(X,n,dim)
transitions(transitions~=1)=0;
freq=filter(ones(1,freq_window)/freq_window,1,transitions,[],2);
freq=movmean(freq,[window-1 0],2,"omitnan"); %window:2+0+1
%%%%cal WS frequency

%%smooth
wsi_s=movmean(wsi,[window-1 0],2,"omitnan"); %window:2+0+1
wsi_s(wsi_s<=0)=0;
wss=wsi_s;
wss(wss <= threshold) = 0; % set non-water scarcity to 0
wss(wss >  threshold) = 1; % set water scarcity to 1
wss_s=wss;
%%%%%%%%%smooth



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:row %
       %%%%%%%%%%%%%%% 
        %%%%%%%%%%%%%%% 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%cal Unprecedented intensity
        %% intensity
         %%%%%%%%%%%%%%%his
        TFE_i(i)=2101;
        TFE_w(i)=2101;
        temp = wsi_s(i,1:hist_end); % only account for historical period
        %%historical maximum 
        [thres_i,loc_i_h]=nanmax(temp,[],2);
        max_i_h(i)=thres_i;
        if thres_i<threshold
           thres_i=single(threshold);%if there is water scarcity (<0.4), set 0
        end
        %%historical maximum base period
        if loc_i_h>window
          base_series=temp(1,loc_i_h-window+1:loc_i_h);
        else 
          base_series=temp(1,1:loc_i_h);
        end
        %% intensity
        
        %%TFE of intensity
        future_time = future_start:future_end; % future time
        
        
        %%%%%%%%%%%conventional WS
        %%Traditional WS
        temp = wss_s(i,future_time); %  future time
        ws=squeeze(nansum(reshape(temp,12,[]),1));%save
        %%
        tmp=ws;
        tmp(tmp<min_month)=0;
        tmp(tmp>=min_month)=1;
        %%  
        %%
        idw=find(ws>=min_month);%have lager than his, only need last for longer
        if idw
          TFE_w(i)=2015+idw(1)-1;
        end
        %%%%%%%%%%%finish conventional WS
        
        %%Unprecedented intensity of WS
        temp = wsi_s(i,future_time); %  future time
        wsi_s_f = wsi_s(i,future_time); %  future time
        max_i_f(i)=nanmax(temp,[],2);
        %%
        temp(temp <= thres_i) = 0; % if there is unprecedented water scarcity, set to 0
        temp(temp >  thres_i) = 1; % set unprecedented water scarcity to 1
        %%
        %%calculate the duration of UIWS
        run_lengths = diff([0, temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
        start_index = find(run_lengths == 1); % 
        end_index = find(run_lengths == -1) - 1; 
        durations = end_index - start_index + 1; 
        %%mask the event less than min_duration
        id_i_dur=find(durations>=min_duration);%at leat 6 months
        widths=durations(id_i_dur);
        %
        if length(widths)>0
        indices=start_index(id_i_dur)+future_start-1;%convert to all series
        %find origional data
        %%%the index of UIWS at leat 6 months
        id_mask = [];
        for k = 1:length(indices)
            id_mask = [id_mask, indices(k):indices(k)+widths(k)-1];
        end
        %%Extract future series  
        temp0=zeros(1,length(future_time));
        for k = 1:length(id_mask)
            loc_i_f=id_mask(k);
            future_series=wsi_s(i,loc_i_f-window+1:loc_i_f);
            %%significant test:  two-sample Kolmogorov-Smirnov test
            [h,p]=my_kstest2(base_series',future_series');%null hypothesis vectors x1 and x2 are from the same continuous, 0:same; 1,significant difference. 
            if h==1 %significant difference
              %%Repalce origional data
              temp0(1,loc_i_f-future_start+1)=1;
            end
        end  %end for k
        
        uiws=squeeze(nansum(reshape(temp0,12,[]),1));
        %%
        tmp=uiws;
        tmp(tmp<min_month)=0;
        tmp(tmp>=min_month)=1;
        %%
        %%UIWS            
        %
        idi=find(uiws>=min_month);%have lager than his, only need last for longer                      
        if idi
          TFE_i(i)=2015+idi(1)-1;
        end  
             
        end %end if length(widths)>0
        %%Unprecedented intensity
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        %%Unprecedented duration   
        TFE_d(i)=2101;
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
        id_min_dur2=find(durations2<min_duration);
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
       %%%%%%%%%%%%%%% 
        %%%%%%
      if length(durations)>0 %if no WS duration, skip...
        %% In whole peried, find if the year 2015 is within the WS duaration and dates
        %% Detect if the 2015 is within the duration of WS events
        id_d=all_time(find(temp==1));% find the dates of WS events
        if temp(hist_end)==1 %if year 2015 is with the duration of WS events, very complex
          %%Detect the start and end date of WS events corss 2015
          id_2015=hist_end;
          id_st=2015; 
          id_ed=2015; 
          c=id_2015;
          while c<=future_end&temp(c)
            id_ed=all_time(c);
            c=c+1;
          end
          %%
          c=id_2015;
          while c>=2&temp(c)
            c=c-1;
            id_st=all_time(c);
          end  
          %%Detect the start and end date of WS events corss 2015    
          max_d2=id_ed-id_st+1;
          max_d21=id_2015-id_st+1; 
          max_d22=id_ed-id_2015;
          max_d_h(i)=max_d2;
          %%%%%%%%%%%%
          %%1：id_st-1,id_st:id_ed,id_ed+1:end 
           temp_h1=temp(1:id_st-1);% only account for historical first period
           run_lengths = diff([0, temp_h1, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
           start_index = find(run_lengths == 1); % 
           end_index = find(run_lengths == -1) - 1; 
           durations = end_index - start_index + 1;     
                  
           if length(durations)>0 %Or wrong if no value
            [max_d1,id_max_dur_h] = nanmax(durations,[],2); % find historical maximum duration
            if max_d1>max_d2
               max_d_h(i)=max_d1;
            end
           else
             max_d1=min_duration;%%very important, or keep the origonal value
           end
           %%%%%%%%%
           temp_f2=temp(1,id_ed+1:end);% only account for future second period
           run_lengths = diff([0, temp_f2, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
           start_index = find(run_lengths == 1); % 
           end_index = find(run_lengths == -1) - 1; 
           durations = end_index - start_index + 1;    
           if length(durations)>0 %Or wrong if no value
            [max_d3,~] = nanmax(durations,[],2); % find historical maximum duration
           else
             max_d3=min_duration;
           end   
           
           %%muti scnerios
           [~,id] = max([max_d1 max_d2 max_d3]);
           if id==1 %No TFE in the future
              TFE_d(i)=2101;%denote before 2015
              max_d_f(i)=max([max_d2,max_d3]);
           elseif id==2 % the maximum duration of WS event include 2015,before 2015, ignore these events
              if max_d21>max_d1
                TFE_d(i)=2014;%occurr before 2015
              else %max_d21<max_d1 %will occurr in the future
                temp_f2=temp(id_ed+1:end);% only account for future second period
                run_lengths = diff([0,temp_f2,0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
                start_index = find(run_lengths == 1); % 
                end_index = find(run_lengths == -1) - 1; 
                durations = end_index - start_index + 1;    
                id_min_dur=find(durations>max_d1);
                widths=durations(id_min_dur);
                indices=start_index(id_min_dur);
                %% Repalce origional data
                temp0=zeros(1,length(temp_f2));
                for k = 1:length(indices)
                   temp0(1,indices(k)+max_d1:indices(k)+widths(k)-1)=1;
                end
                temp0(id_ed+1:end)=0; %%the current is the maximum duration,after this no UDWS
                temp_f1 =temp(id_2015+1:id_ed);
                temp_f1(1:max_d1-max_d21+1)= 0;%%this is different form if !!!!
                udws=squeeze(nansum(reshape([temp_f1 temp0],12,[]),1));
                %%
                tmp=udws;
                tmp(tmp<min_month)=0;
                tmp(tmp>=min_month)=1;
                %%
                max_d_f(i)=max([max_d2,max_d3]);
                %%
                idd=find(udws>=min_month);%have lager than his, only need last for longer
                if idd
                TFE_d(i)=2015+idd(1)-1;
                end
                %%Unprecedented durations             
            end      
           else %id==3 % the maximum WS event occurr after 2015
                temp_f2=temp(id_ed+1:end);% only account for future second period
                run_lengths = diff([0,temp_f2,0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
                start_index = find(run_lengths == 1); % 
                end_index = find(run_lengths == -1) - 1; 
                durations = end_index - start_index + 1;   
             if max_d1<max_d2     
                id_min_dur=find(durations>max_d2);
                widths=durations(id_min_dur);
                indices=start_index(id_min_dur);
                %% Repalce origional data
                temp0=zeros(1,length(temp_f2));
                temp1=zeros(1,length(all_time));
                for k = 1:length(indices)
                   temp0(1,indices(k)+max_d2:indices(k)+widths(k)-1)=1;
                end
                temp_f1 =temp1(id_2015+1:id_ed);
             else %max_d1>max_d2    
                id_min_dur=find(durations>max_d1);
                widths=durations(id_min_dur);
                indices=start_index(id_min_dur);
                %% Repalce origional data
                temp0=zeros(1,length(temp_f2));
                temp1=zeros(1,length(all_time));
                for k = 1:length(indices)
                   temp0(1,indices(k)+max_d1:indices(k)+widths(k)-1)=1;
                end
                temp_f1 =temp1(id_2015+1:id_ed);  
             end
            
                udws=squeeze(nansum(reshape([temp_f1 temp0],12,[]),1));
                %%
                tmp=udws;
                tmp(tmp<min_month)=0;
                tmp(tmp>=min_month)=1;
                %%
                max_d_f(i)=max([max_d2,max_d3]);
                %%
                idd=find(udws>=min_month);%have lager than his, only need last for longer
                if idd
                TFE_d(i)=2015+idd(1)-1;
                end
                %%Unprecedented durations
            end %%%end muti scnerios
           end %end   if temp(hist_end)==1
              %%%%%%%%%%%%
              %%%%%%%%%%%%
        %%%%%%
        else %if not in ,just seperately excute
           temp_h=temp(1,1:hist_end);% only account for historical period
           run_lengths = diff([0, temp_h, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
           start_index = find(run_lengths == 1); % 
           end_index = find(run_lengths == -1) - 1; 
           durations = end_index - start_index + 1;         
           if length(durations)>0 %Or wrong if no value
            [thres,id_max_dur_h] = nanmax(durations,[],2); % find historical maximum duration
            if thres<min_duration
               max_d_h(i)=thres;
               thres=min_duration
            end
           else
             max_d_h(i)=0;
             thres=min_duration;
           end
           %%find future unprecedented duration of WS
           temp_f=temp(1,future_time);% only account for future period
           run_lengths = diff([0,temp_f,0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
           start_index = find(run_lengths == 1); % 
           end_index = find(run_lengths == -1) - 1; 
           durations = end_index - start_index + 1;   
           id_min_dur=find(durations>thres);
           widths=durations(id_min_dur);
           indices=start_index(id_min_dur);
           %%Repalce origional data
           temp0=zeros(1,length(future_time));
           for k = 1:length(indices)
             temp0(1,indices(k)+thres:indices(k)+widths(k)-1)=1;
           end
           udws=squeeze(nansum(reshape(temp0,12,[]),1));
           %%
           tmp=udws;
           tmp(tmp<min_month)=0;
           tmp(tmp>=min_month)=1;
           %%
          
           if length(durations)>0
             max_d_f(i)= nanmax(durations,[],2);
           end
           %%
           %%Unprecedented durations                      
           idd=find(udws>=min_month);%%%very important!!!
           if idd
             TFE_d(i)=2015+idd(1)-1;
           end   
        end %end if length(durations)>0  
        %%%%%%
        %%%%%%           
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%
        %%Unprecedented frenquency
        % Identify the maximum value during the historical baseline period
        TFE_f(i)=2101;
        temp = freq(i,1:hist_end); % only account for historical period
        [thres_f,loc_f_h]=nanmax(temp,[],2);
        max_f_h(i)=thres_f;
        %%historical maximum base period
        if loc_f_h>window
          base_series=temp(1,loc_f_h-window+1:loc_f_h);
        else 
          base_series=temp(1,1:loc_f_h);
        end
        
        %%historical maximum
        if thres_f<threshold_f
          thres_f=single(threshold_f);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%TFE of  frequency   
        temp=freq(i,future_time);
        freq_f=freq(i,future_time);
        max_f_f(i)=nanmax(temp,[],2);
        temp(temp <= thres_f) = 0; % if there is unprecedented water scarcity, set to 0
        temp(temp > thres_f) = 1; % set unprecedented water scarcity to 1
        %%      
        %%
        %%calculate the duration of UFWS
        run_lengths = diff([0, temp, 0]); % % calculate the transition form  non-water scarcity (0) to water scarcity (1)
        start_index = find(run_lengths == 1); % 
        end_index = find(run_lengths == -1) - 1; 
        durations = end_index - start_index + 1; 
        %
        %id_f_dur=find(durations>=min_duration);%at leat 6 months
        %widths=durations(id_f_dur);
        widths=durations;  %not consider the duaration
        %
        if length(widths)>0
        indices=start_index+future_start-1;%convert to all series
        %find origional data
        %%%the index of UIWS at leat 6 months
        id_mask = [];
        for k = 1:length(indices)
            id_mask = [id_mask, indices(k):indices(k)+widths(k)-1];
        end
        %%Extract future series  
        temp0=zeros(1,length(future_time));
        for k = 1:length(id_mask)
            loc_f_f=id_mask(k);
            future_series=freq(i,loc_f_f-window+1:loc_f_f);
            %%significant test:  two-sample Kolmogorov-Smirnov test
            [h,p]=my_kstest2(base_series',future_series');%null hypothesis vectors x1 and x2 are from the same continuous, 0:same; 1,significant difference. 
            if h==1 %significant difference
              temp0(1,loc_f_f-future_start+1)=1;
            end
        end  %end for k
         %%Unprecedented frequency of WS
        %
        ufws=squeeze(nansum(reshape(temp0,12,[]),1));
        %%
        tmp=ufws;
        tmp(tmp<min_month)=0;
        tmp(tmp>=min_month)=1;
        %%
        %%UFWS            
        %
        idf=find(ufws>=min_month);%have lager than his, only need last for longer                      
        if idf
          TFE_f(i)=2015+idf(1)-1;
        end        
        end %end if length(widths)>0
       %%Unprecedented frequency  
%%%%%%%%%%%%%%%  
end %end loops rows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
