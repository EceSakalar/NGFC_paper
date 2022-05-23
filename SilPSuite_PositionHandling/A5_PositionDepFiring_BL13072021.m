clear all
%files to be considered
bs_name='B204b1';
bs_chname='Bsp204b_Ch';
%tr_idnt should be set to 0 if the first trial in the file is partial (for
%example the animal is at the reward zone) and set to 1 if it could be included
tr_idnt=1;
%tr_length (the track length in metres) & tr_bin (the bin size alos in
%metres)
tr_length=7.1;
tr_binsiz=0.05;
%movement threshold (m/s)
mov_thrs=0.04;
%THETA
%filter properties 
flt_ncoeff=512; flt_low=5; flt_high=12;
sh_pyr=[2];
ch_pyr=15;
ch_dwnsmpl=25;                                       %give the downsampling rate
ch_rate=10000; ch_rate=ch_rate/ch_dwnsmpl;

%UNITS DATA ENTRY
%unit spikes channels (in spike2 files)
un.chan=[130:193];
%B207c:[130:178]
%B207b:[130:156]
%B207a:[130:172]
%B206b:[130:175]
%B206a:[130:172]
%B204a:[130:169]
%B200b:[130:164];
%B200a:[130:190];
%B198b:[130:202];
%B198a:[66:100]
%B196a:[65:106];
%B193b:[21:50,17]
%B192d:[37:66,36];
%B192c:[36:63];
%B190c:[21:55,20];
%B188b:[37:59,36];
%B187b:[35:55,58];
%B186d:[35:53,55];
%B186c:[36:51];
%B184c:[36:75,77];
%B183a:[23:50,22]
%B182a:[23:55,22];
%B181g:[20]
%B179b:[23:52,22];
%B178g:[21:39,20];
%B177d:[20:37,19]
%unit numbers (as identifiers from clustering)
un.nums=[4,5,7,8,11,159,197,19,21,24,25,29,33,162,200,34,40,42,44,45,49,167,195,208,219,50,52,55,57,63,64,67,69,170,171,229,232,79,80,90,181,203,233,236,249,250,98,100,106,183,184,185,119,121,123,129,130,131,186,260,140,150,156,264];
%B207c:[2,18,141,143,171,173,177,181,19,32,33,34,147,148,189,195,40,43,44,52,54,218,219,61,62,69,153,234,71,84,159,238,245,255,99,101,103,105,106,107,108,114,121,123,127,128,133,137,162]
%B207b:[6,8,110,15,16,118,17,22,24,29,102,103,108,30,33,36,50,53,54,56,57,132,138,65,68,75,85]
%B207a:[28,249,30,60,250,340,64,65,66,71,72,73,76,81,91,92,105,108,325,135,145,238,239,319,320,161,177,240,243,302,307,181,202,203,269,295,204,207,215,221,222,230,242]
%B206b:[0,5,9,10,14,15,21,127,129,145,155,31,34,36,51,52,55,131,57,68,72,142,147,165,79,82,84,171,86,87,89,91,96,100,137,186,101,108,114,117,119,120,121,123,125,126]
%B206a:[2,5,188,20,23,26,32,34,190,36,37,48,49,55,59,173,186,239,245,248,63,71,86,89,93,97,107,108,109,182,257,110,121,132,133,177,134,138,156,158,166,167,171]
%B193b:[0,4,13,15,18,19,21,22,111,132,135,151,33,53,54,158,62,63,66,68,69,73,75,81,83,84,88,166,172,184,1]
%B204a:[1,4,6,7,8,12,19,20,27,31,121,33,34,37,39,40,44,141,144,46,50,54,62,64,65,66,68,111,71,72,78,80,81,158,89,90,91,102,103,106]
%B200b:[8,18,19,23,27,30,205,47,53,68,70,71,78,85,207,238,95,100,138,139,196,199,211,213,232,142,145,157,158,164,165,179,181,216,218]
%B200a:[0,3,4,25,31,32,33,55,56,169,215,220,299,300,59,61,62,63,67,69,75,78,81,87,89,98,100,101,103,105,186,194,202,207,213,268,269,280,286,307,106,122,128,198,223,267,130,140,141,142,144,156,158,159,162,163,200,237,244,252,254]
%B198b:[10,15,18,20,32,33,34,35,36,40,46,47,185,187,189,230,231,232,233,51,57,61,66,68,70,74,75,79,89,91,93,94,97,101,200,203,205,224,225,236,237,238,248,249,263,277,278,104,119,120,121,145,146,149,211,212,213,246,269,271,151,152,154,169,170,172,177,179,181,214,242,243,244];
%B198a:[3,9,12,14,16,17,18,19,20,23,29,30,83,86,92,94,104,105,44,48,53,54,55,111,62,64,69,73,74,76,77,79,81,90,91];
%B196a:[0,1,3,4,13,16,19,21,129,142,45,47,48,49,51,52,53,144,74,76,77,81,87,90,92,133,141,143,94,95,97,102,104,105,106,107,110,113,116,118,138,146]
%B192d:[0,4,8,9,12,20,23,81,105,293,382,385,40,44,91,106,55,60,86,107,131,132,141,151,168,175,182,186,188,191,1]
%B192c:[0,3,4,14,22,28,141,145,172,181,182,50,51,58,161,164,70,71,72,79,88,92,93,100,101,157,159,1]
%B190c:[0,6,8,10,11,13,15,16,17,19,104,105,113,118,119,160,36,39,40,44,106,47,51,52,53,56,59,61,62,99,138,142,144,151,152,1];
%B188b:[6,11,21,32,38,147,309,47,49,56,57,67,72,126,133,75,81,82,90,106,108,134,135,1];
%B187b:[5,7,8,15,18,101,122,28,29,38,45,47,48,53,57,63,65,67,86,89,222,1];                 
%B186d:[0,5,7,19,23,77,85,104,170,172,202,209,211,50,56,58,73,195,196,1]
%B186c:[1,4,8,108,109,33,54,55,281,287,295,73,79,96,318,1];
%B184c:[0,1,6,7,8,10,16,17,75,21,26,27,34,40,43,85,117,141,152,162,45,50,52,56,57,62,73,81,95,96,108,170,172,245,286,303,304,306,310,322,1]
%B183a:[0,4,5,8,13,20,24,26,27,31,62,72,80,38,44,46,49,55,69,103,109,207,227,228,239,264,280,294,1]
%B182a:[0,1,6,8,9,17,26,78,89,106,110,125,132,134,141,143,34,35,38,41,43,45,47,53,58,70,71,73,80,158,163,166,188,1]
%B181g:[1];
%B179b:[5,6,7,10,22,61,62,79,99,101,121,132,30,31,34,39,45,52,54,55,57,64,71,118,123,147,154,157,158,159,1];
%B178g:[11,14,15,16,20,22,31,40,56,64,80,85,89,92,99,103,105,108,109,1];
%B177d:[2,3,4,7,18,22,24,69,36,38,41,42,50,51,52,54,65,82,1];
%shank identifier (should match the dat file order, 0 stands for glass electrode recording)

un.shnk = [repmat([1],1,7),repmat([2],1,8),repmat([3],1,10),repmat([4],1,12),repmat([5],1,9),repmat([6],1,6),repmat([7],1,8),repmat([8],1,4)];
%B207c:[repmat([1],1,8),repmat([2],1,8),repmat([3],1,7),repmat([4],1,5),repmat([5],1,6),repmat([6],1,6),repmat([7],1,4),repmat([8],1,5)]
%B207b:[repmat([1],1,3),repmat([2],1,3),repmat([3],1,7),repmat([4],1,3),repmat([5],1,7),repmat([6],1,1),repmat([7],1,2),repmat([8],1,1)]
%B207a:[repmat([1],1,2),repmat([2],1,4),repmat([3],1,8),repmat([4],1,5),repmat([5],1,6),repmat([6],1,6),repmat([7],1,5),repmat([8],1,7)]
%B206b:[repmat([1],1,11),repmat([2],1,3),repmat([3],1,4),repmat([4],1,6),repmat([5],1,4),repmat([6],1,8),repmat([7],1,3),repmat([8],1,7)]
%B206a:[repmat([1],1,3),repmat([2],1,6),repmat([3],1,11),repmat([4],1,3),repmat([5],1,8),repmat([6],1,5),repmat([7],1,2),repmat([8],1,5)]
%B193b:[repmat([1],1,12),repmat([2],1,4),repmat([3],1,14),0]
%B204a:[repmat([1],1,6),repmat([2],1,5),repmat([3],1,8),repmat([4],1,3),repmat([5],1,6),repmat([6],1,6),repmat([7],1,3),repmat([8],1,3)];
%B200b:[repmat([1],1,7),repmat([2],1,9),repmat([3],1,9),repmat([4],1,10)]
%B200a:[repmat([1],1,14),repmat([2],1,26),repmat([3],1,6),repmat([4],1,15)]
%B198b:[repmat([1],1,19),repmat([2],1,28),repmat([3],1,13),repmat([4],1,13)]
%B198a:[repmat([1],1,9),repmat([2],1,9),repmat([3],1,6),repmat([4],1,11)]
%B196a:[repmat([1],1,10),repmat([2],1,8),repmat([3],1,10),repmat([4],1,14)]
%B192d:[repmat([1],1,12),repmat([2],1,4),repmat([3],1,14),0]
%B192c:[repmat([1],1,11),repmat([2],1,5),repmat([3],1,11),0];
%B190c:[repmat([1],1,16),repmat([2],1,5),repmat([3],1,14),0];
%B188b:[repmat([1],1,7),repmat([2],1,8),repmat([3],1,8),0];
%B187b:[repmat([1],1,7),repmat([2],1,6),repmat([3],1,8),0];   
%B186d:[repmat([1],1,13),repmat([2],1,2),repmat([3],1,4),0]
%B186c:[repmat([1],1,5),repmat([2],1,6),repmat([3],1,4),0]
%B184c:[repmat([1],1,9),repmat([2],1,11),repmat([3],1,20),repmat([4],1,0),1]
%B183a:[repmat([1],1,13),repmat([2],1,0),repmat([3],1,15),0]
%B182a:[repmat([1],1,16),repmat([2],1,0),repmat([3],1,17),repmat([4],1,0),0]
%B181g:[0]
%B179b:[repmat([1],1,12),repmat([2],1,0),repmat([3],1,18),0]
%B178g:[repmat([1],1,1),repmat([2],1,0),repmat([3],1,18),0];
%B177d:[repmat([1],1,8),repmat([2],1,0),repmat([3],1,10),0];
%boundaries of glass electrode recording:
gl_on=[ ]; gl_off=[ ];

%CALCULUS ON ANIMAL POSITION
%read the time vs. position file (time in seconds, position in metres)
pos_name=strcat(bs_name,'POSrs.txt');
pos_orig=dlmread(pos_name);

%determine the actual speed for each time point
%column
for je=1:(size(pos_orig,1))
    if je==1
        pos_orig(je,3)=(pos_orig(je+1,2)-pos_orig(je,2))/(pos_orig(je+1,1)-pos_orig(je,1));
    elseif je==size(pos_orig,1)
        pos_orig(je,3)=(pos_orig(je,2)-pos_orig(je-1,2))/(pos_orig(je,1)-pos_orig(je-1,1));
    else
        pos_orig(je,3)=(((pos_orig(je,2)+pos_orig(je+1,2))/2)-((pos_orig(je-1,2)+pos_orig(je,2))/2))/(pos_orig(je,1)-pos_orig(je-1,1));
    end
end; clear je

%column4: track identifier
for je=1:(size(pos_orig,1))
    pos_orig(je,4)=tr_idnt;
    if pos_orig(je,3)<((tr_length/(pos_orig(2,1)-pos_orig(3,1)))/4)
        if min(pos_orig(je-10:je-1,3))>((tr_length/(pos_orig(2,1)-pos_orig(3,1)))/4)
            tr_idnt=tr_idnt+1;
        end
    end
end

%column5: smoothed speed track + correcting smaller errors (error
%considered anything above 2m/s or <-1 m/s)
%column6: also every point with speed < mov_thrs (movement threshold)
%replaced by NaN
tempor=pos_orig(:,3);
for je=2:size(tempor,1)
    if tempor(je,1)>2
        tempor(je,1)=tempor(je-1,1);
    elseif tempor(je,1)<-1
        tempor(je,1)=tempor(je-1,1);
    end
end; clear je
tempor_S=smooth(tempor,5);
for je=2:size(tempor,1)
    if pos_orig(je,3)>2
        tempor_S(je,1)=NaN;
    elseif pos_orig(je,3)<-1
        tempor_S(je,1)=NaN;
    end
end; clear je
pos_orig(1:size(pos_orig,1),5)=tempor_S;
pos_orig(1:size(pos_orig,1),6)=tempor_S;
clear tempor tempor_S
%cloumn 6: removing values below mov_thrs
for je=1:size(pos_orig,1)
    if pos_orig(je,5)<=mov_thrs
       pos_orig(je,6)=NaN;
    end
end; clear je
%column7: position with NaN at positions that are not included
for je=1:size(pos_orig,1)
    if isnan(pos_orig(je,6))
       pos_orig(je,7)=NaN;
    else
        pos_orig(je,7)= pos_orig(je,2);
    end
end; clear je
%extract the position sample time interval
pos_stim=mean(pos_orig(2:end,1)-pos_orig(1:end-1,1));

%OCCUPACY maps
for je=1:round(tr_length/tr_binsiz)
    %cycling across spatial bins
    binst=tr_binsiz*(je-1); binen=tr_binsiz*(je); binct=tr_binsiz*(je-0.5);
    %cycling accross time points
    occupancy_all{je}(1)=NaN;
    speedbind_all{je}(1,1)=NaN;
    for ka=2:(size(pos_orig,1))
        if isnan(pos_orig(ka,7))
            occupancy_all{je}(ka)=NaN;
            speedbind_all{je}(ka,1)=NaN;
        elseif pos_orig(ka,7)<=binst
            occupancy_all{je}(ka)=NaN;
            speedbind_all{je}(ka,1)=NaN;
        elseif pos_orig(ka,7)>binen
            occupancy_all{je}(ka)=NaN;
            speedbind_all{je}(ka,1)=NaN;
%         elseif ka==(size(pos_orig,1))
%             occupancy_all{je}(ka)=NaN;
%             speedbind_all{je}(ka,1)=NaN;
        else
            speedbind_all{je}(ka,1)=pos_orig(ka,3);
            if pos_orig(ka-1,7)<=binst %if this is the first point in the bin
                if pos_orig(ka+1,7)>binen            %if this is the last point as well
                    occupancy_all{je}(ka)=(binen-binst)/pos_orig(ka,6);
                else                                 %if it is not the last
                    occupancy_all{je}(ka)=(pos_orig(ka,7)-binst)/(pos_orig(ka,6));
                end
            else                                     %if this is not the first point in the bin
                if pos_orig(ka+1,7)>binen            %if this is the last point as well
                    occupancy_all{je}(ka)=(pos_orig(ka,1)-pos_orig(ka-1,1))+((binen-(pos_orig(ka,7)))/pos_orig(ka,6));
                else                                 %if it is not the last
                    occupancy_all{je}(ka)=(pos_orig(ka,1)-pos_orig(ka-1,1));
                end
            end
        end
        speedbind_all{je}(ka,2)=pos_orig(ka,4);
    end; clear ka
    clear binst binen binct
end; clear je

%SPEED MAPS
%speed binned on the linear map
speedmap_all(:,1)=[tr_binsiz/2:tr_binsiz:tr_length-(tr_binsiz/2)];
if size(speedbind_all,2)==size(speedmap_all,1)
    %ALL TRIALS
    %step throuhg trials
    for tr=1:max(pos_orig(:,4))
        %step through position bins
        for ka=1:size(speedmap_all,1)
            speedmap_all(ka,tr+1)=nanmean(speedbind_all{ka}(speedbind_all{ka}(:,2)==tr,1));
        end; clear ka
    end; clear tr
    speedmap_allmean=[nanmean(speedmap_all(:,2:end),2),(std((speedmap_all(:,2:end))','omitnan'))'];
    %LIMITED TO GLASS ELECTRODE RECORDING
    if ~isempty(gl_on)
        speedmap_lim=speedmap_all;
        speedmap_lim(1:end,2:end)=NaN;
        for tr=min(pos_orig((gl_on/pos_stim):(gl_off/pos_stim),4)):max(pos_orig((gl_on/pos_stim):(gl_off/pos_stim),4))
            %step through position bins
            for ka=1:size(speedmap_lim,1)
                speedmap_lim(ka,tr+1)=nanmean(speedbind_all{ka}(speedbind_all{ka}(:,2)==tr,1));
            end; clear ka
        end; clear tr
        speedmap_limmean=[nanmean(speedmap_lim(:,2:end),2),(std((speedmap_lim(:,2:end))','omitnan'))'];
    end
else
    'size of bins in occupancy map & spike count does not match'
end
%PLOTTING THE SPEED MAPS
figname=strcat(bs_name, '_speedmap');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
%all trials
subplot(2,1,1)
line(pos_orig(:,2),pos_orig(:,6),'Color',[0.7 0.7 0.7])
line(speedmap_all(:,1),speedmap_allmean(:,1),'Color',[0 0 0]);
line(speedmap_all(:,1),speedmap_allmean(:,1)+speedmap_allmean(:,2),'Color',[0 0 0]);
line(speedmap_all(:,1),speedmap_allmean(:,1)-speedmap_allmean(:,2),'Color',[0 0 0]);
axis([-0.1,tr_length+0.1,0,max(pos_orig(:,6))]);
xlabel('track position(m)'); ylabel ('speed(m/s)')
title(strcat('speed profile for all trials (',bs_name,')'));
if ~isempty(gl_on)
    subplot(2,1,2)
    line(pos_orig((gl_on/pos_stim):(gl_off/pos_stim),2),pos_orig((gl_on/pos_stim):(gl_off/pos_stim),6),'Color',[0.7 0.7 0.7]);
    line(speedmap_lim(:,1),speedmap_limmean(:,1),'Color',[0 0 0]);
    line(speedmap_lim(:,1),speedmap_limmean(:,1)+speedmap_limmean(:,2),'Color',[0 0 0]);
    line(speedmap_lim(:,1),speedmap_limmean(:,1)-speedmap_limmean(:,2),'Color',[0 0 0]);
    axis([-0.1,tr_length+0.1,0,max(pos_orig(:,6))]);
    xlabel('track position(m)'); ylabel ('speed(m/s)')
    title(strcat('speed profile for limited trials (',bs_name,')'));
end
%saving the figure
saveas(gcf,figname)
clear figname
        
%THETA ANALYSIS
%reading pyramidal layer LFP
if isempty(sh_pyr)
    lfp_name=strcat(bs_name,'LFP.mat');
else
    lfp_name=strcat(bs_name,'sh',num2str(sh_pyr),'LFP.mat');
end
load(lfp_name);
LFPred=downsample(smooth(LFP(ch_pyr,:),ch_dwnsmpl),ch_dwnsmpl, (ch_dwnsmpl-1));
clear LFP lfp_name;
%creating theta filter
th_filter=dfilt.dffir(fir1(flt_ncoeff, [2.*(flt_low./ch_rate) 2.*(flt_high./ch_rate)],'bandpass', gausswin(flt_ncoeff+1)));
%filtering for theta
LFPflt=filtfilt(th_filter.Numerator,1,LFPred')'; clear LFPred th_filter;

%detecting theta troughs for pyramidal layer LFP
LFPthreshold=0.001*(mean(abs(LFPflt)));
[~,th_troughsLFP]=findpeaks(LFPflt.*-1,'MINPEAKHEIGHT', LFPthreshold, 'MINPEAKDISTANCE', floor(ch_rate/flt_high)-3);
th_troughsLFP(:,2)=(th_troughsLFP./ch_rate);
clear LFPthreshold

%CALULUS ON POSITION DEPENDENT CELL FIRING
%read the spike times (s)
spk_name=strcat(bs_name,'WAW.mat');
for je=1:length(un.chan)
    spk_chan=strcat(bs_chname,num2str(un.chan(je)));
    spk_load=load(spk_name,spk_chan);
    %load spike times
    spks=(spk_load.(spk_chan).times);
    clear spk_load spk_chan
    %extract positions, includedness, speed, trial identifier for each spike 
    %column2: spike positions
    pos_tsrs=timeseries(pos_orig(:,2),pos_orig(:,1));
    pos_spks=resample(pos_tsrs,spks);
    spks(1:size(spks,1),2)=pos_spks.Data;
    %column3: spike trial affiliations
    trl_tsrs=timeseries(pos_orig(:,4),pos_orig(:,1));
    trl_spks=resample(trl_tsrs,spks(:,1));
    spks(1:size(spks,1),3)=trl_spks.Data;
    %column4&5 position&trial with non-include replaced by NaN
    for ka=1:size(spks,1)
        before=size(pos_orig(pos_orig(:,1)<=spks(ka,1),6),1);
        if isnan(pos_orig(before,6))
            spks(ka,4)=NaN;
            spks(ka,5)=NaN;
            spks(ka,6)=NaN;
        elseif isnan(pos_orig(before+1,6))
            spks(ka,4)=NaN;
            spks(ka,5)=NaN;
            spks(ka,6)=NaN;
        else
            spks(ka,4)=spks(ka,2);
            spks(ka,5)=spks(ka,3);
            %determine theta phase
            cyc_bi=size(th_troughsLFP(th_troughsLFP(:,2)<=spks(ka,1),1),1); cyc_ei=cyc_bi+1;
            cyc_bt=th_troughsLFP(cyc_bi,2);cyc_et=th_troughsLFP(cyc_ei,2);
            cyc=LFPflt(th_troughsLFP(cyc_bi,1):th_troughsLFP(cyc_ei,1));
            cyc_ct=((th_troughsLFP(cyc_bi,1)+find(cyc==(max(cyc)),1))./ch_rate);
            if spks(ka,1)<cyc_ct
                spks(ka,6)=((spks(ka,1)-cyc_bt)/(cyc_ct-cyc_bt))*180;
            else
                spks(ka,6)=(((spks(ka,1)-cyc_ct)/(cyc_et-cyc_ct))*180)+180;
            end
            clear cyc cyc_bi cyc_bt cyc_ei cyc_et cyc_ct
        end
        clear before
    end
    spk_posi{je}=spks;
    clear ka spks  pos_tsrs pos_spks trl_tsrs trl_spks
end
clear spk_name je

%CALCULUS FOR SPATIAL MAPS
for je=1:length(un.chan)
    %spatial bins
    spk_map{je}(:,1)=[tr_binsiz/2:tr_binsiz:tr_length-(tr_binsiz/2)];
    %spike count
    spk_map{je}(:,2)=hist(spk_posi{je}(:,4),spk_map{je}(:,1));
    if size(occupancy_all,2)==size(spk_map{je},1)
        if un.shnk(je)~=0   %the unit is not glass unit: all file included
            for ka=1:size(spk_map{je},1)
                spk_map{je}(ka,3)=nansum(occupancy_all{ka});
            end; clear ka
        else                %the unit has a limited recording range
            for ka=1:size(spk_map{je},1)
                spk_map{je}(ka,3)=nansum(occupancy_all{ka}((gl_on/pos_stim):(gl_off/pos_stim)));
            end; clear ka
        end
    else
        'size of bins in occupancy map & spike count does not match'
    end
    %column4: firing rate
    spk_map{je}(:,4)=spk_map{je}(:,2)./spk_map{je}(:,3);
    %removing bins that have no more than 0.4 s total recording
    for ka=1:size(spk_map{je},1)
        if spk_map{je}(ka,3)<0.4
            spk_map{je}(ka,2:4)=NaN;
        end
    end; clear ka
    
    %CALCULUS FOR FIRING STATISTICS statistics
    %item1: overall firing rate
    spk_stat(je,1)=nansum(spk_map{je}(:,2))./nansum(spk_map{je}(:,3));
    %item2: time spent on maze
    spk_stat(je,2)=nansum(spk_map{je}(:,3));
    for ka=1:size(spk_map{je},1)
        %column1: occupancy ratio (P(i))
        calc_spinf(ka,1)=(spk_map{je}(ka,3))/(spk_stat(je,2));
        %column2: firing rate ratio (lambda(i))
        calc_spinf(ka,2)=spk_map{je}(ka,4)/(spk_stat(je,1));
        %column3:sptial info stats
        calc_spinf(ka,3)=(calc_spinf(ka,1))*(spk_map{je}(ka,4))*(log2(calc_spinf(ka,2)));
    end; clear ka
    %item3: spatial info (bits per second)
    spk_stat(je,3)=nansum(calc_spinf(:,3));
    %item4: spatial info (bits per spike)
    spk_stat(je,4)=spk_stat(je,3)/spk_stat(je,1);
    %item5: spike number    
    spk_stat(je,5)=nansum(spk_map{je}(:,2));
    clear calc_spinf
end; clear je
        
%PLOT THE DATA
pltindx=1; figindx=0;
for je=1:length(un.chan)
    if pltindx==1
        figindx=figindx+1;
        figname=strcat(bs_name, '_spatialmod_fg',num2str(figindx));
        figure ('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
    end
    %TEXTUAL
    subplot(4,3,pltindx)
    axis off
    text(0,1,strcat(bs_name, ' sh',num2str(un.shnk(je)),' un',num2str(un.nums(je))))
    %length of recording
    actualtext=strcat('time included: ', num2str(spk_stat(je,2)), ' s');
    text(0,0.8,actualtext);
    %mean AP frequency (included spikes)
    actualtext=strcat('overall FR: ', num2str(spk_stat(je,1)), ' Hz');
    text(0,0.6,actualtext);
    %N of included spikes
    actualtext=strcat('overall N: ', num2str(spk_stat(je,5)));
    text(0,0.4,actualtext);
    %information(bits/sec)
    actualtext=strcat('info: ', num2str(spk_stat(je,3)),'(bit/sec)');
    text(0,0.2,actualtext);
    %information(bits/spike)
    actualtext=strcat('info: ', num2str(spk_stat(je,4)),'(bit/spike)');
    text(0,0,actualtext);
    
    %BY TRIALS
    subplot(4,3,3+pltindx)
    scatter(spk_posi{je}(:,4),spk_posi{je}(:,5),'.','k','SizeData',0.5)
    axis([-0.1,tr_length+0.1,0,max(spk_posi{je}(:,5))+1]);
    xlabel('track position(m)'); ylabel ('trial#')
    
    %FR
    subplot(4,3,6+pltindx)
    line(spk_map{je}(:,1),spk_map{je}(:,4))
    axis([-0.1,tr_length+0.1,0,1.1*(max(spk_map{je}(:,4)))]);
    xlabel('track position(m)'); ylabel ('FR(Hz)')
    
    %THETA PHASE - POSITION
    subplot(4,3,9+pltindx)
    scatter([spk_posi{je}(:,4);spk_posi{je}(:,4)],[spk_posi{je}(:,6);spk_posi{je}(:,6)+360],'.','k','SizeData',0.5)
    axis([-0.1,tr_length+0.1,0,720]);
    xlabel('track position(m)'); ylabel ('theta phase (°)')
    set (gca, 'Ytick', 0:90:720);
    
    %saving figures
    if pltindx==3||je==length(un.chan)
        saveas(gcf,figname)
        clear figname
    end
    
    if pltindx==3
        pltindx=1;
    else
        pltindx=pltindx+1;
    end
    
    
end
clear je figindx pltindx

%saving the relevant data
datum = date;
writename=strcat(bs_name, '_spatialfiring');
if ~isempty(gl_on)
    save(writename, 'occupancy_all', 'pos_orig', 'speedmap_allmean', 'speedmap_limmean','spk_map', 'spk_posi', 'spk_stat', 'tr_length', 'tr_binsiz', 'datum');
else
    save(writename, 'occupancy_all', 'pos_orig', 'speedmap_allmean', 'spk_map', 'spk_posi', 'spk_stat', 'tr_length', 'tr_binsiz', 'datum');
end


