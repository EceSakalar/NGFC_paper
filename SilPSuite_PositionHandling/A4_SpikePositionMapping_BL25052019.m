%This program generates 1D position dependent firing rate maps for units
clear all; close all

%give the experimens file name (.mat and spike2 files)
bs_name='B196a0WAW.mat';
bs_sp2n='Bsx196a';

%give the unit data (shank, 
un_shnk=[repmat([1],1,10), repmat([2],1,8), repmat([3],1,10), repmat([4],1,14)];  
un_numr=[0,1,3,4,13,16,19,21,129,142,45,47,48,49,51,52,53,144,74,76,77,81,87,90,92,133,141,143,94,95,97,102,104,105,106,107,110,113,116,118,138,146];
un_chan=[65:106];

%give the channel for the position along the corridor
%and the starting point of the corridor
%and the speed threshold (in cm/s)
%and the position binning size (in cm) !! the position resolution fro time estimetes is always 1 cm
ch_pos=107;
corridor_start=-360;
speed_thrs=5;
pos_bin=5;

%give the details of the theta filter and the channels where to take the theta from
ch.pyrsh=3;                                         %shank for the pyramidal layer theta to be taken from
ch.pyrch=14;                                        %channel number (ventralmost LFP contact is 1) for the pyramidal layer
smp.rate=10000;                                     %give the sampling rate of the channels (CSD and LFP original)

%THETA properties
%theta filter properties (number of filter coefficients and filter corner frequencies in Hz)
th.ncoeff=512; th.fltlow=5; th.flthigh=12; th.dwnsmpl=25; th.rate=smp.rate/th.dwnsmpl;
%creating theta filter
th_filter=dfilt.dffir(fir1(th.ncoeff, [2.*(th.fltlow./th.rate) 2.*(th.flthigh./th.rate)],'bandpass', gausswin(th.ncoeff+1)));
tic
%Read theta channel do theta filtering and detect peaks and troughs
%read the LFP for theta
lfp_name=strrep(bs_name,'WAW.mat',strcat('sh', num2str(ch.pyrsh),'LFP.mat'));
load(lfp_name);
thLFP = LFP(ch.pyrch,:); clear LFP lfp_name;
%downsample and filter LFP theta
thLFP_red=downsample(smooth(thLFP,th.dwnsmpl),th.dwnsmpl);
thLFP_flt=filtfilt(th_filter.Numerator,1,thLFP_red')';

%detecting theta troughs for pyramidal layer LFP
thLFP_THRSHLD = 0.001*(mean(abs(thLFP_flt)));
[~,thLFP_troughs]=findpeaks(thLFP_flt.*-1,'MINPEAKHEIGHT', thLFP_THRSHLD, 'MINPEAKDISTANCE', floor(th.rate/th.flthigh)-3);
thLFP_troughsT=(thLFP_troughs./th.rate);
clear thLFP_THRSHLD thLFP_red thLFP

'theta filtering done'
toc

%Read available position information
ch_name=strcat(bs_sp2n,'_Ch', num2str(ch_pos)); ch_load=load(num2str(bs_name),num2str(ch_name));
ps_int=ch_load.(ch_name).('interval');
ps_str=ch_load.(ch_name).('start');

%variable ps_seq will contain the following columns
%1. time (in seconds)
%2. position (in cm with the corridor strat being 0)
%3. speed (cm/s)
%4. 1 if speed is above the movement threshold, Nan if not
%5. distance covered (cm)
%6. Trial number

ps_seq=(([0:1:(ch_load.(ch_name).('length'))-1]').*ps_int)+ps_str;
ps_seq(:,2)=(((ch_load.(ch_name).('values'))')./10)-corridor_start;
ps_seq(1,3)=0; ps_seq(end,3)=0;
ps_seq(1,4)=NaN; ps_seq(end,4)=NaN;
ps_seq(1,5)=NaN; ps_seq(end,5)=NaN;
clear ch_name ch_load

for je=2:(size(ps_seq,1)-1)
    ps_seq(je,3)=(ps_seq(je+1,2)-ps_seq(je-1,2))/(ps_seq(je+1,1)-ps_seq(je-1,1));
    ps_seq(je,5)=(ps_seq(je+1,2)-ps_seq(je-1,2))./10;
    if ps_seq(je,3)>=speed_thrs
        ps_seq(je,4)=1;
    else
        ps_seq(je,4)=NaN;
    end
end; clear je

ps_seq(1:4,6)=1;
for je=5:size(ps_seq,1)
    if ps_seq(je,2)-ps_seq(je-4,2)<-500
        ps_seq(je,6)=ps_seq(je-6,6)+1;
    else
         ps_seq(je,6)=ps_seq(je-1,6);
    end
end; clear je

'variable ps seq done'
toc

%variable mz_pos with the following colums
%1. position bin centers
%2. cumulative time spent within
%3. cumulative time with NaN for values <30% of the mean
%4. probability of finding animal in a position (no NaNs in col3 included)
mz_pos=[floor(min(ps_seq(:,2))):1:ceil(max(ps_seq(:,2)))]';
for je=1:max(ps_seq(:,6))
    a_trial=ps_seq(ps_seq(:,6)==je & (ps_seq(:,3))>speed_thrs,:);
    a_posseq=resample(timeseries(a_trial(:,1),a_trial(:,2)),[mz_pos(1,1)-0.5:1:mz_pos(end,1)+0,5]);
    for ka=1:size(mz_pos,1)
        if isnan(a_posseq.data(ka+1))||isnan(a_posseq.data(ka))
            mz_tmp(ka,je)=NaN;
        elseif (1/(a_posseq.data(ka+1)-a_posseq.data(ka)))<speed_thrs
            mz_tmp(ka,je)=NaN;
        else
            mz_tmp(ka,je)=a_posseq.data(ka+1)-a_posseq.data(ka);
        end
    end; clear ka
    clear a_posseq a_trial
end; clear je
mz_pos(:,2)=nansum(mz_tmp');
mz_pos(:,3)=mz_pos(:,2);
mz_pos(mz_pos(:,2)<mean(mz_pos(:,2))*0.3,3)=NaN;
mz_pos(:,4)=mz_pos(:,3)./nansum(mz_pos(:,3));

'variable maze done'
toc

%variable un_place with the following columns (array with each unit separate)
%1. spike times
%2. position at the spike
%3. speed
%4. position for included (above speed threshold spikes
%5. Spike theta phase
 
for je=1:length(un_numr)
    ch_name=strcat(bs_sp2n,'_Ch', num2str(un_chan(je))); un_load=load(num2str(bs_name),num2str(ch_name));
    un_place{je}=un_load.(ch_name).('times');
    un_place{je}(un_place{je}<ps_seq(2,1),:)=[];
    un_place{je}(un_place{je}>ps_seq(size(ps_seq,1)-1,1),:)=[];
    clear ch_name un_load
    for ka=1:size(un_place{je},1)
        ka_tim=un_place{je}(ka,1);
        ka_bef=find(ps_seq(:,1)==max(ps_seq(ps_seq(:,1)<=ka_tim,1)));
        ka_aft=find(ps_seq(:,1)==min(ps_seq(ps_seq(:,1)>ka_tim,1)));
        ka_pos=ps_seq(ka_bef,2)+((ps_seq(ka_aft,2)-ps_seq(ka_bef,2))*((ka_tim-ps_seq(ka_bef,1))/(ps_seq(ka_aft,1)-ps_seq(ka_bef,1))));
        ka_spd=(ps_seq(ka_bef,3)+ps_seq(ka_aft,3))/2;
        un_place{je}(ka,2)=ka_pos;
        un_place{je}(ka,3)=ka_spd;
    end; clear ka ka_tim ka_bef ka_aft ka_pos ka_spd
    un_place{je}(un_place{je}(:,3)>=speed_thrs,4)=un_place{je}(un_place{je}(:,3)>=speed_thrs,2);
    un_place{je}(un_place{je}(:,3)<speed_thrs,4)=NaN;
    %get the theta phase
    for ka=1:size(un_place{je},1)
        if isnan(un_place{je}(ka,4))
            un_place{je}(ka,5)=NaN;
        else
            cyc_beg=max(thLFP_troughsT(thLFP_troughsT<=un_place{je}(ka,1)));
            cyc_end=min(thLFP_troughsT(thLFP_troughsT>un_place{je}(ka,1)));
            cyc_LFP=thLFP_flt((round(cyc_beg*th.rate)):round((cyc_end*th.rate))); [~,cyc_cent]=max(cyc_LFP);
            cyc_cent=(cyc_cent/th.rate)+cyc_beg;
            if  un_place{je}(ka,1)<=cyc_cent
                tht_phase=180*((un_place{je}(ka,1)-cyc_beg)/(cyc_cent-cyc_beg));
            else
                tht_phase=180*((un_place{je}(ka,1)-cyc_cent)/(cyc_end-cyc_cent))+180;
            end
            un_place{je}(ka,5)=tht_phase;
            clear cyc_beg cyc_end cyc_LFP cyc_cent tht_phase
        end
    end; clear ka
    
    strcat('variable un_place for unit', num2str(un_numr(je)), 'done')
    toc
    
end; clear je 
clear thLFP_flt thLFP_troughs thLFP_troughsT
%un_stat variable will be created with the following output
%1.-3. shank, unit, channel
%4. overall firing rate
%5. spatial information

un_stat=[un_shnk' un_numr' un_chan'];
for je=1:length(un_numr)
    un_ratemap{je}=((hist(un_place{je}(:,4),mz_pos(:,1)))')./mz_pos(:,3);
    un_stat(je,4)=sum(((hist(un_place{je}(:,4),mz_pos(:,1)))').*(~isnan(mz_pos(:,3))))/nansum(mz_pos(:,3));
    un_stat(je,5)=nansum(mz_pos(:,4).*un_ratemap{je}.*(log2((un_ratemap{je}./un_stat(je,4)))));
end

%PLOTTING
for ka=1:ceil(length(un_numr)/4)
    figname=strcat(strrep(bs_name, 'WAW.mat','place_'), 'plot', num2str(ka));
    figure ('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
    for el=1:4
        je=((ka-1)*4)+el;
        if je<=length(un_numr)
            subplot(4,4,el)
            axis off
            text(0,0.9,strcat(strrep(bs_name, 'WAW.mat',' sh'), num2str(un_shnk(je)), ' unit', num2str(un_numr(je))))
            text(0,0.6,strcat('spatial info: ', num2str(un_stat(je,5))));
            text(0,0.3,strcat('mean rate: ', num2str(un_stat(je,4)),'Hz'));
            text(0,0,strcat('peak rate: ', num2str(max(un_ratemap{je})),'Hz'));
            subplot(4,4,el+4)
            scatter(un_place{je}(:,4),un_place{je}(:,1),1,'.')
            axis([0,max(mz_pos(:,1)),0,max(un_place{je}(:,1))])
            xlabel('position'); ylabel('exp. time (s)');
            subplot(4,4,el+8)
            plot(mz_pos(:,1),((hist(un_place{je}(:,4),mz_pos(:,1)))')./mz_pos(:,3),'color','r');
            xlabel('position'); ylabel('FR (Hz)');
            axis([0, max(mz_pos(:,1)), 0,(max(((hist(un_place{je}(:,4),mz_pos(:,1)))')./mz_pos(:,3)))*1.1])
            subplot(4,4,el+12)
            scatter([un_place{je}(:,4);un_place{je}(:,4)],[un_place{je}(:,5);un_place{je}(:,5)+360],1,'.')
            xlabel('position'); ylabel('theta phase');
            axis([0,max(mz_pos(:,1)),0,720])
            yticks([0,180,360,540,720])
        end
        clear je
    end; clear el
end;clear ka
            



        
