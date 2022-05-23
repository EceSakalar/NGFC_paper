%This program shuffles theta cycles for spikes of one cell
%assigns spikes in one theta cycles to another randomly chosen theta cycle

%expected input files: ..._gammacomponent_....mat
clear all; %close all

%basic data entry
bs.name='B'; bs.num='193'; bs.exp='b'; bs.typ='sx';                 %unit basic data entry
gm.comp = 'CA1mid'; 

%units data entry
un.chan=[33];                                   %channels containing the spike times for the units to shuffle
%unit numbers;
un.nums=[33];                                    %unit numbers (in kwik file)
%shanks from which unit was isolated;
un.shnk =[2];                       %shanks from which unit was isolated

%THETA & GAMMA properties
th.nbins = 20; th .binsize = 360/th.nbins;                                  %give the number of theta bins 
gm.nbins = 20; gm.binsize = 360/gm.nbins;                                  %give the number of gamma bins 
 
%reading previously generated data from the MATLAB file
readnam=strcat(bs.name, bs.num, bs.exp,'_gammacomponent_', gm.comp, '.mat');
load(readnam); clear readnam
bs.flord=input.base.flord;
ch.pyrch=input.pyrchan.pyrch;
smp.rate=input.gamma.rate;
%extract gamma cycle information from 'OUTPUT' of read file
gmCSD_pts=output.GMcycle;

%---CALCULUS---

tic
%-------read and stich theta 
for fl_indx=1:length(bs.flord) %FILE STEPPER
    
    %read the LFP for theta
    lfp_name=strcat(bs.name, bs.num, bs.exp, num2str(bs.flord(fl_indx)),'LFP.mat');
    load(lfp_name);
    thLFP_file = LFP(ch.pyrch,:); clear LFP lfp_name;
    
    %read theta periods
    th_pname=strcat(bs.name, bs.num, bs.exp, num2str(bs.flord(fl_indx)),'THP.txt');
    if exist(th_pname,'file')==2
        thPER_file=dlmread(th_pname); clear th_pname;
    else
        thPER_file=[]; clear th_pname;
    end
    
    %read the spikes
    fl_name=strcat(bs.name, bs.num, bs.exp, num2str(bs.flord(fl_indx)),'WAW.mat'); %file name
    ch_name=strcat(bs.name, bs.typ, bs.num, bs.exp,'_Ch', num2str(un.chan));
    ch_load=load(num2str(fl_name),num2str(ch_name));
    spks_fl =(ch_load.(ch_name).('times')); clear ch_name ch_load;
    
    
    %concatenate file data (theta trace, theta periods, spikes)
    if fl_indx==1
        flDUR(fl_indx) = length(thLFP_file);
        thLFP_orig = thLFP_file;
        thPER = thPER_file;
        spks=spks_fl;
    else
        flDUR(fl_indx) = flDUR(fl_indx-1)+length(thLFP_file);
        thLFP_orig = [thLFP_orig thLFP_file];
        thePER = [thPER; thPER_file+(flDUR(fl_indx-1)/smp.rate)];
        spks = [spks; spks_fl+(flDUR(fl_indx-1)/smp.rate)];
    end
    clear thLFP_file thPER_file spks_fl fl_name
end; clear fl_indx flDUR

%---CALCULUS----theta cycles
%downsampling & filtering for theta
th.ncoeff=input.theta.ncoeff; th.fltlow=input.theta.fltlow; th.flthigh=input.theta.flthigh; th.dwnsmpl=input.theta.dwnsmpl; th.rate=smp.rate/th.dwnsmpl;
th_filter=dfilt.dffir(fir1(th.ncoeff, [2.*(th.fltlow./th.rate) 2.*(th.flthigh./th.rate)],'bandpass', gausswin(th.ncoeff+1)));
thLFP_red=downsample(smooth(thLFP_orig,th.dwnsmpl),th.dwnsmpl,th.dwnsmpl-1);
thLFP_flt=filtfilt(th_filter.Numerator,1,thLFP_red')';

%detecting theta troughs for pyramidal layer LFP
thLFP_THRSHLD = 0.001*(mean(abs(thLFP_flt)));
[~,thLFP_troughs]=findpeaks(thLFP_flt.*-1,'MINPEAKHEIGHT', thLFP_THRSHLD, 'MINPEAKDISTANCE', floor(th.rate/th.flthigh)-3);
thLFP_troughsT=(thLFP_troughs./th.rate);
clear thLFP_red thLFP_THRSHL thLFP_orig

%Extracting theta cycle data into thCYC file
%theta period stepper
for je=1:size(thPER,1)
    thLFP_pertrough=thLFP_troughsT(thLFP_troughsT>=thPER(je,1)&thLFP_troughsT<=thPER(je,2));
    thCYCper=[(ones([1,(size(thLFP_pertrough,1)-1)])).*je; thLFP_pertrough(1:end-1)';thLFP_pertrough(2:end)'];
    if je==1
        thCYC=thCYCper;
    else
        thCYC=[thCYC thCYCper];
    end
end; clear thLFP_pertrough thCYCper
thCYC(4,:)=[1:size(thCYC,2)];
%order randomizer
[~,razer]=sort(rand(size(thCYC,2),1));
thCYC(5:8,:)=thCYC(1:4,razer); 
clear razer
'theta filtering and spikes done'
toc

%---CALCULUS----extract spike data
%col 1 is spike time in seconds, col 2 is spike time in sample (gamma rate based)
spks(:,2) = spks(:,1).*smp.rate;
%col 3 is 1 if spike is in the theta periods
spks(:,3) = NaN;
for je=1:size(thPER,1)
    spks((spks(:,1)>=thPER(je,1)&spks(:,1)<=thPER(je,2)),3)=1;
end; clear je
%col 4 is theta phase (degrees)
%col 5 is closest gamma cycle time (sec)
for sp_indx=1:size(spks,1)
    if isnan(spks(sp_indx,3))
        spks(sp_indx,4:9)=NaN;
    else
        %4 theta phase
        sp_time=spks(sp_indx,1);
        [cyc_beg, cyc_begi]=max(thLFP_troughsT(thLFP_troughsT<=sp_time));
        cyc_end = thLFP_troughsT(cyc_begi+1);
        cyc_LFP = thLFP_flt(thLFP_troughs(cyc_begi):thLFP_troughs(cyc_begi+1));
        [~,cyc_cnti]=max(cyc_LFP); cyc_cnt=(cyc_cnti/th.rate)+cyc_beg;
        if  sp_time<=cyc_cnt;    spks(sp_indx,4)=180*((sp_time-cyc_beg)/(cyc_cnt-cyc_beg));
        else                     spks(sp_indx,4)=180*((sp_time-cyc_cnt)/(cyc_end-cyc_cnt))+180;
        end
        clear cyc_cnti cyc_begi cyc_cnt cyc_cnti cyc_LFP cyc_end cyc_beg
        
        %5-8 closest gamma cycle properties
        [~,cyc_ind]=min(abs((gmCSD_pts(1,:)./smp.rate)-sp_time));
        spks(sp_indx,5) = gmCSD_pts(1,cyc_ind)/smp.rate;
        spks(sp_indx,6) = gmCSD_pts(3,cyc_ind);
        spks(sp_indx,7) = gmCSD_pts(4,cyc_ind);
        spks(sp_indx,8) = gmCSD_pts(7,cyc_ind);
        
        %getting gamma phase
        if gmCSD_pts(3,cyc_ind)==1
            cyc_cnti=cyc_ind;
        elseif sp_time>=(gmCSD_pts(1,cyc_ind)/smp.rate)
            cyc_cnti=cyc_ind+1;
        else
            cyc_cnti=cyc_ind-1;
        end
        % get cycle begin
        if gmCSD_pts(3,cyc_cnti-1)==-1
            cyc_beg=gmCSD_pts(1,(cyc_cnti-1))/smp.rate;
        else
            'error occured: gamma cycle begins with peak'
        end
        % get cycle center
        if gmCSD_pts(3,cyc_cnti)==1
            cyc_cnt=gmCSD_pts(1,(cyc_cnti))/smp.rate;
        else
            'error occured: gamma cycle center is trough'
        end
        % get cycle end
        if gmCSD_pts(3,cyc_cnti+1)==-1
            cyc_end=gmCSD_pts(1,(cyc_cnti+1))/smp.rate;
        else
            'error occured: gamma cycle begins with peak'
        end
        
        %get gamma phase
        if  sp_time<=cyc_cnt;    gma_phase=180*((sp_time-cyc_beg)/(cyc_cnt-cyc_beg));
        else;                    gma_phase=180*((sp_time-cyc_cnt)/(cyc_end-cyc_cnt))+180;
        end
        spks(sp_indx,9) = gma_phase;
        clear sp_time cyc_ind cyc_cnti cyc_beg cyc_cnt cyc_end gma_phase
    end
end; clear sp_indx
    
strcat('unit ', num2str(un.nums), 'on shank ', num2str(un.shnk), ' is done')
toc

%generate shuffled spike train 
for je=1:size(spks,1)
    cyc_origi=find(thCYC(2,:)<=spks(je,1)&thCYC(3,:)>spks(je,1));
    if isempty(cyc_origi)
        spks_sh(je,1)=spks(je);
    else
        %identify theta peak
        cyc_beg = thCYC(6,cyc_origi);
        cyc_end = thCYC(7,cyc_origi);
        cyc_begi=find(thLFP_troughsT==thCYC(6,cyc_origi));
        cyc_LFP = thLFP_flt(thLFP_troughs(cyc_begi):thLFP_troughs(cyc_begi+1));
        [~,cyc_cnti]=max(cyc_LFP); cyc_cnt=(cyc_cnti/th.rate)+cyc_beg;
        clear cyc_begi cyc_cnti
        %get the shuffled time
        if spks(je,4)<180
            spks_t=cyc_beg+(cyc_cnt-cyc_beg)*(spks(je,4)/180);
        else
            spks_t=cyc_cnt+(cyc_end-cyc_cnt)*((spks(je,4)-180)/180);
        end
        spks_sh(je,1)=spks_t;
        clear spks_t cyc_cnt cyc_end cyc_beg cyc_LFP
    end
    clear cyc_origi
end; clear je

%---CALCULUS----extract data for the shuffled spike
%col 1 is shuffled spike time in seconds, col 2 is time in sample (gamma rate based)
spks_sh(:,2) = spks_sh(:,1).*smp.rate;
%col 3 is 1 if spike is in the theta periods
spks_sh(:,3) = NaN;
for je=1:size(thPER,1)
    spks_sh((spks_sh(:,1)>=thPER(je,1)&spks_sh(:,1)<=thPER(je,2)),3)=1;
end; clear je
%col 4 is theta phase (degrees)
%col 5 is closest gamma cycle time (sec)
for sp_indx=1:size(spks_sh,1)
    if isnan(spks_sh(sp_indx,3))
        spks_sh(sp_indx,4:9)=NaN;
    else
        %4 theta phase
        sp_time=spks_sh(sp_indx,1);
        [cyc_beg, cyc_begi]=max(thLFP_troughsT(thLFP_troughsT<=sp_time));
        cyc_end = thLFP_troughsT(cyc_begi+1);
        cyc_LFP = thLFP_flt(thLFP_troughs(cyc_begi):thLFP_troughs(cyc_begi+1));
        [~,cyc_cnti]=max(cyc_LFP); cyc_cnt=(cyc_cnti/th.rate)+cyc_beg;
        if  sp_time<=cyc_cnt;    spks_sh(sp_indx,4)=180*((sp_time-cyc_beg)/(cyc_cnt-cyc_beg));
        else                     spks_sh(sp_indx,4)=180*((sp_time-cyc_cnt)/(cyc_end-cyc_cnt))+180;
        end
        clear cyc_cnti cyc_begi cyc_cnt cyc_cnti cyc_LFP cyc_end cyc_beg
        
        %5-8 closest gamma cycle properties
        [~,cyc_ind]=min(abs((gmCSD_pts(1,:)./smp.rate)-sp_time));
        spks_sh(sp_indx,5) = gmCSD_pts(1,cyc_ind)/smp.rate;
        spks_sh(sp_indx,6) = gmCSD_pts(3,cyc_ind);
        spks_sh(sp_indx,7) = gmCSD_pts(4,cyc_ind);
        spks_sh(sp_indx,8) = gmCSD_pts(7,cyc_ind);
        
        %getting gamma phase
        if gmCSD_pts(3,cyc_ind)==1
            cyc_cnti=cyc_ind;
        elseif sp_time>=(gmCSD_pts(1,cyc_ind)/smp.rate)
            cyc_cnti=cyc_ind+1;
        else
            cyc_cnti=cyc_ind-1;
        end
        % get cycle begin
        if gmCSD_pts(3,cyc_cnti-1)==-1
            cyc_beg=gmCSD_pts(1,(cyc_cnti-1))/smp.rate;
        else
            'error occured: gamma cycle begins with peak'
        end
        % get cycle center
        if gmCSD_pts(3,cyc_cnti)==1
            cyc_cnt=gmCSD_pts(1,(cyc_cnti))/smp.rate;
        else
            'error occured: gamma cycle center is trough'
        end
        % get cycle end
        if gmCSD_pts(3,cyc_cnti+1)==-1
            cyc_end=gmCSD_pts(1,(cyc_cnti+1))/smp.rate;
        else
            'error occured: gamma cycle ends with peak'
        end
        
        %get gamma phase
        if  sp_time<=cyc_cnt;    gma_phase=180*((sp_time-cyc_beg)/(cyc_cnt-cyc_beg));
        else;                    gma_phase=180*((sp_time-cyc_cnt)/(cyc_end-cyc_cnt))+180;
        end
        spks_sh(sp_indx,9) = gma_phase;
        clear sp_time cyc_ind cyc_cnti cyc_beg cyc_cnt cyc_end gma_phase
    end
end; clear sp_indx
    
strcat(' shuffled unit ', num2str(un.nums), 'on shank ', num2str(un.shnk), ' is done')
toc
clear thLFP_flt thLFP_troughs thLFP_troughsT thPER thLFP_THRSHLD


%%----CALCULUS---statistical summary for original spike
thbins=(th.binsize/2):th.binsize:(360-(th.binsize/2));
gmbins=(gm.binsize/2):gm.binsize:(360-(gm.binsize/2));

%original spike train (position 1)
thtHIST_all{1}=hist(spks(:,4),thbins);
thtHIST_sel{1}=hist(spks(spks(:,8)==1,4),thbins);
gmaHIST_all{1}=hist(spks(:,9),gmbins);
gmaHIST_sel{1}=hist(spks(spks(:,8)==1,9),gmbins);

gm_stats(1,1)=nansum(spks(:,3));
gm_stats(2,1)=nansum(spks(:,8));
gm_stats(3,1)=circ_r(deg2rad((spks(spks(:,8)==1,9))));
gm_stats(4,1)=rad2deg(circ_mean(deg2rad((spks(spks(:,8)==1,9)))));
gm_stats(5,1)=circ_rtest(deg2rad((spks(spks(:,8)==1,9))));
if gm_stats(4,1)<0
    gm_stats(4,1)=gm_stats(4,1)+360;
end
%shuffled spike train (position 2)
thtHIST_all{2}=hist(spks_sh(:,4),thbins);
thtHIST_sel{2}=hist(spks_sh(spks_sh(:,8)==1,4),thbins);
gmaHIST_all{2}=hist(spks_sh(:,9),gmbins);
gmaHIST_sel{2}=hist(spks_sh(spks_sh(:,8)==1,9),gmbins);

gm_stats(1,2)=nansum(spks_sh(:,3));
gm_stats(2,2)=nansum(spks_sh(:,8));
gm_stats(3,2)=circ_r(deg2rad((spks_sh(spks_sh(:,8)==1,9))));
gm_stats(4,2)=rad2deg(circ_mean(deg2rad((spks_sh(spks_sh(:,8)==1,9)))));
gm_stats(5,2)=circ_rtest(deg2rad((spks_sh(spks_sh(:,8)==1,9))));
if gm_stats(4,2)<0
    gm_stats(4,2)=gm_stats(4,2)+360;
end

%%----DISPLAY--------
figname=strcat(bs.name, bs.num, bs.exp, gm.comp,'_shuffling_cellcoupling');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
for un_indx = 1:2
    %textual results
    subplot(2,3,((un_indx-1)*3)+1)
    axis off
    actext=strcat(bs.name, bs.num, bs.exp, ' shank', num2str(un.shnk),' unit', num2str(un.nums));
    text(0,0.8,actext); clear actext;
    if un_indx==1
        actext='original spikes';
    else
        actext='shuffled spikes';
    end
    text(0,0.7,actext); clear actext;
    actext=strcat('coupling to: ', gm.comp, ' gamma');
    text(0,0.6,actext); clear actext;
    actext=strcat('N=', num2str(gm_stats(2,un_indx)), ' of total:' ,num2str(gm_stats(1,un_indx)));
    text(0,0.5,actext); clear actext;
    actext=strcat('R=', num2str(gm_stats(3,un_indx)));
    text(0,0.4,actext); clear actext;
    actext=strcat('MU=', num2str(gm_stats(4,un_indx)),'°');
    text(0,0.3,actext); clear actext;
    actext=strcat('P=', num2str(gm_stats(5,un_indx)));
    text(0,0.2,actext); clear actext;
    text(0,0.1,'Rayleigh test');
    
    %gamma phase histopgram
    subplot(2,3,((un_indx-1)*3)+2)
    bar([gmbins gmbins+360], [gmaHIST_all{un_indx} gmaHIST_all{un_indx}])
    hold on
    bar([gmbins gmbins+360], [gmaHIST_sel{un_indx} gmaHIST_sel{un_indx}])
    hold off
    axis([0 720 0 1.2*(max(gmaHIST_all{un_indx}))]);
    xlabel('CSD gamma phase (deg)'); set(gca, 'XTick',[0:180:720]); ylabel ('spike count');
    
    %theta phase histogram
    subplot(2,3,((un_indx-1)*3)+3)
    bar([thbins thbins+360], [thtHIST_all{un_indx} thtHIST_all{un_indx}])
    hold on
    bar([thbins thbins+360], [thtHIST_sel{un_indx} thtHIST_sel{un_indx}])
    hold off
    axis([0 720 0 1.2*(max(thtHIST_all{un_indx}))]);
    xlabel('theta phase (deg)'); set(gca, 'XTick',[0:180:720]); ylabel ('spike count');
    
end

%saving figures
fg_name=strcat('fig',bs.name, bs.num, bs.exp,'_gm',gm.comp,'_un', num2str(un.nums), 'sh', num2str(un.shnk), '_shuffled');
saveas(gcf,fg_name)
clear fg_name

clear un_indx input output

%saving the relevant data
datum = date;
writename=strcat(bs.name, bs.num, bs.exp, '_un', num2str(un.nums), 'sh', num2str(un.shnk),'_gammaunits_shuffle', gm.comp);
save(writename, 'bs', 'ch', 'un', 'gm', 'th', 'gmaHIST_all', 'gmaHIST_sel', 'thtHIST_all', 'thtHIST_sel', 'spks', 'spks_sh','gm_stats', 'thCYC','datum');
