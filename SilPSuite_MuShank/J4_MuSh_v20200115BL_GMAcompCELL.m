%This program generat es gamma coupling analysis of units to a specific
%gamma component, taking data for this component from the saved variables 
%of previously run G4 multiple shank gamma component characterization
%script.
%It uses band pass filtered CSD traces to take the gamma phase from a
%specified contact for each shank and includes cycles based on an amplitude
%threshold. Only Theta periods are considered

clear all
%close all

%basic data entry
bs.name='B'; bs.num='204'; bs.exp='b'; bs.typ='sp'; 
gm.comp = 'DEGslw'; 

%units data entry
%channels in the spike2 files
un.chan = [130:193];
%unit numbers (in kwik file)
un.nums = [4,5,7,8,11,159,197,19,21,24,25,29,33,162,200,34,40,42,44,45,49,167,195,208,219,50,52,55,57,63,64,67,69,170,171,229,232,79,80,90,181,203,233,236,249,250,98,100,106,183,184,185,119,121,123,129,130,131,186,260,140,150,156,264];                 
%shank identifier (should match the dat file order)
un.shnk = [1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,7,7,8,8,8,8];                       %shanks from which unit was isolated
%shank used for the refernece gamma oscillation extraction
un.rfsh = [3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,3,3,3,3,3,3,3,4,3,5,5,5,5,5,5,5,5,5,5,5,5,6,6,5,5,5,5,6,6,6,5,5,6,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6];

%THETA & GAMMA properties
%number of theta bins
th.nbins = 20; th .binsize = 360/th.nbins;                                  %give the number of theta bins 
gm.nbins = 20; gm.binsize = 360/gm.nbins;                                  %give the number of gamma bins 

%---CALCULUS---
%reading previously generated data from the MATLAB file
readnam=strcat(bs.name, bs.num, bs.exp,'_gammacomponent_', gm.comp, '.mat');
load(readnam); clear readnam
bs.flord=input.base.flord;
ch.pyrsh=input.pyrchan.pyrsh;
ch.pyrch=input.pyrchan.pyrch;
smp.rate=input.gamma.rate;
%extract gamma cycle information from 'OUTPUT' of read file
gmCSD_pts=output.GMcycle;

tic
%-------read and stich theta and spike data for files
for fl_indx=1:length(bs.flord) %FILE STEPPER
    
    %read the LFP for theta
    lfp_name=strcat(bs.name, bs.num, bs.exp, num2str(bs.flord(fl_indx)),'sh',num2str(ch.pyrsh),'LFP.mat');
    load(lfp_name);
    thLFP_file = LFP(ch.pyrch,:); clear LFP lfp_name;
    
    %read theta periods
    th_pname=strcat(bs.name, bs.num, bs.exp, num2str(bs.flord(fl_indx)),'THP.txt');
    if exist(th_pname,'file')==2
        thPER_file=dlmread(th_pname); clear th_pname;
    else
        thPER_file=[]; clear th_pname;
    end
    
    %read the spikes for each cell
    fl_name=strcat(bs.name, bs.num, bs.exp, num2str(bs.flord(fl_indx)),'WAW.mat'); %file name
    for un_indx = 1:length(un.nums)
        ch_name=strcat(bs.name, bs.typ, bs.num, bs.exp,'_Ch', num2str(un.chan(un_indx)));
        ch_load=load(num2str(fl_name),num2str(ch_name));
        spks_fl{un_indx} =(ch_load.(ch_name).('times')); clear ch_name ch_load;
    end
    
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
        for un_indx = 1:length(un.nums)
            spks{un_indx} = [spks{un_indx}; spks_fl{un_indx}+(flDUR(fl_indx-1)/smp.rate)];
        end
    end
    clear thLFP_file thPER_file spks_fl fl_name 
end; clear fl_indx flDUR

%---CALCULUS----theta phase
%downsampling & filtering for theta
th.nbins=input.theta.nbins; th.binsize=360/th.nbins;
th.ncoeff=input.theta.ncoeff; th.fltlow=input.theta.fltlow; th.flthigh=input.theta.flthigh; th.dwnsmpl=input.theta.dwnsmpl; th.rate=smp.rate/th.dwnsmpl;
th_filter=dfilt.dffir(fir1(th.ncoeff, [2.*(th.fltlow./th.rate) 2.*(th.flthigh./th.rate)],'bandpass', gausswin(th.ncoeff+1)));
thLFP_red=downsample(smooth(thLFP_orig,th.dwnsmpl),th.dwnsmpl,th.dwnsmpl-1);
thLFP_flt=filtfilt(th_filter.Numerator,1,thLFP_red')';

%detecting theta troughs for pyramidal layer LFP
thLFP_THRSHLD = 0.001*(mean(abs(thLFP_flt)));
[~,thLFP_troughs]=findpeaks(thLFP_flt.*-1,'MINPEAKHEIGHT', thLFP_THRSHLD, 'MINPEAKDISTANCE', floor(th.rate/th.flthigh)-3);
thLFP_troughsT=(thLFP_troughs./th.rate);
clear thLFP_red thLFP_THRSHL thLFP_orig


'theta filtering and spikes done'
toc


%---CALCULUS----extract spike data
for un_indx = 1:length(un.nums)
    %col 1 is spike time in seconds, col 2 is spike time in sample (gamma rate based)
    spks{un_indx}(:,2) = spks{un_indx}(:,1).*smp.rate;
    %col 3 is 1 if spike is in the theta periods 
    spks{un_indx}(:,3) = NaN;
    for je=1:size(thPER,1)
        spks{un_indx}((spks{un_indx}(:,1)>=thPER(je,1)&spks{un_indx}(:,1)<=thPER(je,2)),3)=1;
    end; clear je
    %col 4 is theta phase (degrees)
    %col 5 is closest gamma cycle time (sec)
    for sp_indx=1:size(spks{un_indx},1)
        if isnan(spks{un_indx}(sp_indx,3))
            spks{un_indx}(sp_indx,4:9)=NaN;
        else
            %4 theta phase
            sp_time=spks{un_indx}(sp_indx,1);
            [cyc_beg, cyc_begi]=max(thLFP_troughsT(thLFP_troughsT<=sp_time));
            cyc_end = thLFP_troughsT(cyc_begi+1);
            cyc_LFP = thLFP_flt(thLFP_troughs(cyc_begi):thLFP_troughs(cyc_begi+1));
            [~,cyc_cnti]=max(cyc_LFP); cyc_cnt=(cyc_cnti/th.rate)+cyc_beg;
            if  sp_time<=cyc_cnt;    spks{un_indx}(sp_indx,4)=180*((sp_time-cyc_beg)/(cyc_cnt-cyc_beg));
            else                     spks{un_indx}(sp_indx,4)=180*((sp_time-cyc_cnt)/(cyc_end-cyc_cnt))+180;
            end
            clear cyc_cnti cyc_begi cyc_cnt cyc_cnti cyc_LFP cyc_end cyc_beg
            
            %5-8 closest gamma cycle properties
            gm_shank=gmCSD_pts{un.rfsh(un_indx)};
            [~,cyc_ind]=min(abs((gm_shank(1,:)./smp.rate)-sp_time));
            spks{un_indx}(sp_indx,5) = gm_shank(1,cyc_ind)/smp.rate;
            spks{un_indx}(sp_indx,6) = gm_shank(3,cyc_ind);
            spks{un_indx}(sp_indx,7) = gm_shank(4,cyc_ind);
            spks{un_indx}(sp_indx,8) = gm_shank(7,cyc_ind);
            
            %getting gamma phase
            if gm_shank(3,cyc_ind)==1
                cyc_cnti=cyc_ind;
            elseif sp_time>=(gm_shank(1,cyc_ind)/smp.rate)
                cyc_cnti=cyc_ind+1;
            else
                cyc_cnti=cyc_ind-1;
            end
            % get cycle begin
            if gm_shank(3,cyc_cnti-1)==-1
                cyc_beg=gm_shank(1,(cyc_cnti-1))/smp.rate;
            else
                'error occured: gamma cycle begins with peak'
            end
            % get cycle center
            if gm_shank(3,cyc_cnti)==1
                cyc_cnt=gm_shank(1,(cyc_cnti))/smp.rate;
            else
                'error occured: gamma cycle center is trough'
            end
            % get cycle end
            if gm_shank(3,cyc_cnti+1)==-1
                cyc_end=gm_shank(1,(cyc_cnti+1))/smp.rate;
            else
                'error occured: gamma cycle begins with peak'
            end
            
            %get gamma phase
            if  sp_time<=cyc_cnt;    gma_phase=180*((sp_time-cyc_beg)/(cyc_cnt-cyc_beg));
            else;                    gma_phase=180*((sp_time-cyc_cnt)/(cyc_end-cyc_cnt))+180;
            end
            spks{un_indx}(sp_indx,9) = gma_phase;
            clear sp_time cyc_ind cyc_cnti cyc_beg cyc_cnt cyc_end gma_phase
        end
    end; clear sp_indx gm_shank
    
    strcat('unit ', num2str(un.nums(un_indx)), 'on shank ', num2str(un.shnk(un_indx)), ' is done') 
    toc
    
end; clear un_indx
clear thLFP_flt thLFP_troughs thLFP_troughsT thPER
clear thLFP_THRSHLD

%%----CALCULUS---statistical summary
thbins=(th.binsize/2):th.binsize:(360-(th.binsize/2));
gmbins=(gm.binsize/2):gm.binsize:(360-(gm.binsize/2));
for un_indx = 1:length(un.nums)
    thtHIST_all{un_indx}=hist(spks{un_indx}(:,4),thbins);
    thtHIST_sel{un_indx}=hist(spks{un_indx}(spks{un_indx}(:,8)==1,4),thbins);
    gmaHIST_all{un_indx}=hist(spks{un_indx}(:,9),gmbins);
    gmaHIST_sel{un_indx}=hist(spks{un_indx}(spks{un_indx}(:,8)==1,9),gmbins);
    
    gm_stats(1,un_indx)=nansum(spks{un_indx}(:,3));
    gm_stats(2,un_indx)=nansum(spks{un_indx}(:,8));
    gm_stats(3,un_indx)=circ_r(deg2rad((spks{un_indx}(spks{un_indx}(:,8)==1,9))));
    gm_stats(4,un_indx)=rad2deg(circ_mean(deg2rad((spks{un_indx}(spks{un_indx}(:,8)==1,9)))));
    gm_stats(5,un_indx)=circ_rtest(deg2rad((spks{un_indx}(spks{un_indx}(:,8)==1,9))));
    if gm_stats(4,un_indx)<0
        gm_stats(4,un_indx)=gm_stats(4,un_indx)+360;
    end
end
clear un_indx datum

%%----DISPLAY--------
for un_indx = 1:length(un.nums)
    figindex=floor((un_indx-1)/6)+1;
    plotrow=un_indx-((figindex-1)*6);
    if plotrow==1
        figname=strcat(bs.name, bs.num, bs.exp, gm.comp,'_cellcoupling_',num2str(figindex));
        figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
    end
    %textual results
    subplot(6,3,((plotrow-1)*3)+1)
    axis off
    actext=strcat(bs.name, bs.num, bs.exp, ' shank', num2str(un.shnk(un_indx)),' unit', num2str(un.nums(un_indx)));
    text(0,1.2,actext); clear actext;
    actext=strcat('coupling to: ', gm.comp, ' gamma');
    text(0,1,actext); clear actext;
    actext=strcat('N=', num2str(gm_stats(2,un_indx)), ' of total:' ,num2str(gm_stats(1,un_indx)));
    text(0,0.8,actext); clear actext;
    actext=strcat('R=', num2str(gm_stats(3,un_indx)));
    text(0,0.6,actext); clear actext;
    actext=strcat('MU=', num2str(gm_stats(4,un_indx)),'°');
    text(0,0.4,actext); clear actext;
    actext=strcat('P=', num2str(gm_stats(5,un_indx)));
    text(0,0.2,actext); clear actext;
    text(0,0,'Rayleigh test');
    
    %gamma phase histopgram
    subplot(6,3,((plotrow-1)*3)+2)
    bar([gmbins gmbins+360], [gmaHIST_all{un_indx} gmaHIST_all{un_indx}])
    hold on
    bar([gmbins gmbins+360], [gmaHIST_sel{un_indx} gmaHIST_sel{un_indx}])
    hold off
    axis([0 720 0 1.2*(max(gmaHIST_all{un_indx}))]);
    xlabel('CSD gamma phase (deg)'); set(gca, 'XTick',[0:180:720]); ylabel ('spike count');
    
    %theta phase histogram
    subplot(6,3,((plotrow-1)*3)+3)
    bar([thbins thbins+360], [thtHIST_all{un_indx} thtHIST_all{un_indx}])
    hold on
    bar([thbins thbins+360], [thtHIST_sel{un_indx} thtHIST_sel{un_indx}])
    hold off
    axis([0 720 0 1.2*(max(thtHIST_all{un_indx}))]);
    xlabel('theta phase (deg)'); set(gca, 'XTick',[0:180:720]); ylabel ('spike count');
    
    %saving figures
    if plotrow==6 || un_indx == length(un.nums)
        fg_name=strcat('fig',bs.name, bs.num, bs.exp,'_gm',gm.comp,'_cellsa_figure', num2str(figindex));
        saveas(gcf,fg_name)
        clear fg_name
    end
end
clear un_indx input output

%saving the relevant data
datum = date;
writename=strcat(bs.name, bs.num, bs.exp, '_gammaunits_', gm.comp);
save(writename, 'bs', 'ch', 'un', 'gm', 'th', 'gmaHIST_all', 'gmaHIST_sel', 'thtHIST_all', 'thtHIST_sel', 'spks', 'gm_stats', 'datum');

            


            
        

        
           
            
           
           
 