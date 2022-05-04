function [ output_args ] = cellGmComp( input_args )
%This program generates gamma coupling analysis of units to a specific
%gamma component, taking data for this component from the saved variables 
%of previously run G4 multiple shank gamma component characterization
%script.
%It uses band pass filtered CSD traces to take the gamma phase from a
%specified contact for each shank and includes cycles based on an amplitude
%threshold. Only Theta periods are considered

%---CALCULUS---
%reading previously generated data from the MATLAB file
readnam=strcat(bs.name, bs.num, '_', bs.exp,'_gammacomponent_', gm.comp, '.mat');
load(readnam); 
clear readnam
bs.flord=input.base.flord;
ch.pyrch=input.pyrchan.pyrch;
smp.rate=input.gamma.rate;
%extract gamma cycle information from 'OUTPUT' of read file
gmCSD_pts=output.GMcycle;

tic
%-------read and stich theta and spike data for files

%read the LFP for theta
lfp_name=strcat(bs.name, bs.num, '_', bs.exp, '_LFP.mat');
load(lfp_name);
thLFP_file = LFP(ch.pyrch,:); 
clear LFP lfp_name;

%read theta periods
th_pname=strcat(bs.name, bs.num, '_', bs.exp, '_THP.txt');
if exist(th_pname,'file')==2
    thPER_file=dlmread(th_pname); 
    clear th_pname;
else
    thPER_file=[]; clear th_pname;
end

%read the spikes for each cell
fl_name=strcat(bs.name, bs.num, '_', bs.exp, '_WAW.mat'); %file name
for un_indx = 1:length(un.nums)
    if un.shnk(un_indx)==0
        ch_name=strcat('Gl_unt_', num2str(un.chan(un_indx).','%02d'));
    elseif un.shnk(un_indx)==1
        ch_name=strcat('Unit_', num2str(un.chan(un_indx).','%03d'));
    end
    ch_load=load(num2str(fl_name),num2str(ch_name));
    spks_fl{un_indx} =(ch_load.(ch_name).('times')); clear ch_name ch_load;
end

flDUR = length(thLFP_file);
thLFP_orig = thLFP_file;
thPER = thPER_file;
spks=spks_fl;

clear thLFP_file thPER_file spks_fl fl_name
clear flDUR

%---CALCULUS----theta phase
%downsampling & filtering for theta
th.nbins=input.theta.nbins; 
th.binsize=360/th.nbins;
th.ncoeff=input.theta.ncoeff; 
th.fltlow=input.theta.fltlow; 
th.flthigh=input.theta.flthigh; 
th.dwnsmpl=input.theta.dwnsmpl; 
th.rate=smp.rate/th.dwnsmpl;
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
            [~,cyc_cnti]=max(cyc_LFP); 
            cyc_cnt=(cyc_cnti/th.rate)+cyc_beg;
            if  sp_time<=cyc_cnt   
                spks{un_indx}(sp_indx,4)=180*((sp_time-cyc_beg)/(cyc_cnt-cyc_beg));
            else
                spks{un_indx}(sp_indx,4)=180*((sp_time-cyc_cnt)/(cyc_end-cyc_cnt))+180;
            end
            clear cyc_cnti cyc_begi cyc_cnt cyc_cnti cyc_LFP cyc_end cyc_beg
            
            %5-8 closest gamma cycle properties
            [~,cyc_ind]=min(abs((gmCSD_pts(1,:)./smp.rate)-sp_time));
            spks{un_indx}(sp_indx,5) = gmCSD_pts(1,cyc_ind)/smp.rate;
            spks{un_indx}(sp_indx,6) = gmCSD_pts(3,cyc_ind);
            spks{un_indx}(sp_indx,7) = gmCSD_pts(4,cyc_ind);
            spks{un_indx}(sp_indx,8) = gmCSD_pts(7,cyc_ind);
            
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
            if  sp_time<=cyc_cnt
                gma_phase=180*((sp_time-cyc_beg)/(cyc_cnt-cyc_beg));
            else                    
                gma_phase=180*((sp_time-cyc_cnt)/(cyc_end-cyc_cnt))+180;
            end
            spks{un_indx}(sp_indx,9) = gma_phase;
            clear sp_time cyc_ind cyc_cnti cyc_beg cyc_cnt cyc_end gma_phase
        end
    end
    clear sp_indx 
    
    strcat('unit ', num2str(un.nums(un_indx)), 'on shank ', num2str(un.shnk(un_indx)), ' is done') 
    toc
end    
clear un_indx
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

end

