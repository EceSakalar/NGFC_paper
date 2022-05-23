%This program generates gamma characterisation in single channels from each SilProbe shank  
%Based on filtered CSD and trough-peak detection
%The following analyses are done
%1. Theta phase modulation of gamma amplitude & wavelength
%2. Threshold identification
%2. Gamma phase difference as a function of theta phase 

clear all
%close all

%basic data entry
bs.name='B'; bs.num='182'; bs.unit='a'; bs.typ='sx'; 
bs.flord=[1];

ch.pyrch=14;                                        %channel number (ventralmost LFP contact is 1) for the pyramidal layer
smp.rate=2000;                                      %give the sampling rate of the channels (CSD and LFP original)

%THETA properties
%number of theta bins
th.nbins=40;th.binsize=360/th.nbins;
%theta filter properties (number of filter coefficients and filter corner frequencies in Hz)
th.ncoeff=512; th.fltlow=5; th.flthigh=12; th.dwnsmpl=5; th.rate=smp.rate/th.dwnsmpl;
%creating theta filter
th_filter=dfilt.dffir(fir1(th.ncoeff, [2.*(th.fltlow./th.rate) 2.*(th.flthigh./th.rate)],'bandpass', gausswin(th.ncoeff+1)));

%GAMMA properties
gm.comp = 'CA1slw_cln';                                                        %name of gamma component
gm.chan = 10;                                 %CSD channel for the actual gamma component for each shank (write NaN for not include)
gm.cfrq = [29.8322, 43.9963];                                              %band pass filter corner frequencies (low and high, included in the boundary)
gm.ncoeff = 1024;
gm.rate = smp.rate;
gm.nbins = 40; gm.binsize = 360/gm.nbins;                                  %give the number of theta bins 
%creating gamma filter
gm_filter=dfilt.dffir(fir1(gm.ncoeff, [2.*(gm.cfrq(1)./gm.rate) 2.*(gm.cfrq(2)./gm.rate)], 'bandpass', gausswin(gm.ncoeff+1)));
tic
%---CALCULUS---concatenating files
%The output: gmCSD_orig; thLFP_orig; thPER; 
for fl_indx=1:length(bs.flord) %FILE STEPPER
    
    %read the LFP for theta
    lfp_name=strcat(bs.name, bs.num, bs.unit, num2str(bs.flord(fl_indx)),'LFP.mat');
    load(lfp_name);
    thLFP_file = LFP(ch.pyrch,:); clear LFP lfp_name;
    
    %read the CSD for gamma; invert the CSD channel so that sink is downwards (-)
    
    csd_name=strcat(bs.name, bs.num, bs.unit, num2str(bs.flord(fl_indx)),'CSDc1.mat');
    load(csd_name);
    gmCSD_file = CSDc(gm.chan,:)*-1; clear CSDc csd_name;
    
    %read theta periods
    th_pname=strcat(bs.name, bs.num, bs.unit, num2str(bs.flord(fl_indx)),'THP.txt');
    if exist(th_pname,'file')==2
        thPER_file=dlmread(th_pname); clear th_pname;
    else
        thPER_file=[]; clear th_pname;
    end
    
    %accumulating CSD and LFP traces, theta segments and file segmentation point 
    if fl_indx==1
        gmCSD_orig = gmCSD_file;
        thLFP_orig = thLFP_file;
        flSEG = length(thLFP_file);
        thPER = thPER_file;
    else
        
        gmCSD_orig = [gmCSD_orig gmCSD_file];
        thLFP_orig = [thLFP_orig thLFP_file];
        flSEG(fl_indx)=flSEG(fl_indx-1)+length(thLFP_file)+1;
        thePER = [thPER; thPER_file+((flSEG(fl_indx-1))/smp.rate)];
    end; clear  gmCSD_file thLFP_file thPER_file
end; clear fl_indx flSEG

%---CALCULUS---analysisng gamma CSD
%gamma filtering output: gmCSD_flt gmCSD_sel gmCSD_pts
gmCSD_flt = filtfilt(gm_filter.Numerator,1,gmCSD_orig);
%generate the empty gamma selector variable
gmCSD_sel=NaN(size(gmCSD_flt));
[gmCSD_pVAL,gmCSD_pPOS] = findpeaks(gmCSD_flt);
[gmCSD_tVAL,gmCSD_tPOS] = findpeaks(gmCSD_flt*-1);
gmCSD_ptTMP = [[gmCSD_pPOS;gmCSD_pVAL;zeros(size(gmCSD_pPOS))+1] [gmCSD_tPOS;gmCSD_tVAL*-1;zeros(size(gmCSD_tPOS))-1]];
[~,sIND] = sort(gmCSD_ptTMP(1,:));
%generating a file with peaks and troughs extracted and sorted along the time
gmCSD_pts = gmCSD_ptTMP(:,sIND);
%calculate the amplitude (as SD) of the cycles (this will be the fourth row)
for je=2:(size(gmCSD_pts,2)-1)
    gmCSD_pts(4,je) = rms(gmCSD_flt(gmCSD_pts(1,je-1):gmCSD_pts(1,je+1)));
end
clear je gmCSD_pVAL gmCSD_tVAL gmCSD_pPOS gmCSD_tPOS sIND gmCSD_ptTMP
clear gmCSD_orig 
toc
'filtering and stuff done'
%generate a gamma selector trace where only the theta periods are included
for je=1:size(thPER,1)
    gmCSD_sel(round(thPER(je,1)*gm.rate):round(thPER(je,2)*gm.rate))=1;
end; clear je thPER

%---CALCULUS----theta phase
%downsampling & filtering for theta
thLFP_red=downsample(smooth(thLFP_orig,th.dwnsmpl),th.dwnsmpl,(th.dwnsmpl-1));
thLFP_flt=filtfilt(th_filter.Numerator,1,thLFP_red')';

%detecting theta troughs for pyramidal layer LFP
thLFP_THRSHLD = 0.001*(mean(abs(thLFP_flt)));
[~,thLFP_troughs]=findpeaks(thLFP_flt.*-1,'MINPEAKHEIGHT', thLFP_THRSHLD, 'MINPEAKDISTANCE', floor(th.rate/th.flthigh)-3);
thLFP_troughsT=(thLFP_troughs./th.rate);
clear thLFP_THRSHLD thLFP_red thLFP_orig 

%Get the theta phase for each gamma cycle in theta periods
%get theta phase if menaingful
for je=1:size(gmCSD_pts,2)
    
    if gmCSD_sel(gmCSD_pts(1,je))==1
        act_gmpkT=gmCSD_pts(1,je)/gm.rate;
        cyc_beg=max(thLFP_troughsT(thLFP_troughsT<=act_gmpkT));
        cyc_end=min(thLFP_troughsT(thLFP_troughsT>act_gmpkT));
        cyc_LFP=thLFP_flt((round(cyc_beg*th.rate)):round((cyc_end*th.rate))); [~,cyc_cent]=max(cyc_LFP);
        cyc_cent=(cyc_cent/th.rate)+cyc_beg;
        if  act_gmpkT<=cyc_cent;    tht_phase=180*((act_gmpkT-cyc_beg)/(cyc_cent-cyc_beg));
        else                        tht_phase=180*((act_gmpkT-cyc_cent)/(cyc_end-cyc_cent))+180;
        end
        gmCSD_pts(5,je)=tht_phase;
        clear cyc_beg cyc_end cyc_LFP cyc_cent act_gmpkT tht_phase
    else
        gmCSD_pts(5,je)=NaN;
    end
end; clear je
        
%get gamma wavelength (in ms)
for je=2:size(gmCSD_pts,2)
    gmCSD_pts(6,je)=((gmCSD_pts(1,je)-gmCSD_pts(1,je-1))./(gm.rate/1000))*2;
end
toc
clear thLFP_troughs thLFP_troughsT thLFP_flt 

%---DATA ARRANGEMENT & SORTING-----------
%sorted by theta phase, descriptive data
gmCSD_sort=gmCSD_pts;
gmCSD_sort(:,isnan(gmCSD_sort(5,:)))=[];
gmCSD_sort=(sortrows(gmCSD_sort',5))';

%---DATA SUMMARIZING-----------------------------
%get maximal amplitudes and wavelengths
gmCSD_ampmx=max(gmCSD_sort(4,:));
gmCSD_wlgmx=max(gmCSD_sort(6,:));
gmCSD_ampbins=0.1:0.2:(((ceil(gmCSD_ampmx*5))/5)-0.1);
gmCSD_wlgbins=0.5:1:ceil(gmCSD_wlgmx)+0.5;
clear gmCSD_ampmx gmCSD_wlgmx
toc
'theta filtering and data arrangement done'
%amplitude and wavelength of all cycles as a function of theta phase  
%CALCULUS---descriptive
for je=1:th.nbins
    %amplitude
    amps_bin=gmCSD_sort(4,(gmCSD_sort(5,:))>=((je-1)*th.binsize)&(gmCSD_sort(5,:))<(je*th.binsize));
    gmCSD_amp(1,je)=mean(amps_bin);
    gmCSD_amp(2,je)=std(amps_bin);
    gmCSD_amp(3,je)=median(amps_bin);
    gmCSD_amp(4,je)=prctile(amps_bin,5);
    gmCSD_amp(5,je)=prctile(amps_bin,95);
    gmCSD_amp(6,je)=prctile(amps_bin,75);
    gmCSD_amphist(:,je)=hist(amps_bin,gmCSD_ampbins);
    %wavelength
    wlgs_bin=gmCSD_sort(6,(gmCSD_sort(5,:))>=((je-1)*th.binsize)&(gmCSD_sort(5,:))<(je*th.binsize));
    gmCSD_wlghist(:,je)=hist(wlgs_bin,gmCSD_wlgbins);
end;clear je amps_bin wlgs_bin
gmCSD_wlgRS=hist((gmCSD_sort(6,gmCSD_sort(3,:)==1)),gmCSD_wlgbins);
gmCSD_wlgFL=hist((gmCSD_sort(6,gmCSD_sort(3,:)==-1)),gmCSD_wlgbins);

%%--------------------IDENTIFYING THE AMPLITUDE THRESHOLD 
%find the threshold for the actual shank
%1-find the theta bin with the smallest median amplitude
%2-fit gaussian with above peak values excluded except zeros
%3-threshold is centroid+2SD
[~,m_indx]=min(gmCSD_amp(3,:));
thrs_h=[gmCSD_ampbins' gmCSD_amphist(:,m_indx)];
thrs_h(:,3)=0;
thrs_h(:,4)=NaN;
[PEAKval,PEAKpos]=max(thrs_h(:,2));
thrs_h(1:max(PEAKpos),3)=1;
thrs_h(thrs_h(:,2)==0,3)=1;
thrs_h(1:max(PEAKpos),4)=thrs_h(1:max(PEAKpos),2);
thrs_h(thrs_h(:,2)==0,4)=0;

fOpt=fitoptions;
fOpt.Weights=thrs_h(:,3);
fRes=fit(thrs_h(:,1),thrs_h(:,2),'gauss1',fOpt);
gmCSD_thrshld=fRes.b1+((fRes.c1)*2);
clear m_indx thrs_h PEAKval PEAKpos fOpt fRes

%decide if above threshold for each gamma cycle (this will be line7 of gmCSD_pts)
for je=1:size(gmCSD_pts,2)
    if isnan(gmCSD_pts(5,je))
        gmCSD_pts(7,je)=NaN;
    elseif gmCSD_pts(4,je)>=gmCSD_thrshld
        gmCSD_pts(7,je)=1;
    else
        gmCSD_pts(7,je)=NaN;
    end
end; clear je
toc
'thresholding and gamma cycles done'

%PLOTTING--descriptive

figname=strcat(bs.name, bs.num, bs.unit,gm.comp,'amplitude');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')

subplot(2,4,1)
surf([-th.binsize/2:th.binsize:360+th.binsize/2], gmCSD_ampbins, [gmCSD_amphist(:,end) gmCSD_amphist gmCSD_amphist(:,1)], 'EdgeColor', 'interp', 'FaceColor', 'interp'); shading flat; axis tight; colormap hot; view(2);
axis([0 360 min(gmCSD_ampbins) 1.2*(max(gmCSD_amp(5,:)))]);
xlabel('theta phase (deg)'); set(gca, 'XTick',[0:180:360]); ylabel ('amplitude (rms)');

subplot(2,4,2)
errorbar ((th.binsize/2):th.binsize:(360-(th.binsize/2)),gmCSD_amp(1,:),gmCSD_amp(2,:));
actualtext=strcat(bs.name,bs.num,bs.unit,' ',gm.comp);
text(0, 0.9*(max(gmCSD_amp(5,:))), actualtext); clear actualtext;
hold on
line([0 360], [gmCSD_thrshld gmCSD_thrshld],'Color','r')
hold off
axis([0 360 min(gmCSD_ampbins) 1.2*(max(gmCSD_amp(5,:)))]);
xlabel('theta phase (deg)'); set(gca, 'XTick',[0:180:360]); ylabel ('amplitude (rms)');

%wavelength
subplot(2,4,3)
surf([-th.binsize/2:th.binsize:360+th.binsize/2], gmCSD_wlgbins, [gmCSD_wlghist(:,end) gmCSD_wlghist gmCSD_wlghist(:,1)], 'EdgeColor', 'interp', 'FaceColor', 'interp'); shading flat; axis tight; colormap hot; view(2);
axis([0 360 ((1/gm.cfrq(2)))*1000 ((1/gm.cfrq(1)))*1000]);
xlabel('theta phase (deg)'); set(gca, 'XTick',[0:180:360]); ylabel ('wavelength (ms)');
subplot(2,4,4);
hold on
line(gmCSD_wlgRS,gmCSD_wlgbins,'Color','r')
line(gmCSD_wlgFL,gmCSD_wlgbins,'Color','b')
axis([0 1.2*(max(gmCSD_wlgRS)) ((1/gm.cfrq(2)))*1000 ((1/gm.cfrq(1)))*1000]);
xlabel('count'); ylabel ('wavelength (ms)');
hold off

%PLOTTING THE SUPPLEMENTARY STUFF
%subsampling
if size(gmCSD_sort,2)>15000
    picker=randperm(size(gmCSD_sort,2),15000);
else
    picker=randperm(size(gmCSD_sort,2));
end
subplot(2,2,3)
gmCSD_resort=gmCSD_sort(:,picker);
scatter(gmCSD_resort(4,:),gmCSD_resort(6,:),5,gmCSD_resort(5,:)); colormap(gca,hsv);
xlabel('amplitude (rms)'); ylabel ('wavelength (ms)');
hold on
line([gmCSD_thrshld gmCSD_thrshld],[min(gmCSD_resort(6,:)) max(gmCSD_resort(6,:))],'Color','r')
hold off
clear gmCSD_resort picker

%PLOTTING FURTHER
subplot(2,2,4)
hist(gmCSD_sort(5,gmCSD_sort(4,:)>gmCSD_thrshld),0:th.binsize:360)
axis([0, 360, 0, max(hist(gmCSD_sort(5,gmCSD_sort(4,:)>gmCSD_thrshld),0:th.binsize:360))*1.2])
xlabel('theta phase (deg)'); ylabel ('cycle count'); set(gca, 'XTick',[-0:90:360]);

%---CLEANUP & SAVE--------------------------------------

%creating variables to be saved
%INPUT
input.base=bs;
input.gamma=gm;
input.theta=th;
input.pyrchan=ch;
input.samplerate=smp.rate;
%clear bs th gm ch smp th_filter gm_filter 

%OUTPUT
output.ampbins=gmCSD_ampbins;
output.GMtrace=gmCSD_flt;
output.GMthrsh=gmCSD_thrshld;
output.GMcycle=gmCSD_pts;
output.GMtrcSL=gmCSD_sel;
output.GMcycST=gmCSD_amp;
output.GMcycHS=gmCSD_amphist;

%clear gmCSD_ampbins gmCSD_flt gmCSD_thrshld gmCSD_pts gmCSD_sel gmCSD_amp gmCSD_amphist gmCSD_int gmCSD_intphH gmCSD_intbins gmCSDintHISTph gmCSD_intSTAT
%clear gmCSD_wlgRS gmCSD_wlgFL gmCSD_wlgbins gmCSD_wlghist gmCSD_sort

%saving the relevant data
datum = date;
writename=strcat(bs.name, bs.num, bs.unit, '_gammacomponent_', gm.comp);
save(writename, 'input', 'output', 'datum','-v7.3');

%2020.01.21 - MINOR CORRECTION - During the downsampling a phase offset is included to have the first sample at the samplin-rate-defined time-point                     
%2020.02.13 - SHORTENING COMPUTATION TIME - introducing parfor 3 lines>1 parpool 2 parfors        
%MAJOR CHANGE
%2020.03.29 - ALL CYCLES DURING THETA  INCLUDED FOR THE CALCULUS OF AMPLITUDE CORRELATIONS (originak amplitude threshold not applied           
            
           
           
 