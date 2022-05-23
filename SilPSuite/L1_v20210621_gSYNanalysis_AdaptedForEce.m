%This script analyses gSYN trace as a function of theta phase
clear all; %close all

bs.name='ES55_d6_';
bs.cell='glass01';

ch.pyrch=13;                                        %channel number (ventralmost LFP contact is 1) for the pyramidal layer
smp.rate=2000;                                      %give the sampling rate of the channels (LFP channel in original LFP.mat file)

%THETA properties
%number of theta bins
th.nbins=40;th.binsize=360/th.nbins;
%theta filter properties (number of filter coefficients and filter corner frequencies in Hz)
th.ncoeff=512; th.fltlow=5; th.flthigh=12; th.dwnsmpl=5; th.rate=smp.rate/th.dwnsmpl;
%creating theta filter
th_filter=dfilt.dffir(fir1(th.ncoeff, [2.*(th.fltlow./th.rate) 2.*(th.flthigh./th.rate)],'bandpass', gausswin(th.ncoeff+1)));

%POSTSYN CONDUCTANCE (Gs) properties
gs.trace = 'gSyn_tr1'; 
gs.cycev = 'gSyn_ev1';

%read the LFP for theta
lfp_name=strcat(bs.name,'LFP.mat');
load(lfp_name);
thLFP_file = LFP(ch.pyrch,:); clear LFP lfp_name;
%time axis of the LFP trace (original)
for je=1:length(thLFP_file); fl_thT(je)=(je*(1/smp.rate)); end

%read the gSyn trace
fl_name=strcat(bs.name, 'Gsyn.mat');
ch_name=gs.trace; ch_load=load(num2str(fl_name),num2str(ch_name));
fl_inval=ch_load.(ch_name).('interval'); fl_start=ch_load.(ch_name).('start'); fl_length=ch_load.(ch_name).('length');
%original time axis of the gSyn trace
for je=1:fl_length; fl_time(je)=((je-1)*fl_inval)+fl_start; end
fl_gsynIN=ch_load.(ch_name).('values');
clear ch_load ch_name je

%resample timeseries of gSyn to match the theta LFP
fl_gsynTS=resample((timeseries(fl_gsynIN,fl_time)),fl_thT);
fl_gsyn=fl_gsynTS.Data;
clear fl_gsynIN fl_gsynTS fl_inval fl_length fl_start fl_thT fl_time

%read the cycle-by cycle gSYN event starts
ch_name=gs.cycev; ch_load=load(num2str(fl_name),num2str(ch_name));
fl_gevn=ch_load.(ch_name).('times');
clear fl_name ch_name ch_load

%read theta periods
th_pname=strcat(bs.name, 'THP.txt');
if exist(th_pname,'file')==2
    thPER_file=dlmread(th_pname); clear th_pname;
else
    thPER_file=[]; clear th_pname;
end

thLFP_orig = thLFP_file;
flSEG = length(thLFP_file);
thPER = thPER_file;
gsTRC = fl_gsyn';
gsEVN = fl_gevn;
clear  thLFP_file thPER_file fl_gsyn fl_gevn

%THETA phase extraction for all theta periods
%filtering the LFP for theta (also downsampling)
thLFP_red=downsample(smooth(thLFP_orig,th.dwnsmpl),th.dwnsmpl,(th.dwnsmpl-1));
thLFP_flt=filtfilt(th_filter.Numerator,1,thLFP_red')';

%detecting theta troughs for pyramidal layer LFP
thLFP_THRSHLD = 0.001*(mean(abs(thLFP_flt)));
[~,thLFP_troughs]=findpeaks(thLFP_flt.*-1,'MINPEAKHEIGHT', thLFP_THRSHLD, 'MINPEAKDISTANCE', floor(th.rate/th.flthigh)-3);
thLFP_troughsT=(thLFP_troughs./th.rate);
[~,thLFP_peaks]=findpeaks(thLFP_flt,'MINPEAKHEIGHT', thLFP_THRSHLD, 'MINPEAKDISTANCE', floor(th.rate/th.flthigh)-3);
thLFP_peaksT=(thLFP_peaks./th.rate);
clear thLFP_THRSHLD thLFP_red %thLFP_orig

%timescale
tmORG=[1:size(thLFP_orig,2)].*(1/smp.rate);

%extracting theta phase for all the points within the 
for je=1:size(thLFP_peaksT,1);   thPHSph_p(je)=180; thPHStm_p(je)=thLFP_peaksT(je);  end; clear je
for je=1:size(thLFP_troughsT,1); thPHSph_t(je)=0;   thPHStm_t(je)=thLFP_troughsT(je);end; clear je
thPHStm_pt_t=[min(tmORG) thPHStm_p thPHStm_t max(tmORG)]; thPHSph_pt_t=[0 thPHSph_p thPHSph_t 0];
[thPHStm_pt,thPHSix_pt] = sort(thPHStm_pt_t); thPHSph_pt=thPHSph_pt_t(thPHSix_pt);
for je=1:size(tmORG,2); thPHSph_all(je)=NaN; end; clear je;

%interpolating theta phases for all points
for je=1:(size(thPHStm_pt,2)-1)
    %if the starting point is a trough

    indx_st=round(thPHStm_pt(je)*smp.rate);
    indx_en=round(thPHStm_pt(je+1)*smp.rate);
    if thPHSph_pt(je)==0
        thPHSph_all(indx_st)=0;
        phse_st=0; phse_en=thPHSph_pt(je+1);
    elseif thPHSph_pt(je)==180
        thPHSph_all(indx_st)=180;
        phse_st=180;
        if thPHSph_pt(je+1)==0
           phse_en=360; 
        else
           phse_en=thPHSph_pt(je+1);
        end
    end
    for ka=(indx_st+1):indx_en
        thPHSph_all(ka)=((ka-indx_st).*(phse_en-phse_st)/(indx_en-indx_st))+phse_st;
    end; clear ka
    clear phse_st phse_en indx_st indx_en
end; clear je
clear thLFP_peaks thLFP_peaksT thLFP_troughs thLFP_troughsT thPHSix_pt thPHSph_p thPHSph_t thPHSph_pt thPHSph_pt_t thPHStm_p thPHStm_t thPHStm_pt thPHStm_pt_t
            

%taking the gsEVN events only for theta periods
for je=1:size(thPER,1)
    if je==1
        gsEVNTtht=gsEVN(gsEVN>=thPER(je,1)&gsEVN<=thPER(je,2));
    else
        gsEVNTtmp=gsEVN(gsEVN>=thPER(je,1)&gsEVN<=thPER(je,2));
        gsEVNTtht=[gsEVNTtht; gsEVNTtmp]; clear gsEVNTtmp
    end
end; clear je

%Figure displaying analysis results
figname=strcat(bs.name, 'gSYN_analysis');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
%Theta figure accessories
th.phasebins=[th.binsize/2:th.binsize:360-(th.binsize/2)];
%ANALYSIS one - theta phase histogram of cycle inhibition onset
gsEVNTtht_phsALL=thPHSph_all(round(gsEVNTtht.*smp.rate));
gsEVNTtht_phsHST=hist(gsEVNTtht_phsALL,th.phasebins);

%theta phase histograms of 
subplot(3,2,1)
bar(th.phasebins,gsEVNTtht_phsHST)
axis([0 360 0 1.2*(max(gsEVNTtht_phsHST))]);
xlabel('theta phase (deg)'); set(gca, 'XTick',[0:90:360]); ylabel ('count (gSYN event onset)');

%extracting theta cycles, postsyn currents and theta phases aligned to the gSYN event onset
%exwin is the extraction halwindow in seconds
exwin=0.15;
for je=1:size(gsEVNTtht,1)
    centro=round(gsEVNTtht(je)*smp.rate);
    trigTHTindi(je,:) = thLFP_orig(centro-(exwin*smp.rate):centro+(exwin*smp.rate))-mean(thLFP_orig(centro-(exwin*smp.rate):centro+(exwin*smp.rate))); %baseline corrected
    trigGSYindi(je,:) = gsTRC(centro-(exwin*smp.rate):centro+(exwin*smp.rate));
    trigTHPindi(je,:) = thPHSph_all(centro-(exwin*smp.rate):centro+(exwin*smp.rate));
end; clear centro je

%theta oscillation
subplot(3,2,2)
line([(-1*exwin):(1/smp.rate):exwin],trigTHTindi, 'Color', [0.8 0.8 0.8]','LineWidth',0.1);
line([(-1*exwin):(1/smp.rate):exwin],mean(trigTHTindi), 'Color', [0 0 0]','LineWidth',1);
axis([-exwin exwin 1.2*min(mean(trigTHTindi)) 1.2*max(mean(trigTHTindi))])
xlabel('time (s)'); ylabel ('LFP (mV)');
clear thLFP_orig thLFP_flt flSEG   

%gSYN trace
subplot(3,2,4)
line([(-1*exwin):(1/smp.rate):exwin],trigGSYindi, 'Color', [0.8 0.8 0.8]','LineWidth',0.1);
line([(-1*exwin):(1/smp.rate):exwin],mean(trigGSYindi), 'Color', [0 0 0]','LineWidth',1);
axis([-exwin exwin 0 max(max(trigGSYindi))])
xlabel('time (s)'); ylabel ('gSYN (AU)');

%THETA phase
subplot(3,2,6)
line([(-1*exwin):(1/smp.rate):exwin],trigTHPindi, 'Color', [0.8 0.8 0.8]','LineWidth',0.1);
for je=1:size(trigTHPindi,2)
    ka=rad2deg(circ_mean(deg2rad(trigTHPindi(:,je))));
    if ka<0; trigTHPindi_mean(je)=ka+360; else; trigTHPindi_mean(je)=ka; end
    clear ka
end; clear je
line([(-1*exwin):(1/smp.rate):exwin],trigTHPindi_mean, 'Color', [0 0 0]','LineWidth',1);
axis([-exwin exwin 0 360])
xlabel('time (s)'); ylabel ('theta phase (°)');  set(gca, 'YTick',[0:90:360]);

%CALCULUS THTphase dependent current gSYN
%THETA periods selected for perids when transients were present
thPERcell=thPER;
for je=1:size(thPERcell,1)
    if thPER(je,1)>max(gsEVN); thPERcell(je,:)=NaN; end
    if thPER(je,2)<min(gsEVN); thPERcell(je,:)=NaN; end
end; clear je
%limiting    
gsTRCtht=zeros(2,size(gsTRC,2)); gsTRCtht(:,:)=NaN;
for je=1:size(thPERcell,1)
    if ~isnan(thPERcell(je,1))
        for ka=round(thPER(je,1)*smp.rate):round(thPER(je,2)*smp.rate)
            gsTRCtht(1,ka)=thPHSph_all(ka);
            gsTRCtht(2,ka)=gsTRC(ka);
        end
    end; clear ka
end; clear je
%binning
for je=1:40
    bingSYNthtphase{1,je}=[((je*9)-9) (je*9)]; %binning boundaries
    bingSYNthtphase{2,je}=find(gsTRCtht(1,:)>=((je*9)-9)&gsTRCtht(1,:)<(je*9));%bin indices 
    bingSYNthtphase{3,je}=gsTRCtht(2,bingSYNthtphase{2,je});% synaptic current amplitude
    bingSYNmean(je)=mean(bingSYNthtphase{3,je});
end; clear je
%hsitograms for the surface plot
% for je=1:40
%     bingSYNhist2D(:,je)=hist(bingSYNthtphase{3,je},[0.05:0.1:2.95]);
% end
%PLOTTING theta phase dependent gSYN current dynamics

subplot(3,2,3)
% %surface plot 
% surf([-4.5:9:364.5],[0.05:0.1:2.95],[bingSYNhist2D(:,end) bingSYNhist2D bingSYNhist2D(:,1)], 'EdgeColor', 'interp', 'FaceColor', 'interp'); shading flat; axis tight; colormap hot; view(2);
% axis([0 360 0 3]); caxis([0 max(max(bingSYNhist2D(2:end,:)))])
% xlabel('theta phase (deg)'); set(gca, 'XTick',[0:180:360]); ylabel ('amplitude (rms)');
%mean current
line(th.phasebins, bingSYNmean)
axis([0 360 0 1.2*max(bingSYNmean')])
xlabel('theta phase (°)'); ylabel('mean gsyn') ;  set(gca, 'XTick',[0:90:360]);

%textual results
gsyn_stats(1)=size(gsEVNTtht,1);
gsyn_stats(2)=circ_r(deg2rad(gsEVNTtht_phsALL));    
mphs=rad2deg(circ_mean(deg2rad(gsEVNTtht_phsALL)));
if mphs<0
    gsyn_stats(3)= mphs+360;
else
    gsyn_stats(3)= mphs;
end; clear mphs
gsyn_stats(4)=circ_rtest(deg2rad((gsEVNTtht_phsALL)));

gsyn_stats(5)=th.phasebins(find(bingSYNmean==min(bingSYNmean)));
gsyn_stats(6)=th.phasebins(find(bingSYNmean==max(bingSYNmean)));
for je=1:40
    if je==1
        bingSYNmeanRT(je)=(bingSYNmean(2)-bingSYNmean(40))/18;
    elseif je==40
        bingSYNmeanRT(je)=(bingSYNmean(1)-bingSYNmean(39))/18;
    else
        bingSYNmeanRT(je)=(bingSYNmean(je+1)-bingSYNmean(je-1))/18;
    end
end; clear je
gsyn_stats(7)=th.phasebins(find(bingSYNmeanRT==max(bingSYNmeanRT)));

subplot(3,2,5)
axis off
actext=strcat(bs.name, ' cell:' ,bs.cell);
text(0,1,actext); clear actext;
actext=strcat('evnt n: ', num2str(gsyn_stats(1)));
text(0,0.8,actext); clear actext;
actext=strcat('evnt r/P: ', num2str(gsyn_stats(2)),' / ',num2str(gsyn_stats(4)));
text(0,0.6,actext); clear actext;
actext=strcat('evnt mu: ', num2str(gsyn_stats(3)));
text(0,0.4,actext); clear actext;
text(0,0.2,'min / max / maxrate (gsyn)')
actext=strcat(num2str(gsyn_stats(5)),' / ',num2str(gsyn_stats(6)),' / ',num2str(gsyn_stats(7)));
text(0,0,actext); clear actext;

%saving figure
fg_name=strcat('fig',bs.name, '_cell',bs.cell,'_gSYNanalys');
saveas(gcf,fg_name); clear fg_name figname 

%generating saveable varables
gsevnt.all=gsEVN;
gsevnt.thtonly=gsEVNTtht;
gsevnt.thtphases=gsEVNTtht_phsALL;
gsevnt.thtphasehisto=gsEVNTtht_phsHST;

thtper.all=thPER;
thtper.cell=thPERcell;

timseri.time=tmORG;
timseri.gsynall=gsTRC;
timseri.gsyntht=gsTRCtht;
timseri.tphall=thPHSph_all;

summar.binGmean=bingSYNmean;
summar.binGmeanrate=bingSYNmeanRT;
summar.binTHTphase=bingSYNthtphase;
summar.trgwindow=exwin;
summar.trgGall=trigGSYindi;
summar.trgTHphase=trigTHPindi;
summar.trgTHtrace=trigTHTindi;
summar.stati=gsyn_stats;


clear gsEVN gsEVNTtht gsEVNTtht_phsALL gsEVNTtht_phsHST thPER thPERcell tmORG gsTRC gsTRCtht thPHSph_all
clear bingSYNmean bingSYNmeanRT bingSYNthtphase exwin trigGSYindi trigTHPindi trigTHTindi trigTHPindi_mean gsyn_stats

%saving the collected data
datum = date;
writename=strcat(bs.name, 'gSYNfor_cell_', bs.cell);
save(writename, 'bs', 'ch', 'datum', 'gs', 'gsevnt', 'smp', 'summar', 'th','thtper','timseri','-v7.3');

