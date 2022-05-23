%This script analyses gSYN trace as a function of gamma phase in the cycle
clear all; %close all

bs.name='B'; bs.num='207'; bs.unit='b'; bs.typ='sp'; bs.cell='s5u050';
bs.flord=[1];
%shank number for gamma (leave empty if not applicable)
bs.shnk =['sh4'];

%POSTSYN CONDUCTANCE (Gs) properties
gs.trace = 166; 
gs.cycev = 168;
gs.spkst = 146; 

%GAMMA properties
gm.comp = 'CA1mid';                                                        %name of gamma component
gm.chan = 7;                                 %CSD channel for the actual gamma component for each shank (write NaN for not include)
gm.cfrq = [55.2359, 89.5522];                                              %band pass filter corner frequencies (low and high, included in the boundary)
gm.ncoeff = 1024;
gm.rate = 10000;
gm.nbins = 40; gm.binsize = 360/gm.nbins;                                  %give the number of theta bins 
%creating gamma filter
gm_filter=dfilt.dffir(fir1(gm.ncoeff, [2.*(gm.cfrq(1)./gm.rate) 2.*(gm.cfrq(2)./gm.rate)], 'bandpass', gausswin(gm.ncoeff+1)));
                                   
%CALCULUIS START - ONLY SINGLE FILE RECORDINGS PERMITTED!
%read the CSD for gamma
csd_name=strcat(bs.name, bs.num, bs.unit, num2str(bs.flord),bs.shnk,'CSD.mat');
load(csd_name);
CSD_file = (CSD(gm.chan,:)).*-1; clear CSD csd_name;
%time axis of the LFP trace (original)
for je=1:length(CSD_file); fl_csdtim(je)=(je*(1/gm.rate)); end; clear je

%read the gSyn trace
fl_name=strcat(bs.name, bs.num, bs.unit, num2str(bs.flord),'Gsyn.mat');
ch_name=strcat(bs.name, bs.typ, bs.num, bs.unit,'_Ch', num2str(gs.trace)); ch_load=load(num2str(fl_name),num2str(ch_name));
fl_inval=ch_load.(ch_name).('interval'); fl_start=ch_load.(ch_name).('start'); fl_length=ch_load.(ch_name).('length');
%original time axis of the gSyn trace
for je=1:fl_length; fl_time(je)=((je-1)*fl_inval)+fl_start; end
fl_gsynIN=ch_load.(ch_name).('values');
clear ch_load ch_name je fl_name

%resample timeseries of gSyn to match the theta LFP
fl_gsynTS=resample((timeseries(fl_gsynIN,fl_time)),fl_csdtim);
fl_gsyn=fl_gsynTS.Data;
clear fl_gsynIN fl_gsynTS fl_inval fl_length fl_start fl_thT fl_time

%read the spikes for the cell
fl_name=strcat(bs.name, bs.num, bs.unit, num2str(bs.flord),'WAW.mat'); %file name
ch_name=strcat(bs.name, bs.typ, bs.num, bs.unit,'_Ch', num2str(gs.spkst));
ch_load=load(num2str(fl_name),num2str(ch_name));
gSyn_spks =(ch_load.(ch_name).('times')); clear ch_name ch_load fl_name;

%read theta periods
th_pname=strcat(bs.name, bs.num, bs.unit, num2str(bs.flord),'THP.txt');
th_pers=dlmread(th_pname); clear th_pname;

%streamline variable names
CSD_org = CSD_file;
CSD_tim = fl_csdtim;
gSyn_trac = fl_gsyn';
clear CSD_file fl_csdtim fl_gsyn

%gamma filtering output: gmCSD_flt gmCSD_sel gmCSD_pts
CSD_gma = filtfilt(gm_filter.Numerator,1,CSD_org);

%spikes for theta periods only!
gSyn_thtspks=[];
for je=1:size(th_pers,1)
    if isempty(gSyn_thtspks)
        gSyn_thtspks = gSyn_spks(gSyn_spks>th_pers(je,1)&gSyn_spks<th_pers(je,2));
    else
        gSyn_thtspks = [gSyn_thtspks; gSyn_spks(gSyn_spks>th_pers(je,1)&gSyn_spks<th_pers(je,2))];
    end
end
clear je

%spike triggered time axis, gma traces, gma phase traces, and gSyn traces
sptrg.time=[-0.05:1/gm.rate:0.1];
for je=1:length(gSyn_thtspks)
    smp_frst=(round((gSyn_thtspks(je)*gm.rate)))-(0.05*gm.rate);
    smp_zero=round((gSyn_thtspks(je)*gm.rate));
    smp_last=(round((gSyn_thtspks(je)*gm.rate)))+(0.1*gm.rate);
    sptrg.gSyn(je,:) = gSyn_trac(smp_frst:smp_last);
    sptrg.gmatrc(je,:) = CSD_gma(smp_frst:smp_last);
    %gamma phase in the spike triggered trace
    [~,p_loc]=findpeaks(sptrg.gmatrc(je,:));
    [~,t_loc]=findpeaks(sptrg.gmatrc(je,:)*-1);
    t_loc(t_loc<min(p_loc))=[]; t_loc(t_loc>max(p_loc))=[];
    phsL=-180-(((length(p_loc(p_loc<(0.05*gm.rate)+1)))-1)*360);
    sptrg.gmaphs(je,1:p_loc(1))=phsL;
    for ka=1:length(p_loc)-1
        for el=((p_loc(ka))+1):t_loc(ka)
            sptrg.gmaphs(je,el)=phsL+((ka-1)*360)+((el-(p_loc(ka)))*(180/(t_loc(ka)-p_loc(ka))));
        end; clear el
        for el=((t_loc(ka))+1):p_loc(ka+1)
            sptrg.gmaphs(je,el)=phsL+180+((ka-1)*360)+((el-(t_loc(ka)))*(180/(p_loc(ka+1)-t_loc(ka))));
        end; clear el
    end; clear ka
    for el=(max(p_loc)+1):length(sptrg.gmatrc(je,:))
        sptrg.gmaphs(je,el)=sptrg.gmaphs(je,max(p_loc));
    end; clear el
    clear smp_frst smp_zero smp_last p_loc t_loc phsL
end; clear je

%selection based on ghSyn not crossing a threshold at the spike time
sptrgS.gmaphs=sptrg.gmaphs; sptrgS.gmatrc=sptrg.gmatrc; sptrgS.gSyn = sptrg.gSyn;
for je=1:size(sptrgS.gSyn,1)
    if sptrgS.gSyn(je,((0.05*gm.rate)+1))>=0.35
       sptrgS.gSyn(je,:)=NaN;
       sptrgS.gmaphs(je,:)=NaN;
       sptrgS.gmatrc(je,:)=NaN;
    end
end; clear je

%picking the gammaM phase bins
out.phsaxis=[-360:18:540];
for je=1:51
    out.gsynALL{je}=sptrg.gSyn(sptrg.gmaphs>=(out.phsaxis(je)-9)&sptrg.gmaphs<(out.phsaxis(je)+9));
    out.gSynM(je)=mean(out.gsynALL{je});
    out.gSynSD(je)=std(out.gsynALL{je});
end; clear je
for je=1:51
    outS.gsynALL{je}=sptrgS.gSyn(sptrgS.gmaphs>=(out.phsaxis(je)-9)&sptrgS.gmaphs<(out.phsaxis(je)+9));
    outS.gSynM(je)=mean(outS.gsynALL{je});
    outS.gSynSD(je)=std(outS.gsynALL{je});
end; clear je

%generating and saving the figure
figname=strcat(bs.name,bs.num,bs.unit,'_gSyn_gammaMcyc_', bs.cell);
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
subplot(3,1,1);
title(strcat(bs.name,bs.num,bs.unit,' gSyn during a gammaMcycle; unit:', bs.cell,'; all cycles with spike'));
line(out.phsaxis,[out.gSynM;out.gSynM-out.gSynSD;out.gSynM+out.gSynSD], 'Color',[0 0 0]);
axis([-360 540 min(out.gSynM-out.gSynSD)-0.05 max(out.gSynM+out.gSynSD)+0.05]);
text(-340,1,strcat('N=',num2str(length(sptrgS.gSyn(:,1)))));
set(gca,'XTick',[-360:90:540]);
xlabel('gammaM CSD phase (°)'); ylabel('mean gsyn') ;

subplot(3,1,2);
title(strcat('gSyn during a gammaMcycles selected for gSyn at spike < 0.35'));
line(out.phsaxis,[outS.gSynM;outS.gSynM-outS.gSynSD;outS.gSynM+outS.gSynSD], 'Color',[0 0 0]);
axis([-360 540 min(outS.gSynM-outS.gSynSD)-0.05 max(outS.gSynM+outS.gSynSD)+0.05]);
text(-340,1,strcat('N=',num2str(length(sptrgS.gSyn(~isnan(sptrgS.gSyn(:,1)),1)))));
set(gca,'XTick',[-360:90:540]);
xlabel('gammaM CSD phase (°)'); ylabel('mean gsyn') ; 

subplot(3,1,3)
title('spike triggered gammaM CSD average')
line(sptrg.time,mean(sptrg.gmatrc));
axis([-0.02 0.03 min(mean(sptrg.gmatrc))*1.1 max(mean(sptrg.gmatrc))*1.1]);
set(gca,'XTick',[-0.02:0.01:0.03]);
xlabel('time (s)'); ylabel('mean CSD') ; 



%saving figure
fg_name=strcat('fig',bs.name, bs.num, bs.unit,'_',bs.cell,'_gSYNgammaMcyc.fig');
saveas(gcf,fg_name)
clear fg_name figname

%saving the collected data
datum = date;
writename=strcat(bs.name, bs.num, bs.unit, '_gSYNgammaM_cell_', bs.cell);
save(writename, 'bs', 'gm', 'datum', 'gs', 'sptrg', 'sptrgS', 'out', 'outS','-v7.3');



