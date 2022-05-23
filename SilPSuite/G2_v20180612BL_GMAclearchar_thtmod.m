%This program analyses the gamma oscillatory coupling relative to CSD traces
%multichannel LFP and CSD traces 
clear all; close all

bs_name='B'; bs_num='182'; bs_unit='a'; bs_typ='sx'; %unit basic data entry

ch_n=16;                                              %number of channels
ch_rate=2000;                                         %give the silicon probe sampling rate (in Hz)
ch_pyr=14;                                            %channel number (ventralmost LFP contact is 1) for the pyramidal layer
ch_def=[];                                            %defective channel (one is standing for the ventralmost, and the linear number of channel is required not the identifier
fl_ord=[1];                                           %file extension numbers belonging to the present cell

%number of theta bins
th_nbins=20;th_binsize=360/th_nbins;
%theta filter properties (number of filter coefficients and filter corner frequencies in Hz)
th_ncoeff=512; th_fltlow=5; th_flthigh=12; th_dwnsmpl=5; th_rate=ch_rate/th_dwnsmpl;
%creating theta filter
th_filter=dfilt.dffir(fir1(th_ncoeff, [2.*(th_fltlow./th_rate) 2.*(th_flthigh./th_rate)],'bandpass', gausswin(th_ncoeff+1)));
    
%gamma wavelet properties
gm_frsta=200; gm_frfin=15; gm_cfnum=80; gm_wavelet{1}='cmor1-1.5'; gm_plimit=0.05;
gm_frcnt=centfrq(gm_wavelet{1}); gm_scsta=gm_frcnt/(gm_frsta/ch_rate); gm_scfin=gm_frcnt/(gm_frfin/ch_rate); gm_scint=(gm_scfin-gm_scsta)/gm_cfnum; gm_fraxis=gm_frcnt./((gm_scsta:gm_scint:gm_scfin)./ch_rate);
gm_nbins=20;gm_binsize=360/gm_nbins;

for je=1:ch_n; %CSD channel stepper
    if isempty(ch_def(ch_def==(je))); %exclude defected channels
        if je==1; th_percount=1; end
        
        for fl_indx=1:length(fl_ord); %file stepper
            
            th_pname=strcat(bs_name, bs_num, bs_unit, num2str(fl_ord(fl_indx)),'THP.txt');  % theta periods
            
            if exist(th_pname,'file')==2; th_pers=dlmread(th_pname); clear th_pname;
                %reading LFP for gamma to theta cpoupling
                if je==1;
                    lfp_name=strcat(bs_name, bs_num, bs_unit, num2str(fl_ord(fl_indx)),'LFP.mat');
                    load(lfp_name);
                    LFPpyr=LFP(ch_pyr,:);clear LFP lfp_name
                    LFPred = downsample(smooth(LFPpyr,th_dwnsmpl),th_dwnsmpl); clear LFPpyr;
                    LFPtht=filtfilt(th_filter.Numerator,1,LFPred'); clear LFPred;
                    [~,LFPtroughs]=findpeaks(LFPtht.*-1,'MINPEAKHEIGHT', 0.01, 'MINPEAKDISTANCE', floor(th_rate/th_flthigh)-3);
                    for ka=1:length(LFPtroughs)-1
                        [~,tempmax]=max(LFPtht(LFPtroughs(ka):LFPtroughs(ka+1)));LFPpeaks(ka)=tempmax+LFPtroughs(ka);clear tempmax;
                    end; clear ka;
                    LFPtroughs=(LFPtroughs.*th_dwnsmpl)-th_dwnsmpl+1;
                    LFPpeaks=(LFPpeaks.*th_dwnsmpl)-th_dwnsmpl+1;
                    clear LFPtht;
                end
                
                %invert the CSD so that the sink is down and the source is up
                csd_name=strcat(bs_name, bs_num, bs_unit, num2str(fl_ord(fl_indx)),'CSDc1.mat');
                load(csd_name);
                if isnan(CSDc(je,1))
                    CSDchan=NaN;
                else
                    CSDchan=CSDc(je,:)*-1;
                end
                clear CSDc csd_name;
                if je==1;
                    LFPtroughpoints=zeros(size(CSDchan));
                    LFPtroughpoints(LFPtroughs)=1;
                    LFPtroughpoints(LFPpeaks)=-1;
                    clear LFPtroughs LFPpeaks
                end
                
                %calculate the theta segment wavelet transforms and theta trough segments
                for ka=1:size(th_pers,1); %theta period stepper
                    sg_start=round(th_pers(ka,1)*ch_rate); sg_end=round(th_pers(ka,2)*ch_rate);
                    
                    if isnan(CSDchan(1))
                        sg_cwt=NaN;
                    else
                        sg_CSD=CSDchan((sg_start-0.5*ch_rate):(sg_end+0.5*ch_rate));
                        %calculating the wavelet transform for the actual channel of the actual segment.
                        sg_cwt=cwt(sg_CSD, gm_scsta:gm_scint:gm_scfin, gm_wavelet{1});
                        sg_cwt=conj(sg_cwt);
                        sg_cwt=sg_cwt(:,(1+0.5*ch_rate):(1+0.5*ch_rate+sg_end-sg_start));
                    end
                    if je==1
                        sg_LFPtroughpoints=LFPtroughpoints(sg_start:sg_end);
                        sg_LFPtroughpoints(sg_LFPtroughpoints==1)=th_percount;
                        sg_LFPtroughpoints(sg_LFPtroughpoints==-1)=-1*th_percount;
                    end
                                        
                    if ~exist('all_cwt','var')
                        all_cwt=sg_cwt;
                        if je==1; all_thttroughs=sg_LFPtroughpoints; end
                        if isnan(sg_cwt(1,1)); all_cwt=NaN; end
                    else
                        if ~isnan(all_cwt(1,1))
                            all_cwt=[all_cwt sg_cwt];
                        end
                        if je==1;
                            all_thttroughs=[all_thttroughs sg_LFPtroughpoints];
                        end
                    end; clear sg_cwt sg_LFPtroughpoints sg_start sg_end sg_CSD
                    if je==1; th_percount=th_percount+1; end;
                end; clear ka th_pers LFP_troughpoints CSDchan
            end;
        end
        if je==1; th_percount=th_percount-1; end
        %calculus for all individual channels for all the different measurements
        %theta dependent gamma modulation calculus
        %zscore (as no licence for statistics toolbox is allways available)
        all_cwt_amp=abs(all_cwt);
        all_cwt_Z=(all_cwt_amp-(repmat(mean(all_cwt_amp,2),1,size(all_cwt,2))))./(repmat(std(all_cwt_amp,0,2),1,size(all_cwt,2)));
        cyc_n=1;
        if isnan(all_cwt(1,1))
            all_ampmatrix=NaN;
            gm_thtmodulation{je}=NaN;
            gm_thtcycnum{je}=NaN;
            gm_oscampspect{je}=NaN;
            gm_oscampmedian{je}=NaN;
            gm_stdvspect{je}=NaN;
        else
            for el=1:th_percount
                sg_troughindex=find(all_thttroughs==el);
                sg_peakindex=find(all_thttroughs==-1*el);
                for em=1:length(sg_troughindex)-1; %cycle stepper
                    cyc_start=sg_troughindex(em); cyc_end=sg_troughindex(em+1);
                    cyc_cent=sg_peakindex(sg_peakindex>cyc_start&sg_peakindex<cyc_end);
                    cyc_step1=(cyc_cent-cyc_start)/(th_nbins/2);
                    cyc_step2=(cyc_end-cyc_cent)/(th_nbins/2);
                    cyc_ampmatrix=zeros(size(all_cwt,1),th_nbins);
                    for en=1:(th_nbins/2); %theta bin stepper
                        cyc_ampmatrix(:,en)=mean(all_cwt_Z(:,round(cyc_start+(en-1)*cyc_step1):round(cyc_start+(en*cyc_step1))),2);
                    end;
                    for en=((th_nbins/2)+1):th_nbins; %theta bin stepper
                        cyc_ampmatrix(:,en)=mean(all_cwt_Z(:,round(cyc_cent+((en-((th_nbins/2)+1))*cyc_step2)):round(cyc_cent+((en-(th_nbins/2))*cyc_step2))),2);
                    end;
                    
                    if cyc_n~=1;
                        all_ampmatrix=all_ampmatrix+cyc_ampmatrix;
                    else all_ampmatrix=cyc_ampmatrix;
                    end
                    cyc_n=cyc_n+1;
                end; clear cyc_start cyc_end cyc_step cyc_ampmatrix sg_troughindex
            end; clear el em en; cyc_n=cyc_n-1;
            all_ampmatrix=all_ampmatrix./cyc_n;
            
            %if channel is functional
            gm_thtmodulation{je}=all_ampmatrix;
            gm_thtcycnum{je}=cyc_n;
            gm_oscampspect{je}=mean(abs(all_cwt),2);
            gm_oscampmedian{je}=median(abs(all_cwt),2);
            gm_stdvspect{je}=std(abs(all_cwt)');
            strcat('CSD_channel_', num2str(je), '_done')
        end
        clear all_ampmatrix all_cwt all_cwt_Z all_cwt_amp cyc_n
    else
        gm_oscampspect{je}=NaN;
        gm_thtmodulation{je}=NaN;
        gm_thtcycnum{je}=NaN;
        strcat('CSD_channel_', num2str(je), '_is_defected')
    end
end

%plotting results
%theta-gamma coordination at the reference electrode
figname=strcat(bs_name, bs_num, bs_unit,' gamma modulation by theta');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
subplot(4,4,1);
axis([0 1 0 4]); axis off;text(0,3,strcat(bs_name,bs_num,bs_unit));text(0,2,'gamma CSD');text(0,1,'by theta')
for je=1:ch_n
    if isnan(gm_thtmodulation{je})
        cmax_t(je)=0;
    else
        cmax_t(je)=max(max(abs(gm_thtmodulation{je})));
    end
end; cmax=max(cmax_t); clear cmax_t;

 
for je=2:ch_n-1
    subplot(4,4,17-je)
    if ~isnan(gm_thtmodulation{je})
        surf([-th_binsize/2:th_binsize:360+th_binsize/2], gm_fraxis, [gm_thtmodulation{je}(:,end) gm_thtmodulation{je} gm_thtmodulation{je}(:,1)], 'EdgeColor', 'interp', 'FaceColor', 'interp'); shading flat; axis tight; colormap jet; caxis([-1*cmax cmax]); view(2);
        axis([0 360 min(gm_fraxis) max(gm_fraxis)]);
        xlabel('theta phase (deg)'); ylabel ('frequency (Hz)'); set(gca, 'YTick',[(min(gm_fraxis)):10:(max(gm_fraxis))]); set(gca, 'XTick',[0:180:360]);
    end
end; clear je cmax

%gamma amplitude spectra
figname=strcat(bs_name, bs_num, bs_unit,' gamma amplitude during theta');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
subplot(4,4,1);
axis([0 1 0 4]); axis off;text(0,3,strcat(bs_name,bs_num,bs_unit));text(0,2,'gamma CSD ampl');text(0,1,'during theta')
for je=2:ch_n-2; smax_t(je)=max(gm_oscampspect{je}); end; smax=max(smax_t);clear smax_t;
 
for je=2:ch_n-1
    subplot(4,4,17-je)
    if ~isnan(gm_oscampspect{je})
        xlabel('CSD amplitude (mv/mm2)'); ylabel ('frequency (Hz)');
        set(gca, 'YTick',[(min(gm_fraxis)):10:(max(gm_fraxis))]);
        plot(gm_oscampspect{je}, gm_fraxis, gm_oscampmedian{je}, gm_fraxis, 'r', gm_stdvspect{je}, gm_fraxis, '--b')
        axis([0 smax*1.1 (min(gm_fraxis)) (max(gm_fraxis))]);
    end
end; clear je cmax

%gamma amplitude modulation amplitude spectra
figname=strcat(bs_name, bs_num, bs_unit,' gamma amplitude modulation');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
subplot(4,4,1);
axis([0 1 0 4]); axis off;text(0,3,strcat(bs_name,bs_num,bs_unit));text(0,2,'CSD Z difference');text(0,1,'during theta')

for je=2:ch_n-1
    if ~isnan(gm_thtmodulation{je})
        gm_modZampspect{je}=max(gm_thtmodulation{je}')-min(gm_thtmodulation{je}');
        smax_t(je)=max(gm_modZampspect{je});
    else
        gm_modZampspect{je}=NaN;
        smax_t(je)=NaN;
    end
end; clear je cmax
smax=max(smax_t);clear je;

for je=2:ch_n-1
    subplot(4,4,17-je)
    if ~isnan(gm_modZampspect{je})
        xlabel('Z difference (max-min)'); ylabel ('frequency (Hz)');
        set(gca, 'YTick',[(min(gm_fraxis)):10:(max(gm_fraxis))]);
        axis([0 smax*1.1 (min(gm_fraxis)) (max(gm_fraxis))]);
        line(gm_modZampspect{je}, gm_fraxis)
    end
end; clear je cmax

%cleanup
clear fl_indx fl_ord ch_def ch_spks fl_name ans figname       

gmtht.gm_fraxis=gm_fraxis;
gmtht.gm_osc_amplspectMEAN=gm_oscampspect;
gmtht.gm_osc_ampldtdvs=gm_stdvspect;
gmtht.gm_osc_thtmod=gm_thtmodulation;
gmtht.gm_osc_thtmodamp=gm_modZampspect;
gmtht.gm_osc_thtcycnum=gm_thtcycnum;
gmtht.gm_osc_amplspectMEDIAN=gm_oscampmedian;
gmtht.ch_pyr=ch_pyr;
gmtht.ch_samplrate=ch_rate;
datum=date;
%saving the relevant data
writename=strcat(bs_name, bs_num, bs_unit, '_a_gmaclean1_thtCHAR');
save(writename, 'datum', 'gmtht');
 