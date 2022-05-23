%This program analyses the gamma oscillatory activity in CSD traces of multiple
%shanks during and in relation to theta oscillation

clear all; close all

%experiment base data
bs.name='IV'; bs.num='004'; bs.unit='a'; bs.typ='sp';
bs.flord=[1];                                         %file extension numbers belonging to the present cell

ch_def={}; 



ch.n=[15 15 15 15 15 15 15 15];                                   %number of channels in each shank
ch.pyrsh=9;                                           %shank for the pyramidal layer theta to be taken from
ch.pyrch=3;                                          %channel number (ventralmost LFP contact is 1) for the pyramidal layer
ch.def={[],[],[],[],[3],[],[14],[11]};                                 %defective channel (one is standing for the ventralmost, and the linear number of channel is required not the identifier

smp.rate=10000;                                       %give the silicon probe sampling rate (in Hz)
smp.target=2000;                                      %give the target sampling rate you want to yuse for this calculus
smp.down=smp.rate/smp.target;

%number of theta bins
th.nbins=20;th.binsize=360/th.nbins;
%theta filter properties (number of filter coefficients and filter corner frequencies in Hz)
th.ncoeff=512; th.fltlow=5; th.flthigh=12; th.dwnsmpl=5; th.rate=smp.target/th.dwnsmpl;
%creating theta filter
th_filter=dfilt.dffir(fir1(th.ncoeff, [2.*(th.fltlow./th.rate) 2.*(th.flthigh./th.rate)],'bandpass', gausswin(th.ncoeff+1)));
    
%gamma wavelet properties
gm.frsta=200; gm.frfin=15; gm.cfnum=80; gm.wavelet{1}='cmor1-1.5'; gm.plimit=0.05;
gm.frcnt=centfrq(gm.wavelet{1}); gm.scsta=gm.frcnt/(gm.frsta/smp.target); gm.scfin=gm.frcnt/(gm.frfin/smp.target); gm.scint=(gm.scfin-gm.scsta)/gm.cfnum; gm.fraxis=gm.frcnt./((gm.scsta:gm.scint:gm.scfin)./smp.target);
gm.nbins=20;gm.binsize=360/gm.nbins;

%shank stepper
for indx_sh=1:length(ch.n)
    %CSD channel stepper
    
    for indx_ch=1:ch.n(indx_sh) 
        if isempty(find((ch.def{indx_sh}==(indx_ch)), 1)) %exclude defected channels
            if indx_ch==1 && indx_sh==1; thr.percount=1; end
            
            for indx_fl=1:length(bs.flord) %file stepper
                
                th_pname=strcat(bs.name, bs.num, bs.unit, num2str(bs.flord(indx_fl)),'THP.txt');  % theta periods
                
                if exist(th_pname,'file')==2; th_pers=dlmread(th_pname); clear th_pname;
                    %reading LFP for gamma to theta cpoupling
                    if indx_ch==1&&indx_sh==1
                        lfp_name=strcat(bs.name, bs.num, bs.unit, num2str(bs.flord(indx_fl)),'sh',num2str(ch.pyrsh),'LFP.mat');
                        load(lfp_name);
                        LFPpyr=downsample(LFP(ch.pyrch,:),smp.down);clear LFP lfp_name
                        LFPred = downsample(smooth(LFPpyr,th.dwnsmpl),th.dwnsmpl); clear LFPpyr;
                        LFPtht=filtfilt(th_filter.Numerator,1,LFPred'); clear LFPred;
                        [~,LFPtroughs]=findpeaks(LFPtht.*-1,'MINPEAKHEIGHT', 0.01, 'MINPEAKDISTANCE', floor(th.rate/th.flthigh)-3);
                        for ka=1:length(LFPtroughs)-1
                            [~,tempmax]=max(LFPtht(LFPtroughs(ka):LFPtroughs(ka+1)));LFPpeaks(ka)=tempmax+LFPtroughs(ka);clear tempmax;
                        end; clear ka;
                        LFPtroughs=(LFPtroughs.*th.dwnsmpl)-th.dwnsmpl+1;
                        LFPpeaks=(LFPpeaks.*th.dwnsmpl)-th.dwnsmpl+1;
                        clear LFPtht;
                    end
                    
                    %load the CSD
                    csd_name=strcat(bs.name, bs.num, bs.unit, num2str(bs.flord(indx_fl)),'sh',num2str(indx_sh),'CSD.mat');
                    load(csd_name); clear csd_name;
                    
                    %extract the LFPtroughpoints
                    if indx_ch==1&&indx_sh==1
                        LFPtroughpoints=zeros(1,size(CSD,2));
                        LFPtroughpoints(LFPtroughs)=1;
                        LFPtroughpoints(LFPpeaks)=-1;
                        clear LFPtroughs LFPpeaks
                    end
                    
                    %invert the CSD so that the sink is down and the source is up (for 1 chan)
                    if isnan(CSD(indx_ch,1))
                        CSDchan=NaN;
                    else
                        CSDchan=downsample((CSD(indx_ch,:)*-1),smp.down);
                    end
                    clear CSD 
                    
                    %calculate the theta segment wavelet transforms and theta trough segments
                    for ka=1:size(th_pers,1) %theta period stepper
                        sg_start=round(th_pers(ka,1)*smp.target); sg_end=round(th_pers(ka,2)*smp.target);
                        
                        if indx_ch==1&&indx_sh==1
                            sg_LFPtroughpoints=LFPtroughpoints(sg_start:sg_end);
                            sg_LFPtroughpoints(sg_LFPtroughpoints==1)=thr.percount;
                            sg_LFPtroughpoints(sg_LFPtroughpoints==-1)=-1*thr.percount;
                        end
                        
                        %perform segment cwt
                        if isnan(CSDchan(1))
                            sg_cwt=NaN;
                        else
                            sg_CSD=CSDchan((sg_start-0.5*smp.target):(sg_end+0.5*smp.target));
                            %calculating the wavelet transform for the actual channel of the actual segment.
                            sg_cwt=cwt(sg_CSD, gm.scsta:gm.scint:gm.scfin, gm.wavelet{1});
                            sg_cwt=conj(sg_cwt);
                            sg_cwt=sg_cwt(:,(1+0.5*smp.target):(1+0.5*smp.target+sg_end-sg_start));
                        end
                        
                        if ~exist('all_cwt','var')
                            all_cwt=sg_cwt;
                            if indx_ch==1&&indx_sh==1; all_thttroughs=sg_LFPtroughpoints; end
                            if isnan(sg_cwt(1,1)); all_cwt=NaN; end
                        else
                            if ~isnan(all_cwt(1,1))
                                all_cwt=[all_cwt sg_cwt];
                            end
                            if indx_ch==1&&indx_sh==1;
                                all_thttroughs=[all_thttroughs sg_LFPtroughpoints];
                            end
                        end; clear sg_cwt sg_LFPtroughpoints sg_start sg_end sg_CSD
                        if indx_ch==1&&indx_sh==1; thr.percount=thr.percount+1; end;
                    end; clear ka th_pers LFP_troughpoints CSDchan
                end;
            end; clear indx_fl ans
            if indx_ch==1&&indx_sh==1; thr.percount=thr.percount-1; end
            %calculus for all individual channels for all the different measurements
            %theta dependent gamma modulation calculus
            %zscore (as no licence for statistics toolbox is allways available)
            
            %from this point 'all_' refers to a variable summarising all the data for a particular silicon probe channel
            all_cwt_amp=abs(all_cwt);
            all_cwt_Z=(all_cwt_amp-(repmat(mean(all_cwt_amp,2),1,size(all_cwt,2))))./(repmat(std(all_cwt_amp,0,2),1,size(all_cwt,2)));
            cyc_n=1;
            
            if isnan(all_cwt(1,1))
                all_ampmatrix=NaN;
                gmr.thtmodulation{indx_ch}=NaN;
                gmr.thtcycnum{indx_ch}=NaN;
                gmr.oscampspect{indx_ch}=NaN;
                gmr.oscampmedian{indx_ch}=NaN;
                gmr.stdvspect{indx_ch}=NaN;
            else
                for el=1:thr.percount
                    sg_troughindex=find(all_thttroughs==el);
                    sg_peakindex=find(all_thttroughs==-1*el);
                    for em=1:length(sg_troughindex)-1; %cycle stepper
                        cyc_start=sg_troughindex(em); cyc_end=sg_troughindex(em+1);
                        cyc_cent=sg_peakindex(sg_peakindex>cyc_start&sg_peakindex<cyc_end);
                        cyc_step1=(cyc_cent-cyc_start)/(th.nbins/2);
                        cyc_step2=(cyc_end-cyc_cent)/(th.nbins/2);
                        cyc_ampmatrix=zeros(size(all_cwt,1),th.nbins);
                        for en=1:(th.nbins/2); %theta bin stepper
                            cyc_ampmatrix(:,en)=mean(all_cwt_Z(:,round(cyc_start+(en-1)*cyc_step1):round(cyc_start+(en*cyc_step1))),2);
                        end;
                        for en=((th.nbins/2)+1):th.nbins; %theta bin stepper
                            cyc_ampmatrix(:,en)=mean(all_cwt_Z(:,round(cyc_cent+((en-((th.nbins/2)+1))*cyc_step2)):round(cyc_cent+((en-(th.nbins/2))*cyc_step2))),2);
                        end;
                        
                        if cyc_n~=1;
                            all_ampmatrix=all_ampmatrix+cyc_ampmatrix;
                        else all_ampmatrix=cyc_ampmatrix;
                        end
                        cyc_n=cyc_n+1;
                    end; clear cyc_start cyc_end cyc_step cyc_cent cyc_step2 cyc_step1 cyc_ampmatrix sg_troughindex sg_peakindex
                end; clear el em en; cyc_n=cyc_n-1;
                all_ampmatrix=all_ampmatrix./cyc_n;
                
                %if channel is functional
                gmr.thtmodulation{indx_ch}=all_ampmatrix;
                gmr.thtcycnum{indx_ch}=cyc_n;
                gmr.oscampspect{indx_ch}=mean(abs(all_cwt),2);
                gmr.oscampmedian{indx_ch}=median(abs(all_cwt),2);
                gmr.stdvspect{indx_ch}=std(abs(all_cwt)');
                strcat('CSD_channel_', num2str(indx_ch), 'shank', num2str(indx_sh), '_done')
            end
            clear all_ampmatrix all_cwt all_cwt_Z all_cwt_amp cyc_n
        else
            gmr.oscampspect{indx_ch}=NaN;
            gmr.thtmodulation{indx_ch}=NaN;
            gmr.thtcycnum{indx_ch}=NaN;
            strcat('CSD_channel_', num2str(indx_ch), '_is_defected')
        end
    end; clear indx_ch
    
    %plotting results
    %theta-gamma coordination at the reference electrode
    figname=strcat(bs.name, bs.num, bs.unit,'sh',num2str(indx_sh),' gamma modulation by theta');
    figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
    subplot(4,4,1);
    axis([0 1 0 4]); axis off;text(0,3.5,strcat(bs.name,bs.num,bs.unit));text(0,2.5,strcat('shank ',num2str(indx_sh)));text(0,1.5,'gamma CSD');text(0,0.5,'by theta')
    for indx_ch=1:ch.n(indx_sh)
        if isnan(gmr.thtmodulation{indx_ch})
            cmax_t(indx_ch)=0;
        else
            cmax_t(indx_ch)=max(max(abs(gmr.thtmodulation{indx_ch})));
        end
    end; cmax=max(cmax_t); clear cmax_t;
    
    for indx_ch=2:ch.n-1
        subplot(4,4,17-indx_ch)
        if ~isnan(gmr.thtmodulation{indx_ch})
            surf([-th.binsize/2:th.binsize:360+th.binsize/2], gm.fraxis, [gmr.thtmodulation{indx_ch}(:,end) gmr.thtmodulation{indx_ch} gmr.thtmodulation{indx_ch}(:,1)], 'EdgeColor', 'interp', 'FaceColor', 'interp'); shading flat; axis tight; colormap jet; caxis([-1*cmax cmax]); view(2);
            axis([0 360 min(gm.fraxis) max(gm.fraxis)]);
            xlabel('theta phase (deg)'); ylabel ('frequency (Hz)'); set(gca, 'YTick',[(min(gm.fraxis)):10:(max(gm.fraxis))]); set(gca, 'XTick',[0:180:360]);
        end
    end; clear indx_ch cmax figname
    
    %gamma amplitude spectra
    figname=strcat(bs.name, bs.num, bs.unit,'sh',num2str(indx_sh),' gamma amplitude during theta');
    figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
    subplot(4,4,1);
    axis([0 1 0 4]); axis off;text(0,3.5,strcat(bs.name,bs.num,bs.unit));text(0,2.5,strcat('shank ',num2str(indx_sh)));text(0,1.5,'gamma CSD ampl');text(0,0.5,'during theta')
    for indx_ch=2:ch.n(indx_sh)-2; smax_t(indx_ch)=max(gmr.oscampspect{indx_ch}); end; smax=max(smax_t);clear smax_t;
    
    for indx_ch=2:ch.n(indx_sh)-1
        subplot(4,4,17-indx_ch)
        if ~isnan(gmr.oscampspect{indx_ch})
            xlabel('CSD amplitude (mv/mm2)'); ylabel ('frequency (Hz)');
            set(gca, 'YTick',[(min(gm.fraxis)):10:(max(gm.fraxis))]);
            plot(gmr.oscampspect{indx_ch}, gm.fraxis, gmr.oscampmedian{indx_ch}, gm.fraxis, 'r', gmr.stdvspect{indx_ch}, gm.fraxis, '--b')
            axis([0 smax*1.1 (min(gm.fraxis)) (max(gm.fraxis))]);
        end
    end; clear indx_ch smax figname
    
    %gamma amplitude modulation amplitude spectra
    figname=strcat(bs.name, bs.num, bs.unit,'sh',num2str(indx_sh),' gamma amplitude modulation');
    figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
    subplot(4,4,1);
    axis([0 1 0 4]); axis off;text(0,3.5,strcat(bs.name,bs.num,bs.unit));text(0,2.5,strcat('shank ',num2str(indx_sh)));text(0,1.5,'CSD Z difference');text(0,0.5,'during theta')
    
    for indx_ch=2:ch.n(indx_sh)-1
        if ~isnan(gmr.thtmodulation{indx_ch})
            gmr.modZampspect{indx_ch}=max(gmr.thtmodulation{indx_ch}')-min(gmr.thtmodulation{indx_ch}');
            smax_t(indx_ch)=max(gmr.modZampspect{indx_ch});
        else
            gmr.modZampspect{indx_ch}=NaN;
            smax_t(indx_ch)=NaN;
        end
    end; clear indx_ch
    smax=max(smax_t);clear smax_t;
    
    for indx_ch=2:ch.n(indx_sh)-1
        subplot(4,4,17-indx_ch)
        if ~isnan(gmr.modZampspect{indx_ch})
            xlabel('Z difference (max-min)'); ylabel ('frequency (Hz)');
            set(gca, 'YTick',[(min(gm.fraxis)):10:(max(gm.fraxis))]);
            axis([0 smax*1.1 (min(gm.fraxis)) (max(gm.fraxis))]);
            line(gmr.modZampspect{indx_ch}, gm.fraxis)
        end
    end; clear indx_ch smax figname
    
    datum=date;
    %saving the relevant data
    writename=strcat(bs.name, bs.num, bs.unit,'_sh',num2str(indx_sh),'_a_gma_thtCHAR');
    save(writename,'datum','gmr','gm','ch','bs','smp','th','thr');
    clear writename gmr
end; clear sh_indx
 