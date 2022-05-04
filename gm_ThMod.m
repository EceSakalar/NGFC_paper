function [gmtht] = gm_ThMod( input_args )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
for je=1:ch.n %CSD channel stepper
    if isempty(ch.def(ch.def==(je))) %exclude defected channels
        if je==1 
            th_percount=1; 
        end
        
        th_pname=strcat(bs.name, bs.num, bs.exp, 'THP.txt');  % theta periods
        
        if exist(th_pname,'file')==2
            th_pers=dlmread(th_pname); 
            clear th_pname;
            %reading LFP for gamma to theta cpoupling
            if je==1
                lfp_name=strcat(bs.name, bs.num, bs.exp, 'LFP.mat');
                load(lfp_name);
                LFPpyr=LFP(ch_pyr,:);
                clear LFP lfp_name
                LFPred = downsample(smooth(LFPpyr,th_dwnsmpl),th_dwnsmpl);
                clear LFPpyr;
                LFPtht=filtfilt(th_filter.Numerator,1,LFPred');
                clear LFPred;
                [~,LFPtroughs]=findpeaks(LFPtht.*-1,'MINPEAKHEIGHT', 0.01, 'MINPEAKDISTANCE', floor(th_rate/th_flthigh)-3);
                for ka=1:length(LFPtroughs)-1
                    [~,tempmax]=max(LFPtht(LFPtroughs(ka):LFPtroughs(ka+1)));
                    LFPpeaks(ka)=tempmax+LFPtroughs(ka);
                    clear tempmax;
                end; clear ka;
                LFPtroughs=(LFPtroughs.*th_dwnsmpl)-th_dwnsmpl+1;
                LFPpeaks=(LFPpeaks.*th_dwnsmpl)-th_dwnsmpl+1;
                clear LFPtht;
            end
            
            %invert the CSD so that the sink is down and the source is up
            csd_name=strcat(bs_name, bs_num, bs_unit, 'CSD.mat');
            load(csd_name);
            if isnan(CSD(je,1))
                CSDchan=NaN;
            else
                CSDchan=CSD(je,:)*-1;
            end
            clear CSD csd_name;
            if je==1
                LFPtroughpoints=zeros(size(CSDchan));
                LFPtroughpoints(LFPtroughs)=1;
                LFPtroughpoints(LFPpeaks)=-1;
                clear LFPtroughs LFPpeaks
            end
            
            %calculate the theta segment wavelet transforms and theta trough segments
            for ka=1:size(th_pers,1) %theta period stepper
                sg_start=round(th_pers(ka,1)*ch_rate); 
                sg_end=round(th_pers(ka,2)*ch_rate);
                
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
                    if je==1 
                        all_thttroughs=sg_LFPtroughpoints; 
                    end
                    if isnan(sg_cwt(1,1)) 
                        all_cwt=NaN; 
                    end
                else
                    if ~isnan(all_cwt(1,1))
                        all_cwt=[all_cwt sg_cwt];
                    end
                    if je==1
                        all_thttroughs=[all_thttroughs sg_LFPtroughpoints];
                    end
                end 
                clear sg_cwt sg_LFPtroughpoints sg_start sg_end sg_CSD
                if je==1 
                    th_percount=th_percount+1; 
                end
            end 
            clear ka th_pers LFP_troughpoints CSDchan
        end
        
        if je==1 
            th_percount=th_percount-1; 
        end
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
                for em=1:length(sg_troughindex)-1 %cycle stepper
                    cyc_start=sg_troughindex(em); 
                    cyc_end=sg_troughindex(em+1);
                    cyc_cent=sg_peakindex(sg_peakindex>cyc_start&sg_peakindex<cyc_end);
                    cyc_step1=(cyc_cent-cyc_start)/(th_nbins/2);
                    cyc_step2=(cyc_end-cyc_cent)/(th_nbins/2);
                    cyc_ampmatrix=zeros(size(all_cwt,1),th_nbins);
                    for en=1:(th_nbins/2) %theta bin stepper
                        cyc_ampmatrix(:,en)=mean(all_cwt_Z(:,round(cyc_start+(en-1)*cyc_step1):round(cyc_start+(en*cyc_step1))),2);
                    end
                    for en=((th_nbins/2)+1):th_nbins %theta bin stepper
                        cyc_ampmatrix(:,en)=mean(all_cwt_Z(:,round(cyc_cent+((en-((th_nbins/2)+1))*cyc_step2)):round(cyc_cent+((en-(th_nbins/2))*cyc_step2))),2);
                    end
                    
                    if cyc_n~=1
                        all_ampmatrix=all_ampmatrix+cyc_ampmatrix;
                    else
                        all_ampmatrix=cyc_ampmatrix;
                    end
                    cyc_n=cyc_n+1;
                end
                clear cyc_start cyc_end cyc_step cyc_ampmatrix sg_troughindex
            end
            clear el em en; 
            cyc_n=cyc_n-1;
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
writename=strcat(bs_name, bs_num, bs_unit, '_a_gma_thtCHAR');
save(writename, 'datum', 'gmtht');
end

