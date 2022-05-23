%This program generates gamma characterisation in single channels from each SilProbe shank  
%Based on filtered CSD and trough-peak detection
%The following analyses are done
%1. Theta phase modulation of gamma amplitude & wavelength
%2. Threshold identification
%2. Gamma phase difference as a function of theta phase 

clear all
%close all

%basic data entry
bs.name='B'; bs.num='192'; bs.unit='d'; bs.typ='sp'; 
bs.flord=[1];

ch.pyrsh=1;                                         %shank for the pyramidal layer theta to be taken from
ch.pyrch=16;                                        %channel number (ventralmost LFP contact is 1) for the pyramidal layer
smp.rate=10000;                                     %give the sampling rate of the channels (CSD and LFP original)

%THETA properties
%number of theta bins
th.nbins=40;th.binsize=360/th.nbins;
%theta filter properties (number of filter coefficients and filter corner frequencies in Hz)
th.ncoeff=512; th.fltlow=5; th.flthigh=12; th.dwnsmpl=25; th.rate=smp.rate/th.dwnsmpl;
%creating theta filter
th_filter=dfilt.dffir(fir1(th.ncoeff, [2.*(th.fltlow./th.rate) 2.*(th.flthigh./th.rate)],'bandpass', gausswin(th.ncoeff+1)));

%GAMMA properties
gm.comp = 'DEGslw';                                                        %name of gamma component
gm.chan = [NaN,2,2,2,NaN,NaN,NaN,NaN];                                 %CSD channel for the actual gamma component for each shank (write NaN for not include)
gm.refr = [NaN,1,1,1,NaN,NaN,NaN,NaN];                                 %channels to be used as a refernce 
gm.cfrq = [29.8322,50.9014];                                              %band pass filter corner frequencies (low and high, included in the boundary)
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
    lfp_name=strcat(bs.name, bs.num, bs.unit, num2str(bs.flord(fl_indx)),'sh',num2str(ch.pyrsh),'LFP.mat');
    load(lfp_name);
    thLFP_file = LFP(ch.pyrch,:); clear LFP lfp_name;
    
    %read the CSD for gamma from each shank; invert the CSD channel so that sink is downwards (-)
    for sh_indx=1:length(gm.chan)
        if ~isnan(gm.chan(sh_indx))
            csd_name=strcat(bs.name, bs.num, bs.unit, num2str(bs.flord(fl_indx)),'sh',num2str(sh_indx),'CSD.mat');
            load(csd_name);
            gmCSD_file{sh_indx} = CSD(gm.chan(sh_indx),:)*-1; clear CSD csd_name;
        else
            gmCSD_file{sh_indx} = NaN;
        end
    end; clear sh_indx
    
    %read theta periods
    th_pname=strcat(bs.name, bs.num, bs.unit, num2str(bs.flord(fl_indx)),'THP.txt');
    if exist(th_pname,'file')==2
        thPER_file=dlmread(th_pname); clear th_pname;
    else
        thPER_file=[]; clear th_pname;
    end
    
    %accumulating CSD and LFP traces, theta segments and file segmentation point 
    if fl_indx==1
        for sh_indx=1:length(gm.chan)
            gmCSD_orig{sh_indx} = gmCSD_file{sh_indx};
        end; clear sh_indx
        thLFP_orig = thLFP_file;
        flSEG = length(thLFP_file);
        thPER = thPER_file;
    else
        for sh_indx=1:length(gm.chan)
            gmCSD_orig{sh_indx} = [gmCSD_orig{sh_indx} gmCSD_file{sh_indx}];
        end; clear sh_indx
        thLFP_orig = [thLFP_orig thLFP_file];
        flSEG(fl_indx)=flSEG(fl_indx-1)+length(thLFP_file)+1;
        thePER = [thPER; thPER_file+((flSEG(fl_indx-1))/smp.rate)];
    end; clear  gmCSD_file thLFP_file thPER_file
end; clear fl_indx flSEG

%---CALCULUS---analysisng gamma CSD
%gamma filtering output: gmCSD_flt gmCSD_sel gmCSD_pts
for sh_indx=1:length(gm.chan)
   if ~isnan(gm.chan(sh_indx))
        gmCSD_flt{sh_indx} = filtfilt(gm_filter.Numerator,1,gmCSD_orig{sh_indx});
        %generate the empty gamma selector variable
        if ~exist ('gmCSD_sel','var')
            gmCSD_sel=NaN(size(gmCSD_flt{sh_indx}));
        end
        [gmCSD_pVAL,gmCSD_pPOS] = findpeaks(gmCSD_flt{sh_indx});
        [gmCSD_tVAL,gmCSD_tPOS] = findpeaks(gmCSD_flt{sh_indx}*-1);
        gmCSD_ptTMP = [[gmCSD_pPOS;gmCSD_pVAL;zeros(size(gmCSD_pPOS))+1] [gmCSD_tPOS;gmCSD_tVAL*-1;zeros(size(gmCSD_tPOS))-1]];
        [~,sIND] = sort(gmCSD_ptTMP(1,:));
        
        %generating a file with peaks and troughs extracted and sorted along the time
        gmCSD_pts{sh_indx} = gmCSD_ptTMP(:,sIND);
        %calculate the amplitude (as SD) of the cycles (this will be the fourth row)
        for je=2:(size(gmCSD_pts{sh_indx},2)-1)
            gmCSD_pts{sh_indx}(4,je) = rms(gmCSD_flt{sh_indx}(gmCSD_pts{sh_indx}(1,je-1):gmCSD_pts{sh_indx}(1,je+1)));
        end
        clear je gmCSD_pVAL gmCSD_tVAL gmCSD_pPOS gmCSD_tPOS sIND gmCSD_ptTMP
    else
        gmCSD_flt{sh_indx}=NaN;
        gmCSD_pts{sh_indx}=NaN;
    end
end; clear sh_indx
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
for sh_indx=1:length(gm.chan)
    if ~isnan(gm.chan(sh_indx))
        %get theta phase if menaingful
        for je=1:size(gmCSD_pts{sh_indx},2)
            
            if gmCSD_sel(gmCSD_pts{sh_indx}(1,je))==1
                act_gmpkT=gmCSD_pts{sh_indx}(1,je)/gm.rate;
                cyc_beg=max(thLFP_troughsT(thLFP_troughsT<=act_gmpkT));
                cyc_end=min(thLFP_troughsT(thLFP_troughsT>act_gmpkT));
                cyc_LFP=thLFP_flt((round(cyc_beg*th.rate)):round((cyc_end*th.rate))); [~,cyc_cent]=max(cyc_LFP);
                cyc_cent=(cyc_cent/th.rate)+cyc_beg;
                if  act_gmpkT<=cyc_cent;    tht_phase=180*((act_gmpkT-cyc_beg)/(cyc_cent-cyc_beg));
                else                        tht_phase=180*((act_gmpkT-cyc_cent)/(cyc_end-cyc_cent))+180;
                end
                gmCSD_pts{sh_indx}(5,je)=tht_phase;
                clear cyc_beg cyc_end cyc_LFP cyc_cent act_gmpkT tht_phase
            else
                gmCSD_pts{sh_indx}(5,je)=NaN;
            end
        end; clear je
        
        %get gamma wavelength (in ms)
        for je=2:size(gmCSD_pts{sh_indx},2)
            gmCSD_pts{sh_indx}(6,je)=((gmCSD_pts{sh_indx}(1,je)-gmCSD_pts{sh_indx}(1,je-1))./(gm.rate/1000))*2;
        end
        toc
    end
end; clear sh_indx
clear thLFP_troughs thLFP_troughsT thLFP_flt 

%---DATA ARRANGEMENT & SORTING-----------
%sorted by theta phase, descriptive data
for sh_indx=1:length(gm.chan)
    if isnan(gmCSD_pts{sh_indx})
        gmCSD_sort{sh_indx}=NaN;
    else
        gmCSD_sort{sh_indx}=gmCSD_pts{sh_indx};
        gmCSD_sort{sh_indx}(:,isnan(gmCSD_sort{sh_indx}(5,:)))=[];
        gmCSD_sort{sh_indx}=(sortrows(gmCSD_sort{sh_indx}',5))';
    end
end; clear sh_indx

%---DATA SUMMARIZING-----------------------------
%get maximal amplitudes and wavelengths
for sh_indx=1:length(gm.chan)
    if isnan(gmCSD_sort{sh_indx})
        gmCSD_ampmx(1,sh_indx)=NaN;
        gmCSD_wlgmx(1,sh_indx)=NaN;
    else
        gmCSD_ampmx(1,sh_indx)=max(gmCSD_sort{sh_indx}(4,:));
        gmCSD_wlgmx(1,sh_indx)=max(gmCSD_sort{sh_indx}(6,:));
    end
end; clear sh_indx
gmCSD_ampbins=0.1:0.2:(((ceil(max(gmCSD_ampmx)*5))/5)-0.1);
gmCSD_wlgbins=0.05:0.2:ceil(max(gmCSD_wlgmx));
clear gmCSD_ampmx gmCSD_wlgmx
toc
'theta filtering and data arrangement done'
%amplitude and wavelength of all cycles as a function of theta phase  
%CALCULUS---descriptive
for sh_indx=1:length(gm.chan)
    if isnan(gmCSD_sort{sh_indx})
        gmCSD_amp{sh_indx}=NaN;
        gmCSD_wlghist{sh_indx}=NaN;
        gmCSD_wlgRS{sh_indx}=NaN;
        gmCSD_wlgFL{sh_indx}=NaN;
    else
        for je=1:th.nbins
            %amplitude
            amps_bin=gmCSD_sort{sh_indx}(4,(gmCSD_sort{sh_indx}(5,:))>=((je-1)*th.binsize)&(gmCSD_sort{sh_indx}(5,:))<(je*th.binsize));
            gmCSD_amp{sh_indx}(1,je)=mean(amps_bin);
            gmCSD_amp{sh_indx}(2,je)=std(amps_bin);
            gmCSD_amp{sh_indx}(3,je)=median(amps_bin);
            gmCSD_amp{sh_indx}(4,je)=prctile(amps_bin,5);
            gmCSD_amp{sh_indx}(5,je)=prctile(amps_bin,95);
            gmCSD_amp{sh_indx}(6,je)=prctile(amps_bin,75);
            gmCSD_amphist{sh_indx}(:,je)=hist(amps_bin,gmCSD_ampbins);
            %wavelength
            wlgs_bin=gmCSD_sort{sh_indx}(6,(gmCSD_sort{sh_indx}(5,:))>=((je-1)*th.binsize)&(gmCSD_sort{sh_indx}(5,:))<(je*th.binsize));
            gmCSD_wlghist{sh_indx}(:,je)=hist(wlgs_bin,gmCSD_wlgbins);
        end;clear je amps_bin wlgs_bin
        gmCSD_wlgRS{sh_indx}=hist((gmCSD_sort{sh_indx}(6,gmCSD_sort{sh_indx}(3,:)==1)),gmCSD_wlgbins);
        gmCSD_wlgFL{sh_indx}=hist((gmCSD_sort{sh_indx}(6,gmCSD_sort{sh_indx}(3,:)==-1)),gmCSD_wlgbins);
    end
end; clear sh_indx

%%--------------------IDENTIFYING THE AMPLITUDE THRESHOLD FOR EACH CHANNEL
for sh_indx=1:length(gm.chan)
    if isnan(gmCSD_sort{sh_indx})
        gmCSD_thrshld{sh_indx}=NaN;
    else
        %find the threshold for the actual shank
        %1-find the theta bin with the smallest median amplitude
        %2-fit gaussian with above peak values excluded except zeros
        %3-threshold is centroid+2SD
        [~,m_indx]=min(gmCSD_amp{sh_indx}(3,:));
        thrs_h=[gmCSD_ampbins' gmCSD_amphist{sh_indx}(:,m_indx)];
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
        gmCSD_thrshld{sh_indx}=fRes.b1+((fRes.c1)*2);
        
        clear m_indx thrs_h PEAKval PEAKpos fOpt fRes
        
        %decide if above threshold for each gamma cycle (this will be line7 of gmCSD_pts)
        for je=1:size(gmCSD_pts{sh_indx},2)
            if isnan(gmCSD_pts{sh_indx}(5,je))
                 gmCSD_pts{sh_indx}(7,je)=NaN;
            elseif gmCSD_pts{sh_indx}(4,je)>=gmCSD_thrshld{sh_indx}
                gmCSD_pts{sh_indx}(7,je)=1;
            else
                gmCSD_pts{sh_indx}(7,je)=NaN;
            end
        end; clear je
    end
end; clear sh_indx
toc
'thresholding and gamma cycles done'
%PLOTTING--descriptive
for fg_indx=1:ceil(length(gm.chan)/4)
    figname=strcat(bs.name, bs.num, bs.unit,gm.comp,'amplitude',num2str(fg_indx));
    figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
    %shank 1 is on the top
    for sh_indx=(((fg_indx-1)*4)+1):(((fg_indx-1)*4)+4)
        if sh_indx<=length(gm.chan)
            if ~isnan(gmCSD_amp{sh_indx})
                subplot(4,4,(((sh_indx-((fg_indx-1)*4))-1)*4)+1)
                surf([-th.binsize/2:th.binsize:360+th.binsize/2], gmCSD_ampbins, [gmCSD_amphist{sh_indx}(:,end) gmCSD_amphist{sh_indx} gmCSD_amphist{sh_indx}(:,1)], 'EdgeColor', 'interp', 'FaceColor', 'interp'); shading flat; axis tight; colormap hot; view(2);
                axis([0 360 min(gmCSD_ampbins) 1.2*(max(gmCSD_amp{sh_indx}(5,:)))]);
                xlabel('theta phase (deg)'); set(gca, 'XTick',[0:180:360]); ylabel ('amplitude (rms)');
                subplot(4,4,(((sh_indx-((fg_indx-1)*4))-1)*4)+2)
                errorbar ((th.binsize/2):th.binsize:(360-(th.binsize/2)),gmCSD_amp{sh_indx}(1,:),gmCSD_amp{sh_indx}(2,:));
                actualtext=strcat(bs.name,bs.num,bs.unit,' ',gm.comp,' shank',num2str(sh_indx));
                text(0, 0.9*(max(gmCSD_amp{sh_indx}(5,:))), actualtext); clear actualtext;
                hold on
                line([0 360], [gmCSD_thrshld{sh_indx} gmCSD_thrshld{sh_indx}],'Color','r')
                hold off
                axis([0 360 min(gmCSD_ampbins) 1.2*(max(gmCSD_amp{sh_indx}(5,:)))]);
                xlabel('theta phase (deg)'); set(gca, 'XTick',[0:180:360]); ylabel ('amplitude (rms)');
                %wavelength
                subplot(4,4,(((sh_indx-((fg_indx-1)*4))-1)*4)+3)
                surf([-th.binsize/2:th.binsize:360+th.binsize/2], gmCSD_wlgbins, [gmCSD_wlghist{sh_indx}(:,end) gmCSD_wlghist{sh_indx} gmCSD_wlghist{sh_indx}(:,1)], 'EdgeColor', 'interp', 'FaceColor', 'interp'); shading flat; axis tight; colormap hot; view(2);
                axis([0 360 ((1/gm.cfrq(2)))*1000 ((1/gm.cfrq(1)))*1000]);
                xlabel('theta phase (deg)'); set(gca, 'XTick',[0:180:360]); ylabel ('wavelength (ms)');
                subplot(4,4,(((sh_indx-((fg_indx-1)*4))-1)*4)+4);
                hold on
                line(gmCSD_wlgRS{sh_indx},gmCSD_wlgbins,'Color','r')
                line(gmCSD_wlgFL{sh_indx},gmCSD_wlgbins,'Color','b')
                axis([0 1.2*(max(gmCSD_wlgRS{sh_indx})) ((1/gm.cfrq(2)))*1000 ((1/gm.cfrq(1)))*1000]);
                xlabel('count'); ylabel ('wavelength (ms)');
                hold off
            end
        end
    end; clear sh_indx
end; clear fg_indx

%PLOTTING THE SUPPLEMENTARY STUFF
for fg_indx=1:ceil(length(gm.chan)/4)
    figname=strcat(bs.name, bs.num, bs.unit,gm.comp,'wlg_amp',num2str(fg_indx));
    figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
    %shank 1 is on the top
    for sh_indx=(((fg_indx-1)*4)+1):(((fg_indx-1)*4)+4)
        if sh_indx<=length(gm.chan)
            if ~isnan(gmCSD_amp{sh_indx})
                %subsampling
                if size(gmCSD_sort{sh_indx},2)>15000
                    picker=randperm(size(gmCSD_sort{sh_indx},2),15000);
                else
                    picker=randperm(size(gmCSD_sort{sh_indx},2));
                end
                subplot(4,1,(sh_indx-((fg_indx-1)*4)))
                gmCSD_resort=gmCSD_sort{sh_indx}(:,picker);
                scatter(gmCSD_resort(4,:),gmCSD_resort(6,:),5,gmCSD_resort(5,:)); colormap HSV;
                xlabel('amplitude (rms)'); ylabel ('wavelength (ms)');
                hold on
                line([gmCSD_thrshld{sh_indx} gmCSD_thrshld{sh_indx}],[min(gmCSD_resort(6,:)) max(gmCSD_resort(6,:))],'Color','r')
                hold off
                clear gmCSD_resort picker
            end
        end
    end; clear sh_indx
end; clear fg_indx

%DISPLAY OF SELECTED GAMMA CYCLES (OPTIONAL)
% for sh_indx=1:length(gm.chan)
%     figure
%     for je=1:(size(gmCSD_pts{sh_indx},2))-1
%         if gmCSD_pts{sh_indx}(7,je)==1|gmCSD_pts{sh_indx}(7,je+1)==1
%             gmCSD_sig(gmCSD_pts{sh_indx}(1,je):gmCSD_pts{sh_indx}(1,je+1))=gmCSD_sel(gmCSD_pts{sh_indx}(1,je):gmCSD_pts{sh_indx}(1,je+1));
%         end
%     end
% end

%---CALCULUS----gamma phase and amplitude accross shanks
%Generate an interaction variable for each gamma cycle

parpool
for sh_indx=1:length(gm.chan)
    if gm.refr(sh_indx)==1
        gmCSD_ref=gmCSD_pts{sh_indx};
        gmCSD_ref(11:((length(gm.chan)-1)*4)+14,1:size(gmCSD_ref,2))=NaN;
        for sh_xshi=1:length(gm.chan)
            rowindx=((sh_xshi-1)*4)+11;
            gmCSD_red=gmCSD_pts{sh_xshi};
            if sh_indx==sh_xshi
                gmCSD_ref(rowindx:rowindx+3,:)=NaN;
            elseif isnan(gm.chan(sh_xshi))
                gmCSD_ref(rowindx:rowindx+3,:)=NaN;
            else
                %find the closest neighbour (omit the first and the last two points (label as NaN))
                %before and after & before_i and after_i are the troughs
                %and peaks of the referred gamma around 
                before(1:size(gmCSD_ref,2))=NaN;
                before_i(1:size(gmCSD_ref,2))=NaN;
                after(1:size(gmCSD_ref,2))=NaN;
                after_i(1:size(gmCSD_ref,2))=NaN;
                
                gmCSD_redcomp=gmCSD_red(1,:);
                gmCSD_refcomp=gmCSD_ref(1,:);
                parfor je=3:(size(gmCSD_ref,2)-2)
                %for je=3:(size(gmCSD_ref,2)-2)
                    before_i(je)=size(gmCSD_redcomp(gmCSD_redcomp<=gmCSD_refcomp(je)),2);
                end
                after_i=before_i+1;
                parfor je=3:(size(gmCSD_ref,2)-2)
                %for je=3:(size(gmCSD_ref,2)-2)
                    after(je)=gmCSD_red(1,after_i(je));
                    before(je)=gmCSD_red(1,before_i(je));
                end
                clear gmCSD_redcomp gmCSD_refcomp je    
                
                %position and amplitude of closest gamma in the other shank
                %take if that cycle is above threshold amplitude
                for je=3:(size(gmCSD_ref,2)-2)
                    if abs(gmCSD_ref(1,je)-before(je))<=abs(gmCSD_ref(1,je)-after(je))
                        gmCSD_ref(rowindx,je)=before(je);
                        gmCSD_ref(rowindx+1,je)=gmCSD_red(4,before_i(je));
                        gmCSD_ref(rowindx+3,je)=gmCSD_red(7,before_i(je));
                        %determine the phase offset
                        %the dependant relative to reference: positive if the dependant lags negative if the dependant leads
                        if gmCSD_ref(3,je)==gmCSD_red(3,before_i(je))
                            phase_rel=((before(je)-gmCSD_ref(1,je))/(gmCSD_ref(1,je)-(gmCSD_ref(1,je-1))))*180;
                        else
                            phase_rel=(((before(je)-gmCSD_ref(1,je))/(gmCSD_ref(1,je)-(gmCSD_ref(1,je-1))))*180)+180;
                        end
                    else
                        gmCSD_ref(rowindx,je)=after(je);
                        gmCSD_ref(rowindx+1,je)=gmCSD_red(4,after_i(je));
                        gmCSD_ref(rowindx+3,je)=gmCSD_red(7,after_i(je));
                        %determine the phase offset
                        %the dependant relative to reference: positive if the dependant lags negative if the dependant leads
                        if gmCSD_ref(3,je)==gmCSD_red(3,after_i(je))
                            phase_rel=((after(je)-gmCSD_ref(1,je))/((gmCSD_ref(1,je+1))-gmCSD_ref(1,je)))*180;
                        else
                            phase_rel=(((after(je)-gmCSD_ref(1,je))/((gmCSD_ref(1,je+1))-gmCSD_ref(1,je)))*180)-180;
                        end
                    end
                    gmCSD_ref(rowindx+2,je)=phase_rel;
                    clear phase_rel 
                end; clear je
                clear gmCSD_red before before_i after after_i
            end
        end; clear sh_xshi rowindx
        strcat('interaction on shank ', num2str(sh_indx), ' done')
        toc
        %replace the variable
        gmCSD_pts{sh_indx}=gmCSD_ref; clear gmCSD_ref
    end
end; clear sh_indx 
toc
'interaction first half done'

%---DATA ARRANGEMENT & SORTING AGAIN -----------
%sorted by theta phase, descriptive data
for sh_indx=1:length(gm.chan)
    if isnan(gmCSD_pts{sh_indx})
        gmCSD_sort{sh_indx}=NaN;
    else
        gmCSD_sort{sh_indx}=gmCSD_pts{sh_indx};
        gmCSD_sort{sh_indx}(:,isnan(gmCSD_sort{sh_indx}(5,:)))=[];
        gmCSD_sort{sh_indx}=(sortrows(gmCSD_sort{sh_indx}',5))';
    end
end; clear sh_indx

%CALCULUS---interaction theta phase dependent
gmCSD_intbins=[-180:(360/gm.nbins):180];
for sh_indx=1:length(gm.chan)
    if gm.refr(sh_indx)==1
        for je=1:th.nbins
            %interaction
            ints_bin=gmCSD_sort{sh_indx}(:,(gmCSD_sort{sh_indx}(5,:))>=((je-1)*th.binsize)&(gmCSD_sort{sh_indx}(5,:))<(je*th.binsize));
            for dp_indx=1:length(gm.chan)
                if sh_indx~=dp_indx
                    if isnan(gmCSD_amp{dp_indx})
                        gmCSDint{sh_indx,dp_indx}=NaN;
                        gmCSD_intphH{sh_indx,dp_indx}=NaN;
                    else
                        [er,pe]=corr(ints_bin(4,:)', ints_bin((((dp_indx-1)*4)+12),:)');
                        gmCSD_int{sh_indx,dp_indx}(1,je)=er;
                        gmCSD_int{sh_indx,dp_indx}(2,je)=pe;
                        clear er pe
                        [er,pe]=corr(ints_bin(4,:)', ints_bin((((dp_indx-1)*4)+12),:)','type', 'Spearman');
                        gmCSD_int{sh_indx,dp_indx}(3,je)=er;
                        gmCSD_int{sh_indx,dp_indx}(4,je)=pe;
                        clear er pe
                        gmCSD_int{sh_indx,dp_indx}(5,je)=rad2deg(circ_mean((deg2rad(ints_bin((((dp_indx-1)*4)+13),:)))'));
                        gmCSD_int{sh_indx,dp_indx}(6,je)=circ_r((deg2rad(ints_bin((((dp_indx-1)*4)+13),:)))');
                        gmCSD_int{sh_indx,dp_indx}(7,je)=circ_rtest((deg2rad(ints_bin((((dp_indx-1)*4)+13),:)))');
                        gmCSD_intphH{sh_indx,dp_indx}(1:(gm.nbins+1),je)=hist(ints_bin((((dp_indx-1)*4)+13),:),gmCSD_intbins);
                        gmCSD_intphH{sh_indx,dp_indx}(1,je)=(gmCSD_intphH{sh_indx,dp_indx}(1,je)+gmCSD_intphH{sh_indx,dp_indx}((gm.nbins+1),je));
                        gmCSD_intphH{sh_indx,dp_indx}((gm.nbins+1),je)=gmCSD_intphH{sh_indx,dp_indx}(1,je);
                    end
                end
            end
            clear dp_indx ints_bin
        end; clear je
    end
end; clear sh_indx
toc
'intraction completely done'
%PLOT---interaction
for sh_indx=1:length(gm.chan)
    if gm.refr(sh_indx)==1
        figname=strcat(bs.name, bs.num, bs.unit, gm.comp,'_shank',num2str(sh_indx),'_interaction_1');
        figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
        for dp_indx=1:length(gm.chan)
            if ~isnan(gmCSD_amp{dp_indx})
                %amplitude correlations
                subplot(length(gm.chan),4,((dp_indx-1)*4)+1)
                if sh_indx==dp_indx
                    text(0,0.2,strcat(bs.name, bs.num, bs.unit, gm.comp))
                    text(0,0,strcat('interactions of shank ', num2str(sh_indx)))
                    axis off
                else
                    hold on
                    line([(th.binsize/2):th.binsize:(360-(th.binsize/2))],gmCSD_int{sh_indx,dp_indx}(1,:),'Color','r');
                    line([(th.binsize/2):th.binsize:(360-(th.binsize/2))],gmCSD_int{sh_indx,dp_indx}(3,:),'Color','b')
                    axis([0 360 0 0.8]); set(gca, 'XTick',[0:180:360]); set(gca, 'YTick',[0:0.2:0.8])
                    xlabel('theta phase (deg)'); ylabel ('correlation (r:Pr b:Rk)');
                    hold off
                end
                %phase difference histogram
                if sh_indx~=dp_indx
                    subplot(length(gm.chan),4,((dp_indx-1)*4)+2)
                    surf([-th.binsize/2:th.binsize:360+th.binsize/2], gmCSD_intbins, [gmCSD_intphH{sh_indx,dp_indx}(:,end) gmCSD_intphH{sh_indx,dp_indx} gmCSD_intphH{sh_indx,dp_indx}(:,1)], 'EdgeColor', 'interp', 'FaceColor', 'interp'); shading flat; axis tight; colormap bone; view(2);
                    axis([0 360 -90 90]);
                    xlabel('theta phase (deg)'); set(gca, 'XTick',[0:180:360]); ylabel ('gamma phase diff (deg)'); set(gca, 'YTick',[-90:30:90]);
                end
                %mean phase difference
                if sh_indx~=dp_indx
                    subplot(length(gm.chan),4,((dp_indx-1)*4)+3)
                    line([(th.binsize/2):th.binsize:(360-(th.binsize/2))],gmCSD_int{sh_indx,dp_indx}(5,:),'Color','r');
                    axis([0 360 -90 90]);
                    xlabel('theta phase (deg)'); set(gca, 'XTick',[0:180:360]); ylabel ('gamma phase diff (deg)'); set(gca, 'YTick',[-90:30:90]);
                end
                %mean vector length of phase difference
                if sh_indx~=dp_indx
                    subplot(length(gm.chan),4,((dp_indx-1)*4)+4)
                    line([(th.binsize/2):th.binsize:(360-(th.binsize/2))],gmCSD_int{sh_indx,dp_indx}(6,:),'Color','b');
                    axis([0 360 0 1]);
                    xlabel('theta phase (deg)'); set(gca, 'XTick',([0:180:360])); ylabel ('vector length (diff)'); set(gca, 'YTick',[0:0.2:1]);
                end
            end
        end; clear dp_indx
    end
end; clear sh_indx

%PLOTTING THE SUPPLEMENTARY STUFF
for sh_indx=1:length(gm.chan)
    if gm.refr(sh_indx)==1
        figname=strcat(bs.name, bs.num, bs.unit, gm.comp,'_shank',num2str(sh_indx),'_interaction_2');
        figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
        for dp_indx=1:length(gm.chan)
            if ~isnan(gmCSD_amp{dp_indx})
                %amplitude correlations
                subplot(length(gm.chan),1,dp_indx)
                if sh_indx==dp_indx
                    text(0,0.2,strcat(bs.name, bs.num, bs.unit, gm.comp))
                    text(0,0,strcat('interactions of shank ', num2str(sh_indx)))
                    axis off
                else
                    if size(gmCSD_sort{sh_indx},2)>15000
                        picker=randperm(size(gmCSD_sort{sh_indx},2),15000);
                    else
                        picker=randperm(size(gmCSD_sort{sh_indx},2));
                    end
                    
                    gmCSD_intsort=gmCSD_sort{sh_indx}(:,picker);
                    scatter((gmCSD_intsort(4,:)+gmCSD_intsort((((dp_indx-1)*4)+12),:)),gmCSD_intsort((((dp_indx-1)*4)+13),:),5,gmCSD_intsort(5,:)); colormap HSV;
                    axis([0 max(gmCSD_intsort(4,:)+gmCSD_intsort((((dp_indx-1)*4)+12),:)) -180 180]);
                    xlabel('amplitude sum (rms)'); ylabel ('gamma phase diff (deg)'); set(gca, 'YTick',[-180:180:180]);
                    clear gmCSD_intsort picker
                end
            end
        end; clear dp_indx
    end
end; clear sh_indx

%PHASE relationship for 'significant (strong)' gamma cycles - CALCULUS & PLOT
for sh_indx=1:length(gm.chan)
    if gm.refr(sh_indx)==1
        figname=strcat(bs.name, bs.num, bs.unit, gm.comp,'_shank',num2str(sh_indx),'_interaction_SEL1');
        figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
        for dp_indx=1:length(gm.chan)
            if sh_indx==dp_indx
                subplot(length(gm.chan),3,((dp_indx-1)*3)+1)
                axis off
                text(0,0.6,strcat(bs.name, bs.num, bs.unit, gm.comp))
                text(0,0.2,strcat('actual reference shank is',num2str(sh_indx)))
            elseif ~isnan(gm.chan(sh_indx))&&~isnan(gm.chan(dp_indx))
                gmCSD_intSELLOOSE{sh_indx,dp_indx}=gmCSD_pts{sh_indx}(:,gmCSD_pts{sh_indx}((((dp_indx-1)*4)+14),:)==1|gmCSD_pts{sh_indx}(7,:)==1);
                gmCSD_intSELSTRNG{sh_indx,dp_indx}=gmCSD_pts{sh_indx}(:,gmCSD_pts{sh_indx}((((dp_indx-1)*4)+14),:)==1&gmCSD_pts{sh_indx}(7,:)==1);
                
                %calculate and plot phase difference histogram
                gmCSDintHISTph{sh_indx,dp_indx}=hist(gmCSD_intSELSTRNG{sh_indx,dp_indx}((((dp_indx-1)*4)+13),:),gmCSD_intbins);
                subplot(length(gm.chan),3,((dp_indx-1)*3)+2)
                bar(gmCSD_intbins,gmCSDintHISTph{sh_indx,dp_indx})
                axis([-180 180 0 max((gmCSDintHISTph{sh_indx,dp_indx})')])
                xlabel('phase difference (deg)'); ylabel ('cycle count'); set(gca, 'XTick',[-180:45:180]);
                
                subplot(length(gm.chan),3,((dp_indx-1)*3)+3)
                hist(gmCSD_intSELSTRNG{sh_indx,dp_indx}(5,:),0:th.binsize:360)
                axis([0, 360, 0, max(hist(gmCSD_intSELSTRNG{sh_indx,dp_indx}(5,:),0:th.binsize:360))*1.2])
                xlabel('theta phase (deg)'); ylabel ('cycle count'); set(gca, 'XTick',[-0:90:360]);
                
                %calculate and plot statistical data
                gmCSD_intSTAT{sh_indx,dp_indx}(1)=rad2deg(circ_mean((deg2rad(gmCSD_intSELSTRNG{sh_indx,dp_indx}((((dp_indx-1)*4)+13),:)))'));
                gmCSD_intSTAT{sh_indx,dp_indx}(2)=circ_r((deg2rad(gmCSD_intSELSTRNG{sh_indx,dp_indx}((((dp_indx-1)*4)+13),:)))');
                gmCSD_intSTAT{sh_indx,dp_indx}(3)=circ_rtest((deg2rad(gmCSD_intSELSTRNG{sh_indx,dp_indx}((((dp_indx-1)*4)+13),:)))');
                gmCSD_intSTAT{sh_indx,dp_indx}(4)=size(gmCSD_intSELSTRNG{sh_indx,dp_indx}((((dp_indx-1)*4)+13),:),2);
                subplot(length(gm.chan),3,((dp_indx-1)*3)+1)
                axis off
                text (0,0.9,strcat('dep. shank is:',num2str(dp_indx)))
                text (0,0.7,strcat('mean phase:',num2str(gmCSD_intSTAT{sh_indx,dp_indx}(1))))
                text (0,0.5,strcat('r:',num2str(gmCSD_intSTAT{sh_indx,dp_indx}(2))))
                text (0,0.3,strcat('P:',num2str(gmCSD_intSTAT{sh_indx,dp_indx}(3))))
                text (0,0.1,strcat('n pairs:',num2str(gmCSD_intSTAT{sh_indx,dp_indx}(4))))
                
                %calculate amplitude correlations
                %[er,pe]=corr(gmCSD_intSELLOOSE{sh_indx,dp_indx}(4,:)',gmCSD_intSELLOOSE{sh_indx,dp_indx}((((dp_indx-1)*4)+12),:)');
                [er,pe]=corr(gmCSD_sort{sh_indx}(4,:)',gmCSD_sort{sh_indx}((((dp_indx-1)*4)+12),:)');
                gmCSD_intSTAT{sh_indx,dp_indx}(5)=er;
                gmCSD_intSTAT{sh_indx,dp_indx}(6)=pe;
                clear er pe
                %[er,pe]=corr(gmCSD_intSELLOOSE{sh_indx,dp_indx}(4,:)',gmCSD_intSELLOOSE{sh_indx,dp_indx}((((dp_indx-1)*4)+12),:)','type', 'Spearman');
                [er,pe]=corr(gmCSD_sort{sh_indx}(4,:)',gmCSD_sort{sh_indx}((((dp_indx-1)*4)+12),:)','type', 'Spearman');
                gmCSD_intSTAT{sh_indx,dp_indx}(7)=er;
                gmCSD_intSTAT{sh_indx,dp_indx}(8)=pe;
                clear er pe
            end
        end; clear figname dp_indx
    end
end; clear sh_indx

%AMPLITUDE correlations for gamma cycles - PLOT
for sh_indx=1:length(gm.chan)
    if gm.refr(sh_indx)==1
        figname=strcat(bs.name, bs.num, bs.unit, gm.comp,'_shank',num2str(sh_indx),'_interaction_SEL2');
        figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
        for dp_indx=1:length(gm.chan)
            if sh_indx==dp_indx
                %stats on gamma cycles of reference channel
                subplot(length(gm.chan),3,((dp_indx-1)*3)+1)
                axis off
                text(0,0.9,strcat(bs.name, bs.num, bs.unit, gm.comp))
                text(0,0.6,strcat('actual reference shank is',num2str(sh_indx)))
                text(0,0.3,strcat('overall gm cycle number:',num2str(size(gmCSD_pts{sh_indx}(~isnan(gmCSD_pts{sh_indx}(5,:))),2))));
                text(0,0.1,strcat('significant gm cycle number:',num2str(nansum(gmCSD_pts{sh_indx}(7,:)))));
            elseif ~isnan(gm.chan(sh_indx))&&~isnan(gm.chan(dp_indx))
                %stats on gamma cycles of dependant channel 
                subplot(length(gm.chan),3,((dp_indx-1)*3)+1)
                axis off
                text(0,0.9,strcat('dependant shank is ', num2str(dp_indx)))
                text(0,0.7,strcat('overall gm cycle: ',num2str(size(gmCSD_pts{dp_indx}(~isnan(gmCSD_pts{dp_indx}(5,:))),2))));
                text(0,0.5,strcat('significant gm cycle: ',num2str(nansum(gmCSD_pts{dp_indx}(7,:)))));
                text(0,0.3,strcat('in at least one shank: ',num2str(size(gmCSD_intSELLOOSE{sh_indx,dp_indx},2))));
                text(0,0.1,strcat('in both shanks: ',num2str(size(gmCSD_intSELSTRNG{sh_indx,dp_indx},2))));
                
                %scatter plot for the amplitudes
                subplot(length(gm.chan),3,((dp_indx-1)*3)+2)
                scatter (gmCSD_sort{sh_indx}(4,:),gmCSD_sort{sh_indx}((((dp_indx-1)*4)+12),:),5,gmCSD_sort{sh_indx}(5,:)); colormap HSV;
                caxis([0 360]); xlabel('amplitude reference (rms)'); ylabel ('amplitude dependant (rms)');

                %correlations
                subplot(length(gm.chan),3,((dp_indx-1)*3)+3)
                axis off
                text(0,0.9,strcat('r (Pearson)', num2str(gmCSD_intSTAT{sh_indx,dp_indx}(5))))
                text(0,0.6,strcat('P (Pearson)', num2str(gmCSD_intSTAT{sh_indx,dp_indx}(6))))
                text(0,0.3,strcat('r (Spearman)', num2str(gmCSD_intSTAT{sh_indx,dp_indx}(7))))
                text(0,0,strcat('P (Spearman)', num2str(gmCSD_intSTAT{sh_indx,dp_indx}(8))))
            end
        end; clear dp_indx figname
    end
end; clear sh_indx

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
output.intbins=gmCSD_intbins;
output.GMtrace=gmCSD_flt;
output.GMthrsh=gmCSD_thrshld;
output.GMcycle=gmCSD_pts;
output.GMtrcSL=gmCSD_sel;
output.GMcycST=gmCSD_amp;
output.GMcycHS=gmCSD_amphist;
if length(gm.chan(~isnan(gm.chan)))>1
    output.GMintST=gmCSD_int;
    output.GMintHS=gmCSD_intphH;
    output.GMselcycHS=gmCSDintHISTph;
    output.GMselcycST=gmCSD_intSTAT;
end

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
            
           
           
 