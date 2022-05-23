%This program analyses the SWR coupling and ripple oscillatory coupling relative to LFP and CSD traces
clear all; %close all

bs_name='B'; bs_num='192'; bs_unit='d'; bs_typ='sx'; %unit basic data entry

%basic parameters
ch_n=16;                                              %number of channels
ch_rate=2000;                                         %give the silicon probe sampling rate (in Hz)
ch_def=[ ];                                           %defective channel (one is standing for the ventralmost, and the linear number of channel is required not the identifier
fl_ord=[1];                                           %file extension numbers belonging to the present cell
ch_step=1;                                            %indicate the step size of the CSD used for the calculus in A1LFPtoCSDtransf  

%bins for time normalisation:before, 1st half, 2nd half, after
n_tbins=[15 20 15];             %wavelet ripple characterisation

%ripple wavelet properties
rp_frsta=250; rp_frfin=100; rp_cfnum=54 ; rp_wavelet{1}='cmor1-1.5'; rp_plimit=0.05;
rp_frcnt=centfrq(rp_wavelet{1});
rp_scsta=rp_frcnt/(rp_frsta/ch_rate);
rp_scfin=rp_frcnt/(rp_frfin/ch_rate);
rp_scint=(rp_scfin-rp_scsta)/rp_cfnum;
rp_fraxis=rp_frcnt./((rp_scsta:rp_scint:rp_scfin)./ch_rate);

%SWR counter reset
SWR_no=0;

%%LFP%%
%preallocate variables
for je=1:ch_n
    LFP_rips{je}=zeros(length(rp_fraxis),sum(n_tbins));
end; LFP_ripsN=LFP_rips; clear je
LFP_segavg=zeros(ch_n,sum(n_tbins));
   
for je=1:ch_n  %LFP channel stepper
    if isempty(ch_def(ch_def==(je))); %exclude defected channels
        
        for fl_indx=1:length(fl_ord); %file stepper
            swr_pname=strcat(bs_name, bs_num, bs_unit, num2str(fl_ord(fl_indx)),'SWP.txt');     %SWR periods
            
            if exist(swr_pname,'file')==2;
                swr_pers=dlmread(swr_pname); clear swr_pname;
                %reading LFP for ripple LFP characterisation
                lfp_name=strcat(bs_name, bs_num, bs_unit, num2str(fl_ord(fl_indx)),'LFP.mat');
                load(lfp_name); LFPact=LFP(je,:); clear LFP lfp_name LFP;
                
                for ka=1:size(swr_pers,1)
                    if je==1; SWR_no=SWR_no+1; end;
                    sta_rip=round(ch_rate*swr_pers(ka,1)); end_rip=round(ch_rate*swr_pers(ka,2));
                    bin_siz=(end_rip-sta_rip)/n_tbins(2);
                    sta_wlt=floor(sta_rip-(bin_siz*(n_tbins(1)*1.4))); end_wlt=ceil(end_rip+(bin_siz*(n_tbins(3)*1.4)));
                    lfp_seg=LFPact(sta_wlt:end_wlt);
                    wlt_rip=conj(cwt(lfp_seg,rp_scsta:rp_scint:rp_scfin, rp_wavelet{1}));
                    sta_rip=sta_rip-sta_wlt+1; end_rip=end_rip-sta_wlt+1;
                    bin_1st=sta_rip-(bin_siz*n_tbins(1));
                    amp_rip=zeros(length(rp_fraxis),sum(n_tbins));
                    for el=0:sum(n_tbins)-1
                        wlt_el=abs(wlt_rip(:,(ceil(bin_1st+(el*bin_siz))):(floor(bin_1st+((el+1)*bin_siz)))));
                        lfp_el=lfp_seg((ceil(bin_1st+(el*bin_siz))):(floor(bin_1st+((el+1)*bin_siz))));
                        amp_rip(:,el+1)=mean(wlt_el,2); clear wlt_el
                        amp_lfp(el+1)=mean(lfp_el); clear lfp_el
                    end; clear el
                    LFP_rips{je}=LFP_rips{je}+amp_rip; clear amp_rip
                    LFP_segavg(je,:)=LFP_segavg(je,:)+amp_lfp; clear amp_lfp
                end; clear sta_rip end_rip bin_siz sta_wlt end_wlt wlt_rip bin_1st ka lfp_seg LFPact
            end; clear swr_pers
        end; clear fl_indx
    else LFP_rips{je}=NaN;LFP_ripsN{je}=NaN;
    end
end; clear je

%%CSD%%
%preallocate variable
for je=1:ch_n
    CSD_rips{je}=zeros(length(rp_fraxis),sum(n_tbins));
end; CSD_ripsN=CSD_rips; clear je
CSD_segavg=zeros(ch_n,sum(n_tbins));
   
for je=(ch_step+1):(ch_n-ch_step)  %CSD channel stepper
    if isempty(ch_def(ch_def==(je))); %exclude defected channels
        
        for fl_indx=1:length(fl_ord); %file stepper
            swr_pname=strcat(bs_name, bs_num, bs_unit, num2str(fl_ord(fl_indx)),'SWP.txt');     %SWR periods
            
            if exist(swr_pname,'file')==2;
                swr_pers=dlmread(swr_pname); clear swr_pname;
                %reading CSD for ripple CSD characterisation
                csd_name=strcat(bs_name, bs_num, bs_unit, num2str(fl_ord(fl_indx)),'CSD.mat');
                %CSD traces are inverted to have the source + and the sink -
                load(csd_name); CSDact=CSD(je,:)*-1; clear CSD csd_name;
                
                if isnan(CSDact)
                  CSD_rips{je}=NaN;
                  CSD_segavg(je,:)=NaN;  
                else
                    for ka=1:size(swr_pers,1)
                        sta_rip=round(ch_rate*swr_pers(ka,1)); end_rip=round(ch_rate*swr_pers(ka,2));
                        bin_siz=(end_rip-sta_rip)/n_tbins(2);
                        sta_wlt=floor(sta_rip-(bin_siz*(n_tbins(1)*1.4))); end_wlt=ceil(end_rip+(bin_siz*(n_tbins(3)*1.4)));
                        csd_seg=CSDact(sta_wlt:end_wlt);
                        wlt_rip=conj(cwt(csd_seg, rp_scsta:rp_scint:rp_scfin, rp_wavelet{1}));
                        sta_rip=sta_rip-sta_wlt+1; end_rip=end_rip-sta_wlt+1;
                        bin_1st=sta_rip-(bin_siz*n_tbins(1));
                        amp_rip=zeros(length(rp_fraxis),sum(n_tbins));
                        for el=0:sum(n_tbins)-1
                            wlt_el=abs(wlt_rip(:,(ceil(bin_1st+(el*bin_siz))):(floor(bin_1st+((el+1)*bin_siz)))));
                            csd_el=csd_seg((ceil(bin_1st+(el*bin_siz))):(floor(bin_1st+((el+1)*bin_siz))));
                            amp_rip(:,el+1)=mean(wlt_el,2); clear wlt_el
                            amp_csd(el+1)=mean(csd_el); clear csd_el
                        end
                        CSD_rips{je}=CSD_rips{je}+amp_rip; clear amp_rip
                        CSD_segavg(je,:)=CSD_segavg(je,:)+amp_csd; clear amp_csd
                    end; clear sta_rip end_rip bin_siz sta_wlt end_wlt wlt_rip bin_1st ka csd_seg CSDact
                end
            end; clear swr_pers
        end; clear fl_indx
    else CSD_rips{je}=NaN;CSD_ripsN{je}=NaN;
    end
end; clear je

%%plotting results
%%LFP ripple frequency profiles for 
figname=strcat(bs_name, bs_num, bs_unit,' LFP ripple amplitude');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
cmax=1; cmin=1;
for je=1:ch_n
    subplot(4,ch_n/4,ch_n-je+1)
    if ~isnan(LFP_rips{je})
        LFP_rips{je}=LFP_rips{je}./SWR_no;
        LFP_calib=repmat(mean(LFP_rips{je}(:,1:n_tbins(1)),2),1,sum(n_tbins));
        LFP_ripsN{je}=LFP_rips{je}./LFP_calib;
        cmin_A=min(min(LFP_ripsN{je}));cmax_A=max(max(LFP_ripsN{je}));
        if cmin_A<cmin; cmin=cmin_A; end;
        if cmax_A>cmax; cmax=cmax_A; end;
    end
end; 
cmin=cmin*0.9; cmax=cmax*1.1; clear cmax_A cmin_A;


for je=1:ch_n
    subplot(4,ch_n/4,ch_n-je+1)
    if ~isnan(LFP_ripsN{je})
        surf([(-1*(sum(n_tbins)/2))+0.5:1:((sum(n_tbins)/2)-0.5)], rp_fraxis, LFP_ripsN{je}, 'EdgeColor', 'none', 'FaceColor', 'flat'); shading flat; axis tight; colormap jet; caxis([cmin cmax]); view(2);
        xlabel('SWR phase (normalised time)'); ylabel ('frequency (Hz)'); set(gca, 'YTick',[(min(rp_fraxis)):20:(max(rp_fraxis))]);
    end
end;
LFP_cmin=cmin; LFP_cmax=cmax; clear cmin cmax

%%CSD ripple frequency profiles for 
figname=strcat(bs_name, bs_num, bs_unit,' CSD ripple amplitude');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
cmax=1; cmin=1;
for je=ch_step+1:ch_n-ch_step
    subplot(4,ch_n/4,ch_n-je+1)
    if ~isnan(CSD_rips{je})
        CSD_rips{je}=CSD_rips{je}./SWR_no;
        CSD_calib=repmat(mean(CSD_rips{je}(:,1:n_tbins(1)),2),1,sum(n_tbins));
        CSD_ripsN{je}=CSD_rips{je}./CSD_calib;
        cmin_A=min(min(CSD_ripsN{je}));cmax_A=max(max(CSD_ripsN{je}));
        if cmin_A<cmin; cmin=cmin_A; end;
        if cmax_A>cmax; cmax=cmax_A; end;
    end
end;
cmin=cmin*0.9; cmax=cmax*1.1; clear cmax_A cmin_A


for je=ch_step+1:ch_n-ch_step
    subplot(4,ch_n/4,ch_n-je+1)
    if ~isnan(CSD_ripsN{je})
        surf([(-1*(sum(n_tbins)/2))+0.5:1:((sum(n_tbins)/2)-0.5)], rp_fraxis, CSD_ripsN{je}, 'EdgeColor', 'none', 'FaceColor', 'flat'); shading flat; axis tight; colormap jet; caxis([cmin cmax]); view(2);
        xlabel('SWR phase (normalised time)'); ylabel ('frequency (Hz)'); set(gca, 'YTick',[(min(rp_fraxis)):20:(max(rp_fraxis))]);
    end
end;
CSD_cmin=cmin; CSD_cmax=cmax; clear cmin cmax

%%SWR figure
figname=strcat(bs_name, bs_num, bs_unit,' SW averages');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')

%LFP SW
LFP_segavg=LFP_segavg./SWR_no;
cmax=max(abs([1.1*min(min(LFP_segavg)) 1.1*max(max(LFP_segavg))]));
subplot(1,2,1)
surf([(-1*(sum(n_tbins)/2))+0.5:1:((sum(n_tbins)/2)+0.5)],[0.5:1:ch_n+0.5], [[LFP_segavg zeros(ch_n,1)];zeros(1,sum(n_tbins)+1)], 'EdgeColor', 'none', 'FaceColor', 'flat'); shading flat; axis tight; colormap jet; caxis([-cmax cmax]); view(2);
xlabel('SWR phase (normalised time)'); ylabel ('contact#'); set(gca, 'YTick',[1:1:ch_n]);
LFPavg_cmax=cmax; clear cmax

%CSD SW
CSD_segavg=CSD_segavg./SWR_no;
cmax=max(abs([1.1*min(min(CSD_segavg)) 1.1*max(max(CSD_segavg))]));
subplot(1,2,2)
surf([(-1*(sum(n_tbins)/2))+0.5:1:((sum(n_tbins)/2)+0.5)],[0.5:1:ch_n+0.5], [[CSD_segavg zeros(ch_n,1)];zeros(1,sum(n_tbins)+1)], 'EdgeColor', 'none', 'FaceColor', 'flat'); shading flat; axis tight; colormap jet; caxis([-cmax cmax]); view(2);
xlabel('SWR phase (normalised time)'); ylabel ('contact#'); set(gca, 'YTick',[1:1:ch_n]);
CSDavg_cmax=cmax; clear cmax

%saving the relevant data
results.CSDramp=CSD_rips;
results.CSDrampN=CSD_ripsN;
results.CSDrCMAX=CSD_cmax;
results.CSDrCMIN=CSD_cmin;
results.CSDswamp=CSD_segavg;
results.CSDswCMAX=CSDavg_cmax;
results.LFPramp=LFP_rips;
results.LFPrampN=LFP_ripsN;
results.LFPrCMAX=LFP_cmax;
results.LFPrCMIN=LFP_cmin;
results.LFPswamp=CSD_segavg;
results.LFPswCMAX=LFPavg_cmax;
results.swrn=SWR_no;
properties.chdef=ch_def;
properties.chn=ch_n;
properties.smplrate=ch_rate;
properties.files=fl_ord;
properties.swrbinning=n_tbins;
properties.rfraxis=rp_fraxis;
properties.rwlt=rp_wavelet;

datum = date;
%saving the relevant data
writename=strcat(bs_name, bs_num, bs_unit, '_SWRchar');
save(writename,'results','properties','datum');



                    
                    
                               