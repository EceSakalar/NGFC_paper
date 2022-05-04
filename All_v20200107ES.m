%give the experimenter identifier, the experiment number,m and the experiment type
bs_name='ES'; bs_num='49'; bs_exp='d5'; % also on 1335

%This program generates CSD from LFP array data

%give the experiment details
dz=0.05;                                            %contact site spacing (in mm)
ch_n=16;                                            %number of channels
ch_ord=[1:16];    %order of channels, separated by space, starting from the ventralmost to the dorsalmost
ch_def=[];                                          %identifier of deffective channels, separated by space; if extreme channels are deffected that decreases the probe site and does not include here
fl_ord=[1];                                         %file extension numbers belonging to the present cell

%give analysis parameters to be used
ch_inv=1;                                           %channel inversion; write 1 if not needed, write -1 if needed (or a scaling factor if necessary) %0.001 for intan recorded
sam_rate=2000;                                      %requested sampling rate in Hz (LFP resampling will be performed if non-matching)
ch_step = 1;                                        %requested step size for the CSD calculus
plotornot = 1;                                	    %should I plot or not the results

%file index stepper
for fl_indx=1:length(fl_ord);
    
    % name of the actually processed file
    fl_name=strcat(bs_name, bs_num, '_', bs_exp,'_WAW.mat');
    
    %get the length (seconds) of all channels and take the minimum
    fl_length=zeros(size(ch_ord));
    for je=1:ch_n
        ch_name=strcat('Lin1_', num2str(ch_ord(je).','%02d')); ch_load=load(num2str(fl_name),num2str(ch_name));
        fl_length(je)=(((ch_load.(ch_name).('length'))-1)*(ch_load.(ch_name).('interval')))+(ch_load.(ch_name).('start'));
    end; clear je;
    fl_length=min(fl_length);
    sam_length=floor(fl_length*sam_rate);
    
    %LFP import channel by channel; resample if necessary
    LFP = zeros(ch_n,sam_length);
    LFP_tbase = [(1/sam_rate):(1/sam_rate):sam_length*(1/sam_rate)];
    for je=1:ch_n
        %load data for the actual channel
        ch_name=strcat('Lin1_', num2str(ch_ord(je).','%02d')); ch_load=load(num2str(fl_name),num2str(ch_name));
        LFP_in = ch_inv.*(ch_load.(ch_name).('values'));
        LFP_samint = ch_load.(ch_name).('interval');
        LFP_samst = ch_load.(ch_name).('start');
        
        if isempty(ch_def==ch_ord(je))
            
            %check if the time offset and the sampling interval is the desired
            if LFP_samint == 1/sam_rate && LFP_samst == 1/sam_rate
               %put LFP channels into the final variable
                LFP(je,:) = LFP_in(1:sam_length);
            else
                %resample timeseries and transform it back
                LFP_ts = resample((timeseries(LFP_in, [LFP_samst:LFP_samint:((LFP_samint*(length(LFP_in)-1))+LFP_samst)])),LFP_tbase);
                LFP(je,:) = LFP_ts.Data;
                clear LFP_ts
            end
            
            %check the amplitude unit, end value should be in mV
            if ch_load.(ch_name).('units')=='uV'
                LFP(je,:)=LFP(je,:)/1000;
            end
            
            clear LFP_samint LFP_samst LFP_in ch_name 
        end
    end;
    clear je ch_load LFP_tbase
    
    %interpolating the deffective channels (linear interpolation from the adjacent channels)
    if ~isempty(ch_def);
        for ka=1:length(ch_def);
            je=find(ch_ord==(ch_def(ka)));
            LFP(je,:)=(LFP(je-1,:)+LFP(je+1,:))./2;
        end; clear ka je;
    end
    clear fl_length
    
    %generate the CSD operator matrix
    CSD_matr = zeros(ch_n-(2*ch_step), ch_n);
    for je=1:(ch_n-(2*ch_step))
        CSD_matr(je,je)=1;
        CSD_matr(je,je+ch_step)=-2;
        CSD_matr(je,je+(2*ch_step))=1;
    end; clear je
    
    %CSD calculus
    CSD = (CSD_matr*LFP)./((dz*ch_step)^2);
    CSD_dum(1:ch_step,1:size(CSD,2)) = NaN;
    CSD = [CSD_dum; CSD; CSD_dum];
    
    %writing files
    fl_write=strcat(bs_name, bs_num, '_', bs_exp, '_LFP');
    save(fl_write, 'LFP');
    fl_write=strcat(bs_name, bs_num, '_', bs_exp, '_CSD');
    save(fl_write, 'CSD');
    
    if plotornot==1
        for je=1:ch_n
            LFPplot(je,:)=LFP(je,1:5000)+(0.7*je);
        end; clear je;
        for je=1:(ch_n-2)
            CSDplot(je,:)=CSD(je,1:5000)+(100*je);
        end; clear je;
        figure
        plot(LFPplot');xlim([1 2000]);title('LFP');
        plot_write=strcat(bs_name, bs_num, '_', bs_exp, '_LFP.fig');
        savefig(plot_write)
        figure
        plot(CSDplot');xlim([1 2000]);title('CSD');
        plot_write=strcat(bs_name, bs_num, '_', bs_exp, '_CSD.fig');
        savefig(plot_write)
    end
    clear CSD LFP CSD_dum CSDplot LFPplot CSD_matr fl_length fl_name fl_write LFP_tbase sam_length
end

   
%This program analyses the SWR coupling and ripple oscillatory coupling relative to LFP and CSD traces

%basic parameters
ch_n=16;                                              %number of channels
ch_rate=2000;                                         %give the silicon probe sampling rate (in Hz)
ch_def=[];                                           %defective channel (one is standing for the ventralmost, and the linear number of channel is required not the identifier
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
            swr_pname=strcat(bs_name, bs_num, '_', bs_exp,'_SWP.txt');     %SWR periods
            
            if exist(swr_pname,'file')==2;
                swr_pers=dlmread(swr_pname); clear swr_pname;
                %reading LFP for ripple LFP characterisation
                lfp_name=strcat(bs_name, bs_num, '_', bs_exp,'_LFP.mat');
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
            swr_pname=strcat(bs_name, bs_num, '_', bs_exp,'_SWP.txt');     %SWR periods
            
            if exist(swr_pname,'file')==2;
                swr_pers=dlmread(swr_pname); clear swr_pname;
                %reading CSD for ripple CSD characterisation
                csd_name=strcat(bs_name, bs_num, '_', bs_exp,'_CSD.mat');
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
figname=strcat(bs_name, bs_num, '_', bs_exp, '_LFP ripple amplitude');
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
plot_write=strcat(bs_name, bs_num, '_', bs_exp, '_LFP_SWP.fig');
savefig(plot_write)
%%CSD ripple frequency profiles for 
figname=strcat(bs_name, bs_num, '_', bs_exp, '_CSD ripple amplitude');
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
plot_write=strcat(bs_name, bs_num, '_', bs_exp, '_CSD_SWP.fig');
savefig(plot_write)

%%SWR figure
figname=strcat(bs_name, bs_num, '_', bs_exp, ' SW averages');
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

plot_write=strcat(bs_name, bs_num, '_', bs_exp, '_SWaver.fig');
savefig(plot_write)

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
writename=strcat(bs_name, bs_num, '_', bs_exp, '_SWRchar');
save(writename,'results','properties','datum');


%This code analyses SWR dependent firing: It does histograms for SWR dependent firing 

%unit firing information content
fl_ord=[1];                                       %file extension numbers belonging to the present cell

un_chans=[1];                                                          %channels containing the spike times for the units
%unit numbers;
un_nums=[1];                                    %unit numbers (in kwik file)
%shanks from which unit was isolated;
un_shnk =[repmat([1],1,1) repmat([2],1,0) ];                              %shanks from which unit was isolated

%analysis information
%give the number of bins for the center (middle) and surroundings (side) of SWRs histogram
middle=8; side=12;
%give the binsize (ms) and the time before and after (sammple number) for the startaligned and endaligned histograms
binsize=20; binsize=binsize/1000; before=20; after=20;

%script starts
SWR_totaltime=0;
    
for un_indx=1:length(un_chans)
    SWR_spikecount{un_indx}=[];
end
for un_indx=1:length(un_chans)
    
    for fl_indx=1:length(fl_ord); %file stepper
        fl_name=strcat(bs_name, bs_num, '_', bs_exp, '_WAW.mat');   %file name
        sw_pname=strcat(bs_name, bs_num, '_', bs_exp, '_SWP.txt');  %swr periods
        
        %reading the spike sequence for the cell
        ch_name=strcat('Gl_unt_', num2str(un_chans(un_indx).','%02d'));
        ch_load=load(num2str(fl_name),num2str(ch_name));
        fl_spks=(ch_load.(ch_name).('times')); clear ch_name ch_load fl_name;
        
        if exist(sw_pname,'file')==2; sw_pers=dlmread(sw_pname); clear sw_pname;
            
            fl_swrspks=periodcut(fl_spks,sw_pers(:,1:2),1);
            
            %I calculate the SWR time spent in this file
            if un_indx==1
                SWR_totaltime_file=sum(((sw_pers(:,2))-(sw_pers(:,1))),1);
                SWR_totaltime=SWR_totaltime+SWR_totaltime_file; clear SWR_totaltime_file;
            end
            
            %I calculate the bin sizes
            bin_start=(sw_pers(:,3)-sw_pers(:,1))/(middle/2); bin_end=(sw_pers(:,2)-sw_pers(:,3))/(middle/2); bin_side=(bin_start+bin_end)/2;
            
            %I calculate binmatrix for the four periods (before the SWR, first and second SWR halves and after the SWR)
            SWR_binstarts_k=repmat((sw_pers(:,1)-(bin_side*side)),1,(middle+2*side));
            SWR_binstarts_a=cumsum([zeros(size(bin_side,1),1) repmat(bin_side,1,side) repmat(bin_start,1,middle/2) repmat(bin_end,1,middle/2) repmat(bin_side,1,(side-1))],2);
            SWR_binends_a=cumsum([repmat(bin_side,1,side) repmat(bin_start,1,middle/2) repmat(bin_end,1,middle/2) repmat(bin_side,1,side)],2);
            SWR_binstarts_f=SWR_binstarts_k+SWR_binstarts_a;
            SWR_binends_f=SWR_binstarts_k+SWR_binends_a;
            clear bin_start bin_end bin_side SWR_binstarts_k SWR_binstarts_a SWR_binends_a;
            
            %I calculate the spike numbers for all the bins
            for j=1:size(SWR_binstarts_f,1)
                for k=1:size(SWR_binstarts_f,2)
                    SWR_spikecount_f(j,k)=length(find(fl_spks>=SWR_binstarts_f(j,k)&fl_spks<SWR_binends_f(j,k)));
                end
            end; clear j k
            SWR_bintime_f=SWR_binends_f-SWR_binstarts_f;
            
            SWR_spikecounttoblank_f=zeros(size(SWR_spikecount_f));
            SWR_bintimetoblank_f=zeros(size(SWR_bintime_f));
            
            %I'm blanking the spikes included in neighbouring SWRs
            for j=1:size(SWR_binstarts_f,1)
                for k=1:side
                    SWR_spikecounttoblank_f(j,k)=length(find(fl_swrspks>=SWR_binstarts_f(j,k)&fl_swrspks<SWR_binends_f(j,k)));
                end; clear k
                for k=(side+middle+1):(side+middle+side)
                    SWR_spikecounttoblank_f(j,k)=length(find(fl_swrspks>SWR_binstarts_f(j,k)&fl_swrspks<=SWR_binends_f(j,k)));
                end; clear k
            end; clear j
            
            %I'm blanking the time included in neighbouring SWRs
            for j=1:size(SWR_binstarts_f,1)
                if j>1
                    for k=1:side
                        ind=find(sw_pers(1:j-1,1)<SWR_binends_f(j,k)&sw_pers(1:j-1,2)>SWR_binstarts_f(j,k));
                        if size(ind,1)>0
                            for ind2=1:size(ind,1)
                                if SWR_binstarts_f(j,k)>sw_pers(ind(ind2),1)&&SWR_binends_f(j,k)<sw_pers(ind(ind2),2);
                                    SWR_bintimetoblank_f(j,k)=SWR_binends_f(j,k)-SWR_binstarts_f(j,k);
                                elseif SWR_binstarts_f(j,k)<=sw_pers(ind(ind2),1)
                                    SWR_bintimetoblank_f(j,k)=SWR_bintimetoblank_f(j,k)+(SWR_binends_f(j,k)-sw_pers(ind(ind2),1));
                                elseif SWR_binends_f(j,k)>=sw_pers(ind(ind2),2)
                                    SWR_bintimetoblank_f(j,k)=SWR_bintimetoblank_f(j,k)+(sw_pers(ind(ind2),2)-SWR_binstarts_f(j,k));
                                end
                            end; clear ind2;
                        end
                    end; clear ind k;
                end
                
                if j<size(SWR_binstarts_f, 1);
                    for k=(middle+side+1):(middle+(2*side))
                        ind=find(sw_pers(:,1)<SWR_binends_f(j,k)&sw_pers(:,2)>SWR_binstarts_f(j,k));
                        if size(ind,1)>0
                            for ind2=1:size(ind,1)
                                if SWR_binstarts_f(j,k)>sw_pers(ind(ind2),1)&&SWR_binends_f(j,k)<sw_pers(ind(ind2),2);
                                    SWR_bintimetoblank_f(j,k)=SWR_binends_f(j,k)-SWR_binstarts_f(j,k);
                                elseif SWR_binstarts_f(j,k)<=sw_pers(ind(ind2),1)
                                    SWR_bintimetoblank_f(j,k)=SWR_bintimetoblank_f(j,k)+(SWR_binends_f(j,k)-sw_pers(ind(ind2),1));
                                elseif SWR_binends_f(j,k)>=sw_pers(ind(ind2),2)
                                    SWR_bintimetoblank_f(j,k)=SWR_bintimetoblank_f(j,k)+(sw_pers(ind(ind2),2)-SWR_binstarts_f(j,k));
                                end
                            end; clear ind2
                        end
                    end; clear ind k
                end
            end; clear j;
            
            %calculating the SWR-instant spike-rate plots
            for j=1:size(sw_pers,1)
                sw_actspks=periodcut(fl_spks,sw_pers(j,:),1);
                if length(sw_actspks)>1
                    sw_actspks_pos=(sw_actspks-sw_pers(j,1))./(sw_pers(j,2)-sw_pers(j,1));
                    for k=1:length(sw_actspks)-1
                        sw_ratedat_per(k,1)=(sw_actspks_pos(k)+sw_actspks_pos(k+1))/2;
                        sw_ratedat_per(k,2)=1/(sw_actspks(k+1)-sw_actspks(k));
                    end; clear sw_actspks_pos k
                    
                    if exist('sw_ratedat','var')
                        sw_ratedat=[sw_ratedat; sw_ratedat_per];
                    else
                        sw_ratedat=sw_ratedat_per;
                    end; clear sw_ratedat_per;
                end; clear sw_actspks;
            end; clear j;
            
            %generating the endaligned and the startaligned binned histograms
            aligned_binner=repmat((-1*before*binsize)-(binsize/2):binsize:(after*binsize)+(binsize/2),size(sw_pers,1),1);
            start_binner=repmat(sw_pers(:,1),1,size(aligned_binner,2));
            end_binner=repmat(sw_pers(:,2),1,size(aligned_binner,2));
            startaligned_binner=start_binner+aligned_binner;
            endaligned_binner=end_binner+aligned_binner;
            startaligned_f=zeros(size(startaligned_binner)); endaligned_f=zeros(size(endaligned_binner));
            clear aligned_binner start_binner end_binner
            %actual binning
            for j=1:size(sw_pers,1)
                startaligned_f(j,:)=hist(fl_spks,startaligned_binner(j,:));
                endaligned_f(j,:)=hist(fl_spks,endaligned_binner(j,:));
            end; clear j startaligned_binner endaligned_binner
            
            %I'm gathering the data for the actual unit SWR
            if isempty(SWR_spikecount{un_indx})
                SWR_startaligned{un_indx}=startaligned_f;
                SWR_endaligned{un_indx}=endaligned_f;
                SWR_spikecount{un_indx}=SWR_spikecount_f;
                SWR_bintime{un_indx}=SWR_bintime_f;
                SWR_spikecounttoblank=SWR_spikecounttoblank_f;
                SWR_bintimetoblank=SWR_bintimetoblank_f;
            else
                SWR_startaligned{un_indx}=[SWR_startaligned{un_indx}; startaligned_f];
                SWR_endaligned{un_indx}=[SWR_endaligned{un_indx}; endaligned_f];
                SWR_spikecount{un_indx}=[SWR_spikecount{un_indx}; SWR_spikecount_f];
                SWR_bintime{un_indx}=[SWR_bintime{un_indx}; SWR_bintime_f];
                SWR_spikecounttoblank=[SWR_spikecounttoblank; SWR_spikecounttoblank_f];
                SWR_bintimetoblank=[SWR_bintimetoblank; SWR_bintimetoblank_f];
            end
            clear SWR_spikecount_f SWR_bintime_f SWR_spikecounttoblank_f SWR_bintimetoblank_f SWR_binstarts_f SWR_binends_f
            clear sw_pers fl_swrspks startaligned_f endaligned_f
        end; clear fl_spks
    end; clear fl_indx
    
    SWR_spikecountblanked{un_indx}=SWR_spikecount{un_indx}-SWR_spikecounttoblank;
    SWR_bintimeblanked{un_indx}=SWR_bintime{un_indx}-SWR_bintimetoblank;
    clear SWR_bintimetoblank SWR_spikecounttoblank
end

SWR.Xaxis=[1:1:middle+(2*side)];
SWR.totaltime=SWR_totaltime;
SWR.aligned_xaxis=[((-1*(binsize*(before-0.5)))*1000):(binsize*1000):(((after-0.5)*binsize)*1000)];
%plotting the results
figindex=0;
for un_indx=1:length(un_chans)
    figindex_T=ceil(un_indx/6);
    if figindex_T~=figindex
        figindex=figindex+1;
        figname=strcat(bs_name, bs_num, '_', bs_exp,'_autocorrelogram');
        figure ('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
        plotindex=1;
    else plotindex=((un_indx-((figindex-1)*6))*3)-2;
    end
    subplot (6,3,plotindex);
    %textual results
    if max(max(SWR_spikecountblanked{un_indx}))>0
        axis off
        text(0,1,strcat(bs_name, bs_num, bs_exp, ' t', num2str(un_shnk(un_indx)),'u',num2str(un_nums(un_indx))))
        %calculus
        SWR.spikenumberH{un_indx}=sum(SWR_spikecountblanked{un_indx},1);
        SWR.number{un_indx}=size(SWR_spikecount{un_indx},1);
        SWR.activenumber{un_indx}=size(find((sum((SWR_spikecountblanked{un_indx}(:,(side+1):(side+middle))),2))>0),1);
        SWR.spikenumber{un_indx}=sum(SWR.spikenumberH{un_indx}((side+1):(side+middle)));
        SWR.spikefrequency{un_indx}=SWR.spikenumber{un_indx}/SWR.totaltime;
        %number of SWRs
        actualtext= strcat('Number of SWRs analysed: ', num2str (SWR.number{un_indx})); text (0, 0.8, actualtext);
        %number of active SWRs
        actualtext= strcat('Number of active SWRs : ', num2str (SWR.activenumber{un_indx})); text (0, 0.6, actualtext);
        %mean number of APs in SWRs
        actualtext= strcat('Mean number of APs in an SWR:', num2str((SWR.spikenumber{un_indx})/(SWR.activenumber{un_indx}))); text (0, 0.4, actualtext);
        %mean AP frequency
        actualtext=strcat('Mean AP frequency under SWR:', num2str (SWR.spikefrequency{un_indx})); text (0, 0.2, actualtext);
        %overall time spent in SWR
        actualtext=strcat('overall SWR time:', num2str (SWR.totaltime), ' s'); text (0, 0, actualtext);
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
        %SWR spike frq histrogram
        subplot(6,3,plotindex+1)
        SWR.spikefrqH{un_indx}=(sum(SWR_spikecountblanked{un_indx},1))./(sum(SWR_bintimeblanked{un_indx},1));
        bar (SWR.Xaxis, SWR.spikefrqH{un_indx});
        set (gca, 'Xtick', 0:1:middle+(2*side));
        axis ([0 (middle+(2*side)+1) 0 max(SWR.spikefrqH{un_indx})*1.1]);
        xlabel('phase (of SWR)'); ylabel('spike/s');
        line([side+0.5 side+0.5],[0 max(SWR.spikefrqH{un_indx})*1.1], 'Linewidth',2,'Color','red');
        line([side+middle+0.5 side+middle+0.5],[0 max(SWR.spikefrqH{un_indx})*1.1], 'Linewidth',2,'Color','red');
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
        %histograms aligned to the start or to the end of SWRs
        subplot(6,3,plotindex+2)
        SWR_startaligned{un_indx}(:,1)=[];SWR_startaligned{un_indx}(:,end)=[];
        SWR_endaligned{un_indx}(:,1)=[]; SWR_endaligned{un_indx}(:,end)=[];
        SWR.startaligned{un_indx}=sum(SWR_startaligned{un_indx},1)./(binsize*size(SWR_startaligned{un_indx},1));
        SWR.endaligned{un_indx}=sum(SWR_endaligned{un_indx},1)./(binsize*size(SWR_endaligned{un_indx},1));
        line(SWR.aligned_xaxis,SWR.startaligned{un_indx},'Linewidth',2,'Color','red');line(SWR.aligned_xaxis,SWR.endaligned{un_indx},'Linewidth',2,'Color','blue');
        set(gca, 'Xtick', (before*binsize*-1000):binsize*2000:(after*binsize*1000));
        axis([min(SWR.aligned_xaxis) max(SWR.aligned_xaxis) 0 (max([max(SWR.startaligned{un_indx}) max(SWR.endaligned{un_indx})]))*1.1]);
        xlabel ('time (ms)'); ylabel ('spike rate(s-1)');  grid on;
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
    end
end

plot_write=strcat(bs_name, bs_num, '_', bs_exp, '_SWP_Gl_unt.fig');
savefig(plot_write)

clear SWR_startaligned SWR_endaligned
clear SWR_spikecountblanked SWR_bintimeblanked SWR_spikecounttoblank SWR_bintimetoblank SWR_bintime SWR_spikecount
clear after before binsize middle side actualtext

datum=date;

%saving the relevant data
writename=strcat(bs_name, bs_num, '_', bs_exp, '_a_swr_eventcoupling');
save(writename, 'SWR', 'datum');

clear writename figname
        
        
%This code analyses SWR dependent firing: It does a shuffling for the SWR
%window to determine whether there are significantly more or less spikes
%found than expected in the surrounding time period

fl_ord=[1];                                                 %file extension numbers belonging to the present cell
un_chans=[1];                       %channels containing the spike times for the units
un_num=[1];                                                %unit numbers (in kwik file)
un_shnk =[repmat([1],1,1) repmat([2],1,0)];                  %tetrodes from which unit was isolated


%USER INTERVENTION
%set the numer of random placements to be performed for one SWR and the size of the random window to be included in the random placements (in sec).
Rnum = 10;
Rwin = 2;
%indicate one lfp channel for the file length extraction
ch_lfp=15;

%get the file lengths (starting points in seconds in one vector)
Fstarts=0;
for fl_indx=1:length(fl_ord)
    fl_name=strcat(bs_name, bs_num, '_', bs_exp,'_WAW.mat');   %file name
    ch_name=strcat('Lin1_', num2str(ch_lfp.','%02d'));
    ch_load=load(num2str(fl_name),num2str(ch_name));
    Fdata=(ch_load.(ch_name).('length'));
    Finter=(ch_load.(ch_name).('interval'));
    Flength=Fdata*Finter; clear Fdata Finter h_load ch_name fl_name;
    Fstarts=[Fstarts Fstarts(fl_indx)+Flength]; clear Flength;
end; clear fl_indx ch_lfp

%get the SWR periods together for the whole recording
sw_pers=[];
for fl_indx=1:length(fl_ord)
    fl_name=strcat(bs_name, bs_num, '_', bs_exp, '_WAW.mat');   %file name
    sw_pname=strcat(bs_name, bs_num, '_', bs_exp, '_SWP.txt');  %swr periods
    if exist(sw_pname,'file')==2;
        sw_persF=dlmread(sw_pname)+Fstarts(fl_indx); 
        if isempty(sw_pers)
            sw_pers=sw_persF;
            clear sw_persF
        else
            sw_pers=[sw_pers; sw_persF];
            clear sw_persF
        end
    end;
    clear sw_pname fl_name;
end; clear fl_indx;
sw_lengths=sw_pers(:,2)-sw_pers(:,1);
sw_cumlength=cumsum(sw_lengths);
sw_gaps=sw_pers(:,1)-[0; sw_cumlength(1:(end-1))];
%get the unit spike times together for the whole recording
spks{length(un_chans)}=[];
for un_indx=1:length(un_chans)
    for fl_indx=1:length(fl_ord); %file stepper
        fl_name=strcat(bs_name, bs_num, '_', bs_exp,'_WAW.mat');   %file name
        %reading the spike sequence for the cell file
        ch_name=strcat('Gl_unt_', num2str(un_indx.','%02d'));
        ch_load=load(num2str(fl_name),num2str(ch_name));
        fl_spks=(ch_load.(ch_name).('times')); clear ch_name ch_load fl_name;
        fl_spks=fl_spks+Fstarts(fl_indx);
        if isempty(spks{un_indx})
            spks{un_indx}=fl_spks;
            clear fl_spks
        else
            spks{un_indx}=[spks{un_indx}; fl_spks];
            clear fl_spks
        end
    end; clear fl_indx
end;clear un_indx

%making of modified variables
for un_indx=1:length(un_chans)
    spks_all=spks{un_indx};
    spks_out=periodcut(spks_all, sw_pers(:,1:2), 0);
    spks_adj=zeros(size(spks_out));
    for je=1:size(sw_pers,1)
        spks_adj(spks_out>sw_pers(je,1))=sw_cumlength(je);
    end; clear je
    spks_rem=spks_out-spks_adj; clear spks_adj
    
    %extract data spk nums and shuffled numbers
    num_in{un_indx} = zeros(size(sw_pers,1),1);
    num_shf{un_indx}= zeros(size(sw_pers,1),Rnum);
    
    %stepper for swrs
    for je=1:size(sw_pers,1)
        %number of spikes in swr
        num_in{un_indx}(je)=length(spks_all(spks_all>sw_pers(je,1)&spks_all<=sw_pers(je,2)));
        
        %shuffling
        if sw_gaps(je,1)<Rwin/2
            for ka=1:Rnum
                sh_sta = (rand(1)*(sw_gaps(je)+(Rwin/2)-sw_lengths(je)));
                sh_end = sh_sta+sw_lengths(je);
                num_shf{un_indx}(je,ka)=length(spks_rem(spks_rem>sh_sta&spks_rem<=sh_end));
                clear sh_sta sh_end
            end; clear ka
        elseif sw_gaps(je,1)+sw_lengths(je,1)>(Fstarts(end)-sum(sw_lengths)-(Rwin/2))
            for ka=1:Rnum
                sh_sta = (Fstarts(end)-sum(sw_lengths))-(rand(1)*((Rwin/2)+(Fstarts(end)-sum(sw_lengths))-sw_gaps(je)));
                sh_end = sh_sta+sw_lengths(je);
                num_shf{un_indx}(je,ka)=length(spks_rem(spks_rem>sh_sta&spks_rem<=sh_end));
                clear sh_sta sh_end
            end; clear ka
        else
            for ka=1:Rnum
                sh_sta=(sw_gaps(je)-(rand(1)*(Rwin-sw_lengths(je))))+(Rwin/2);
                sh_end = sh_sta+sw_lengths(je);
                num_shf{un_indx}(je,ka)=length(spks_rem(spks_rem>sh_sta&spks_rem<=sh_end));
                clear sh_sta sh_end
            end; clear ka
        end
    end; clear je
end; clear un_indx


%----PLOTS
figindex=0;
for un_indx=1:length(un_chans)
    figindex_T=ceil(un_indx/6);
    if figindex_T~=figindex
        figindex=figindex+1;
        figname=strcat(bs_name, bs_num, '_', bs_exp,'_SWRstats');
        figure ('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
        plotindex=1;
    else plotindex=((un_indx-((figindex-1)*6))*3)-2;
    end
    
    if (max(num_in{un_indx})+max(max(num_shf{un_indx})))>0
        %1st histogram - spike numbers in SWRs
        subplot (6,3,plotindex)
        Xaxis = [0:1:max([max(num_in{un_indx}) max(max(num_shf{un_indx}))])];
        hist(num_in{un_indx}, Xaxis);
        title ('distribution in SWRs'); xlabel ('spike count'); ylabel ('SWR count');
        axis ([-1 max([max(num_in{un_indx}) max(max(num_shf{un_indx}))])+1 0 1.1*(max(hist(num_in{un_indx}, Xaxis)))]);
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
        %2nd histogram - spike numbers in shuffling windows
        subplot(6,3,plotindex+1)
        Xaxis = [0:1:max([max(num_in{un_indx}) max(max(num_shf{un_indx}))])];
        bar(Xaxis,sum(hist(num_shf{un_indx}, Xaxis)'));
        title ('distribution in shuffled windows'); xlabel ('spike count'); ylabel ('window count');
        axis ([-1 max([max(num_in{un_indx}) max(max(num_shf{un_indx}))])+1 0 1.1*max(sum(hist(num_shf{un_indx}, Xaxis)'))]);
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
    end
    %textual results - statistic
    subplot(6,3,plotindex+2)
    axis off;
    actualtext= strcat('unit:', bs_name, bs_num, bs_exp, 't', num2str(un_shnk(un_indx)),'u',num2str(un_num(un_indx))); text(-0.2, 1, actualtext);
    actualtext= strcat('mean SWR duration:', num2str(mean(sw_lengths)*1000), 'ms'); text(-0.2, 0.88, actualtext);
    actualtext= strcat('in/out N:', num2str(sum(num_in{un_indx})),'/',num2str(mean(sum(num_shf{un_indx})))); text(-0.2, 0.76, actualtext);
    actualtext= strcat('mean in SWR frq:', num2str(sum(num_in{un_indx})/sum(sw_lengths)), 'Hz'); text(-0.2, 0.64, actualtext);
    actualtext= strcat('mean out SWR frq:', num2str((sum(sum(num_shf{un_indx})))/(sum(sw_lengths)*Rnum)), 'Hz'); text(-0.2, 0.52, actualtext);
    [p,h,stati]=ranksum(num_in{un_indx}, reshape(num_shf{un_indx},size(num_shf{un_indx},1)*size(num_shf{un_indx},2),1));
    actualtext=strcat('P/ranksum(Mann-Whitney U):', num2str(p),'/',num2str(stati.ranksum)); text(-0.2, 0.40, actualtext);
    [h,p]=kstest2(num_in{un_indx}, reshape(num_shf{un_indx},size(num_shf{un_indx},1)*size(num_shf{un_indx},2),1));
    actualtext=strcat('P(Kolmogorov-Smirnov):', num2str(p)); text(-0.2, 0.28, actualtext);
    actualtext=strcat('shuffling number:', num2str(Rnum)); text(-0.2, 0.16, actualtext);
    actualtext=strcat('shuffling window:', num2str(Rwin), 's'); text(-0.2, 0.04, actualtext);
    %expand the graph
    pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
    set(gca, 'Clipping', 'off')
end

plot_write=strcat(bs_name, bs_num, '_', bs_exp, '_SWP_Gl_unt_stats.fig');
savefig(plot_write)

clear p h spks_all spks_out spks_rem sw_cumlength Xaxis actualtext figindex figindex_T figname plotindex stati

datum=date;

%saving the relevant data
writename=strcat(bs_name, bs_num, '_', bs_exp, '_a_swr_eventsignificance');
save(writename, 'Fstarts', 'Rnum','Rwin','num_in','num_shf','spks','sw_pers','sw_gaps','datum');
 
clear writename figname
         
        
%2018.04.13 bug corrected in line 82:
%spks_rem=spks_out+spks_adj; clear spks_adj changed to spks_rem=spks_out-spks_adj; clear spks_adj

       
 %This program generates theta modulation of unit firing relative to 'pyramidal layer' LFP 

ch_dwnsmpl=5;                                       %give the downsampling rate
ch_rate=2000; ch_rate=ch_rate/ch_dwnsmpl;           %give the silicon probe sampling rate (in Hz)

un_chans=[1];                                                           %channels containing the spike times for the units
un_nums=[1];          %unit numbers (in kwik file)
%shanks from which unit was isolated; shank 0 is glass electrode and allows shortening of the recording time
un_shnk =[repmat([1],1,1) repmat([2],1,0)];                              

ch_pyr=12;                                       %channel number (ventralmost LFP contact is 1) for the pyramidal layer1
fl_ord=[0];                                     %file extension numbers belonging to the present cell

%here you can give the limits of inclusion for each file (rows) for each unit (colums) on glass electrode
fl_limS = [0]; %starts
fl_limE = [10000000]; %ends

%number of theta bins
th_nbins = 20;                                      
th_binsize = 360/th_nbins;

%filter properties 
flt_ncoeff=512; flt_low=5; flt_high=12;

%generating empty variables
for un_indx=1:length(un_chans)
    th_coupvect_LFP{un_indx}=[];
    cycdur{un_indx}=[];
end
th_time=zeros(length(un_nums),1)';
th_cyccount=zeros(length(un_nums),1)';

%CALCULUS
for fl_indx=1:length(fl_ord); %FILE STEPPER
    fl_name=strcat(bs_name, bs_num, '_', bs_exp, '_WAW.mat');   %file name
    
    %THETA ANALYSIS
    %reading pyramidal layer LFP
    lfp_name=strcat(bs_name, bs_num, '_', bs_exp, '_LFP.mat');
    load(lfp_name);
    LFPred=downsample(smooth(LFP(ch_pyr,:),ch_dwnsmpl),ch_dwnsmpl);
    clear LFP;
    
    %creating theta filter
    th_filter=dfilt.dffir(fir1(flt_ncoeff, [2.*(flt_low./ch_rate) 2.*(flt_high./ch_rate)],'bandpass', gausswin(flt_ncoeff+1)));
    
    %filtering for theta
    LFPflt=filtfilt(th_filter.Numerator,1,LFPred')'; clear LFPred;
    
    %reading theta periods
    th_pname=strcat(bs_name, bs_num, '_', bs_exp, '_THP.txt');
    
    if exist(th_pname,'file')==2;
        
        %detecting theta troughs for pyramidal layer LFP
        LFPthreshold=0.001*(mean(abs(LFPflt)));
        [~,th_troughsLFP]=findpeaks(LFPflt.*-1,'MINPEAKHEIGHT', LFPthreshold, 'MINPEAKDISTANCE', floor(ch_rate/flt_high)-3);
        th_troughsLFP=(th_troughsLFP./ch_rate);%-(1/ch_rate)+(1/(ch_rate*ch_dwnsmpl));
        
        for un_indx=1:length(un_chans)
            %reading theta periods
            th_pers=dlmread(th_pname);
            
            %limiting theta periods if necessary
            if un_shnk(un_indx)==0;
                for je=1:size(th_pers,1)
                    if th_pers(je,2)<fl_limS(fl_indx, un_nums(un_indx))
                        th_pers(je,:)=0;
                    elseif th_pers(je,1)<fl_limS(fl_indx, un_nums(un_indx))
                        th_pers(je,1)=fl_limS(fl_indx, un_nums(un_indx));
                    end
                    if th_pers(je,1)>fl_limE(fl_indx, un_nums(un_indx))
                        th_pers(je,:)=0;
                    elseif th_pers(je,2)>fl_limE(fl_indx, un_nums(un_indx))
                        th_pers(je,2)=fl_limE(fl_indx, un_nums(un_indx));
                    end
                end; clear je
                th_persN=[th_pers(th_pers(:,1)>0,1) th_pers(th_pers(:,2)>0,2)];
                th_pers=th_persN; clear th_persN
            end
            
            %calculating theta time
            fl_th_time(un_indx)=sum(th_pers(:,2)-th_pers(:,1));
            th_time(un_indx)=th_time(un_indx)+fl_th_time(un_indx);
            
            %reading the spike sequence for the cell
            ch_name=strcat('Gl_unt_', num2str(un_nums(un_indx).','%02d'));
            ch_load=load(num2str(fl_name),num2str(ch_name));
            fl_spks=(ch_load.(ch_name).('times')); clear ch_name ch_load;
            
            fl_thspks=periodcut(fl_spks,th_pers,1);
            
            if isempty(fl_thspks)
                th_coupvect_LFPact=[];
            else
                for el=1:length(fl_thspks)
                    spk_act=fl_thspks(el);
                    cyc_beg=max(th_troughsLFP(th_troughsLFP<=spk_act));cyc_end=min(th_troughsLFP(th_troughsLFP>spk_act));
                    cyc_LFP=LFPflt((cyc_beg*ch_rate):(cyc_end*ch_rate)); [~,cyc_cent]=max(cyc_LFP);cyc_cent=(cyc_cent/ch_rate)+cyc_beg;
                    if  spk_act<=cyc_cent;      th_coupvect_LFPact(el)=180*((spk_act-cyc_beg)/(cyc_cent-cyc_beg));
                    else                        th_coupvect_LFPact(el)=180*((spk_act-cyc_cent)/(cyc_end-cyc_cent))+180;
                    end; clear spk_act cyc_beg cyc_cent cyc_end cyc_LFP;
                    
                end; clear el;
            end
            
            %cycle statistics - spike count
            fl_th_cycspkcount=histc(fl_thspks,th_troughsLFP);
            fl_th_cycspkcount(fl_th_cycspkcount==0)=[];
            %cycle statistics - cycle count
            fl_th_cyccount=length(periodcut(th_troughsLFP,th_pers,1))+size(th_pers,1);
            th_cyccount(un_indx)=th_cyccount(un_indx)+fl_th_cyccount;
            clear fl_th_cyccount
            
            %sum coupling vectors accross files
            if isempty(th_coupvect_LFP{un_indx})
                th_coupvect_LFP{un_indx}=th_coupvect_LFPact;
                th_cycspkcount{un_indx}=fl_th_cycspkcount;
            else
                th_coupvect_LFP{un_indx}=[th_coupvect_LFP{un_indx}  th_coupvect_LFPact];
                th_cycspkcount{un_indx}=[th_cycspkcount{un_indx}; fl_th_cycspkcount];
            end; clear th_coupvect_LFPact fl_th_cycspkcount
            
            %bin time calculus
            for el=1:size(th_pers,1)
                perdur=th_pers(el,2)-th_pers(el,1);
                th_acttroughs=periodcut(th_troughsLFP',th_pers(el,:),1);
                
                %calculates cycle boundaries for each cycle in this THperiod
                for em=1:(length(th_acttroughs)+1);
                    if em==1
                        cyc_end=th_acttroughs(em);
                        cyc_beg=max(th_troughsLFP(th_troughsLFP<th_acttroughs(em)));
                    elseif em<=length(th_acttroughs)
                        cyc_end=th_acttroughs(em);
                        cyc_beg=th_acttroughs(em-1);
                    else
                        cyc_end=min(th_troughsLFP(th_troughsLFP>th_acttroughs(em-1)));
                        cyc_beg=th_acttroughs(em-1);
                    end
                    cyc_LFP=LFPflt((cyc_beg*ch_rate):(cyc_end*ch_rate));
                    [~,cyc_cent]=max(cyc_LFP);cyc_cent=(cyc_cent/ch_rate)+cyc_beg;
                    cycs(em,1)=cyc_beg; cycs(em,2)=cyc_cent;cycs(em,3)=cyc_end;
                    clear cyc_beg cyc_end cyc_cent;
                end; clear em;
                
                cycs=cycs-th_pers(el,1);
                cycdur_p=zeros(size(cycs,1),2);
                
                %cycle calculates cycle half lengths within a period
                for em=1:size(cycs,1)
                    if em==1
                        if cycs(em,2)<=0
                            cycdur_p(em,2)=cycs(em,3);
                        else
                            cycdur_p(em,1)=cycs(em,2);
                            cycdur_p(em,2)=cycs(em,3)-cycs(em,2);
                        end
                    elseif em==size(cycs,1)
                        if cycs(em,2)>perdur
                            cycdur_p(em,1)=perdur-cycs(em,1);
                        else
                            cycdur_p(em,1)=cycs(em,2)-cycs(em,1);
                            cycdur_p(em,2)=perdur-cycs(em,2);
                        end
                    else
                        cycdur_p(em,1)=cycs(em,2)-cycs(em,1);
                        cycdur_p(em,2)=cycs(em,3)-cycs(em,2);
                    end
                end;
                
                %summing the period data
                if isempty (cycdur{un_indx})
                    cycdur{un_indx}=cycdur_p;
                else cycdur{un_indx}=[cycdur{un_indx}; cycdur_p];
                end; clear cycdur_p em cycs perdur
            end; clear el
        end; clear un_indx
    end;
end; clear th_pname

%%%%%%%%%%%%

%Calculating coupling for all channels for the cell
%LFP
binscale=(th_binsize/2):(th_binsize):(360-(th_binsize/2));
for un_indx=1:length(un_chans)
    if isempty(th_coupvect_LFP{un_indx})
        th_coupvect_LFPm(un_indx)=NaN;
        th_coupvect_LFPr(un_indx)=NaN;
        th_coupvect_LFPh(un_indx,:)=NaN;
        th_coupvect_LFPp(un_indx)=NaN;
    else
        th_coupvect_LFPm(un_indx)=rad2deg(circ_mean(deg2rad(th_coupvect_LFP{un_indx}')));
        if th_coupvect_LFPm(un_indx)<0
            th_coupvect_LFPm(un_indx)=th_coupvect_LFPm(un_indx)+360;
        end;
        th_coupvect_LFPr(un_indx)=circ_r(deg2rad(th_coupvect_LFP{un_indx}'));
        th_coupvect_LFPh(un_indx,:)=hist(rad2deg(deg2rad(th_coupvect_LFP{un_indx}')),binscale);
        n=length(th_coupvect_LFP{un_indx});
        r=th_coupvect_LFPr(un_indx);
        eR=n*r;
        z=(eR.^2)/n; P=exp((sqrt(1+(4*n)+(4*(n.^2-eR.^2))))-(1+(2*n)));
        th_coupvect_LFPp(un_indx)=P; clear n P r z eR;
    end
end; clear un_indx;


%PLOTTING STATISTICS
%plotting figures
figindex=0;

for un_indx=1:length(un_chans)
    figindex_T=ceil(un_indx/6);
    if figindex_T~=figindex
        figindex=figindex+1;
        figname=strcat(bs_name, bs_num, '_', bs_exp,'_thetamodulation');
        figure ('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
        plotindex=1;
    else plotindex=((un_indx-((figindex-1)*6))*3)-2;
    end
    if isempty(th_coupvect_LFP{un_indx})
        subplot(6,3,plotindex);
        axis off;
        text (0, 1, 'no data for this cell');
    else
        %textual results
        subplot (6,3,plotindex);
        axis off;
        actualtext= strcat(bs_name, bs_num, bs_exp, 'sh', num2str(un_shnk(un_indx)),'u',num2str(un_nums(un_indx))); text (0, 1, actualtext);
        %number of theta cycles
        actualtext= strcat('Number of theta cycles analysed: ', num2str (th_cyccount(un_indx))); text (0, 0.9, actualtext);
        %number of active theta cycles
        th_acyccount=length(th_cycspkcount{un_indx});th_acycperc=(th_acyccount/th_cyccount(un_indx))*100; actualtext= strcat('of these active: ', num2str (th_acyccount), ' (', num2str(th_acycperc), '%)'); text (0, 0.8, actualtext);
        %%mean number of APs in active theta cycles
        actualtext= strcat('# of spikes/cycle (SD):', num2str(mean(th_cycspkcount{un_indx})), ' (', num2str(std(th_cycspkcount{un_indx})), ')'); text (0, 0.7, actualtext);
        %mean AP frequency
        th_spkrate=length(th_coupvect_LFP{un_indx})/th_time(un_indx); actualtext=strcat('spike frequency under theta:', num2str (th_spkrate), ' s-1'); text (0, 0.6, actualtext);
        %mean vector legth
        actualtext= strcat('mean vector length:', num2str (th_coupvect_LFPr(un_indx))); text (0, 0.5, actualtext);
        %mean phase angle
        actualtext= strcat('mean vector angle:', num2str (th_coupvect_LFPm(un_indx))); text (0, 0.4, actualtext);
        % angular deviation and circular SD
        s=(180/pi)*(sqrt(2*(1-th_coupvect_LFPr(un_indx)))); so=(180/pi)*(sqrt(-2*(log(th_coupvect_LFPr(un_indx)))));
        actualtext=strcat('Angular deviation (s): ', num2str(s)); text (0, 0.3, actualtext);
        actualtext=strcat('Circular standard deviation (so): ', num2str(so)); text (0, 0.2, actualtext);
        %Rayleigh statistics
        n=length(th_coupvect_LFP{un_indx});
        R=n*th_coupvect_LFPr(un_indx); z=R.^2/n; P=exp((sqrt(1+(4*n)+(4*(n.^2-R.^2))))-(1+(2*n)));
        text (0, 0, 'Rayleigh test for the uniformity of distribution');actualtext=strcat ('z=', num2str(z),' P=', num2str(P),'n=', num2str(n)); text (0, 0.1, actualtext);
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
        subplot (6,3,plotindex+1);
        cycbintimes=[repmat(sum(cycdur{un_indx}(:,1),1)./10, 1, th_nbins/2) repmat(sum(cycdur{un_indx}(:,2),1)./10, 1, th_nbins/2)];
        line ([binscale binscale+360], repmat(th_coupvect_LFPh(un_indx,:)./(cycbintimes),1,2));
        x=[0:1:720]; y=(cosd(x+180)+1)*(max(th_coupvect_LFPh(un_indx,:)./max(cycbintimes))).*0.5; line (x,y, 'Color', 'red');
        axis ([0 720 0 max(th_coupvect_LFPh(un_indx,:)./(th_time(un_indx)/20))*1.2]);
        set (gca, 'Xtick', 0:90:720); xlabel ('theta phase (deg)'); ylabel ('spike rate (s-1)')
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
        %phase plot for all spikes
        subplot (6,3,plotindex+2);
        plot (th_coupvect_LFP{un_indx}, [1:1:length(th_coupvect_LFP{un_indx})], 'b. ');
        axis ([0 360 0 length(th_coupvect_LFP{un_indx})]); xlabel ('phase (deg)'); ylabel ('spike #');
         %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.2*pozi(3) 1.2*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
    end
end
clear th_nbins;

plot_write=strcat(bs_name, bs_num, '_', bs_exp, '_THP_Gl_unt.fig');
savefig(plot_write)

%creating saving variables
thcell.units=[un_nums; un_shnk];

thcell.cyccount=th_cyccount;
thcell.spikespercycles=th_cycspkcount;
thcell.pyrlayer=ch_pyr;
thcell.smplrate=ch_rate;
thcell.cycbintime=cycbintimes;
 

thcell.coupvect_lfp=th_coupvect_LFP;
thcell.coupvect_lfp_r=th_coupvect_LFPr;
thcell.coupvect_lfp_m=th_coupvect_LFPm;
thcell.coupvect_lfp_hr=th_coupvect_LFPh./repmat(cycbintimes,(size(th_coupvect_LFPh,1)),1);
thcell.coupvect_lfp_h=th_coupvect_LFPh;
thcell.coupvect_lfp_p=th_coupvect_LFPp;
 
%saving the relevant data
datum = date;
writename=strcat(bs_name, bs_num, '_', bs_exp, '_a_tht_cellcoupling');
save(writename, 'thcell','datum');
            
                     
%This program analyses the gamma oscillatory coupling relative to CSD traces
%multichannel LFP and CSD traces 

ch_n=16;                                              %number of channels
ch_rate=2000;                                         %give the silicon probe sampling rate (in Hz)
ch_pyr=12;                                            %channel number (ventralmost LFP contact is 1) for the pyramidal layer
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
            
            th_pname=strcat(bs_name, bs_num, '_', bs_exp, '_THP.txt');  % theta periods
            
            if exist(th_pname,'file')==2; th_pers=dlmread(th_pname); clear th_pname;
                %reading LFP for gamma to theta cpoupling
                if je==1;
                    lfp_name=strcat(bs_name, bs_num, '_', bs_exp, '_LFP.mat');
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
                csd_name=strcat(bs_name, bs_num, '_', bs_exp, '_CSD.mat');
                load(csd_name);
                if isnan(CSD(je,1))
                    CSDchan=NaN;
                else
                    CSDchan=CSD(je,:)*-1;
                end
                clear CSD csd_name;
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
figname=strcat(bs_name, bs_num, '_', bs_exp, ' gamma modulation by theta');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
subplot(4,4,1);
axis([0 1 0 4]); axis off;text(0,3,strcat(bs_name,bs_num,bs_exp));text(0,2,'gamma CSD');text(0,1,'by theta')
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
        surf([-th_binsize/2:th_binsize:360+th_binsize/2], gm_fraxis, [gm_thtmodulation{je}(:,end) gm_thtmodulation{je} gm_thtmodulation{je}(:,1)], 'EdgeColor', 'interp', 'FaceColor', 'interp'); shading flat; axis tight; colormap hot; caxis([-1*cmax cmax]); view(2);
        axis([0 360 min(gm_fraxis) max(gm_fraxis)]);
        xlabel('theta phase (deg)'); ylabel ('frequency (Hz)'); set(gca, 'YTick',[(min(gm_fraxis)):10:(max(gm_fraxis))]); set(gca, 'XTick',[0:180:360]);
    end
end; clear je cmax

plot_write=strcat(bs_name, bs_num, '_', bs_exp, '_THP_Gamma.fig');
savefig(plot_write)

%gamma amplitude spectra
figname=strcat(bs_name, bs_num, '_', bs_exp, ' gamma amplitude during theta');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
subplot(4,4,1);
axis([0 1 0 4]); axis off;text(0,3,strcat(bs_name,bs_num,bs_exp));text(0,2,'gamma CSD ampl');text(0,1,'during theta')
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

plot_write=strcat(bs_name, bs_num, '_', bs_exp, '_THP_Gamma_ampl.fig');
savefig(plot_write)

%gamma amplitude modulation amplitude spectra
figname=strcat(bs_name, bs_num, bs_exp,' gamma amplitude modulation');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
subplot(4,4,1);
axis([0 1 0 4]); axis off;text(0,3,strcat(bs_name,bs_num,bs_exp));text(0,2,'CSD Z difference');text(0,1,'during theta')

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

plot_write=strcat(bs_name, bs_num, '_', bs_exp, '_THP_Gamma_ampl_mod.fig');
savefig(plot_write)

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
writename=strcat(bs_name, bs_num, '_', bs_exp, '_a_gma_thtCHAR');
save(writename, 'datum', 'gmtht');
 
%This program analyses the gamma oscillatory coupling relative to CSD
%traces (spectral aproach; multichannel LFP and CSD traces) 
clear all; close all;
bs_name='ES'; bs_num='49'; bs_exp='d5';

ch_n=16;                                                            %number of channels
ch_rate=2000;                                                       %give the silicon probe sampling rate (in Hz)

un_chans=[1];                                                      %channels containing the spike times for the units
un_nums=[1];                                                        %unit numbers (in kwik file)
un_shank =[repmat([1],1,1) repmat([2],1,0)];                       %shanks from which unit was isolated

ch_def=[];                                                          %defective channel (one is standing for the ventralmost, and the linear number of channel is required not the identifier
fl_ord=[1];                                                         %file extension numbers belonging to the present cell

%gamma wavelet properties
gm_frsta=200; gm_frfin=15; gm_cfnum=80; gm_wavelet{1}='cmor1-1.5'; gm_plimit=0.05;
gm_frcnt=centfrq(gm_wavelet{1}); gm_scsta=gm_frcnt/(gm_frsta/ch_rate); gm_scfin=gm_frcnt/(gm_frfin/ch_rate); gm_scint=(gm_scfin-gm_scsta)/gm_cfnum; gm_fraxis=gm_frcnt./((gm_scsta:gm_scint:gm_scfin)./ch_rate);
gm_nbins=20;gm_binsize=360/gm_nbins;

for je=1:ch_n
    for un_indx=1:length(un_chans)
        gm_cellcwtspect{je,un_indx}=[];
    end; clear un_indx
end; clear je

for je=1:ch_n;                                                      %CSD channel stepper
    for fl_indx=1:length(fl_ord);                                   %file stepper
        fl_name=strcat(bs_name, bs_num, '_', bs_exp, '_WAW.mat'); %file name
        th_pname=strcat(bs_name, bs_num, '_', bs_exp, '_THP.txt');%theta periods
        
        if exist(th_pname,'file')==2; th_pers=dlmread(th_pname); clear th_pname;
            
            %load the csd traces
            csd_name=strcat(bs_name, bs_num, '_', bs_exp, '_CSD.mat'); load(csd_name);
            
            %check if it is edge and if it is defected channel
            if ~isnan(CSD(je,1))&&isempty(ch_def(ch_def==(je)))
                %invert the CSD so that the sink is down and the source is up
                CSDchan=CSD(je,:)*-1; clear CSD csd_name;
                
                %reading the spike sequence for the cells, transforming it to sample number instead of time
                for un_indx=1:length(un_chans)
                    ch_name=strcat('Gl_unt_', num2str(un_chans(un_indx).','%02d'));
                    ch_load=load(num2str(fl_name),num2str(ch_name));
                    fl_spks{un_indx}=((ch_load.(ch_name).('times')).*ch_rate); clear ch_name ch_load;
                end
                
                %calculate the theta segment wavelet transforms
                for ka=1:size(th_pers,1); %theta period stepper
                    sg_start=(th_pers(ka,1)*ch_rate); sg_end=(th_pers(ka,2)*ch_rate);
                    sg_CSD=CSDchan((round(sg_start)-0.5*ch_rate):(round(sg_end)+0.5*ch_rate));
                    %calculating the wavelet transform for the actual channel of the actual segment.
                    sg_cwt=cwt(sg_CSD, gm_scsta:gm_scint:gm_scfin, gm_wavelet{1});
                    sg_cwt=conj(sg_cwt);
                    sg_cwt=sg_cwt(:,(1+0.5*ch_rate):(1+0.5*ch_rate+ceil(sg_end)-floor(sg_start)));
                    
                    %extract spike triggered CWT for the particular segment for all units and add to a 2D cell array
                    for un_indx=1:length(un_chans)
                        sg_spks=round((fl_spks{un_indx}(fl_spks{un_indx}>sg_start&fl_spks{un_indx}<sg_end)-sg_start+1)');
                        sg_gspect=sg_cwt(:,sg_spks);
                        if isempty(gm_cellcwtspect{je,un_indx})
                            gm_cellcwtspect{je,un_indx}=sg_gspect;
                        else
                            gm_cellcwtspect{je,un_indx}=[gm_cellcwtspect{je,un_indx} sg_gspect];
                        end;
                        clear sg_spks sg_gspect
                    end; clear un_indx
                    clear sg_cwt sg_start sg_end sg_CSD
                end; clear ka th_pers CSDchan
                strcat('CSD_channel_', num2str(je), '_done for all units')
                clear fl_spks
            else
                strcat('CSD_channel_', num2str(je), '_is_not included')
            end
        end
    end; clear fl_indx
end; clear je


%plotting results gamma coupling (during theta)
%calculus
    for un_indx=1:length(un_chans)
        gm_cellcwtspect_r{un_indx}(1:ch_n,1:length(gm_fraxis))=NaN;
        gm_cellcwtspect_m{un_indx}(1:ch_n,1:length(gm_fraxis))=NaN;
        gm_cellcwtspect_p{un_indx}(1:ch_n,1:length(gm_fraxis))=NaN;
        gm_cellcwtspect_z{un_indx}(1:ch_n,1:length(gm_fraxis))=NaN;
        gm_cellcwtspect_sign{un_indx}(1:ch_n,1:length(gm_fraxis))=NaN;
    end

for un_indx=1:length(un_chans)
    %histograms
    figname=strcat(bs_name, bs_num, bs_exp, ' shk', num2str(un_shank(un_indx)),'u',num2str(un_nums(un_indx)));
    figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
    subplot(15,3,1);
    axis off;
    text(0,1,figname);text(0,0.5,'gamma coupling to CSD during theta');text(0,0,'X gamma phase in deg; Y gamma frq in Hz')
    th_spnum=size(gm_cellcwtspect{2,un_indx},2);
    cmax=(th_spnum/20);
    
    for je=2:ch_n-1
        subplot(15,3,49-(je*3))
        set(gca,'FontSize',5)
        if ~isnan(gm_cellcwtspect{je,un_indx})
            for ka=1:size(gm_cellcwtspect{je,un_indx},1)
                gm_phasehist{je,un_indx}(ka,:)=repmat((hist(rad2deg(angle(gm_cellcwtspect{je,un_indx}(ka,:)))+180,[gm_binsize/2:gm_binsize:360-gm_binsize/2])),1,2);
            end; clear ka
            surf([-gm_binsize/2:gm_binsize:720+gm_binsize/2], gm_fraxis, [gm_phasehist{je,un_indx}(:,end) gm_phasehist{je,un_indx} gm_phasehist{je,un_indx}(:,1)],'EdgeColor', 'interp', 'FaceColor', 'interp'); shading flat; colormap bone; caxis([0 2*cmax]);view(2);
            axis([0 720 min(gm_fraxis) max(gm_fraxis)]);
            set(gca, 'YTick',[(min(gm_fraxis)):10:(max(gm_fraxis))]); set(gca, 'XTick',[0:180:720]);
        else gm_phasehist{je}=NaN;
        end
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
    end; clear je cmax
    
    %calculus
    for je=2:ch_n-1
        if ~isnan(gm_cellcwtspect{je,un_indx});
            gm_cellcwtspect_r{un_indx}(je,:)=circ_r(angle(gm_cellcwtspect{je,un_indx}),[],[],2);
            gm_cellcwtspect_m{un_indx}(je,:)=(rad2deg(circ_mean(angle(gm_cellcwtspect{je,un_indx}),[],2)))+180;
            for ka=1:size(gm_cellcwtspect{je,un_indx},1)
                [pe,ze]=circ_rtest(angle(gm_cellcwtspect{je,un_indx}(ka,:)));
                gm_cellcwtspect_p{un_indx}(je,ka)=pe; gm_cellcwtspect_z{un_indx}(je,ka)=ze;
                if pe<=gm_plimit
                    gm_cellcwtspect_sign{un_indx}(je,ka)=1;
                end
                clear pe ze;
            end; clear ka;
        end
    end
    
    %gamma  coupling (during theta); coupling strength spectra
    for je=2:ch_n-1
        subplot(15,3,(50-(je*3)))
        set(gca,'FontSize',5)
        if ~isnan(gm_cellcwtspect{je,un_indx});
            line(gm_cellcwtspect_r{un_indx}(je,:),gm_fraxis);
            axis([0 max(max(gm_cellcwtspect_r{un_indx})*1.1) (min(gm_fraxis)) (max(gm_fraxis)) ]);
            line((gm_cellcwtspect_r{un_indx}(je,:).*gm_cellcwtspect_sign{un_indx}(je,:)), gm_fraxis,'color','red');
            set(gca, 'YTick',[(min(gm_fraxis)):10:(max(gm_fraxis))]);
        end
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
    end; clear je;
    
    %gamma  coupling (during theta); coupling phase spectra
    for je=2:ch_n-1
        subplot(15,3,(51-(je*3)))
        set(gca,'FontSize',5)
        if ~isnan(gm_cellcwtspect{je,un_indx});
            line(gm_cellcwtspect_m{un_indx}(je,:), gm_fraxis);
            line((gm_cellcwtspect_m{un_indx}(je,:))+360, gm_fraxis);
            line((gm_cellcwtspect_m{un_indx}(je,:).*gm_cellcwtspect_sign{un_indx}(je,:)), gm_fraxis,'color','red');
            line((((gm_cellcwtspect_m{un_indx}(je,:))+360).*gm_cellcwtspect_sign{un_indx}(je,:)), gm_fraxis,'color','red');
            axis([0 720 (min(gm_fraxis)) (max(gm_fraxis))]);
            set(gca, 'YTick',[(min(gm_fraxis)):10:(max(gm_fraxis))]); set(gca, 'XTick',[0:90:720]);
        end; clear je;
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
    end
end

plot_write=strcat(bs_name, bs_num, '_', bs_exp, '_Gamma_Gl_unt.fig');
savefig(plot_write)

%cleanup
clear fl_indx fl_ord ch_def ch_spks fl_name ans figname un_indx       

gmtht.gm_fraxis=gm_fraxis;

gmtht.gm_sptrig_spikenum=th_spnum;
gmtht.gm_sptrig_rawcwt=gm_cellcwtspect;
gmtht.gm_sptrig_phasespect=gm_cellcwtspect_m;
gmtht.gm_sptrig_rspect=gm_cellcwtspect_r;
gmtht.gm_sptrig_zspect=gm_cellcwtspect_z;
gmtht.gm_sptrig_Pspect=gm_cellcwtspect_p;
gmtht.gm_sptrig_signif=gm_cellcwtspect_sign;
gmtht.gm_sptrig_histograms=gm_phasehist;
gmtht.gm_sptrig_plimit=gm_plimit;
gmtht.ch_samplrate=ch_rate;

datum=date

%saving the relevant data
writename=strcat(bs_name, bs_num, '_', bs_exp, '_a_gma_thtCOUP');
save(writename,'gmtht', 'datum','-v7.3');


%---------Changes---------

%2018-01-26
%action: corrected in line 57 to include the last points of the cwt segment
%reason: error (script terminated)
%'round' changed to 'floor' and 'ceil' to expend the cwt segment's duration
 
    
           
 
      
        
       
       
        
    
