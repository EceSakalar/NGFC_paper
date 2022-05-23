%This program removes a selected gamma component from the CSD traces based on correlations
%of significant (high amplitude) gamma cycles from filtered CSD 

clear all
%close all

%basic data entry
bs_name='B'; bs_num='182'; bs_exp='a'; bs_typ='sx'; 
fl_ord=[1];
ch_rate=2000;                                       %give the sampling rate of the channels (CSD and LFP original)

%gamma to be removed properties 
gm_comp = 'CA1mid';                                                        %name of gamma component
gm_chan = 5;                                                               %CSD channel for the actual gamma component
gm_cfrq = [43.9963 136.752];                                                 %band pass filter corner frequencies (low and high, included in the boundary)
gm_ncoeff=1024;
gm_rate = ch_rate;

%creating gamma filter
gm_filter=dfilt.dffir(fir1(gm_ncoeff, [2.*(gm_cfrq(1)./gm_rate) 2.*(gm_cfrq(2)./gm_rate)],'bandpass', gausswin(gm_ncoeff+1)));

%---CALCULUS---concatenating files
for fl_indx=1:length(fl_ord); %FILE STEPPER
    
    %read the CSD for gamma; invert the CSD channel so that sink is downwards (-)
    csd_name=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'CSD.mat');
    load(csd_name);
    CSD_file = CSD(:,:)*-1; clear CSD csd_name;
    
    %read theta periods
    th_pname=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)),'THP.txt');
    if exist(th_pname,'file')==2;
        thPER_file=dlmread(th_pname); clear th_pname;
    else
        thPER_file=[]; clear th_pname;
    end
    
    %accumulating CSD and LFP traces, theta segments and file segmentation point 
    if fl_indx==1;
        CSD_cum = CSD_file;
        flSEG = length(CSD_file);
        thPER = thPER_file;
    else
        CSD_cum = [CSD_cum CSD_file];
        thLFP_orig = [thLFP_orig thLFP_file];
        flSEG(fl_indx)=flSEG(fl_indx-1)+length(CSD_file);
        thePER = [thPER; thPER_file+(flSEG(fl_indx-1)/ch_rate)];
    end; clear  CSD_file thLFP_file thPER_file
end; clear fl_indx

%---CALCULUS----------------------------filtering CSD for gamma
CSD_filt=CSD_cum;
for je=1:size(CSD_cum,1)
    if ~isnan(CSD_cum(je,1))
        CSD_filt(je,:)= filtfilt(gm_filter.Numerator,1,CSD_cum(je,:));
    end
end; clear je

%---CALCULUS----------------------------gamma cycle inclusion

%generating a file with peaks and troughs extracted and sorted along the time
[gmCSD_pVAL,gmCSD_pPOS] = findpeaks(CSD_filt(gm_chan,:));
[gmCSD_tVAL,gmCSD_tPOS] = findpeaks(CSD_filt(gm_chan,:)*-1);
gmCSD_ptTMP = [[gmCSD_pPOS;gmCSD_pVAL;zeros(size(gmCSD_pPOS))+1] [gmCSD_tPOS;gmCSD_tVAL*-1;zeros(size(gmCSD_tPOS))-1]];
[~,sIND] = sort(gmCSD_ptTMP(1,:));
CSD_gmPT = gmCSD_ptTMP(:,sIND);
%calculate the amplitude (as SD) of the cycles (this will be the fourth row)
for je=2:(size(CSD_gmPT,2)-1)
    CSD_gmPT(4,je) = std(CSD_filt(gm_chan,CSD_gmPT(1,je-1):CSD_gmPT(1,je+1)));
end; clear je
clear gmCSD_pVAL gmCSD_tVAL gmCSD_pPOS gmCSD_tPOS sIND gmCSD_ptTMP

%finding threshold
[binV,binP]=hist(CSD_gmPT(4,:),100);
gmCSD_THRSHLD=min(binP(binV==max(binV)));

%generate a gamma selector trace where only cycles above threshold are included
CSD_sel=NaN(1,size(CSD_filt,2));
for je=1:size(CSD_gmPT,2)
    if CSD_gmPT(4,je)>gmCSD_THRSHLD;
        CSD_sel(CSD_gmPT(1,je-1):CSD_gmPT(1,je+1))=1;
    end;
end; clear je
dwn_samplerat=floor((nansum(CSD_sel)/100000));

%calculate the gamma projections on all channels
CSD_gmproj=NaN(size(CSD_filt));

figname=strcat(bs_name, bs_num, bs_exp, gm_comp, 'correlations');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'portrait', 'PaperType', 'A4')
plotnum=sum(~isnan(CSD_cum(:,1)));
ka=1;

for je=1:size(CSD_filt,1)
    if isnan(CSD_filt(je,1))
        CSD_load(je)=NaN;
        CSD_pval(je)=NaN;
    else
        LModell=fitlm(downsample(((CSD_filt(gm_chan,:)).*CSD_sel),dwn_samplerat),downsample((CSD_filt(je,:).*CSD_sel),dwn_samplerat));
        subplot(plotnum,2,((plotnum-ka)*2)+1)
        plot(LModell)
        legend off
        subplot(plotnum,2,((plotnum-ka)*2)+2)
        axis off
        text(0,0.75,strcat('slope:',num2str(LModell.Coefficients.Estimate(2))));
        text(0,0.25,strcat('p-value:',num2str(LModell.Coefficients.pValue(2))));
        CSD_pval(je)=LModell.Coefficients.pValue(2)>0.05;
        
        ka=ka+1;
        if LModell.Coefficients.pValue(2)>0.05;
            CSD_load(je)=0;
        else
            CSD_load(je)=LModell.Coefficients.Estimate(2);
        end
        %calculate projection
        CSD_gmproj(je,:)=CSD_filt(gm_chan,:).*(CSD_load(je));
        %expand the graph
        pozi=get(gca,'Position'); pozi=[pozi(1) pozi(2) 1.1*pozi(3) 1.1*pozi(4)]; set(gca,'Position',pozi)
        set(gca, 'Clipping', 'off')
        
    end
end; clear je

%correction (gamma removal)
CSD_clear = CSD_cum - CSD_gmproj;

%writing the cleared CSD file
for fl_indx=1:length(fl_ord); %FILE STEPPER
    if fl_indx==1
        CSDc=(CSD_clear(:,1:flSEG(fl_indx))).*-1;
    else
        CSDc=(CSD_clear(:,flSEG(fl_indx-1):flSEG(fl_indx))).*-1;
    end
    fl_write=strcat(bs_name, bs_num, bs_exp, num2str(fl_ord(fl_indx)), 'CSDc1');
    save(fl_write, 'CSDc');
    clear CSDc fl_write
end; clear fl_indx

%saving data
clean.files=fl_ord;
clean.gmrate=gm_rate;
clean.gmcomp=gm_comp;
clean.gmchan=gm_chan;
clean.gmfreq=gm_cfrq;
clean.gmncoeff=gm_ncoeff;
clean.CSDslopes=CSD_load;
clean.CSDpvalues=CSD_pval;
datum=date;
%saving the relevant data
writename=strcat(bs_name, bs_num, bs_exp, '_a_gmacleaning1');
save(writename, 'datum', 'clean');

    
    

