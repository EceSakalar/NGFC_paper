function [results, properties] = LFP_CSD_SWR_plt(LFP_rips, CSD_rips, rp, ch)
%%plotting results
%%LFP ripple frequency profiles for 
figname=strcat(bs.name, bs.num, bs.exp,' LFP ripple amplitude');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
cmax=1; 
cmin=1;
for je=1:ch.n
    subplot(4,ch.n/4,ch.n-je+1)
    if ~isnan(LFP_rips{je})
        LFP_rips{je}=LFP_rips{je}./SWR_no;
        LFP_calib=repmat(mean(LFP_rips{je}(:,1:n_tbins(1)),2),1,sum(n_tbins));
        LFP_ripsN{je}=LFP_rips{je}./LFP_calib;
        cmin_A=min(min(LFP_ripsN{je}));
        cmax_A=max(max(LFP_ripsN{je}));
        if cmin_A<cmin 
            cmin=cmin_A; 
        end
        if cmax_A>cmax 
            cmax=cmax_A; 
        end
    end
end 
cmin=cmin*0.9; 
cmax=cmax*1.1; 
clear cmax_A cmin_A;


for je=1:ch.n
    subplot(4,ch.n/4,ch.n-je+1)
    if ~isnan(LFP_ripsN{je})
        surf([(-1*(sum(n_tbins)/2))+0.5:1:((sum(n_tbins)/2)-0.5)], rp.fraxis, LFP_ripsN{je}, 'EdgeColor', 'none', 'FaceColor', 'flat'); 
        shading flat; 
        axis tight; 
        colormap jet; 
        caxis([cmin cmax]); 
        view(2);
        xlabel('SWR phase (normalised time)'); 
        ylabel ('frequency (Hz)'); 
        set(gca, 'YTick',[(min(rp.fraxis)):20:(max(rp.fraxis))]);
    end
end
LFP_cmin=cmin; 
LFP_cmax=cmax; 
clear cmin cmax

%%CSD ripple frequency profiles for 
figname=strcat(bs.name, bs.num, bs.exp,' CSD ripple amplitude');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
cmax=1; 
cmin=1;
for je=ch.step+1:ch.n-ch.step
    subplot(4,ch.n/4,ch.n-je+1)
    if ~isnan(CSD_rips{je})
        CSD_rips{je}=CSD_rips{je}./SWR_no;
        CSD_calib=repmat(mean(CSD_rips{je}(:,1:n_tbins(1)),2),1,sum(n_tbins));
        CSD_ripsN{je}=CSD_rips{je}./CSD_calib;
        cmin_A=min(min(CSD_ripsN{je}));
        cmax_A=max(max(CSD_ripsN{je}));
        if cmin_A<cmin 
            cmin=cmin_A; 
        end
        if cmax_A>cmax 
            cmax=cmax_A; 
        end
    end
end
cmin=cmin*0.9; 
cmax=cmax*1.1; 
clear cmax_A cmin_A


for je=ch.step+1:ch.n-ch.step
    subplot(4,ch_n/4,ch_n-je+1)
    if ~isnan(CSD_ripsN{je})
        surf([(-1*(sum(n_tbins)/2))+0.5:1:((sum(n_tbins)/2)-0.5)], rp.fraxis, CSD_ripsN{je}, 'EdgeColor', 'none', 'FaceColor', 'flat'); 
        shading flat; 
        axis tight; 
        colormap jet; 
        caxis([cmin cmax]); 
        view(2);
        xlabel('SWR phase (normalised time)'); 
        ylabel ('frequency (Hz)'); 
        set(gca, 'YTick',[(min(rp.fraxis)):20:(max(rp.fraxis))]);
    end
end
CSD_cmin=cmin; 
CSD_cmax=cmax; 
clear cmin cmax

%%SWR figure
figname=strcat(bs.name, bs.num, bs.exp,' SW averages');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')

%LFP SW
LFP_segavg=LFP_segavg./SWR_no;
cmax=max(abs([1.1*min(min(LFP_segavg)) 1.1*max(max(LFP_segavg))]));
subplot(1,2,1)
surf([(-1*(sum(n_tbins)/2))+0.5:1:((sum(n_tbins)/2)+0.5)],[0.5:1:ch.n+0.5], [[LFP_segavg zeros(ch.n,1)];zeros(1,sum(n_tbins)+1)], 'EdgeColor', 'none', 'FaceColor', 'flat'); 
shading flat; 
axis tight; 
colormap jet; 
caxis([-cmax cmax]); 
view(2);
xlabel('SWR phase (normalised time)'); 
ylabel ('contact#'); 
set(gca, 'YTick',[1:1:ch.n]);
LFPavg_cmax=cmax; 
clear cmax

%CSD SW
CSD_segavg=CSD_segavg./SWR_no;
cmax=max(abs([1.1*min(min(CSD_segavg)) 1.1*max(max(CSD_segavg))]));
subplot(1,2,2)
surf([(-1*(sum(n_tbins)/2))+0.5:1:((sum(n_tbins)/2)+0.5)],[0.5:1:ch_n+0.5], [[CSD_segavg zeros(ch.n,1)];zeros(1,sum(n_tbins)+1)], 'EdgeColor', 'none', 'FaceColor', 'flat'); 
shading flat; 
axis tight; 
colormap jet; 
caxis([-cmax cmax]); 
view(2);
xlabel('SWR phase (normalised time)'); 
ylabel ('contact#'); 
set(gca, 'YTick',[1:1:ch_n]);
CSDavg_cmax=cmax; 
clear cmax

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
end

