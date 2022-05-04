function [] = gmAmplMod_Plt(bs, ch_n, gmtht)
%gamma amplitude modulation amplitude spectra
figname=strcat(bs.name, bs.num, bs.unit,' gamma amplitude modulation');
figure('Name', figname, 'Numbertitle', 'off', 'PaperOrientation', 'landscape', 'PaperType', 'A4')
subplot(4,4,1);
axis([0 1 0 4]); 
axis off;
text(0,3,strcat(bs.name,bs.num,bs.unit));
text(0,2,'CSD Z difference');
text(0,1,'during theta')

for je=2:ch_n-1
    if ~isnan(gmtht.gm_osc_thtmod{je})
        gmtht.gm_osc_thtmodamp{je}=max(gmtht.gm_osc_thtmod{je}')-min(gmtht.gm_osc_thtmod{je}');
        smax_t(je)=max(gmtht.gm_osc_thtmodamp{je});
    else
        gmtht.gm_osc_thtmodamp{je}=NaN;
        smax_t(je)=NaN;
    end
end; clear je cmax
smax=max(smax_t);clear je;

for je=2:ch_n-1
    subplot(4,4,17-je)
    if ~isnan(gmtht.gm_osc_thtmodamp{je})
        xlabel('Z difference (max-min)'); 
        ylabel ('frequency (Hz)');
        set(gca, 'YTick',[(min(gmtht.gm_fraxis)):10:(max(gmtht.gm_fraxis))]);
        axis([0 smax*1.1 (min(gmtht.gm_fraxis)) (max(gmtht.gm_fraxis))]);
        line(gmtht.gm_osc_thtmodamp{je}, gmtht.gm_fraxis)
    end
end 
clear je cmax
end